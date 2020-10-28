#include "DensitySubGridCreator.hpp"
#include "GridWriter.hpp"
#include "HydroDensitySubGrid.hpp"
#include "Scheduler.hpp"

/*******************************************************************************

This test program illustrates the usage of the task-based algorithm on a
distributed grid.

The idea is conceptually very simple: each cell in the grid is initialized with
a single meaningful value: its number density. During the program, we will
perform the following two tasks on each cell (in order):
 - we first copy the number density as it is into a new variable, "test density"
 - we then sum the values of "test density" in all cells that share a face with
   a cell and store the result in another variable "test neighbour density sum".

The first task is trivial and can be achieved by simply looping over all cells.
The second task is more complex, as it requires a careful loop over pairs of
cells, whereby a pair can span the boundary between two parts of the grid, or
the boundary of the simulation box. During this loop, we need to make sure to
only execute the second task once for each pair. Note that to have a consistent
result, task 2 can only be executed if task 1 has been executed for the cell and
its neighbour. This defines a causal dependency. Additionally, task 2 cannot be
executed in parallel for more than one cell pair that contains the same cell,
since then race conditions could arise when updating the test neighbour density
sum. This defines an explicit task dependency.

The code below executes the tasks in two very distinct ways. The first way is a
serial, naive way, whereby first task 1 and then task 2 is executed for all
subgrids by simply looping over all subgrids:
 - during the first loop, we simply call HydroDensitySubGrid::set_test_density()
   for each subgrid.
 - during the second loop, we first call
   HydroDensitySubGrid::inner_test_neighbour_density_sum_sweep() to get the
   contributions for all cell pairs that are entirely contained within the
   subgrid. We then in turn determine what the neighbouring subgrids are along
   the x, y and z axis (in the positive direction), and call
   HydroDensitySubGrid::outer_test_neighbour_density_sum_sweep() to account for
   the cell pairs that cross the boundary between these two subgrids. If the
   neighbouring subgrid does not exist, we instead call a special function,
   HydroDensitySubGrid::outer_ghost_test_neighbour_density_sum_sweep() to deal
   with the box boundary. This last function is also called if the subgrid
   neighbour in the negative direction does not exist (if it exists, task 2 is
   handled from the perspective of that subgrid and cannot be done again).

The second way uses a task-based parallel algorithm. This algorithm still uses
the same functions as above, but now these functions are executed in parallel
by multiple threads. The order of the tasks is no longer predictable, but we
still impose the causal and task dependency for each cell by making sure that
 - set_test_density() is always called on a subgrid before that subgrid is used
   by any of the other tasks.
 - only one of the sweep() tasks is executed concurrently that involves the same
   subgrid. In other words, if subgrid A is doing its inner_sweep(), then no
   outer_sweep() task can be executed that involves subgrid A, not even if
   subgrid A is only a neighbour participating in the task.
The only thing required to make this work, is some code that creates the
corresponding tasks and sets the dependencies correctly. The parallel execution
(respecting these dependencies) is completely handled by a generic algorithm
that does not care about any of the details of the subgrids or the tasks.

*******************************************************************************/

/**
 * @brief Example initial condition.
 */
class TestDensityFunction : public DensityFunction {
public:
  /**
   * @brief This function is called exactly once for each cell in the simulation
   * at the start of the simulation.
   *
   * We use it to set the number density in each cell to a constant value.
   */
  virtual DensityValues operator()(const Cell &cell) {
    DensityValues values;
    values.set_number_density(1.);
    return values;
  }
};

/**
 * @brief Main program routine.
 *
 * This is the code that gets executed when the program is called.
 */
int main() {

  /// METHOD 1
  {
    // create the initial condition object
    TestDensityFunction test_function;
    // create the distributed grid
    // the simulation box that spans the cube with corners [0,0,0] and [1,1,1]
    // we use a 64x64x64 grid subdivided into 8x8x8 subgrids
    // the x and y boundaries are assumed to be period, the z boundary is not
    DensitySubGridCreator<HydroDensitySubGrid> grid(
        Box<>(0., 1.), 64, 8, CoordinateVector<bool>(true, true, false));
    // set up and initialize the cells and subgrids in the distributed grid
    grid.initialize(test_function);

    // task 1: loop over all subgrids
    for (auto gridit = grid.begin(); gridit != grid.all_end(); ++gridit) {
      // execute task 1 for all cells in this subgrid
      (*gridit).set_test_density();
    }

    // we initialize a boundary condition to also illustrate how this is used
    // note that the boundary condition is not actually used here, since all
    // box boundaries are periodic
    InflowHydroBoundary boundary;
    // task 2: loop over all subgrids
    for (auto gridit = grid.begin(); gridit != grid.all_end(); ++gridit) {
      // execute task 2 for all cell pairs within this subgrid
      (*gridit).inner_test_neighbour_density_sum_sweep();
      // execute task 2 for all cell pairs across the boundary between this
      // and the next subgrid in the positive x direction
      uint_fast32_t ngbx = (*gridit).get_neighbour(TRAVELDIRECTION_FACE_X_P);
      if (ngbx == NEIGHBOUR_OUTSIDE) {
        // if the positive x direction is the box boundary, use the boundary
        // condition instead
        (*gridit).outer_ghost_test_neighbour_density_sum_sweep(
            TRAVELDIRECTION_FACE_X_P, boundary);
      } else {
        (*gridit).outer_test_neighbour_density_sum_sweep(
            TRAVELDIRECTION_FACE_X_P, *grid.get_subgrid(ngbx));
      }
      // also check if the negative x direction is a box boundary
      ngbx = (*gridit).get_neighbour(TRAVELDIRECTION_FACE_X_N);
      if (ngbx == NEIGHBOUR_OUTSIDE) {
        // yes: use the boundary condition
        (*gridit).outer_ghost_test_neighbour_density_sum_sweep(
            TRAVELDIRECTION_FACE_X_N, boundary);
      }
      // now repeat for y
      uint_fast32_t ngby = (*gridit).get_neighbour(TRAVELDIRECTION_FACE_Y_P);
      if (ngby == NEIGHBOUR_OUTSIDE) {
        (*gridit).outer_ghost_test_neighbour_density_sum_sweep(
            TRAVELDIRECTION_FACE_Y_P, boundary);
      } else {
        (*gridit).outer_test_neighbour_density_sum_sweep(
            TRAVELDIRECTION_FACE_Y_P, *grid.get_subgrid(ngby));
      }
      ngby = (*gridit).get_neighbour(TRAVELDIRECTION_FACE_Y_N);
      if (ngby == NEIGHBOUR_OUTSIDE) {
        (*gridit).outer_ghost_test_neighbour_density_sum_sweep(
            TRAVELDIRECTION_FACE_Y_N, boundary);
      }
      // and repeat again for z
      uint_fast32_t ngbz = (*gridit).get_neighbour(TRAVELDIRECTION_FACE_Z_P);
      if (ngbz == NEIGHBOUR_OUTSIDE) {
        (*gridit).outer_ghost_test_neighbour_density_sum_sweep(
            TRAVELDIRECTION_FACE_Z_P, boundary);
      } else {
        (*gridit).outer_test_neighbour_density_sum_sweep(
            TRAVELDIRECTION_FACE_Z_P, *grid.get_subgrid(ngbz));
      }
      ngbz = (*gridit).get_neighbour(TRAVELDIRECTION_FACE_Z_N);
      if (ngbz == NEIGHBOUR_OUTSIDE) {
        (*gridit).outer_ghost_test_neighbour_density_sum_sweep(
            TRAVELDIRECTION_FACE_Z_N, boundary);
      }
    }

    // create a grid output object that writes output files of the form
    //  testXXX.hdf5, where XXX is a three digit counter value
    GridWriter writer("test");
    // write output file test000.hdf5 with the current contents of the grid
    writer.write(grid, 0);
  }

  /// METHOD 2
  {

    // create the initial condition object
    TestDensityFunction test_function;
    // create the distributed grid
    // the simulation box that spans the cube with corners [0,0,0] and [1,1,1]
    // we use a 64x64x64 grid subdivided into 8x8x8 subgrids
    // the x and y boundaries are assumed to be period, the z boundary is not
    DensitySubGridCreator<HydroDensitySubGrid> grid(
        Box<>(0., 1.), 64, 8, CoordinateVector<bool>(true, true, false));
    // set up and initialize the cells and subgrids in the distributed grid
    grid.initialize(test_function);

    // create a scheduler to run the tasks in parallel
    Scheduler sched(
        // number of tasks:
        //  - 1 (test density) + 4 (test neighbour density sum) per subgrid
        //  - (potentially) 1 additional task for each subgrid that borders a
        //    box boundary in a negative direction
        5 * 8 * 8 * 8 + 3 * 8 * 8,
        // number of queues, set to the same number as threads we will use
        4,
        // number of tasks that can be stored in the queue for a single thread
        // if we provide enough space to store all density tasks, then even the
        // single thread case will have sufficient space
        8 * 8 * 8,
        // number of tasks that can be stored in the queue that is shared among
        // the threads. Not used in this case.
        0);

    // add the tasks
    // just like in method 1, we do this in two stages
    //  - first, we create all test density tasks
    //  - second, we create the test neighbour density sum tasks
    // we need the indices of the first tasks to link them as dependencies
    // for the second task. That is why we need to do two loops.

    // we set up an atomic counter for the tasks that are queued
    // this counter is used during parallel execution to decide when to stop
    // looking for tasks
    // atomic means that only one thread can access the counter at any given
    // time, so that all thread always agree on the value of the counter
    AtomicValue<uint_fast32_t> number_of_tasks;
    // loop over the subgrids
    for (auto gridit = grid.begin(); gridit != grid.all_end(); ++gridit) {
      // store the index of the subgrid for convenience
      const uint_fast32_t igrid = gridit.get_index();
      // create a new task and store it in two variables:
      //  - task contains a reference to the actual Task object
      //  - test_density_task contains the index of the task. This index is
      //    used internally to link dependencies
      uint_fast32_t test_density_task;
      Task &task = sched.get_free_task(test_density_task);
      // set the type of the task
      // this type determines what action is performed during parallel
      // execution
      task.set_type(TASKTYPE_TEST_DENSITY);
      // set the dependency for this task: the subgrid
      // the task-based algorithm will make sure that the task gets exclusive
      // access to the dependency during parallel execution
      task.set_dependency((*gridit).get_dependency());
      // store the subgrid information in the task
      // this is used by the task to figure out which subgrid to act on
      task.set_subgrid(igrid);
      // store the index of the task in the subgrid so that we can access it
      // during the second loop (this is the only reason we store the index)
      (*gridit).set_hydro_task(0, test_density_task);
      // test density tasks have no task dependencies, so we can immediately
      // put them in the queue
      // we put the task in the queue for the thread that "owns" the subgrid
      // this is the last thread that had exclusive access to this subgrid
      // this is done for efficiency reasons: a thread that reuses the same
      // subgrid consecutively better exploits memory caches
      // initially, the subgrids are distributed uniformly accross the
      // available threads, so that each thread roughly owns the same number of
      // subgrids
      // during parallel execution, subgrid ownership changes to optimally
      // use the available resources
      sched.enqueue(number_of_tasks, test_density_task,
                    (*gridit).get_owning_thread());
    }

    // loop again over all subgrids
    // note that in this case, tasks are not queued yet, since all test
    // neighbour density sum tasks depend on the test density tasks created
    // in the first loop
    for (auto gridit = grid.begin(); gridit != grid.all_end(); ++gridit) {
      // the first bit is similar to loop 1:
      // we create a single task to compute the sums for cell pairs that are
      // entirely contained within a single subgrid
      const uint_fast32_t igrid = gridit.get_index();
      uint_fast32_t test_neighbour_density_sum_task;
      Task &internal_task =
          sched.get_free_task(test_neighbour_density_sum_task);
      // we use a different task type for a different action
      internal_task.set_type(TASKTYPE_TEST_NEIGHBOUR_SWEEP_INTERNAL);
      internal_task.set_dependency((*gridit).get_dependency());
      internal_task.set_subgrid(igrid);
      // this bit is new: add a dependency link between the test density task
      // and the test neighbour density sum task for this subgrid
      // the latter can now no longer be executed before the former has
      // completed
      sched.add_link((*gridit).get_hydro_task(0),
                     test_neighbour_density_sum_task);

      // now set up additional tasks for cell pairs that span two subgrids
      // the code below is similar in structure to the equivalent part of
      // method 1, except that now we need to create a new task for each pair
      // there is some unnecessary code duplication here, so it looks more
      // scary than it actually is

      // create a variable to store the index of new tasks in (we will reuse
      // this variable a few times below)
      uint_fast32_t ngb_task;

      // get the neighbouring subgrid in the positive x direction
      // this is exactly the same as for method 1
      uint_fast32_t ngbx = (*gridit).get_neighbour(TRAVELDIRECTION_FACE_X_P);
      // decide whether the neighbour exists or not
      if (ngbx == NEIGHBOUR_OUTSIDE) {
        // no: we are dealing with a box boundary
        // create a TEST_GHOST_SWEEP task that only depends on this subgrid and
        // this subgrid's test density task
        Task &external_task = sched.get_free_task(ngb_task);
        external_task.set_type(TASKTYPE_TEST_GHOST_SWEEP);
        external_task.set_dependency((*gridit).get_dependency());
        external_task.set_subgrid(igrid);
        external_task.set_interaction_direction(TRAVELDIRECTION_FACE_X_P);
        sched.add_link((*gridit).get_hydro_task(0), ngb_task);
      } else {
        // yes: we are dealing with an actual subgrid pair
        // create a task that involves two subgrids
        // creating the task is similar to before
        Task &external_task = sched.get_free_task(ngb_task);
        external_task.set_type(TASKTYPE_TEST_NEIGHBOUR_SWEEP_EXTERNAL);
        // now set the dependencies: this task requires exclusive access to
        // both subgrids that are involved
        // we have to be careful here to avoid a concurrency issue known as the
        // dining philosphers' problem: if multiple tasks share the same two
        // dependencies and these dependencies are ordered in an arbitrary way,
        // then the algorithm can reach a state where two threads keep locking
        // and unlocking the same variables over and over again
        // in practice, the algorithm will then deadlock and never finish
        // we can avoid this by always logically ordering dependencies, in this
        // case using the unique subgrid index
        if (igrid < ngbx) {
          external_task.set_dependency((*gridit).get_dependency());
          external_task.set_extra_dependency(
              (*grid.get_subgrid(ngbx)).get_dependency());
        } else {
          external_task.set_dependency(
              (*grid.get_subgrid(ngbx)).get_dependency());
          external_task.set_extra_dependency((*gridit).get_dependency());
        }
        // we now need to store information about two subgrids
        // within CMacIonize, there is however only one field for subgrid
        // information. Luckily, there is another field (buffer) that is not
        // used in this case and that can be used instead.
        external_task.set_subgrid(igrid);
        external_task.set_buffer(ngbx);
        // lastly, we need to store the direction in which the two subgrids
        // touch. CMacIonize also has a field that is not exactly meant for
        // this, but can easily be repurposed for it
        external_task.set_interaction_direction(TRAVELDIRECTION_FACE_X_P);
        // finally, we add the task dependencies
        // there are two in this case: one for this subgrid's test density
        // task and one for that of the neighbour subgrid
        sched.add_link((*gridit).get_hydro_task(0), ngb_task);
        sched.add_link((*grid.get_subgrid(ngbx)).get_hydro_task(0), ngb_task);
      }
      // now repeat for the negative x direction in case that is a box boundary
      ngbx = (*gridit).get_neighbour(TRAVELDIRECTION_FACE_X_N);
      if (ngbx == NEIGHBOUR_OUTSIDE) {
        Task &external_task = sched.get_free_task(ngb_task);
        external_task.set_type(TASKTYPE_TEST_GHOST_SWEEP);
        external_task.set_dependency((*gridit).get_dependency());
        external_task.set_subgrid(igrid);
        external_task.set_interaction_direction(TRAVELDIRECTION_FACE_X_N);
        sched.add_link((*gridit).get_hydro_task(0), ngb_task);
      }

      // and repeat for y
      uint_fast32_t ngby = (*gridit).get_neighbour(TRAVELDIRECTION_FACE_Y_P);
      if (ngby == NEIGHBOUR_OUTSIDE) {
        Task &external_task = sched.get_free_task(ngb_task);
        external_task.set_type(TASKTYPE_TEST_GHOST_SWEEP);
        external_task.set_dependency((*gridit).get_dependency());
        external_task.set_subgrid(igrid);
        external_task.set_interaction_direction(TRAVELDIRECTION_FACE_Y_P);
        sched.add_link((*gridit).get_hydro_task(0), ngb_task);
      } else {
        Task &external_task = sched.get_free_task(ngb_task);
        external_task.set_type(TASKTYPE_TEST_NEIGHBOUR_SWEEP_EXTERNAL);
        if (igrid < ngby) {
          external_task.set_dependency((*gridit).get_dependency());
          external_task.set_extra_dependency(
              (*grid.get_subgrid(ngby)).get_dependency());
        } else {
          external_task.set_dependency(
              (*grid.get_subgrid(ngby)).get_dependency());
          external_task.set_extra_dependency((*gridit).get_dependency());
        }
        external_task.set_subgrid(igrid);
        external_task.set_buffer(ngby);
        external_task.set_interaction_direction(TRAVELDIRECTION_FACE_Y_P);
        sched.add_link((*gridit).get_hydro_task(0), ngb_task);
        sched.add_link((*grid.get_subgrid(ngby)).get_hydro_task(0), ngb_task);
      }
      ngby = (*gridit).get_neighbour(TRAVELDIRECTION_FACE_Y_N);
      if (ngby == NEIGHBOUR_OUTSIDE) {
        Task &external_task = sched.get_free_task(ngb_task);
        external_task.set_type(TASKTYPE_TEST_GHOST_SWEEP);
        external_task.set_dependency((*gridit).get_dependency());
        external_task.set_subgrid(igrid);
        external_task.set_interaction_direction(TRAVELDIRECTION_FACE_Y_N);
        sched.add_link((*gridit).get_hydro_task(0), ngb_task);
      }

      // and z
      uint_fast32_t ngbz = (*gridit).get_neighbour(TRAVELDIRECTION_FACE_Z_P);
      if (ngbz == NEIGHBOUR_OUTSIDE) {
        Task &external_task = sched.get_free_task(ngb_task);
        external_task.set_type(TASKTYPE_TEST_GHOST_SWEEP);
        external_task.set_dependency((*gridit).get_dependency());
        external_task.set_subgrid(igrid);
        external_task.set_interaction_direction(TRAVELDIRECTION_FACE_Z_P);
        sched.add_link((*gridit).get_hydro_task(0), ngb_task);
      } else {
        Task &external_task = sched.get_free_task(ngb_task);
        external_task.set_type(TASKTYPE_TEST_NEIGHBOUR_SWEEP_EXTERNAL);
        if (igrid < ngbz) {
          external_task.set_dependency((*gridit).get_dependency());
          external_task.set_extra_dependency(
              (*grid.get_subgrid(ngbz)).get_dependency());
        } else {
          external_task.set_dependency(
              (*grid.get_subgrid(ngbz)).get_dependency());
          external_task.set_extra_dependency((*gridit).get_dependency());
        }
        external_task.set_subgrid(igrid);
        external_task.set_buffer(ngbz);
        external_task.set_interaction_direction(TRAVELDIRECTION_FACE_Z_P);
        sched.add_link((*gridit).get_hydro_task(0), ngb_task);
        sched.add_link((*grid.get_subgrid(ngbz)).get_hydro_task(0), ngb_task);
      }
      ngbz = (*gridit).get_neighbour(TRAVELDIRECTION_FACE_Z_N);
      if (ngbz == NEIGHBOUR_OUTSIDE) {
        Task &external_task = sched.get_free_task(ngb_task);
        external_task.set_type(TASKTYPE_TEST_GHOST_SWEEP);
        external_task.set_dependency((*gridit).get_dependency());
        external_task.set_subgrid(igrid);
        external_task.set_interaction_direction(TRAVELDIRECTION_FACE_Z_N);
        sched.add_link((*gridit).get_hydro_task(0), ngb_task);
      }
    }

    // the code above contains all the complexity required to deal with the
    // task dependencies
    // the code below actually executes the tasks (in parallel)
    // most of this code does not require any changing and is simply there
    // as an illustration of the algorithm

    // we create a dummy boundary condition that is used for the non-periodic
    // z boundaries of the box
    InflowHydroBoundary boundary;
    // the '#pragma' statement tells the compiler to compile the block below
    // as parallel code. This is basic OpenMP syntax.
#pragma omp parallel default(shared) num_threads(4)
    {
      // everything inside this block is executed by 4 threads simultaneously
      // each variable declared within the block is only accessible to one
      // thread and can be used as you would normally do
      // all other variables should be treated with great care:
      //  - only reading them is fine (so none of the threads ever writes to
      //    them)
      //  - writing to them is never safe, and can only be done if some sort
      //    of locking mechanism is used
      //  - reading a variable that is also used for writing does not
      //    necessarily require a locking mechanism, but almost always does

      // get the index of this thread from the OpenMP library
      // OpenMP assigns an index from [0, number of threads[ to all threads
      // we use this index to record the progress of the task-based algorithm
      // and to know which queue to preferentially get tasks from
      const int_fast32_t thread_id = get_thread_index();
      // the single line below is "critical"
      // OpenMP will make sure that only one thread can execute it at a time
      // we do this to avoid garbling the output messages
      // the only reason we display these messages is to prove that the code
      // is effectively executed in parallel
#pragma omp critical
      std::cout << "Thread " << thread_id << " active" << std::endl;

      // this is the essential part of the task-based algorithm: we enter a
      // loop that only finishes when all tasks have been executed
      // remember from above that number_of_tasks is an atomic variable, so that
      // access to the variable is per definition thread-safe
      // number_of_tasks.value() employs an atomic locking mechanism to
      // guarantee thread-safe reading of its current value
      // other operations in the loop will atomically increment number_of_tasks
      // when new tasks are queued, and decrement it when tasks have finished
      while (number_of_tasks.value() > 0) {
        // get a task from one of the queues
        // the Scheduler object first tries to obtain a task from this thread's
        // queue. If that fails, it tries to steal a task from another thread's
        // queue.
        // the get_task() function takes care of locking here: the thread gets
        // exclusive access to both the returned task and its dependencies
        uint_fast32_t task_index;
        Task &current_task = *sched.get_task(thread_id, task_index);
        // since get_task() can be called in parallel by multiple threads, there
        // is no guarantee that there will still be a task left for this thread
        // we need to check explicitly that this is the case
        if (task_index != NO_TASK) {
          // we got a task!
          // this thread now has guaranteed exclusive read and write access to
          // the Task object
          // if the task has one or two subgrid dependencies, it also has
          // exclusive access to these
          // this means that the test density and test neighbour density sum
          // variables can only be changed by this thread

          // start the task
          // this command simply records the start of the task
          // the start time and executing thread are later output so that we
          // can analyse the performance
          current_task.start(thread_id);

          // now actually execute the task
          // this is the part of the task-based algorithm that needs to be
          // changed when new tasks are created
          // we check the type of the task and call one of the subgrid functions
          // accordingly
          // the function calls are the same as for method 1, except that now
          // we need to retrieve the task arguments from the Task object
          if (current_task.get_type() == TASKTYPE_TEST_DENSITY) {
            // test density task, call set_test_density()
            (*grid.get_subgrid(current_task.get_subgrid())).set_test_density();
          } else if (current_task.get_type() ==
                     TASKTYPE_TEST_NEIGHBOUR_SWEEP_INTERNAL) {
            // internal test neighbour density sum task, call
            // inner_test_neighbour_density_sum_sweep()
            (*grid.get_subgrid(current_task.get_subgrid()))
                .inner_test_neighbour_density_sum_sweep();
          } else if (current_task.get_type() ==
                     TASKTYPE_TEST_NEIGHBOUR_SWEEP_EXTERNAL) {
            // external test neighbour density sum task, call
            // outer_test_neighbour_density_sum_sweep()
            (*grid.get_subgrid(current_task.get_subgrid()))
                .outer_test_neighbour_density_sum_sweep(
                    current_task.get_interaction_direction(),
                    *grid.get_subgrid(current_task.get_buffer()));
          } else if (current_task.get_type() == TASKTYPE_TEST_GHOST_SWEEP) {
            // boundary test neighbour density task, call
            // outer_ghost_test_neighbour_density_sum_sweep()
            (*grid.get_subgrid(current_task.get_subgrid()))
                .outer_ghost_test_neighbour_density_sum_sweep(
                    current_task.get_interaction_direction(), boundary);
          }

          // record the end of the task
          current_task.stop();

          // the task finished
          // before updating number_of_tasks, we unlock its dependencies:
          //  - we unlock the locks that protect the subgrid(s)
          //  - we add tasks that depended on this task being finished to the
          //    appropriate queues if possible
          sched.unlock_dependencies(task_index, number_of_tasks, grid);
          // now decrement number_of_tasks
          // note that we cannot do this before we call unlock_dependencies(),
          // since number_of_tasks is the control variable for the main loop
          // assume for example that this task was the last task that was
          // queued, but that it queues 3 new tasks
          // if we decrement number_of_tasks before queueing the new tasks, then
          // number_of_tasks will temporarily be 0, and other threads will
          // exit the while loop and be unavailable to execute the 3 new tasks
          number_of_tasks.pre_decrement();
        }
      } // end of main loop
    }   // end of parallel block

    // the code below is executed by a single thread
    // all variable access is safe again

    // dump diagnostic information about the task execution to "tasks.txt"
    sched.dump_tasks();

    // create a grid output object that writes output files of the form
    //  testXXX.hdf5, where XXX is a three digit counter value
    GridWriter writer("test");
    // write output file test001.hdf5 with the current contents of the grid
    writer.write(grid, 1);
  }

  // if all went well, test000.hdf5 and test001.hdf5 will be the same, except
  // for the iteration number stored in that file
  // this can be checked using "h5diff test000.hdf5 test001.hdf5"
  // example output:
  // attribute: <Iteration of </RuntimePars>> and <Iteration of </RuntimePars>>
  //  1 differences found


  // we are done
  // UNIX protocol dictates that we tell UNIX that our program exited
  // successfully by returing a status code 0
  return 0;
}
