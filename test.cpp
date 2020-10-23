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

  // create a scheduler (needs proper usage)
  Scheduler sched(10, 4, 10, 10);

  // create the initial condition object
  TestDensityFunction test_function;
  // create the distributed grid
  // the simulation box that spans the cube with corners [0,0,0] and [1,1,1]
  // we use a 64x64x64 grid subdivided into 8x8x8 subgrids
  // all boundaries of the simulation box are assumed to be periodic
  DensitySubGridCreator<HydroDensitySubGrid> grid(Box<>(0., 1.), 64, 8, true);
  // set up and initialize the cells and subgrids in the distributed grid
  grid.initialize(test_function);

  /// METHOD 1

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
      (*gridit).outer_test_neighbour_density_sum_sweep(TRAVELDIRECTION_FACE_X_P,
                                                       *grid.get_subgrid(ngbx));
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
      (*gridit).outer_test_neighbour_density_sum_sweep(TRAVELDIRECTION_FACE_Y_P,
                                                       *grid.get_subgrid(ngby));
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
      (*gridit).outer_test_neighbour_density_sum_sweep(TRAVELDIRECTION_FACE_Z_P,
                                                       *grid.get_subgrid(ngbz));
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

  // we are done
  // UNIX protocol dictates that we tell UNIX that our program exited
  // successfully by returing a status code 0
  return 0;
}
