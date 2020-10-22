#include "DensitySubGridCreator.hpp"
#include "GridWriter.hpp"
#include "HydroDensitySubGrid.hpp"
#include "Scheduler.hpp"

class TestDensityFunction : public DensityFunction {
public:
  virtual DensityValues operator()(const Cell &cell) {
    DensityValues values;
    values.set_number_density(1.);
    return values;
  }
};

int main(int argc, char **argv) {

  Scheduler sched(10, 4, 10, 10);

  TestDensityFunction test_function;
  DensitySubGridCreator<HydroDensitySubGrid> grid(Box<>(0., 1.), 64, 8, true);
  grid.initialize(test_function);

  GridWriter writer("test");
  writer.write(grid, 0);

  return 0;
}
