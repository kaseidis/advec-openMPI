#include "mpi_advec_driver.hpp"

int main(int argc, char* argv[])
{

  MPI_Init(&argc, &argv);

  double L = 1.0;
  int nx = 200;
  int output_interval = 5;
  double u = 5e-7;
  double v = 1.3e-6;
  int nt = 1000;

  // This is numerically stable
  double dt = L / (nx + 1) / std::sqrt(2 * (u * u + v * v));

  MpiAdvecDriver<double> driver(L, nx, dt, u, v);
  driver.run(nt, output_interval);

  MPI_Finalize();
}