#include "math.h"
#include "mpi_advec_driver.h"
#include "simple_mat.h"

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  double L = 1.0;
  int nx = 200;
  int output_interval = 5;
  // double u = 1e-6;
  // double v = 1e-6;
  double u = 5e-7;
  double v = 1.3e-6;
  int nt = 2000;

  // This is numerically stable
  double dt = L / (nx + 1) / sqrt(2 * (u * u + v * v));

  MpiDriverData d = mpi_advec_init(L, nx, dt, u, v);
  mpi_advec_run(&d, nt, output_interval);
  mpi_advec_free(&d);

  MPI_Finalize();
}
