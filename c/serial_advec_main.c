//
// Created by Ronald Rahaman on 3/7/21.
//

#include "serial_advec_driver.h"
#include "simple_mat.h"
#include <math.h>

int main()
{
  double L = 1.0;
  int nx = 200;
  int output_interval = 5;
  double u = 5e-7;
  double v = 1.3e-6;
  int nt = 2000;

  // This is numerically stable
  double dt = L / (nx + 1) / sqrt(2 * (u * u + v * v));

  SerialDriverData d = serial_advec_init(L, nx, dt, u, v);
  // serial_output(&d, 0);
  serial_advec_run(&d, nt, output_interval);
  serial_advec_free(&d);
}
