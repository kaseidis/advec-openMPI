//
// Created by Ronald Rahaman on 3/7/21.
//

#include "serial_advec_driver.h"
#include "bitmap.h"
#include <limits.h>
#include <math.h>
#include <stdio.h>

SerialDriverData serial_advec_init(double L, int nx, double dt, double u, double v)
{
  SerialDriverData d;
  d.L = L;
  d.nx = nx;
  d.ny = nx;
  d.dt = dt;
  d.u = u;
  d.v = u;
  d.dx = L / nx;
  d.dy = d.dx;
  d.xmin = -L / 2.0;
  d.ymin = d.xmin;
  d.xmax = L / 2.0;
  d.ymax = d.xmax;
  d.c = init_simple_mat(d.nx + 2, d.ny + 2);
  d.c_nxt = init_simple_mat(d.nx + 2, d.ny + 2);

  // Assert stability
  if (d.dt > d.dx / sqrt(2.0 * (d.u * d.u + d.v * d.v))) {
    fprintf(stderr, "Input did not meet stability condition");
    exit(EXIT_FAILURE);
  }

  // Initial condtion
  for (int i = 0; i < d.nx + 2; ++i) {
    for (int j = 0; j < d.ny + 2; ++j) {
      double x = d.xmin + i * d.dx;
      double y = d.ymin + j * d.dy;
      d.c.at[i][j] = advec_init_cond(x, y, d.L);
    }
  }

  d.c_max = max_simple_mat(&d.c);

  return d;
}

// Initial condition at point x, y
double advec_init_cond(double x, double y, double L)
{
  return exp(-(8.0 * x * x + 8.0 * y * y) / (L * L));
}

void serial_advec_update_ghost_cells(SerialDriverData* d)
{
  // For readability
  const int nx = d->nx;
  const int ny = d->ny;

  // Update the columns on the edge
  for (int i = 1; i < nx + 1; ++i) {
    d->c.at[i][0] = d->c.at[i][ny];
    d->c.at[i][ny + 1] = d->c.at[i][1];
  }
  // Update rows on edge
  for (int j = 1; j < ny + 1; ++j) {
    d->c.at[0][j] = d->c.at[nx][j];
    d->c.at[nx + 1][j] = d->c.at[1][j];
  }
}

void serial_advec_run(SerialDriverData* d, int nt, int output_interval)
{

  for (int t = 0; t < nt; ++t) {
    serial_advec_advance(d);
    if (output_interval > 0 && t % output_interval == 0) {
      serial_advec_output(d, t);
    }
  }
}

void serial_advec_output(SerialDriverData* d, int t)
{
  RgbTriple* pixmap = malloc(d->nx * d->ny * sizeof(RgbTriple));
  for (int i = 0; i < d->nx; ++i) {
    for (int j = 0; j < d->ny; ++j) {
      unsigned char color = d->c.at[i][j] * 255 / d->c_max;
      RgbTriple* p = &pixmap[i * d->nx + j];
      p->red = color;
      p->green = color;
      p->blue = color;
    }
  }
  char filename[PATH_MAX];
  sprintf(filename, "serial_advec_%04d.bmp", t);
  save_bitmap(pixmap, d->nx, d->ny, filename);
  free(pixmap);
}

void serial_advec_advance(SerialDriverData* d)
{
  // For readability
  const int nx = d->nx;
  const int ny = d->ny;
  const double dt = d->dt;
  const double dx = d->dx;
  const double u = d->u;
  const double v = d->v;

  for (int i = 1; i < nx + 1; ++i) {
    for (int j = 1; j < ny + 1; ++j) {
      d->c_nxt.at[i][j] = 1.0 / 4.0 *
                            (d->c.at[i - 1][j] + d->c.at[i + 1][j] + d->c.at[i][j - 1] +
                             d->c.at[i][j + 1]) -
                          dt / (2.0 * dx) *
                            (u * (d->c.at[i - 1][j] - d->c.at[i + 1][j]) +
                             v * (d->c.at[i][j - 1] - d->c.at[i][j + 1]));
    }
  }
  copy_simple_mat(&d->c, &d->c_nxt);
  serial_advec_update_ghost_cells(d);
}

void serial_advec_free(SerialDriverData* d)
{
  free_simple_mat(&d->c);
  free_simple_mat(&d->c_nxt);
}
