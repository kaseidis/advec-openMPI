#include "mpi_advec_driver.h"

//
// Created by Ronald Rahaman on 3/7/21.
//

#include "bitmap.h"
#include "mpi_advec_driver.h"
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

MpiDriverData mpi_advec_init(double L, int nx, double dt, double u, double v)
{
  MpiDriverData d;
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

  MPI_Comm_size(MPI_COMM_WORLD, &d.world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &d.world_rank);

  if (d.dt > d.dx / sqrt(2.0 * (d.u * d.u + d.v * d.v))) {
    fprintf(stderr, "Input did not meet stability condition");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if (d.world_size % 2 != 0) {
    fprintf(stderr, "Error: algorithm requires and even number of MPI ranks\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if (nx % d.world_size != 0) {
    fprintf(stderr, "Error: domain is not evenly divisible by the number of MPI ranks\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // Desired dimension of cart comm along each dim
  // Here, we have a 1D cart comm with length word size
  d.cart_dims[0] = d.world_size;
  d.cart_periodic[0] = 1;

  // Create cartesian communicator
  MPI_Cart_create(MPI_COMM_WORLD,  // Super communicator
                  1,               // Number of dims in new cart comm
                  d.cart_dims,     // Length of new comm in each dim
                  d.cart_periodic, // Is domain perioidic along each dim?
                  0,               // Allow reording of comms from old comm
                  &d.cart_comm);   // New comm
  MPI_Comm_size(d.cart_comm, &d.cart_size);
  MPI_Comm_rank(d.cart_comm, &d.cart_rank);

  // Find this rank's position in cart comm
  MPI_Cart_get(d.cart_comm,     // Cartesian communicator
               1,               // Number of dimensions in cart comm
               d.cart_dims,     // Length of cart comm along each dim,
               d.cart_periodic, // Is domain periodic along each dim?
               d.cart_coords);  // The coordinates of this rank in the cart comm (result)

  // Find this rank's neighbors on top and bottom
  MPI_Cart_shift(d.cart_comm, 0, 1, &d.nbr_down, &d.nbr_up);

  // Get bounds
  d.nx_local = d.nx / d.cart_dims[0];
  d.xmin_local = d.xmin + d.cart_coords[0] * d.nx_local * d.dx;
  d.xmax_local = d.xmin + (d.cart_coords[0] + 1) * d.nx_local * d.dx;

  d.ny_local = d.ny;
  d.ymin_local = d.ymin;
  d.ymax_local = d.ymax;

  // Useful for debugging
  for (int i = 0; i < d.cart_size; ++i) {
    if (d.cart_rank == i) {
      printf("Card rank: %d:\n", d.cart_rank);
      printf("    Card coord: %d\n", d.cart_coords[0]);
      printf("    Neighbors (down, up): %d, %d\n", d.nbr_down, d.nbr_up);
      printf("    x-limits: %0.4f, %0.4f\n", d.xmin_local, d.xmax_local);
      printf("    y-limits: %0.4f, %0.4f\n", d.ymin_local, d.ymax_local);
    }
    MPI_Barrier(d.cart_comm);
  }

  d.c = init_simple_mat(d.nx_local + 2, d.ny_local + 2);
  d.c_nxt = init_simple_mat(d.nx_local + 2, d.ny_local + 2);

  for (int i = 0; i < d.nx_local + 2; ++i) {
    for (int j = 0; j < d.ny_local + 2; ++j) {
      double x = d.xmin_local + i * d.dx;
      double y = d.ymin_local + j * d.dy;
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

void mpi_advec_update_ghost_cells(MpiDriverData* d)
{
  // For readability
  const int nx_local = d->nx_local;
  const int ny_local = d->ny_local;
  const int nbr_up = d->nbr_up;
  const int nbr_down = d->nbr_down;

  // Ghost columns can be updated without any communication
  for (int i = 1; i < nx_local + 1; ++i) {
    d->c.at[i][0] = d->c.at[i][ny_local];
    d->c.at[i][ny_local + 1] = d->c.at[i][1];
  }

  bool is_even = (d->cart_coords[0] % 2 == 0);

  if (is_even) {
    // 1. Even ranks send top interior row to top neighbor
    MPI_Send(&d->c.at[nx_local][1], ny_local, MPI_DOUBLE, nbr_up, 314159, d->cart_comm);
  } else {
    // 1. Odd ranks recv bottom ghost row from bottom neighbor
    MPI_Recv(&d->c.at[0][1],
             ny_local,
             MPI_DOUBLE,
             nbr_down,
             MPI_ANY_TAG,
             d->cart_comm,
             MPI_STATUS_IGNORE);
  }

  if (is_even) {
    // 2. Even ranks recv bottom ghost row from bottom neighbor
    MPI_Recv(&d->c.at[0][1],
             ny_local,
             MPI_DOUBLE,
             nbr_down,
             MPI_ANY_TAG,
             d->cart_comm,
             MPI_STATUS_IGNORE);
  } else {
    // 2. Odd ranks send top interior row to top neighbor
    MPI_Send(&d->c.at[nx_local][1], ny_local, MPI_DOUBLE, nbr_up, 314159, d->cart_comm);
  }

  if (is_even) {
    // 3. Even ranks send bottom interior row to bottom neighbor
    MPI_Send(&d->c.at[1][1], ny_local, MPI_DOUBLE, nbr_down, 314159, d->cart_comm);
  } else {
    // 3. Odd ranks recv top ghost row from top  neighbor
    MPI_Recv(&d->c.at[nx_local + 1][1],
             ny_local,
             MPI_DOUBLE,
             nbr_up,
             MPI_ANY_TAG,
             d->cart_comm,
             MPI_STATUS_IGNORE);
  }

  if (is_even) {
    // 4. Even ranks recv top ghost row from top  neighbor
    MPI_Recv(&d->c.at[nx_local + 1][1],
             ny_local,
             MPI_DOUBLE,
             nbr_up,
             MPI_ANY_TAG,
             d->cart_comm,
             MPI_STATUS_IGNORE);
  } else {
    // 4. Odd ranks send bottom interior row to bottom neighbor
    MPI_Send(&d->c.at[1][1], ny_local, MPI_DOUBLE, nbr_down, 314159, d->cart_comm);
  }
}

void mpi_advec_run(MpiDriverData* d, int nt, int output_interval)
{

  for (int t = 0; t < nt; ++t) {
    mpi_advec_advance(d);
    if (output_interval > 0 && t % output_interval == 0) {
      mpi_advec_output(d, t);
    }
  }
}

void mpi_advec_output(MpiDriverData* d, int t)
{
  return;
}

void mpi_advec_advance(MpiDriverData* d)
{
  // For readability
  const int nx_local = d->nx_local;
  const int ny_local = d->ny_local;
  const double dt = d->dt;
  const double dx = d->dx;
  const double u = d->u;
  const double v = d->v;

  for (int i = 1; i < nx_local + 1; ++i) {
    for (int j = 1; j < ny_local + 1; ++j) {
      d->c_nxt.at[i][j] = 1.0 / 4.0 *
                            (d->c.at[i - 1][j] + d->c.at[i + 1][j] + d->c.at[i][j - 1] +
                             d->c.at[i][j + 1]) -
                          dt / (2.0 * dx) *
                            (u * (d->c.at[i - 1][j] - d->c.at[i + 1][j]) +
                             v * (d->c.at[i][j - 1] - d->c.at[i][j + 1]));
    }
  }
  copy_simple_mat(&d->c, &d->c_nxt);
  mpi_advec_update_ghost_cells(d);
}

void mpi_advec_free(MpiDriverData* d)
{
  free_simple_mat(&d->c);
  free_simple_mat(&d->c_nxt);
}
