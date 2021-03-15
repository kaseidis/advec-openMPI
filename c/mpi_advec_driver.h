#ifndef MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_C_MPI_ADVEC_DRIVER_H
#define MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_C_MPI_ADVEC_DRIVER_H

#include "simple_mat.h"
#include <mpi.h>
#include <stdlib.h>

typedef struct {
  double L;          //! Length of global domain
  int nx, ny;        //! Number of gridpoints in global domain
  double dt;         //! Length of timestep
  double u, v;       //! Components of velocity
  double dx, dy;     //! Spacing between gridpoints
  double xmin, xmax; //! Boundaries
  double ymin, ymax; //! Boundaries
  double c_max;

  int nx_local, ny_local;        //! Number of gridpoints in this rank's subdomain
  double xmin_local, xmax_local; //! Min and max physical bounds of this rank's subdomain
  double ymin_local, ymax_local; //! Min and max physical bounds of this rank's subdomain
  SimpleMat c;                   //! This rank's subdomain at this timestep
  SimpleMat c_nxt;               //! This rank's subdomain at next timestep

  MPI_Comm cart_comm; //! The cartesian communicator

  int world_size, world_rank; //! Size and rank in MPI_COMM_WORLD
  int cart_size, cart_rank;   //! Size and rank in cart_comm

  int nbr_up, nbr_down; //! The ranks above and below this rank in cart_comm
  int cart_coords[1];   //! The coordinates of this rank in cart_comm
  int cart_dims[1];     //! The dimensions of cart_comm
  int cart_periodic[1]; //! Specifies which boundaries are periodic

  MPI_Datatype mpi_rgb_type;

} MpiDriverData;

MpiDriverData mpi_advec_init(double L, int nx, double dt, double u, double v);

double advec_init_cond(double x, double y, double L);

void mpi_advec_update_ghost_cells(MpiDriverData* d);

void mpi_advec_run(MpiDriverData* d, int nt, int output_interval);

void mpi_advec_output(MpiDriverData* d, int t);

void mpi_advec_advance(MpiDriverData* d);

void mpi_advec_free(MpiDriverData* d);

#endif // MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_C_MPI_ADVEC_DRIVER_H
