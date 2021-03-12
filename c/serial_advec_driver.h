#ifndef MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_C_SERIAL_ADVEC_DRIVER_H
#define MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_C_SERIAL_ADVEC_DRIVER_H

#include "simple_mat.h"
#include <stdlib.h>

typedef struct {
  double L;          //! length of domain
  int nx, ny;        //! Number of gridpoints
  double dt;         //! Length of timestep
  double u, v;       //! Components of velocity
  double dx, dy;     //! Spacing between gridpoints
  double xmin, xmax; //! Boundaries
  double ymin, ymax; //! Boundaries
  SimpleMat c;       //! Domain at this timestep
  SimpleMat c_nxt;   //! Domain at next timestep
  double c_max;      //! Maximum value in c
} SerialDriverData;

SerialDriverData serial_advec_init(double L, int nx, double dt, double u, double v);

double advec_init_cond(double x, double y, double L);

void serial_advec_run(SerialDriverData* d, int nt, int output_interval);

void serial_advec_advance(SerialDriverData* d);

void serial_advec_update_ghost_cells(SerialDriverData* d);

void serial_advec_output(SerialDriverData* d, int t);

void serial_advec_free(SerialDriverData* d);

#endif // MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_C_SERIAL_ADVEC_DRIVER_H
