# C Advection

## Matrix Data Structure 

A datatype for a matrix of doubles is implemented in `simple_mat.c` and `simple_mat.h`.  This acts 
like the matrix data structures we've seen before.  Refer to the comments in `simple_mat.h` for usage

## Serial Advection

A serial advection program is fully implemented in `serial_advec_driver.c`,  `serial_advec_main.c`, 
and the corresponding headers.  It is compiled by the makefile target `serial_advec_c`.

### `serial_advec_main.c`

This is the main function for running the simulation.  The simulation's parameters are hardcoded,

### `serial_advec_driver.c` and `serial_advec_driver.h`

* `SerialDriverData`: A struct that packages the data for the problem.  Refer to the comments in the 
  header for info about the struct members
* `serial_advec_init`:  Sets up problem (including memory allocation) and returns an instance of a `SerialDriverData` struct.  Use this instance when you're calling the rest of the functions.
* `advec_init_cond`:  Returns the initial condition at point (x, y)
* `serial_advec_run`:  Runs the simulation for the given number of timesteps `nt`.  
   A bitmap is written at each timestep interval specified by `output_interval`.  If `output_interval = 0`, 
   there are no images produced.
* `serial_advec_advance`:  Runs a single timestep of the simulation, including ghost cell updates
* `serial_advec_update_ghost_cells`:  Updates all the ghost cells
* `serial_advec_output`:  Writes a bitmap image of the domain.  The file's name includes the timestep `t`
* `serial_advec_free`:  Frees all data structures.

## MPI Advection

A serial advection program is partially implemented in `mpi_advec_driver.c`,  `mpi_advec_main.c`, 
and the corresponding headers.  It is compiled by the makefile target `mpi_advec_c`.

### `mpi_advec_main.c`

This is the main function for running the simulation.  The simulation's parameters are hardcoded,

### `mpi_advec_driver.c` and `mpi_advec_driver.h`

* `MpiDriverData`: A struct that packages the data for the problem.  Refer to the comments in the 
  header for info about the struct members.  
* The other functions serve the same purposes as the corresponding serial functions, but they are implemented
  in parallel. 
* The only function that is **not** implemented is `mpi_advec_output`.  Implementing it is the purpose of Hw 8.  
