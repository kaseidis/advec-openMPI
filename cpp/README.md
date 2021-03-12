# C++ Advection

## Matrix Data Structure 

A datatype for a matrix of doubles is implemented in `simple_mat.hpp`. This acts like the matrix data 
structures we've seen before.  Refer to the comments in `simple_mat.hpp` for usage.

## Serial Advection

A serial advection program is fully implemented in `serial_advec_driver.hpp`,  `serial_advec_main.cpp`, 
and the corresponding headers.

### `serial_advec_main.cpp`

This is the main function for running the simulation.  The simulation's parameters are hardcoded.  
It is compiled by the makefile target `serial_advec_cpp`.

### `serial_advec_driver.hpp`

`SerialAdvecDriver`: A class that packages the data for the problem.  Refer to the comments in the 
header for info about the data members (at bottom of class definition).  The member functions serve this puprpose
  
* Constructor: Sets up problem (including memory allocation)
* `init_cond`:  Returns the initial condition at point (x, y)
* `run`:  Runs the simulation for the given number of timesteps `nt`.  
   A bitmap is written at each timestep interval specified by `output_interval`.  If `output_interval = 0`, 
   there are no images produced.
* `advance`:  Runs a single timestep of the simulation, including ghost cell updates
* `update_ghost_cells`:  Updates all the ghost cells
* `output`:  Writes a bitmap image of the domain.  The file's name includes the timestep `t`

## MPI Advection

A serial advection program is partially implemented in `mpi_advec_driver.c`,  `mpi_advec_main.c`, 
and the corresponding headers.    

### `mpi_advec_main.c`

This is the main function for running the simulation.  The simulation's parameters are hardcoded.  
It is compiled by the makefile target `mpi_advec_cpp`

### `mpi_advec_driver.c` and `mpi_advec_driver.h`

* `MpiAdvec`: A class that packages the data for the problem.  Refer to the comments in the 
  header for info about data members (at bottom of class defintions).
* The member functions serve the same purposes as the corresponding member functions in `SerialAdvecDriver` but they are implemented in parallel. 
* The only member function that is **not** implemented is `output`.  Implementing it is the purpose of Hw 8.  
