# Advection solver

Advection solver with periodic boundaries

## Author Name
Shenghua Chen

## Comments

### Files
Here is file list:
- ``bitmap/*`` Code for generate bitmap
- ``c/serial_advec_driver.c[h]`` Code do advection without MPI.
- ``c/serial_advec_main.c`` Main function for above code.
- ``c/mpi_advec_driver.c[h]`` Code do advection with MPI.
- ``c/mpi_advec_main.c`` Main function for above code.
- ``c/simple_mat.c[h]`` Simple matrix code

## Instruction for compile
Navigate to problem folder, run ``cmake .`` then ``make``.

Usage:
```
./serial_advec_c
mpirun -np [number of process] ./mpi_advec_c
```
