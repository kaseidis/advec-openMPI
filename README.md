# Homework 8 - MPCS HPC Winter 2021

See assignment at:  https://www.dropbox.com/s/jvsdzmsnnb2p2rx/homework%208%20-%20mpcs%20hpc%20winter%202021.pdf?dl=0 

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
- ``c/simple_mat.c[h]`` Simple matrix code provided.

## Instruction for compile
Navigate to problem folder, run ``cmake .`` then ``make``.

Usage:
```
./serial_advec_c
mpirun -np [number of process] ./mpi_advec_c
```
