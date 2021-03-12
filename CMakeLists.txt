cmake_minimum_required(VERSION 3.13)
project(mpcs-hpc-winter-2021-hw8 C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)

find_package(MPI)

# Serial advection, C++
add_executable(serial_advec_cpp cpp/serial_advec_main.cpp bitmap/bitmap.c)
target_include_directories(serial_advec_cpp PUBLIC cpp bitmap)

# Serial advection, C
add_executable(serial_advec_c c/serial_advec_main.c c/serial_advec_driver.c c/simple_mat.c bitmap/bitmap.c)
target_include_directories(serial_advec_c PUBLIC c bitmap)
target_link_libraries(serial_advec_c PUBLIC m)

if (MPI_FOUND)
    # MPI advection, C++
    add_executable(mpi_advec_cpp cpp/mpi_advec_main.cpp bitmap/bitmap.c)
    target_include_directories(mpi_advec_cpp PUBLIC cpp bitmap)
    target_link_libraries(mpi_advec_cpp PUBLIC MPI::MPI_C)

    # MPI advection, C
    add_executable(mpi_advec_c c/mpi_advec_main.c c/mpi_advec_driver.c c/simple_mat.c bitmap/bitmap.c)
    target_include_directories(mpi_advec_c PUBLIC c bitmap)
    target_link_libraries(mpi_advec_c PUBLIC MPI::MPI_C m)
endif ()
