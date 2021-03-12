#ifndef MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_GET_MPI_TYPES_HPP
#define MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_GET_MPI_TYPES_HPP

/* =============================================================

Templates for getting MPI datatypes.  Useful in MPI function calls.
For example:

template<typename T>
void my_recv(std::vector<T> x, int src_rank) {
  MPI_Recv(x.data(),
           x.size(),
           get_mpi_type<T>(),
           MPI_ANY_SOURCE,
           MPI_ANY_TAG,
           MPI_COMM_WORLD,
           src_rank,
           MPI_STATUS_IGNORE);
}

============================================================= */

#include "mpi.h"

// Getting datatypes
template<typename T>
MPI_Datatype get_mpi_type();

// Template specializations
template<>
MPI_Datatype get_mpi_type<int>()
{
  return MPI_INT;
}

template<>
MPI_Datatype get_mpi_type<float>()
{
  return MPI_FLOAT;
}

template<>
MPI_Datatype get_mpi_type<double>()
{
  return MPI_DOUBLE;
}

#endif // MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_GET_MPI_TYPES_HPP
