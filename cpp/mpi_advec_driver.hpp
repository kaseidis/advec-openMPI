#ifndef MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_MPI_ADVEC_DRIVER_HPP
#define MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_MPI_ADVEC_DRIVER_HPP

#include "get_mpi_types.hpp"
#include "mpi.h"
#include "simple_mat.hpp"
extern "C" {
#include "bitmap.h"
};
#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

template<typename T>
class MpiAdvecDriver {
public:
  // Using ghost cells in c and d_tp1.
  // Not necessary for serial impmlementation, but helpful for prototyping MPI version
  MpiAdvecDriver(T L, int nx, T dt, T u, T v)
    : L(L)
    , nx(nx)
    , ny(nx)
    , dt(dt)
    , u(u)
    , v(v)
    , dx(L / nx)
    , dy(dx)
    , xmin(-L / 2.0)
    , ymin(xmin)
    , xmax(L / 2.0)
    , ymax(xmax)
  {
    if (this->dt > this->dx / sqrt(2.0 * (this->u * this->u + this->v * this->v))) {
      std::cerr << "Error: input did not meet stability condition" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    init_comms();
  }

  void init_comms()
  {
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_size % 2 != 0) {
      std::cerr << "Error: Algorithm requires an even number of MPI ranks" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (nx % world_size != 0) {
      std::cerr << "Error: Domain is not evenly divisible by the number of MPI ranks"
                << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Desired dimension of cart comm along each dim
    // Here, we have a 1D cart comm with length word size
    cart_dims[0] = world_size;

    // Create cartesian communicator
    MPI_Cart_create(MPI_COMM_WORLD, // Super communicator
                    1,              // Number of dims in new cart comm
                    cart_dims,      // Length of new comm in each dim
                    cart_periodic,  // Is domain perioidic along each dim?
                    0,              // Allow reording of comms from old comm
                    &cart_comm);    // New comm
    MPI_Comm_size(cart_comm, &cart_size);
    MPI_Comm_rank(cart_comm, &cart_rank);

    // Find this rank's position in cart comm
    MPI_Cart_get(cart_comm,     // Cartesian communicator
                 1,             // Number of dimensions in cart comm
                 cart_dims,     // Length of cart comm along each dim,
                 cart_periodic, // Is domain periodic along each dim?
                 cart_coords);  // The coordinates of this rank in the cart comm (result)

    // Find this rank's neighbors on top and bottom
    MPI_Cart_shift(cart_comm, 0, 1, &nbr_down, &nbr_up);

    // Get bounds
    nx_local = nx / cart_dims[0];
    xmin_local = xmin + cart_coords[0] * nx_local * dx;
    xmax_local = xmin + (cart_coords[0] + 1) * nx_local * dx;

    ny_local = ny;
    ymin_local = ymin;
    ymax_local = ymax;

    // Useful for debugging
    for (int i = 0; i < cart_size; ++i) {
      if (cart_rank == i) {
        std::cout << "Cart Rank: " << cart_rank << std::endl;
        std::cout << "    Cart Coord: " << cart_coords[0] << std::endl;
        std::cout << "    Neighbors (down, up): " << nbr_down << ", " << nbr_up
                  << std::endl;
        std::cout << "    x-limits: " << xmin_local << ", " << xmax_local << std::endl;
        std::cout << "    y-limits: " << ymin_local << ", " << ymax_local << std::endl;
      }
      MPI_Barrier(cart_comm);
    }

    c = SimpleMat<T>(nx_local + 2, ny_local + 2);
    c_nxt = SimpleMat<T>(nx_local + 2, ny_local + 2);

    for (int i = 0; i < nx_local + 2; ++i) {
      for (int j = 0; j < ny_local + 2; ++j) {
        auto x = xmin_local + i * dx;
        auto y = ymin_local + j * dy;
        c.at[i][j] = init_cond(x, y);
      }
    }
  }

  // Initial condition at point x, y
  T init_cond(T x, T y) { return exp(-(8.0 * x * x + 8.0 * y * y) / (L * L)); }

  void update_ghost_cells()
  {

    // Ghost columns can be updated without any communication
    for (int i = 1; i < nx_local + 1; ++i) {
      c.at[i][0] = c.at[i][ny_local];
      c.at[i][ny_local + 1] = c.at[i][1];
    }

    auto is_even = (cart_coords[0] % 2 == 0);

    if (is_even) {
      // 1. Even ranks send top interior row to top neighbor
      MPI_Send(
        &c.at[nx_local][1], ny_local, get_mpi_type<T>(), nbr_up, 314159, cart_comm);
    } else {
      // 1. Odd ranks recv bottom ghost row from bottom neighbor
      MPI_Recv(&c.at[0][1],
               ny_local,
               get_mpi_type<T>(),
               nbr_down,
               MPI_ANY_TAG,
               cart_comm,
               MPI_STATUS_IGNORE);
    }

    if (is_even) {
      // 2. Even ranks recv bottom ghost row from bottom neighbor
      MPI_Recv(&c.at[0][1],
               ny_local,
               get_mpi_type<T>(),
               nbr_down,
               MPI_ANY_TAG,
               cart_comm,
               MPI_STATUS_IGNORE);
    } else {
      // 2. Odd ranks send top interior row to top neighbor
      MPI_Send(
        &c.at[nx_local][1], ny_local, get_mpi_type<T>(), nbr_up, 314159, cart_comm);
    }

    if (is_even) {
      // 3. Even ranks send bottom interior row to bottom neighbor
      MPI_Send(&c.at[1][1], ny_local, get_mpi_type<T>(), nbr_down, 314159, cart_comm);
    } else {
      // 3. Odd ranks recv top ghost row from top  neighbor
      MPI_Recv(&c.at[nx_local + 1][1],
               ny_local,
               get_mpi_type<T>(),
               nbr_up,
               MPI_ANY_TAG,
               cart_comm,
               MPI_STATUS_IGNORE);
    }

    if (is_even) {
      // 4. Even ranks recv top ghost row from top  neighbor
      MPI_Recv(&c.at[nx_local + 1][1],
               ny_local,
               get_mpi_type<T>(),
               nbr_up,
               MPI_ANY_TAG,
               cart_comm,
               MPI_STATUS_IGNORE);
    } else {
      // 4. Odd ranks send bottom interior row to bottom neighbor
      MPI_Send(&c.at[1][1], ny_local, get_mpi_type<T>(), nbr_down, 314159, cart_comm);
    }
  }

  // Run nt timesteps
  void run(int nt, int output_interval)
  {
    for (int t = 0; t < nt; ++t) {
      advance();
      if (output_interval > 0 && t % output_interval == 0) {
        output(t);
      }
    }
  }

  // This is for you!
  void output(int timestep)
  {
    return;
  }

  // Advance one timestep
  void advance()
  {
    for (int i = 1; i < nx_local + 1; ++i) {
      for (int j = 1; j < ny_local + 1; ++j) {
        c_nxt.at[i][j] =
          1.0 / 4.0 *
            (c.at[i - 1][j] + c.at[i + 1][j] + c.at[i][j - 1] + c.at[i][j + 1]) -
          dt / (2.0 * dx) *
            (u * (c.at[i - 1][j] - c.at[i + 1][j]) +
             v * (c.at[i][j - 1] - c.at[i][j + 1]));
      }
    }
    c = c_nxt;
    update_ghost_cells();
  }

  const T L;          //! length of domain
  const int nx, ny;   //! Number of gridpoints
  const T dt;         //! Length of timestep
  const T u, v;       //! Components of velocity
  const T dx, dy;     //! Spacing between gridpoints
  const T xmin, xmax; //! Boundaries
  const T ymin, ymax; //! Boundaries

  int nx_local, ny_local;        //! Number of gridpoints in this rank's subdomain
  double xmin_local, xmax_local; //! Min and max physical bounds of this rank's subdomain
  double ymin_local, ymax_local; //! Min and max physical bounds of this rank's subdomain
  SimpleMat<T> c;                //! This rank's subdomain at this timestep
  SimpleMat<T> c_nxt;            //! This rank's subomain at next timestep

  int world_size, world_rank; //! Size and rank in MPI_COMM_WORLD
  int cart_size, cart_rank;   //! Size and rank in cart_comm

  MPI_Comm cart_comm; //! The cartesian communicator

  int nbr_up, nbr_down;       //! The ranks above and below this rank in cart_comm
  int cart_coords[1];         //! The coordinates of this rank in cart_comm
  int cart_dims[1];           //! The dimensions of cart_comm
  int cart_periodic[1] = {1}; //! Specifies which boundaries are periodic
};

#endif // MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_MPI_ADVEC_DRIVER_HPP
