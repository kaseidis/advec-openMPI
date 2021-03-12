#ifndef MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_SERIAL_ADVEC_DRIVER_HPP
#define MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_SERIAL_ADVEC_DRIVER_HPP

#include "simple_mat.hpp"
extern "C" {
#include "bitmap.h"
};
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

template<typename T>
class SerialAdvecDriver {
public:
  // Using ghost cells in c and d_tp1.
  // Not necessary for serial impmlementation, but helpful for prototyping MPI version
  SerialAdvecDriver(T L, int nx, T dt, T u, T v)
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
    , c(nx + 2, ny + 2)
    , c_nxt(nx + 2, ny + 2)
  {

    // Assert stability
    if (dt > dx / std::sqrt(2 * (u * u + v * v))) {
      throw std::runtime_error("Input did not meet stability condition");
    }

    // Initial condtion
    for (int i = 1; i < nx + 1; ++i) {
      for (int j = 1; j < ny + 1; ++j) {
        T x = xmin + i * dx;
        T y = ymin + j * dy;
        c.at[i][j] = init_cond(x, y);
      }
    }
    c_max = *std::max_element(c.data, c.data + nx * ny);
    c_min = *std::min_element(c.data, c.data + nx * ny);

    // Update the ghost cells
    update_ghost_cells();
  }

  // Initial condition at point x, y
  T init_cond(T x, T y) { return exp(-(8.0 * x * x + 8.0 * y * y) / (L * L)); }

  // Update ghost cells
  void update_ghost_cells()
  {
    // Update the columns on the edge
    for (int i = 1; i < nx + 1; ++i) {
      c.at[i][0] = c.at[i][ny];
      c.at[i][ny + 1] = c.at[i][1];
    }
    // Update rows on edge
    for (int j = 1; j < ny + 1; ++j) {
      c.at[0][j] = c.at[nx][j];
      c.at[nx + 1][j] = c.at[1][j];
    }
  }

  // Run nt timesteps
  void run(int nt, int output_interval)
  {
    std::vector<RgbTriple> pixmap;
    if (output_interval > 0) {
      pixmap.resize(nx * ny);
    }

    for (int t = 0; t < nt; ++t) {
      advance();

      if (output_interval > 0 && t % output_interval == 0) {
        for (int i = 0; i < nx; ++i) {
          for (int j = 0; j < ny; ++j) {
            unsigned char color = c_min + c.at[i][j] * 255 / (c_max - c_min);
            auto& p = pixmap.at(i * nx + j);
            p.red = color;
            p.green = color;
            p.blue = color;
          }
        }
        std::stringstream filename;
        filename << "serial_advec_" << std::setfill('0') << std::setw(4) << t << ".bmp";
        save_bitmap(pixmap.data(), nx, ny, filename.str().c_str());
      }
    }
  }

  // Advance one timestep
  void advance()
  {
    for (int i = 1; i < nx + 1; ++i) {
      for (int j = 1; j < ny + 1; ++j) {
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
  SimpleMat<T> c;     //! Domain at this timestep
  SimpleMat<T> c_nxt; //! Domain at next timestep
  T c_max, c_min;     //! Taken at initial condition
};

#endif // MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_SERIAL_ADVEC_DRIVER_HPP
