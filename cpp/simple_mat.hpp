#ifndef MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_SIMPLE_MAT_HPP
#define MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_SIMPLE_MAT_HPP

/* =============================================================
Data structure for a matrix

Allocate a new m by n matrix of doubles by:
    SimpleMat<double> mat(m,n);

Refer to element i, j by:
    mat.at[i][j]

Get a pointer to the continous data buffer:
    mat.data

The class has a working copy constructor, copy-assignment
operator, and destructor

============================================================= */

#include "math.h"
#include "stdio.h"
#include <cassert>
#include <cstring>

template<typename T>
class SimpleMat {
public:
  SimpleMat(int m = 1, int n = 1)
    : m(m)
    , n(n)
  {
    data = reinterpret_cast<T*>(malloc(m * n * sizeof(T)));
    at = reinterpret_cast<T**>(malloc(m * sizeof(T*)));
    for (int i = 0; i < m; ++i) {
      at[i] = &(data[i * n]);
    }
  }

  SimpleMat(const SimpleMat<T>& other)
    : m(other.m)
    , n(other.n)
  {
    data = reinterpret_cast<T*>(malloc(other.m * other.n * sizeof(T)));
    at = reinterpret_cast<T**>(malloc(other.m * sizeof(T*)));
    for (int i = 0; i < m; ++i) {
      at[i] = &(data[i * n]);
    }
    memcpy(data, other.data, other.m * other.n * sizeof(T));
  }

  SimpleMat& operator=(const SimpleMat<T>& rhs)
  {
    m = rhs.m;
    n = rhs.n;

    free(data);
    free(at);

    data = reinterpret_cast<T*>(malloc(m * n * sizeof(T)));
    memcpy(data, rhs.data, m * n * sizeof(T));

    at = reinterpret_cast<T**>(malloc(m * sizeof(T*)));
    for (int i = 0; i < m; ++i) {
      at[i] = &(data[i * n]);
    }
    return *this;
  }

  ~SimpleMat()
  {
    free(data);
    free(at);
  }

  int m;
  int n;
  T* data;
  T** at;
};

#endif // MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_SIMPLE_MAT_HPP
