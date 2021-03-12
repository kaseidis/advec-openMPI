#ifndef MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_C_SIMPLE_MAT_H
#define MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_C_SIMPLE_MAT_H

#include <stdlib.h>

/* =============================================================
Data structure for a matrix of doubles

Allocate a new m by n matrix by:
    SimpleMat mat = init_simple_mat(m, n);

Refer to element i, j by:
    mat.at[i][j]

Get a pointer to the continous data buffer:
    mat.data

Allocate a new matrix and initialize with the values of another:
    SimpleMat mat1(m, n);
    // ...
    // Initialize values of mat1
    // ...
    SimpleMat mat2 = copy_init_simple_mat(&mat1)

Copy the values from mat3 to mat4:
    copy_simple_mat(&mat4, &mat3);

Get the maximum element in a matrix:
    double max = max_simple_mat(&mat);

Free a matrix:
    free_simple_mat(&mat);

============================================================= */

typedef struct {
  int m;
  int n;
  double* data;
  double** at;
} SimpleMat;

SimpleMat init_simple_mat(int m, int n);

void free_simple_mat(SimpleMat* mat);

SimpleMat copy_init_simple_mat(const SimpleMat* other);

void copy_simple_mat(SimpleMat* dest, const SimpleMat* src);

double max_simple_mat(SimpleMat* mat);

double min_simple_mat(SimpleMat* mat);

#endif // MPCS_HPC_WINTER_2021_WEEK08_03_TOPOLOGIES_C_SIMPLE_MAT_H
