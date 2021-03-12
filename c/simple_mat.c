#include "simple_mat.h"
#include <float.h>
#include <string.h>

SimpleMat init_simple_mat(int m, int n)
{
  SimpleMat mat;
  mat.m = m;
  mat.n = n;
  mat.data = (double*)malloc(m * n * sizeof(double));
  mat.at = (double**)malloc(m * sizeof(double*));
  for (int i = 0; i < m; ++i) {
    mat.at[i] = &(mat.data[i * n]);
  }
  return mat;
}

void free_simple_mat(SimpleMat* mat)
{
  free(mat->data);
  free(mat->at);
  mat->m = 0;
  mat->n = 0;
}

SimpleMat copy_init_simple_mat(const SimpleMat* other)
{
  SimpleMat res = init_simple_mat(other->m, other->n);
  memcpy(res.data, other->data, other->m * other->n * sizeof(double));
  return res;
}

void copy_simple_mat(SimpleMat* dest, const SimpleMat* src)
{
  free_simple_mat(dest);

  dest->m = src->m;
  dest->n = src->n;
  dest->data = (double*)malloc(src->m * src->n * sizeof(double));
  dest->at = (double**)malloc(src->m * sizeof(double*));

  memcpy(dest->data, src->data, src->m * src->n * sizeof(double));

  for (int i = 0; i < dest->m; ++i) {
    dest->at[i] = &(dest->data[i * dest->n]);
  }
}

double max_simple_mat(SimpleMat* mat)
{
  double max = DBL_MIN;
  for (int i = 0; i < mat->m * mat->n; ++i) {
    double val = mat->data[i];
    max = val > max ? val : max;
  }
  return max;
}

double min_simple_mat(SimpleMat* mat)
{
  double min = DBL_MAX;
  for (int i = 0; i < mat->m * mat->n; ++i) {
    double val = mat->data[i];
    min = val < min ? val : min;
  }
  return min;
}
