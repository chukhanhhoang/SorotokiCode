/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mtimes.c
 *
 * Code generation for function 'mtimes'
 *
 */

/* Include files */
#include "mtimes.h"
#include "computeLagrangianFast_emxutil.h"
#include "computeLagrangianFast_types.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void b_mtimes(const emxArray_real_T *A, const emxArray_real_T *B,
              emxArray_real_T *C)
{
  const double *A_data;
  const double *B_data;
  double *C_data;
  int coffset;
  int i;
  int j;
  int k;
  int m;
  int n;
  B_data = B->data;
  A_data = A->data;
  m = A->size[0];
  n = B->size[1];
  coffset = C->size[0] * C->size[1];
  C->size[0] = A->size[0];
  C->size[1] = B->size[1];
  emxEnsureCapacity_real_T(C, coffset);
  C_data = C->data;
  for (j = 0; j < n; j++) {
    int boffset;
    coffset = j * m;
    boffset = j * 6;
    for (i = 0; i < m; i++) {
      double s;
      s = 0.0;
      for (k = 0; k < 6; k++) {
        s += A_data[k * A->size[0] + i] * B_data[boffset + k];
      }
      C_data[coffset + i] = s;
    }
  }
}

void mtimes(const emxArray_real_T *A, const double B[36], emxArray_real_T *C)
{
  const double *A_data;
  double *C_data;
  int coffset;
  int i;
  int j;
  int k;
  int m;
  A_data = A->data;
  m = A->size[1];
  coffset = C->size[0] * C->size[1];
  C->size[0] = A->size[1];
  C->size[1] = 6;
  emxEnsureCapacity_real_T(C, coffset);
  C_data = C->data;
  for (j = 0; j < 6; j++) {
    int boffset;
    coffset = j * m;
    boffset = j * 6;
    for (i = 0; i < m; i++) {
      double s;
      int aoffset;
      aoffset = i * 6;
      s = 0.0;
      for (k = 0; k < 6; k++) {
        s += A_data[aoffset + k] * B[boffset + k];
      }
      C_data[coffset + i] = s;
    }
  }
}

/* End of code generation (mtimes.c) */
