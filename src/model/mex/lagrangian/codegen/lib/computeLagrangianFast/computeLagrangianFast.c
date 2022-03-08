/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * computeLagrangianFast.c
 *
 * Code generation for function 'computeLagrangianFast'
 *
 */

/* Include files */
#include "computeLagrangianFast.h"
#include "computeLagrangianFast_emxutil.h"
#include "computeLagrangianFast_types.h"
#include "mtimes.h"
#include <string.h>

/* Function Declarations */
static void LagrangianODEX(const emxArray_real_T *x, const emxArray_real_T *dx,
                           const emxArray_real_T *Z1,
                           const emxArray_real_T *Theta, const double xia0[6],
                           const double gVec[3], const emxArray_real_T *Ba,
                           const double Mtt[36], const double Ktt[36],
                           emxArray_real_T *dZ1, emxArray_real_T *dZ2);

/* Function Definitions */
static void LagrangianODEX(const emxArray_real_T *x, const emxArray_real_T *dx,
                           const emxArray_real_T *Z1,
                           const emxArray_real_T *Theta, const double xia0[6],
                           const double gVec[3], const emxArray_real_T *Ba,
                           const double Mtt[36], const double Ktt[36],
                           emxArray_real_T *dZ1, emxArray_real_T *dZ2)
{
  emxArray_real_T *C;
  emxArray_real_T *Jg;
  emxArray_real_T *Jgt;
  emxArray_real_T *b_C;
  emxArray_real_T *c_C;
  emxArray_real_T *dM_tmp;
  emxArray_real_T *d_C;
  emxArray_real_T *r;
  emxArray_real_T *r1;
  emxArray_real_T *y_tmp;
  double A[36];
  double Ai[36];
  double adV[36];
  double b_A[36];
  double Rt[9];
  double dv[9];
  double V[6];
  double XI[6];
  double bkj;
  double d;
  double d1;
  int aoffset;
  int b_i;
  int boffset;
  int coffset;
  int i;
  int i1;
  int i2;
  int i3;
  int inner;
  int j;
  int k;
  int n;
  int nc;
  unsigned int u;
  if (5U > x->size[0] + 4U) {
    i = 0;
    i1 = -1;
  } else {
    i = 4;
    i1 = x->size[0] + 3;
  }
  bkj = 2.0 * ((double)x->size[0] - 1.0) + 6.0;
  if (((double)x->size[0] + 6.0) - 1.0 > bkj) {
    i2 = 0;
    i3 = 0;
  } else {
    i2 = x->size[0] + 4;
    i3 = (int)bkj;
  }
  emxInit_real_T(&y_tmp, 2);
  /* Theta_ = ThetaEval;%Model.ShpFnc(s); */
  inner = Ba->size[1];
  nc = Theta->size[1];
  n = y_tmp->size[0] * y_tmp->size[1];
  y_tmp->size[0] = 6;
  y_tmp->size[1] = Theta->size[1];
  emxEnsureCapacity_real_T(y_tmp, n);
  for (j = 0; j < nc; j++) {
    coffset = j * 6;
    boffset = j * Theta->size[0];
    for (b_i = 0; b_i < 6; b_i++) {
      y_tmp->data[coffset + b_i] = 0.0;
    }
    for (k = 0; k < inner; k++) {
      aoffset = k * 6;
      bkj = Theta->data[boffset + k];
      for (b_i = 0; b_i < 6; b_i++) {
        n = coffset + b_i;
        y_tmp->data[n] += Ba->data[aoffset + b_i] * bkj;
      }
    }
  }
  inner = y_tmp->size[1];
  for (b_i = 0; b_i < 6; b_i++) {
    XI[b_i] = 0.0;
  }
  for (k = 0; k < inner; k++) {
    aoffset = k * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      XI[b_i] += y_tmp->data[aoffset + b_i] * x->data[k];
    }
  }
  for (n = 0; n < 6; n++) {
    XI[n] += xia0[n];
  }
  /*  build forward kin - position */
  /* --------------------------------------------------------------------------
   */
  /* --------------------------------------------------------------------------
   */
  /* --------------------------------------------------------------------------
   */
  memset(&A[0], 0, 36U * sizeof(double));
  for (n = 0; n < 3; n++) {
    bkj = Z1->data[6 * n];
    A[6 * n] = bkj;
    nc = 6 * (n + 3);
    A[nc + 3] = bkj;
    inner = 6 * n + 1;
    bkj = Z1->data[inner];
    A[inner] = bkj;
    A[nc + 4] = bkj;
    inner = 6 * n + 2;
    bkj = Z1->data[inner];
    A[inner] = bkj;
    A[nc + 5] = bkj;
  }
  dv[0] = 0.0;
  dv[3] = -Z1->data[20];
  dv[6] = Z1->data[19];
  dv[1] = Z1->data[20];
  dv[4] = 0.0;
  dv[7] = -Z1->data[18];
  dv[2] = -Z1->data[19];
  dv[5] = Z1->data[18];
  dv[8] = 0.0;
  /* --------------------------------------------------------------------------
   */
  for (n = 0; n < 3; n++) {
    bkj = dv[n];
    d = dv[n + 3];
    d1 = dv[n + 6];
    for (inner = 0; inner < 3; inner++) {
      A[(n + 6 * inner) + 3] =
          (bkj * Z1->data[6 * inner] + d * Z1->data[6 * inner + 1]) +
          d1 * Z1->data[6 * inner + 2];
    }
    Rt[3 * n] = Z1->data[n];
    Rt[3 * n + 1] = Z1->data[n + 6];
    Rt[3 * n + 2] = Z1->data[n + 12];
  }
  /* --------------------------------------------------------------------------
   */
  memset(&Ai[0], 0, 36U * sizeof(double));
  for (n = 0; n < 3; n++) {
    bkj = Rt[3 * n];
    Ai[6 * n] = bkj;
    inner = 6 * (n + 3);
    Ai[inner + 3] = bkj;
    bkj = Rt[3 * n + 1];
    Ai[6 * n + 1] = bkj;
    Ai[inner + 4] = bkj;
    bkj = Rt[3 * n + 2];
    Ai[6 * n + 2] = bkj;
    Ai[inner + 5] = bkj;
  }
  dv[0] = 0.0;
  dv[1] = -Z1->data[20];
  dv[2] = Z1->data[19];
  dv[3] = Z1->data[20];
  dv[4] = 0.0;
  dv[5] = -Z1->data[18];
  dv[6] = -Z1->data[19];
  dv[7] = Z1->data[18];
  dv[8] = 0.0;
  for (n = 0; n < 3; n++) {
    bkj = Rt[n];
    d = Rt[n + 3];
    d1 = Rt[n + 6];
    for (inner = 0; inner < 3; inner++) {
      Ai[(n + 6 * inner) + 3] = (bkj * dv[3 * inner] + d * dv[3 * inner + 1]) +
                                d1 * dv[3 * inner + 2];
    }
  }
  emxInit_real_T(&Jg, 2);
  /*  build jacobian */
  inner = i1 - i;
  i1 = Jg->size[0] * Jg->size[1];
  Jg->size[0] = 6;
  Jg->size[1] = inner + 1;
  emxEnsureCapacity_real_T(Jg, i1);
  for (j = 0; j <= inner; j++) {
    nc = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (k = 0; k < 6; k++) {
        i1 = nc + k;
        bkj += Ai[k * 6 + b_i] * Z1->data[i1 % 6 + 6 * (i + i1 / 6)];
      }
      Jg->data[nc + b_i] = bkj;
    }
  }
  emxInit_real_T(&Jgt, 2);
  inner = i3 - i2;
  n = inner - 1;
  i = Jgt->size[0] * Jgt->size[1];
  Jgt->size[0] = 6;
  Jgt->size[1] = inner;
  emxEnsureCapacity_real_T(Jgt, i);
  for (j = 0; j <= n; j++) {
    nc = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (k = 0; k < 6; k++) {
        i = nc + k;
        bkj += Ai[k * 6 + b_i] * Z1->data[i % 6 + 6 * (i2 + i / 6)];
      }
      Jgt->data[nc + b_i] = bkj;
    }
  }
  inner = Jg->size[1];
  for (b_i = 0; b_i < 6; b_i++) {
    V[b_i] = 0.0;
  }
  for (k = 0; k < inner; k++) {
    aoffset = k * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      V[b_i] += Jg->data[aoffset + b_i] * dx->data[k];
    }
  }
  /* --------------------------------------------------------------------------
   */
  memset(&adV[0], 0, 36U * sizeof(double));
  /* --------------------------------------------------------------------------
   */
  Rt[0] = 0.0;
  Rt[3] = -V[2];
  Rt[6] = V[1];
  Rt[1] = V[2];
  Rt[4] = 0.0;
  Rt[7] = -V[0];
  Rt[2] = -V[1];
  Rt[5] = V[0];
  Rt[8] = 0.0;
  /* --------------------------------------------------------------------------
   */
  for (i = 0; i < 3; i++) {
    bkj = Rt[3 * i];
    adV[6 * i] = bkj;
    inner = 6 * (i + 3);
    adV[inner + 3] = bkj;
    bkj = Rt[3 * i + 1];
    adV[6 * i + 1] = bkj;
    adV[inner + 4] = bkj;
    bkj = Rt[3 * i + 2];
    adV[6 * i + 2] = bkj;
    adV[inner + 5] = bkj;
  }
  emxInit_real_T(&dM_tmp, 2);
  adV[3] = 0.0;
  adV[9] = -V[5];
  adV[15] = V[4];
  adV[4] = V[5];
  adV[10] = 0.0;
  adV[16] = -V[3];
  adV[5] = -V[4];
  adV[11] = V[3];
  adV[17] = 0.0;
  /*  compute inertia, coriolis, gravity */
  mtimes(Jg, Mtt, dM_tmp);
  /*  compute (nonlinear stiffness) */
  /*  compute grav. potential energy */
  i = dZ1->size[0] * dZ1->size[1];
  dZ1->size[0] = 6;
  dZ1->size[1] = (int)(2.0 * ((double)x->size[0] - 1.0) + 6.0);
  emxEnsureCapacity_real_T(dZ1, i);
  inner = 6 * (int)(2.0 * ((double)x->size[0] - 1.0) + 6.0);
  for (i = 0; i < inner; i++) {
    dZ1->data[i] = 0.0;
  }
  dv[0] = 0.0;
  dv[3] = -XI[2];
  dv[6] = XI[1];
  dv[1] = XI[2];
  dv[4] = 0.0;
  dv[7] = -XI[0];
  dv[2] = -XI[1];
  dv[5] = XI[0];
  dv[8] = 0.0;
  for (i = 0; i < 3; i++) {
    dZ1->data[i + 18] = 0.0;
    for (i1 = 0; i1 < 3; i1++) {
      i2 = i + 6 * i1;
      dZ1->data[i2] =
          (Z1->data[i] * dv[3 * i1] + Z1->data[i + 6] * dv[3 * i1 + 1]) +
          Z1->data[i + 12] * dv[3 * i1 + 2];
      dZ1->data[i + 18] += Z1->data[i2] * XI[i1 + 3];
    }
  }
  if (5U > x->size[0] + 4U) {
    i = 0;
  } else {
    i = 4;
  }
  emxInit_real_T(&C, 2);
  n = y_tmp->size[1];
  i1 = C->size[0] * C->size[1];
  C->size[0] = 6;
  C->size[1] = y_tmp->size[1];
  emxEnsureCapacity_real_T(C, i1);
  for (j = 0; j < n; j++) {
    nc = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (k = 0; k < 6; k++) {
        bkj += A[k * 6 + b_i] * y_tmp->data[nc + k];
      }
      C->data[nc + b_i] = bkj;
    }
  }
  inner = C->size[1];
  for (i1 = 0; i1 < inner; i1++) {
    for (i2 = 0; i2 < 6; i2++) {
      dZ1->data[i2 + 6 * (i + i1)] = C->data[i2 + 6 * i1];
    }
  }
  if (((double)x->size[0] + 6.0) - 1.0 >
      2.0 * ((double)x->size[0] - 1.0) + 6.0) {
    i = -4;
  } else {
    i = x->size[0];
  }
  for (i1 = 0; i1 < 6; i1++) {
    for (i2 = 0; i2 < 6; i2++) {
      bkj = 0.0;
      for (i3 = 0; i3 < 6; i3++) {
        bkj += A[i1 + 6 * i3] * adV[i3 + 6 * i2];
      }
      b_A[i1 + 6 * i2] = bkj;
    }
  }
  n = y_tmp->size[1];
  i1 = C->size[0] * C->size[1];
  C->size[0] = 6;
  C->size[1] = y_tmp->size[1];
  emxEnsureCapacity_real_T(C, i1);
  for (j = 0; j < n; j++) {
    nc = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (k = 0; k < 6; k++) {
        bkj += b_A[k * 6 + b_i] * y_tmp->data[nc + k];
      }
      C->data[nc + b_i] = bkj;
    }
  }
  inner = C->size[1];
  for (i1 = 0; i1 < inner; i1++) {
    for (i2 = 0; i2 < 6; i2++) {
      dZ1->data[i2 + 6 * ((i + i1) + 4)] = C->data[i2 + 6 * i1];
    }
  }
  dZ1->data[22] =
      (Mtt[21] * Z1->data[18] * -gVec[0] + Mtt[21] * Z1->data[19] * -gVec[1]) +
      Mtt[21] * Z1->data[20] * -gVec[2];
  bkj = 0.0;
  for (i = 0; i < 6; i++) {
    d = 0.0;
    for (i1 = 0; i1 < 6; i1++) {
      d += 0.5 * V[i1] * Mtt[i1 + 6 * i];
    }
    bkj += d * V[i];
  }
  dZ1->data[23] = bkj;
  i = dZ2->size[0] * dZ2->size[1];
  dZ2->size[0] = x->size[0];
  dZ2->size[1] = (int)(4.0 * (double)x->size[0] + 1.0);
  emxEnsureCapacity_real_T(dZ2, i);
  inner = x->size[0] * (int)(4.0 * (double)x->size[0] + 1.0);
  for (i = 0; i < inner; i++) {
    dZ2->data[i] = 0.0;
  }
  emxInit_real_T(&r, 2);
  b_mtimes(dM_tmp, Jg, r);
  inner = r->size[1];
  for (i = 0; i < inner; i++) {
    nc = r->size[0];
    for (i1 = 0; i1 < nc; i1++) {
      dZ2->data[i1 + dZ2->size[0] * i] = r->data[i1 + r->size[0] * i];
    }
  }
  if (x->size[0] + 1U > ((unsigned int)x->size[0] << 1)) {
    i = 0;
  } else {
    i = x->size[0];
  }
  for (i1 = 0; i1 < 6; i1++) {
    for (i2 = 0; i2 < 6; i2++) {
      bkj = 0.0;
      d = 0.0;
      for (i3 = 0; i3 < 6; i3++) {
        n = i3 + 6 * i2;
        bkj += Mtt[i1 + 6 * i3] * adV[n];
        d += adV[i3 + 6 * i1] * Mtt[n];
      }
      nc = i1 + 6 * i2;
      b_A[nc] = d;
      A[nc] = bkj;
    }
  }
  for (i1 = 0; i1 < 36; i1++) {
    A[i1] -= b_A[i1];
  }
  n = Jg->size[1];
  i1 = C->size[0] * C->size[1];
  C->size[0] = 6;
  C->size[1] = Jg->size[1];
  emxEnsureCapacity_real_T(C, i1);
  for (j = 0; j < n; j++) {
    nc = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (k = 0; k < 6; k++) {
        bkj += A[k * 6 + b_i] * Jg->data[nc + k];
      }
      C->data[nc + b_i] = bkj;
    }
  }
  emxInit_real_T(&b_C, 2);
  n = Jgt->size[1];
  i1 = b_C->size[0] * b_C->size[1];
  b_C->size[0] = 6;
  b_C->size[1] = Jgt->size[1];
  emxEnsureCapacity_real_T(b_C, i1);
  for (j = 0; j < n; j++) {
    nc = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (k = 0; k < 6; k++) {
        bkj += Mtt[k * 6 + b_i] * Jgt->data[nc + k];
      }
      b_C->data[nc + b_i] = bkj;
    }
  }
  inner = 6 * C->size[1];
  i1 = C->size[0] * C->size[1];
  C->size[0] = 6;
  emxEnsureCapacity_real_T(C, i1);
  for (i1 = 0; i1 < inner; i1++) {
    C->data[i1] += b_C->data[i1];
  }
  emxFree_real_T(&b_C);
  emxInit_real_T(&c_C, 2);
  inner = Jg->size[1];
  n = C->size[1];
  i1 = c_C->size[0] * c_C->size[1];
  c_C->size[0] = Jg->size[1];
  c_C->size[1] = C->size[1];
  emxEnsureCapacity_real_T(c_C, i1);
  for (j = 0; j < n; j++) {
    coffset = j * inner;
    boffset = j * 6;
    for (b_i = 0; b_i < inner; b_i++) {
      aoffset = b_i * 6;
      bkj = 0.0;
      for (k = 0; k < 6; k++) {
        bkj += Jg->data[aoffset + k] * C->data[boffset + k];
      }
      c_C->data[coffset + b_i] = bkj;
    }
  }
  emxFree_real_T(&C);
  inner = c_C->size[1];
  for (i1 = 0; i1 < inner; i1++) {
    nc = c_C->size[0];
    for (i2 = 0; i2 < nc; i2++) {
      dZ2->data[i2 + dZ2->size[0] * (i + i1)] =
          c_C->data[i2 + c_C->size[0] * i1];
    }
  }
  u = ((unsigned int)x->size[0] << 1) + 1U;
  if (u > 3.0 * (double)x->size[0]) {
    i = 1;
  } else {
    i = (int)u;
  }
  emxInit_real_T(&r1, 2);
  mtimes(y_tmp, Ktt, r1);
  b_mtimes(r1, y_tmp, r);
  inner = r->size[1];
  emxFree_real_T(&y_tmp);
  for (i1 = 0; i1 < inner; i1++) {
    nc = r->size[0];
    for (i2 = 0; i2 < nc; i2++) {
      dZ2->data[i2 + dZ2->size[0] * ((i + i1) - 1)] =
          r->data[i2 + r->size[0] * i1];
    }
  }
  bkj = 3.0 * (double)x->size[0] + 1.0;
  if (bkj > 4.0 * (double)x->size[0]) {
    i = 1;
  } else {
    i = (int)bkj;
  }
  mtimes(Jgt, Mtt, r1);
  b_mtimes(r1, Jg, c_C);
  b_mtimes(dM_tmp, Jgt, r);
  inner = c_C->size[1];
  emxFree_real_T(&r1);
  emxFree_real_T(&dM_tmp);
  emxFree_real_T(&Jgt);
  for (i1 = 0; i1 < inner; i1++) {
    nc = c_C->size[0];
    for (i2 = 0; i2 < nc; i2++) {
      dZ2->data[i2 + dZ2->size[0] * ((i + i1) - 1)] =
          c_C->data[i2 + c_C->size[0] * i1] + r->data[i2 + r->size[0] * i1];
    }
  }
  emxFree_real_T(&r);
  emxFree_real_T(&c_C);
  i = (int)(4.0 * (double)x->size[0] + 1.0) - 1;
  V[0] = 0.0;
  V[1] = 0.0;
  V[2] = 0.0;
  V[3] = -gVec[0];
  V[4] = -gVec[1];
  V[5] = -gVec[2];
  for (i1 = 0; i1 < 6; i1++) {
    bkj = 0.0;
    for (i2 = 0; i2 < 6; i2++) {
      d = 0.0;
      for (i3 = 0; i3 < 6; i3++) {
        d += Ai[i1 + 6 * i3] * Mtt[i3 + 6 * i2];
      }
      bkj += d * V[i2];
    }
    XI[i1] = bkj;
  }
  emxInit_real_T(&d_C, 1);
  inner = Jg->size[1];
  i1 = d_C->size[0];
  d_C->size[0] = Jg->size[1];
  emxEnsureCapacity_real_T(d_C, i1);
  for (b_i = 0; b_i < inner; b_i++) {
    aoffset = b_i * 6;
    bkj = 0.0;
    for (k = 0; k < 6; k++) {
      bkj += Jg->data[aoffset + k] * XI[k];
    }
    d_C->data[b_i] = bkj;
  }
  emxFree_real_T(&Jg);
  inner = d_C->size[0];
  for (i1 = 0; i1 < inner; i1++) {
    dZ2->data[i1 + dZ2->size[0] * i] = d_C->data[i1];
  }
  emxFree_real_T(&d_C);
}

void computeLagrangianFast(
    const emxArray_real_T *x, const emxArray_real_T *dx, double ds,
    const double p0[3], const double Phi0[9], const emxArray_real_T *xia0,
    const emxArray_real_T *Th, const emxArray_real_T *Ba, const double gVec[3],
    const double Ktt[36], const double Mtt[36], double Zeta, emxArray_real_T *M,
    emxArray_real_T *C, emxArray_real_T *K, emxArray_real_T *R,
    emxArray_real_T *G, double p[3], double Phi[9], emxArray_real_T *J,
    double *Vg, double *Kin, emxArray_real_T *Jt, emxArray_real_T *Mt)
{
  emxArray_real_T *K1Z1;
  emxArray_real_T *K1Z2;
  emxArray_real_T *K2Z1;
  emxArray_real_T *K2Z2;
  emxArray_real_T *Z1;
  emxArray_real_T *Z2;
  emxArray_real_T *b_Th;
  emxArray_real_T *b_Z1;
  double tmp[36];
  double Rt[9];
  double dv[9];
  double d;
  double d1;
  double s;
  int b_i;
  int b_loop_ub;
  int i;
  int i1;
  int i2;
  int i3;
  int ii;
  int k;
  int loop_ub;
  int n_tmp;
  unsigned int u;
  emxInit_real_T(&Z1, 2);
  /*  % states */
  /*       % spatial steps */
  /*       % position zero */
  /*     % phi zero */
  /*     % intrinsic strain vector */
  /*       % evaluated Theta matrix */
  /*       % state to strain matrix */
  /*     % gravitational acceleration (in mili-g) */
  /*      % geometric stiffness */
  /*      % geometric inertia */
  /*  compute total length */
  i = Z1->size[0] * Z1->size[1];
  Z1->size[0] = 6;
  Z1->size[1] = (int)(2.0 * ((double)x->size[0] - 1.0) + 6.0);
  emxEnsureCapacity_real_T(Z1, i);
  loop_ub = 6 * (int)(2.0 * ((double)x->size[0] - 1.0) + 6.0);
  for (i = 0; i < loop_ub; i++) {
    Z1->data[i] = 0.0;
  }
  emxInit_real_T(&Z2, 2);
  i = Z2->size[0] * Z2->size[1];
  Z2->size[0] = x->size[0];
  Z2->size[1] = (int)(4.0 * (double)x->size[0] + 1.0);
  emxEnsureCapacity_real_T(Z2, i);
  loop_ub = x->size[0] * (int)(4.0 * (double)x->size[0] + 1.0);
  for (i = 0; i < loop_ub; i++) {
    Z2->data[i] = 0.0;
  }
  for (i = 0; i < 3; i++) {
    Z1->data[6 * i] = Phi0[3 * i];
    Z1->data[6 * i + 1] = Phi0[3 * i + 1];
    Z1->data[6 * i + 2] = Phi0[3 * i + 2];
    Z1->data[i + 18] = p0[i];
  }
  /* NLStiff = false;  */
  i = (int)((double)Th->size[2] / 2.0);
  emxInit_real_T(&K1Z1, 2);
  emxInit_real_T(&K1Z2, 2);
  emxInit_real_T(&K2Z1, 2);
  emxInit_real_T(&K2Z2, 2);
  emxInit_real_T(&b_Th, 2);
  emxInit_real_T(&b_Z1, 2);
  for (ii = 0; ii < i; ii++) {
    /*  first EL-diff eval */
    i1 = (int)((unsigned int)(ii + 1) << 1);
    loop_ub = Th->size[0];
    b_loop_ub = Th->size[1];
    i2 = b_Th->size[0] * b_Th->size[1];
    b_Th->size[0] = Th->size[0];
    b_Th->size[1] = Th->size[1];
    emxEnsureCapacity_real_T(b_Th, i2);
    for (i2 = 0; i2 < b_loop_ub; i2++) {
      for (i3 = 0; i3 < loop_ub; i3++) {
        b_Th->data[i3 + b_Th->size[0] * i2] =
            Th->data[(i3 + Th->size[0] * i2) +
                     Th->size[0] * Th->size[1] * (i1 - 2)];
      }
    }
    LagrangianODEX(x, dx, Z1, b_Th,
                   *(double(*)[6]) &
                       xia0->data[6 * ((int)((unsigned int)(ii + 1) << 1) - 2)],
                   gVec, Ba, Mtt, Ktt, K1Z1, K1Z2);
    /*  second EL-diff eval */
    s = 0.66666666666666663 * ds;
    i2 = b_Z1->size[0] * b_Z1->size[1];
    b_Z1->size[0] = 6;
    b_Z1->size[1] = Z1->size[1];
    emxEnsureCapacity_real_T(b_Z1, i2);
    loop_ub = 6 * Z1->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      b_Z1->data[i2] = Z1->data[i2] + s * K1Z1->data[i2];
    }
    loop_ub = Th->size[0];
    b_loop_ub = Th->size[1];
    i2 = b_Th->size[0] * b_Th->size[1];
    b_Th->size[0] = Th->size[0];
    b_Th->size[1] = Th->size[1];
    emxEnsureCapacity_real_T(b_Th, i2);
    for (i2 = 0; i2 < b_loop_ub; i2++) {
      for (i3 = 0; i3 < loop_ub; i3++) {
        b_Th->data[i3 + b_Th->size[0] * i2] =
            Th->data[(i3 + Th->size[0] * i2) +
                     Th->size[0] * Th->size[1] * (i1 - 1)];
      }
    }
    LagrangianODEX(x, dx, b_Z1, b_Th,
                   *(double(*)[6]) &
                       xia0->data[6 * ((int)((unsigned int)(ii + 1) << 1) - 1)],
                   gVec, Ba, Mtt, Ktt, K2Z1, K2Z2);
    /*  update integrands */
    s = 0.25 * ds;
    loop_ub = 6 * Z1->size[1];
    i1 = Z1->size[0] * Z1->size[1];
    Z1->size[0] = 6;
    emxEnsureCapacity_real_T(Z1, i1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      Z1->data[i1] += s * (K1Z1->data[i1] + 3.0 * K2Z1->data[i1]);
    }
    loop_ub = Z2->size[0] * Z2->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      Z2->data[i1] += s * (K1Z2->data[i1] + 3.0 * K2Z2->data[i1]);
    }
  }
  emxFree_real_T(&b_Z1);
  emxFree_real_T(&b_Th);
  emxFree_real_T(&K2Z2);
  emxFree_real_T(&K2Z1);
  emxFree_real_T(&K1Z2);
  emxFree_real_T(&K1Z1);
  /*  recover the kinematics entities */
  /* --------------------------------------------------------------------------
   */
  for (i = 0; i < 3; i++) {
    p[i] = Z1->data[i + 18];
    Phi[3 * i] = Z1->data[6 * i];
    Rt[3 * i] = Z1->data[i];
    ii = 3 * i + 1;
    Phi[ii] = Z1->data[6 * i + 1];
    Rt[ii] = Z1->data[i + 6];
    ii = 3 * i + 2;
    Phi[ii] = Z1->data[6 * i + 2];
    Rt[ii] = Z1->data[i + 12];
  }
  /* --------------------------------------------------------------------------
   */
  memset(&tmp[0], 0, 36U * sizeof(double));
  for (i = 0; i < 3; i++) {
    s = Rt[3 * i];
    tmp[6 * i] = s;
    ii = 6 * (i + 3);
    tmp[ii + 3] = s;
    s = Rt[3 * i + 1];
    tmp[6 * i + 1] = s;
    tmp[ii + 4] = s;
    s = Rt[3 * i + 2];
    tmp[6 * i + 2] = s;
    tmp[ii + 5] = s;
  }
  dv[0] = 0.0;
  dv[1] = -Z1->data[20];
  dv[2] = Z1->data[19];
  dv[3] = Z1->data[20];
  dv[4] = 0.0;
  dv[5] = -Z1->data[18];
  dv[6] = -Z1->data[19];
  dv[7] = Z1->data[18];
  dv[8] = 0.0;
  for (i = 0; i < 3; i++) {
    s = Rt[i];
    d = Rt[i + 3];
    d1 = Rt[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      tmp[(i + 6 * i1) + 3] =
          (s * dv[3 * i1] + d * dv[3 * i1 + 1]) + d1 * dv[3 * i1 + 2];
    }
  }
  if (5U > x->size[0] + 4U) {
    i = 0;
    i1 = -1;
  } else {
    i = 4;
    i1 = x->size[0] + 3;
  }
  s = 2.0 * ((double)x->size[0] - 1.0) + 6.0;
  if (((double)x->size[0] + 6.0) - 1.0 > s) {
    i2 = 0;
    i3 = 0;
  } else {
    i2 = x->size[0] + 4;
    i3 = (int)s;
  }
  n_tmp = i1 - i;
  i1 = J->size[0] * J->size[1];
  J->size[0] = 6;
  J->size[1] = n_tmp + 1;
  emxEnsureCapacity_real_T(J, i1);
  for (b_loop_ub = 0; b_loop_ub <= n_tmp; b_loop_ub++) {
    loop_ub = b_loop_ub * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      s = 0.0;
      for (k = 0; k < 6; k++) {
        i1 = loop_ub + k;
        s += tmp[k * 6 + b_i] * Z1->data[i1 % 6 + 6 * (i + i1 / 6)];
      }
      J->data[loop_ub + b_i] = s;
    }
  }
  n_tmp = i3 - i2;
  ii = n_tmp - 1;
  i = Jt->size[0] * Jt->size[1];
  Jt->size[0] = 6;
  Jt->size[1] = n_tmp;
  emxEnsureCapacity_real_T(Jt, i);
  for (b_loop_ub = 0; b_loop_ub <= ii; b_loop_ub++) {
    loop_ub = b_loop_ub * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      s = 0.0;
      for (k = 0; k < 6; k++) {
        i = loop_ub + k;
        s += tmp[k * 6 + b_i] * Z1->data[i % 6 + 6 * (i2 + i / 6)];
      }
      Jt->data[loop_ub + b_i] = s;
    }
  }
  /*  recover the dynamics entities */
  if (1 > x->size[0]) {
    loop_ub = 0;
    b_loop_ub = 0;
  } else {
    loop_ub = x->size[0];
    b_loop_ub = x->size[0];
  }
  i = M->size[0] * M->size[1];
  M->size[0] = loop_ub;
  M->size[1] = b_loop_ub;
  emxEnsureCapacity_real_T(M, i);
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      M->data[i1 + M->size[0] * i] = Z2->data[i1 + Z2->size[0] * i];
    }
  }
  if (1 > x->size[0]) {
    loop_ub = 0;
  } else {
    loop_ub = x->size[0];
  }
  u = (unsigned int)x->size[0] << 1;
  if (x->size[0] + 1U > u) {
    i = 0;
    i1 = 0;
  } else {
    i = x->size[0];
    i1 = (int)u;
  }
  i2 = C->size[0] * C->size[1];
  C->size[0] = loop_ub;
  b_loop_ub = i1 - i;
  C->size[1] = b_loop_ub;
  emxEnsureCapacity_real_T(C, i2);
  for (i1 = 0; i1 < b_loop_ub; i1++) {
    for (i2 = 0; i2 < loop_ub; i2++) {
      C->data[i2 + C->size[0] * i1] = Z2->data[i2 + Z2->size[0] * (i + i1)];
    }
  }
  if (1 > x->size[0]) {
    loop_ub = 0;
  } else {
    loop_ub = x->size[0];
  }
  u = ((unsigned int)x->size[0] << 1) + 1U;
  s = 3.0 * (double)x->size[0];
  if (u > s) {
    i = 0;
    i1 = 0;
  } else {
    i = (int)u - 1;
    i1 = (int)s;
  }
  i2 = K->size[0] * K->size[1];
  K->size[0] = loop_ub;
  b_loop_ub = i1 - i;
  K->size[1] = b_loop_ub;
  emxEnsureCapacity_real_T(K, i2);
  for (i1 = 0; i1 < b_loop_ub; i1++) {
    for (i2 = 0; i2 < loop_ub; i2++) {
      K->data[i2 + K->size[0] * i1] = Z2->data[i2 + Z2->size[0] * (i + i1)];
    }
  }
  if (1 > x->size[0]) {
    ii = 0;
  } else {
    ii = x->size[0];
  }
  s = 3.0 * (double)x->size[0] + 1.0;
  d = 4.0 * (double)x->size[0];
  if (s > d) {
    i1 = 0;
    i2 = 0;
  } else {
    i1 = (int)s - 1;
    i2 = (int)d;
  }
  i3 = Mt->size[0] * Mt->size[1];
  Mt->size[0] = ii;
  n_tmp = i2 - i1;
  Mt->size[1] = n_tmp;
  emxEnsureCapacity_real_T(Mt, i3);
  for (i2 = 0; i2 < n_tmp; i2++) {
    for (i3 = 0; i3 < ii; i3++) {
      Mt->data[i3 + Mt->size[0] * i2] = Z2->data[i3 + Z2->size[0] * (i1 + i2)];
    }
  }
  if (1 > x->size[0]) {
    ii = 0;
  } else {
    ii = x->size[0];
  }
  i1 = (int)(4.0 * (double)x->size[0] + 1.0);
  i2 = G->size[0];
  G->size[0] = ii;
  emxEnsureCapacity_real_T(G, i2);
  for (i2 = 0; i2 < ii; i2++) {
    G->data[i2] = Z2->data[i2 + Z2->size[0] * (i1 - 1)];
  }
  *Vg = Z1->data[22];
  *Kin = Z1->data[23];
  i1 = R->size[0] * R->size[1];
  R->size[0] = loop_ub;
  R->size[1] = b_loop_ub;
  emxEnsureCapacity_real_T(R, i1);
  emxFree_real_T(&Z1);
  for (i1 = 0; i1 < b_loop_ub; i1++) {
    for (i2 = 0; i2 < loop_ub; i2++) {
      R->data[i2 + R->size[0] * i1] =
          Zeta * Z2->data[i2 + Z2->size[0] * (i + i1)];
    }
  }
  emxFree_real_T(&Z2);
}

/* End of code generation (computeLagrangianFast.c) */
