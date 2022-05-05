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
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static void LagrangianODEX(const emxArray_real_T *x, const emxArray_real_T *dx,
                           const emxArray_real_T *Z1,
                           const emxArray_real_T *Theta, const double xia0[6],
                           const emxArray_real_T *Ba, const double Mtt[36],
                           const double Ktt[36], const double Gvec[3], double k,
                           double npe, const double r0[3], double rs,
                           emxArray_real_T *dZ1, emxArray_real_T *dZ2);

static void b_binary_expand_op(emxArray_real_T *in1, double in2,
                               const emxArray_real_T *in3,
                               const emxArray_real_T *in4);

static void binary_expand_op(emxArray_real_T *in1, double in2,
                             const emxArray_real_T *in3,
                             const emxArray_real_T *in4);

static void c_binary_expand_op(
    const emxArray_real_T *in1, const emxArray_real_T *in2,
    const emxArray_real_T *in3, double in4, const emxArray_real_T *in5,
    const emxArray_real_T *in6, int in7, const emxArray_real_T *in8,
    const emxArray_real_T *in9, const double in10[36], const double in11[36],
    const double in12[3], double in13, double in14, const double in15[3],
    double in16, emxArray_real_T *in17, emxArray_real_T *in18);

static void plus(emxArray_real_T *in1, const emxArray_real_T *in2);

static double rt_powd_snf(double u0, double u1);

/* Function Definitions */
static void LagrangianODEX(const emxArray_real_T *x, const emxArray_real_T *dx,
                           const emxArray_real_T *Z1,
                           const emxArray_real_T *Theta, const double xia0[6],
                           const emxArray_real_T *Ba, const double Mtt[36],
                           const double Ktt[36], const double Gvec[3], double k,
                           double npe, const double r0[3], double rs,
                           emxArray_real_T *dZ1, emxArray_real_T *dZ2)
{
  emxArray_real_T *C;
  emxArray_real_T *Jg;
  emxArray_real_T *b_C;
  emxArray_real_T *c_C;
  emxArray_real_T *dG;
  emxArray_real_T *d_C;
  emxArray_real_T *r;
  emxArray_real_T *y_tmp;
  double A[36];
  double Ai[36];
  double adV[36];
  double b_A[36];
  double Phi_[9];
  double Rt_tmp[9];
  double Wh[9];
  double V[6];
  double XI[6];
  double b[6];
  double y[6];
  const double *Ba_data;
  const double *Theta_data;
  const double *Z1_data;
  const double *dx_data;
  const double *x_data;
  double a_tmp;
  double abs_rd;
  double absxk;
  double bkj;
  double d;
  double d1;
  double rd_idx_0;
  double rd_idx_1;
  double t;
  double *C_data;
  double *Jg_data;
  double *dG_data;
  double *dZ1_data;
  double *y_tmp_data;
  int aoffset;
  int b_i;
  int b_k;
  int boffset;
  int coffset;
  int i;
  int i1;
  int i2;
  int i3;
  int inner;
  int j;
  int n;
  int nc;
  unsigned int u;
  Ba_data = Ba->data;
  Theta_data = Theta->data;
  Z1_data = Z1->data;
  dx_data = dx->data;
  x_data = x->data;
  for (i = 0; i < 3; i++) {
    Phi_[3 * i] = Z1_data[6 * i];
    Phi_[3 * i + 1] = Z1_data[6 * i + 1];
    Phi_[3 * i + 2] = Z1_data[6 * i + 2];
  }
  if (x->size[0] + 4U < 5U) {
    i = 0;
    i1 = -1;
  } else {
    i = 4;
    i1 = x->size[0] + 3;
  }
  d = 2.0 * ((double)x->size[0] - 1.0) + 6.0;
  if (((double)x->size[0] + 6.0) - 1.0 > d) {
    i2 = 0;
    i3 = 0;
  } else {
    i2 = x->size[0] + 4;
    i3 = (int)d;
  }
  emxInit_real_T(&y_tmp, 2);
  /* Theta_ = ThetaEval;%Model.ShpFnc(s); */
  inner = Ba->size[1];
  nc = Theta->size[1];
  n = y_tmp->size[0] * y_tmp->size[1];
  y_tmp->size[0] = 6;
  y_tmp->size[1] = Theta->size[1];
  emxEnsureCapacity_real_T(y_tmp, n);
  y_tmp_data = y_tmp->data;
  for (j = 0; j < nc; j++) {
    coffset = j * 6;
    boffset = j * Theta->size[0];
    for (b_i = 0; b_i < 6; b_i++) {
      y_tmp_data[coffset + b_i] = 0.0;
    }
    for (b_k = 0; b_k < inner; b_k++) {
      aoffset = b_k * 6;
      bkj = Theta_data[boffset + b_k];
      for (b_i = 0; b_i < 6; b_i++) {
        n = coffset + b_i;
        y_tmp_data[n] += Ba_data[aoffset + b_i] * bkj;
      }
    }
  }
  inner = y_tmp->size[1];
  for (b_i = 0; b_i < 6; b_i++) {
    XI[b_i] = 0.0;
  }
  for (b_k = 0; b_k < inner; b_k++) {
    aoffset = b_k * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      XI[b_i] += y_tmp_data[aoffset + b_i] * x_data[b_k];
    }
  }
  for (n = 0; n < 6; n++) {
    XI[n] += xia0[n];
  }
  /*  build forward kin - position */
  /* --------------------------------------------------------------------------
   */
  memset(&A[0], 0, 36U * sizeof(double));
  for (n = 0; n < 3; n++) {
    d = Phi_[3 * n];
    A[6 * n] = d;
    inner = 6 * (n + 3);
    A[inner + 3] = d;
    d = Phi_[3 * n + 1];
    A[6 * n + 1] = d;
    A[inner + 4] = d;
    d = Phi_[3 * n + 2];
    A[6 * n + 2] = d;
    A[inner + 5] = d;
  }
  Wh[0] = 0.0;
  Wh[3] = -Z1_data[20];
  Wh[6] = Z1_data[19];
  Wh[1] = Z1_data[20];
  Wh[4] = 0.0;
  Wh[7] = -Z1_data[18];
  Wh[2] = -Z1_data[19];
  Wh[5] = Z1_data[18];
  Wh[8] = 0.0;
  for (n = 0; n < 3; n++) {
    d = Wh[n];
    bkj = Wh[n + 3];
    absxk = Wh[n + 6];
    for (nc = 0; nc < 3; nc++) {
      A[(n + 6 * nc) + 3] = (d * Phi_[3 * nc] + bkj * Phi_[3 * nc + 1]) +
                            absxk * Phi_[3 * nc + 2];
    }
  }
  /* --------------------------------------------------------------------------
   */
  for (n = 0; n < 3; n++) {
    Rt_tmp[3 * n] = Phi_[n];
    Rt_tmp[3 * n + 1] = Phi_[n + 3];
    Rt_tmp[3 * n + 2] = Phi_[n + 6];
  }
  /* --------------------------------------------------------------------------
   */
  memset(&Ai[0], 0, 36U * sizeof(double));
  for (n = 0; n < 3; n++) {
    d = Rt_tmp[3 * n];
    Ai[6 * n] = d;
    inner = 6 * (n + 3);
    Ai[inner + 3] = d;
    d = Rt_tmp[3 * n + 1];
    Ai[6 * n + 1] = d;
    Ai[inner + 4] = d;
    d = Rt_tmp[3 * n + 2];
    Ai[6 * n + 2] = d;
    Ai[inner + 5] = d;
  }
  Wh[0] = 0.0;
  Wh[1] = -Z1_data[20];
  Wh[2] = Z1_data[19];
  Wh[3] = Z1_data[20];
  Wh[4] = 0.0;
  Wh[5] = -Z1_data[18];
  Wh[6] = -Z1_data[19];
  Wh[7] = Z1_data[18];
  Wh[8] = 0.0;
  for (n = 0; n < 3; n++) {
    d = Rt_tmp[n];
    bkj = Rt_tmp[n + 3];
    absxk = Rt_tmp[n + 6];
    for (nc = 0; nc < 3; nc++) {
      Ai[(n + 6 * nc) + 3] =
          (d * Wh[3 * nc] + bkj * Wh[3 * nc + 1]) + absxk * Wh[3 * nc + 2];
    }
  }
  emxInit_real_T(&Jg, 2);
  /*  build jacobian */
  inner = i1 - i;
  i1 = Jg->size[0] * Jg->size[1];
  Jg->size[0] = 6;
  Jg->size[1] = inner + 1;
  emxEnsureCapacity_real_T(Jg, i1);
  Jg_data = Jg->data;
  for (j = 0; j <= inner; j++) {
    nc = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (b_k = 0; b_k < 6; b_k++) {
        i1 = nc + b_k;
        bkj += Ai[b_k * 6 + b_i] * Z1_data[i1 % 6 + 6 * (i + i1 / 6)];
      }
      Jg_data[nc + b_i] = bkj;
    }
  }
  inner = Jg->size[1];
  for (b_i = 0; b_i < 6; b_i++) {
    V[b_i] = 0.0;
  }
  for (b_k = 0; b_k < inner; b_k++) {
    aoffset = b_k * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      V[b_i] += Jg_data[aoffset + b_i] * dx_data[b_k];
    }
  }
  /* --------------------------------------------------------------------------
   */
  memset(&adV[0], 0, 36U * sizeof(double));
  /* --------------------------------------------------------------------------
   */
  Wh[0] = 0.0;
  Wh[3] = -V[2];
  Wh[6] = V[1];
  Wh[1] = V[2];
  Wh[4] = 0.0;
  Wh[7] = -V[0];
  Wh[2] = -V[1];
  Wh[5] = V[0];
  Wh[8] = 0.0;
  /* --------------------------------------------------------------------------
   */
  for (i = 0; i < 3; i++) {
    d = Wh[3 * i];
    adV[6 * i] = d;
    inner = 6 * (i + 3);
    adV[inner + 3] = d;
    d = Wh[3 * i + 1];
    adV[6 * i + 1] = d;
    adV[inner + 4] = d;
    d = Wh[3 * i + 2];
    adV[6 * i + 2] = d;
    adV[inner + 5] = d;
  }
  adV[3] = 0.0;
  adV[9] = -V[5];
  adV[15] = V[4];
  adV[4] = V[5];
  adV[10] = 0.0;
  adV[16] = -V[3];
  adV[5] = -V[4];
  adV[11] = V[3];
  adV[17] = 0.0;
  /*  compute inertia, coriolis, conservative force */
  bkj = 3.3121686421112381E-170;
  d = Z1_data[18] - r0[0];
  rd_idx_0 = d;
  absxk = fabs(d);
  if (absxk > 3.3121686421112381E-170) {
    abs_rd = 1.0;
    bkj = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    abs_rd = t * t;
  }
  d = Z1_data[19] - r0[1];
  rd_idx_1 = d;
  absxk = fabs(d);
  if (absxk > bkj) {
    t = bkj / absxk;
    abs_rd = abs_rd * t * t + 1.0;
    bkj = absxk;
  } else {
    t = absxk / bkj;
    abs_rd += t * t;
  }
  d = Z1_data[20] - r0[2];
  absxk = fabs(d);
  if (absxk > bkj) {
    t = bkj / absxk;
    abs_rd = abs_rd * t * t + 1.0;
    bkj = absxk;
  } else {
    t = absxk / bkj;
    abs_rd += t * t;
  }
  abs_rd = bkj * sqrt(abs_rd);
  a_tmp = abs_rd - rs;
  t = npe * k / rt_powd_snf(a_tmp, npe + 1.0);
  b[0] = 0.0;
  b[1] = 0.0;
  b[2] = 0.0;
  bkj = Gvec[0];
  absxk = Gvec[1];
  d1 = Gvec[2];
  for (i = 0; i < 3; i++) {
    b[i + 3] = (Rt_tmp[i] * bkj + Rt_tmp[i + 3] * absxk) + Rt_tmp[i + 6] * d1;
  }
  for (i = 0; i < 6; i++) {
    bkj = 0.0;
    for (i1 = 0; i1 < 6; i1++) {
      bkj += Mtt[i + 6 * i1] * b[i1];
    }
    y[i] = bkj;
  }
  rd_idx_0 = t * rd_idx_0 / abs_rd;
  rd_idx_1 = t * rd_idx_1 / abs_rd;
  bkj = t * d / abs_rd;
  b[0] = y[0];
  b[1] = y[1];
  b[2] = y[2];
  for (i = 0; i < 3; i++) {
    b[i + 3] = y[i + 3] + ((Rt_tmp[i] * rd_idx_0 + Rt_tmp[i + 3] * rd_idx_1) +
                           Rt_tmp[i + 6] * bkj);
  }
  emxInit_real_T(&dG, 1);
  inner = Jg->size[1];
  i = dG->size[0];
  dG->size[0] = Jg->size[1];
  emxEnsureCapacity_real_T(dG, i);
  dG_data = dG->data;
  for (b_i = 0; b_i < inner; b_i++) {
    aoffset = b_i * 6;
    bkj = 0.0;
    for (b_k = 0; b_k < 6; b_k++) {
      bkj += Jg_data[aoffset + b_k] * b[b_k];
    }
    dG_data[b_i] = bkj;
  }
  /*  compute (nonlinear stiffness) */
  /*  compute potential energy */
  /*  assemble matrices */
  i = dZ1->size[0] * dZ1->size[1];
  dZ1->size[0] = 6;
  dZ1->size[1] = (int)(2.0 * ((double)x->size[0] - 1.0) + 6.0);
  emxEnsureCapacity_real_T(dZ1, i);
  dZ1_data = dZ1->data;
  nc = 6 * (int)(2.0 * ((double)x->size[0] - 1.0) + 6.0);
  for (i = 0; i < nc; i++) {
    dZ1_data[i] = 0.0;
  }
  Wh[0] = 0.0;
  Wh[3] = -XI[2];
  Wh[6] = XI[1];
  Wh[1] = XI[2];
  Wh[4] = 0.0;
  Wh[7] = -XI[0];
  Wh[2] = -XI[1];
  Wh[5] = XI[0];
  Wh[8] = 0.0;
  for (i = 0; i < 3; i++) {
    dZ1_data[i + 18] = 0.0;
    for (i1 = 0; i1 < 3; i1++) {
      dZ1_data[i + 6 * i1] =
          (Phi_[i] * Wh[3 * i1] + Phi_[i + 3] * Wh[3 * i1 + 1]) +
          Phi_[i + 6] * Wh[3 * i1 + 2];
      dZ1_data[i + 18] += Phi_[i + 3 * i1] * XI[i1 + 3];
    }
  }
  if (x->size[0] + 4U < 5U) {
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
  C_data = C->data;
  for (j = 0; j < n; j++) {
    nc = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (b_k = 0; b_k < 6; b_k++) {
        bkj += A[b_k * 6 + b_i] * y_tmp_data[nc + b_k];
      }
      C_data[nc + b_i] = bkj;
    }
  }
  nc = C->size[1];
  for (i1 = 0; i1 < nc; i1++) {
    for (n = 0; n < 6; n++) {
      dZ1_data[n + 6 * (i + i1)] = C_data[n + 6 * i1];
    }
  }
  if (((double)x->size[0] + 6.0) - 1.0 >
      2.0 * ((double)x->size[0] - 1.0) + 6.0) {
    i = -4;
  } else {
    i = x->size[0];
  }
  for (i1 = 0; i1 < 6; i1++) {
    for (n = 0; n < 6; n++) {
      d = 0.0;
      for (nc = 0; nc < 6; nc++) {
        d += A[i1 + 6 * nc] * adV[nc + 6 * n];
      }
      b_A[i1 + 6 * n] = d;
    }
  }
  n = y_tmp->size[1];
  i1 = C->size[0] * C->size[1];
  C->size[0] = 6;
  C->size[1] = y_tmp->size[1];
  emxEnsureCapacity_real_T(C, i1);
  C_data = C->data;
  for (j = 0; j < n; j++) {
    nc = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (b_k = 0; b_k < 6; b_k++) {
        bkj += b_A[b_k * 6 + b_i] * y_tmp_data[nc + b_k];
      }
      C_data[nc + b_i] = bkj;
    }
  }
  nc = C->size[1];
  for (i1 = 0; i1 < nc; i1++) {
    for (n = 0; n < 6; n++) {
      dZ1_data[n + 6 * ((i + i1) + 4)] = C_data[n + 6 * i1];
    }
  }
  dZ1_data[22] =
      ((-Mtt[21] * Z1_data[18] * Gvec[0] + -Mtt[21] * Z1_data[19] * Gvec[1]) +
       -Mtt[21] * Z1_data[20] * Gvec[2]) +
      k / rt_powd_snf(a_tmp, npe);
  d = 0.0;
  for (i = 0; i < 6; i++) {
    bkj = 0.0;
    for (i1 = 0; i1 < 6; i1++) {
      bkj += 0.5 * V[i1] * Mtt[i1 + 6 * i];
    }
    d += bkj * V[i];
  }
  dZ1_data[23] = d;
  i = dZ2->size[0] * dZ2->size[1];
  dZ2->size[0] = x->size[0];
  dZ2->size[1] = (int)(3.0 * (double)x->size[0] + 1.0);
  emxEnsureCapacity_real_T(dZ2, i);
  dZ1_data = dZ2->data;
  nc = x->size[0] * (int)(3.0 * (double)x->size[0] + 1.0);
  for (i = 0; i < nc; i++) {
    dZ1_data[i] = 0.0;
  }
  emxInit_real_T(&b_C, 2);
  emxInit_real_T(&r, 2);
  mtimes(Jg, Mtt, r);
  b_mtimes(r, Jg, b_C);
  y_tmp_data = b_C->data;
  nc = b_C->size[1];
  for (i = 0; i < nc; i++) {
    inner = b_C->size[0];
    for (i1 = 0; i1 < inner; i1++) {
      dZ1_data[i1 + dZ2->size[0] * i] = y_tmp_data[i1 + b_C->size[0] * i];
    }
  }
  if (x->size[0] + 1U > ((unsigned int)x->size[0] << 1)) {
    i = 0;
  } else {
    i = x->size[0];
  }
  for (i1 = 0; i1 < 6; i1++) {
    for (n = 0; n < 6; n++) {
      d = 0.0;
      bkj = 0.0;
      for (nc = 0; nc < 6; nc++) {
        inner = nc + 6 * n;
        d += Mtt[i1 + 6 * nc] * adV[inner];
        bkj += adV[nc + 6 * i1] * Mtt[inner];
      }
      inner = i1 + 6 * n;
      b_A[inner] = bkj;
      A[inner] = d;
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
  C_data = C->data;
  for (j = 0; j < n; j++) {
    nc = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (b_k = 0; b_k < 6; b_k++) {
        bkj += A[b_k * 6 + b_i] * Jg_data[nc + b_k];
      }
      C_data[nc + b_i] = bkj;
    }
  }
  emxInit_real_T(&c_C, 2);
  inner = i3 - i2;
  n = inner - 1;
  i1 = c_C->size[0] * c_C->size[1];
  c_C->size[0] = 6;
  c_C->size[1] = inner;
  emxEnsureCapacity_real_T(c_C, i1);
  C_data = c_C->data;
  for (j = 0; j <= n; j++) {
    nc = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (b_k = 0; b_k < 6; b_k++) {
        i1 = nc + b_k;
        bkj += Ai[b_k * 6 + b_i] * Z1_data[i1 % 6 + 6 * (i2 + i1 / 6)];
      }
      C_data[nc + b_i] = bkj;
    }
  }
  emxInit_real_T(&d_C, 2);
  n = c_C->size[1];
  i1 = d_C->size[0] * d_C->size[1];
  d_C->size[0] = 6;
  d_C->size[1] = c_C->size[1];
  emxEnsureCapacity_real_T(d_C, i1);
  y_tmp_data = d_C->data;
  for (j = 0; j < n; j++) {
    nc = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (b_k = 0; b_k < 6; b_k++) {
        bkj += Mtt[b_k * 6 + b_i] * C_data[nc + b_k];
      }
      y_tmp_data[nc + b_i] = bkj;
    }
  }
  emxFree_real_T(&c_C);
  if (C->size[1] == d_C->size[1]) {
    nc = 6 * C->size[1];
    i1 = C->size[0] * C->size[1];
    C->size[0] = 6;
    emxEnsureCapacity_real_T(C, i1);
    C_data = C->data;
    for (i1 = 0; i1 < nc; i1++) {
      C_data[i1] += y_tmp_data[i1];
    }
  } else {
    plus(C, d_C);
    C_data = C->data;
  }
  emxFree_real_T(&d_C);
  inner = Jg->size[1];
  n = C->size[1];
  i1 = b_C->size[0] * b_C->size[1];
  b_C->size[0] = Jg->size[1];
  b_C->size[1] = C->size[1];
  emxEnsureCapacity_real_T(b_C, i1);
  y_tmp_data = b_C->data;
  for (j = 0; j < n; j++) {
    coffset = j * inner;
    boffset = j * 6;
    for (b_i = 0; b_i < inner; b_i++) {
      aoffset = b_i * 6;
      bkj = 0.0;
      for (b_k = 0; b_k < 6; b_k++) {
        bkj += Jg_data[aoffset + b_k] * C_data[boffset + b_k];
      }
      y_tmp_data[coffset + b_i] = bkj;
    }
  }
  emxFree_real_T(&C);
  emxFree_real_T(&Jg);
  nc = b_C->size[1];
  for (i1 = 0; i1 < nc; i1++) {
    inner = b_C->size[0];
    for (i2 = 0; i2 < inner; i2++) {
      dZ1_data[i2 + dZ2->size[0] * (i + i1)] =
          y_tmp_data[i2 + b_C->size[0] * i1];
    }
  }
  u = ((unsigned int)x->size[0] << 1) + 1U;
  if (u > 3.0 * (double)x->size[0]) {
    i = 1;
  } else {
    i = (int)u;
  }
  mtimes(y_tmp, Ktt, r);
  b_mtimes(r, y_tmp, b_C);
  y_tmp_data = b_C->data;
  nc = b_C->size[1];
  emxFree_real_T(&r);
  emxFree_real_T(&y_tmp);
  for (i1 = 0; i1 < nc; i1++) {
    inner = b_C->size[0];
    for (i2 = 0; i2 < inner; i2++) {
      dZ1_data[i2 + dZ2->size[0] * ((i + i1) - 1)] =
          y_tmp_data[i2 + b_C->size[0] * i1];
    }
  }
  emxFree_real_T(&b_C);
  i = (int)(3.0 * (double)x->size[0] + 1.0) - 1;
  nc = dG->size[0];
  for (i1 = 0; i1 < nc; i1++) {
    dZ1_data[i1 + dZ2->size[0] * i] = dG_data[i1];
  }
  emxFree_real_T(&dG);
}

static void b_binary_expand_op(emxArray_real_T *in1, double in2,
                               const emxArray_real_T *in3,
                               const emxArray_real_T *in4)
{
  emxArray_real_T *b_in1;
  const double *in3_data;
  const double *in4_data;
  double *b_in1_data;
  double *in1_data;
  int aux_0_1;
  int aux_1_1;
  int i;
  int i1;
  int i2;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  in4_data = in4->data;
  in3_data = in3->data;
  in1_data = in1->data;
  emxInit_real_T(&b_in1, 2);
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = 6;
  if (in3->size[1] == 1) {
    b_in1->size[1] = in1->size[1];
  } else {
    b_in1->size[1] = in3->size[1];
  }
  emxEnsureCapacity_real_T(b_in1, i);
  b_in1_data = b_in1->data;
  stride_0_1 = (in1->size[1] != 1);
  stride_1_1 = (in3->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (in3->size[1] == 1) {
    loop_ub = in1->size[1];
  } else {
    loop_ub = in3->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < 6; i1++) {
      i2 = i1 + 6 * aux_1_1;
      b_in1_data[i1 + 6 * i] = in1_data[i1 + 6 * aux_0_1] +
                               in2 * (in3_data[i2] + 3.0 * in4_data[i2]);
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = 6;
  in1->size[1] = b_in1->size[1];
  emxEnsureCapacity_real_T(in1, i);
  in1_data = in1->data;
  loop_ub = b_in1->size[1];
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < 6; i1++) {
      i2 = i1 + 6 * i;
      in1_data[i2] = b_in1_data[i2];
    }
  }
  emxFree_real_T(&b_in1);
}

static void binary_expand_op(emxArray_real_T *in1, double in2,
                             const emxArray_real_T *in3,
                             const emxArray_real_T *in4)
{
  emxArray_real_T *b_in1;
  const double *in3_data;
  const double *in4_data;
  double *b_in1_data;
  double *in1_data;
  int aux_0_1;
  int aux_1_1;
  int aux_2_1;
  int b_in4;
  int i;
  int i1;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  int stride_2_0;
  int stride_2_1;
  in4_data = in4->data;
  in3_data = in3->data;
  in1_data = in1->data;
  emxInit_real_T(&b_in1, 2);
  i = b_in1->size[0] * b_in1->size[1];
  if (in4->size[0] == 1) {
    b_in4 = in3->size[0];
  } else {
    b_in4 = in4->size[0];
  }
  if (b_in4 == 1) {
    b_in1->size[0] = in1->size[0];
  } else if (in4->size[0] == 1) {
    b_in1->size[0] = in3->size[0];
  } else {
    b_in1->size[0] = in4->size[0];
  }
  if (in4->size[1] == 1) {
    b_in4 = in3->size[1];
  } else {
    b_in4 = in4->size[1];
  }
  if (b_in4 == 1) {
    b_in1->size[1] = in1->size[1];
  } else if (in4->size[1] == 1) {
    b_in1->size[1] = in3->size[1];
  } else {
    b_in1->size[1] = in4->size[1];
  }
  emxEnsureCapacity_real_T(b_in1, i);
  b_in1_data = b_in1->data;
  stride_0_0 = (in1->size[0] != 1);
  stride_0_1 = (in1->size[1] != 1);
  stride_1_0 = (in3->size[0] != 1);
  stride_1_1 = (in3->size[1] != 1);
  stride_2_0 = (in4->size[0] != 1);
  stride_2_1 = (in4->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  aux_2_1 = 0;
  if (in4->size[1] == 1) {
    b_in4 = in3->size[1];
  } else {
    b_in4 = in4->size[1];
  }
  if (b_in4 == 1) {
    loop_ub = in1->size[1];
  } else if (in4->size[1] == 1) {
    loop_ub = in3->size[1];
  } else {
    loop_ub = in4->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    if (in4->size[0] == 1) {
      b_in4 = in3->size[0];
    } else {
      b_in4 = in4->size[0];
    }
    if (b_in4 == 1) {
      b_in4 = in1->size[0];
    } else if (in4->size[0] == 1) {
      b_in4 = in3->size[0];
    } else {
      b_in4 = in4->size[0];
    }
    for (i1 = 0; i1 < b_in4; i1++) {
      b_in1_data[i1 + b_in1->size[0] * i] =
          in1_data[i1 * stride_0_0 + in1->size[0] * aux_0_1] +
          in2 * (in3_data[i1 * stride_1_0 + in3->size[0] * aux_1_1] +
                 3.0 * in4_data[i1 * stride_2_0 + in4->size[0] * aux_2_1]);
    }
    aux_2_1 += stride_2_1;
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = b_in1->size[0];
  in1->size[1] = b_in1->size[1];
  emxEnsureCapacity_real_T(in1, i);
  in1_data = in1->data;
  loop_ub = b_in1->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_in4 = b_in1->size[0];
    for (i1 = 0; i1 < b_in4; i1++) {
      in1_data[i1 + in1->size[0] * i] = b_in1_data[i1 + b_in1->size[0] * i];
    }
  }
  emxFree_real_T(&b_in1);
}

static void c_binary_expand_op(
    const emxArray_real_T *in1, const emxArray_real_T *in2,
    const emxArray_real_T *in3, double in4, const emxArray_real_T *in5,
    const emxArray_real_T *in6, int in7, const emxArray_real_T *in8,
    const emxArray_real_T *in9, const double in10[36], const double in11[36],
    const double in12[3], double in13, double in14, const double in15[3],
    double in16, emxArray_real_T *in17, emxArray_real_T *in18)
{
  emxArray_real_T *b_in3;
  emxArray_real_T *b_in6;
  const double *in3_data;
  const double *in5_data;
  const double *in6_data;
  const double *in8_data;
  double *b_in3_data;
  int aux_0_1;
  int i;
  int i1;
  int i2;
  int loop_ub;
  int stride_0_1;
  in8_data = in8->data;
  in6_data = in6->data;
  in5_data = in5->data;
  in3_data = in3->data;
  emxInit_real_T(&b_in3, 2);
  i = (int)((unsigned int)(in7 + 1) << 1);
  i1 = b_in3->size[0] * b_in3->size[1];
  b_in3->size[0] = 6;
  b_in3->size[1] = in5->size[1];
  emxEnsureCapacity_real_T(b_in3, i1);
  b_in3_data = b_in3->data;
  stride_0_1 = (in3->size[1] != 1);
  aux_0_1 = 0;
  loop_ub = in5->size[1];
  for (i1 = 0; i1 < loop_ub; i1++) {
    for (i2 = 0; i2 < 6; i2++) {
      int i3;
      i3 = i2 + 6 * i1;
      b_in3_data[i3] = in3_data[i2 + 6 * aux_0_1] + in4 * in5_data[i3];
    }
    aux_0_1 += stride_0_1;
  }
  emxInit_real_T(&b_in6, 2);
  loop_ub = in6->size[0];
  stride_0_1 = in6->size[1];
  i1 = b_in6->size[0] * b_in6->size[1];
  b_in6->size[0] = loop_ub;
  b_in6->size[1] = stride_0_1;
  emxEnsureCapacity_real_T(b_in6, i1);
  b_in3_data = b_in6->data;
  for (i1 = 0; i1 < stride_0_1; i1++) {
    for (i2 = 0; i2 < loop_ub; i2++) {
      b_in3_data[i2 + b_in6->size[0] * i1] =
          in6_data[(i2 + in6->size[0] * i1) +
                   in6->size[0] * in6->size[1] * (i - 1)];
    }
  }
  LagrangianODEX(in1, in2, b_in3, b_in6,
                 *(double(*)[6]) &
                     in8_data[6 * ((int)((unsigned int)(in7 + 1) << 1) - 1)],
                 in9, in10, in11, in12, in13, in14, in15, in16, in17, in18);
  emxFree_real_T(&b_in6);
  emxFree_real_T(&b_in3);
}

static void plus(emxArray_real_T *in1, const emxArray_real_T *in2)
{
  emxArray_real_T *b_in1;
  const double *in2_data;
  double *b_in1_data;
  double *in1_data;
  int aux_0_1;
  int aux_1_1;
  int i;
  int i1;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  in2_data = in2->data;
  in1_data = in1->data;
  emxInit_real_T(&b_in1, 2);
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = 6;
  if (in2->size[1] == 1) {
    b_in1->size[1] = in1->size[1];
  } else {
    b_in1->size[1] = in2->size[1];
  }
  emxEnsureCapacity_real_T(b_in1, i);
  b_in1_data = b_in1->data;
  stride_0_1 = (in1->size[1] != 1);
  stride_1_1 = (in2->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (in2->size[1] == 1) {
    loop_ub = in1->size[1];
  } else {
    loop_ub = in2->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < 6; i1++) {
      b_in1_data[i1 + 6 * i] =
          in1_data[i1 + 6 * aux_0_1] + in2_data[i1 + 6 * aux_1_1];
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = 6;
  in1->size[1] = b_in1->size[1];
  emxEnsureCapacity_real_T(in1, i);
  in1_data = in1->data;
  loop_ub = b_in1->size[1];
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < 6; i1++) {
      stride_0_1 = i1 + 6 * i;
      in1_data[stride_0_1] = b_in1_data[stride_0_1];
    }
  }
  emxFree_real_T(&b_in1);
}

static double rt_powd_snf(double u0, double u1)
{
  double y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    double d;
    double d1;
    d = fabs(u0);
    d1 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d == 1.0) {
        y = 1.0;
      } else if (d > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }
  return y;
}

void computeLagrangianFast(
    const emxArray_real_T *x, const emxArray_real_T *dx, double ds,
    const double p0[3], const double Phi0[9], const emxArray_real_T *xia0,
    const emxArray_real_T *Th, const emxArray_real_T *Ba, const double Ktt[36],
    const double Mtt[36], double Zeta, const double Gvec[3], double k,
    double npe, const double r0[3], double rs, emxArray_real_T *M,
    emxArray_real_T *C, emxArray_real_T *K, emxArray_real_T *R,
    emxArray_real_T *G, emxArray_real_T *p, emxArray_real_T *Phi,
    emxArray_real_T *J, emxArray_real_T *Jt, double *Vg, double *Kin)
{
  emxArray_real_T *K1Z1;
  emxArray_real_T *K1Z2;
  emxArray_real_T *K2Z1;
  emxArray_real_T *K2Z2;
  emxArray_real_T *Z1;
  emxArray_real_T *Z2;
  emxArray_real_T *b_Th;
  emxArray_real_T *b_Z1;
  double Ai[36];
  double dv[9];
  const double *Th_data;
  const double *xia0_data;
  double a;
  double a_tmp;
  double s;
  double *J_data;
  double *Jt_data;
  double *K1Z1_data;
  double *K1Z2_data;
  double *K2Z1_data;
  double *M_data;
  double *Phi_data;
  double *Z1_data;
  double *Z2_data;
  double *p_data;
  int Rt_tmp;
  int b_K1Z2;
  int b_k;
  int b_loop_ub;
  int b_n;
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  int i5;
  int i6;
  int ii;
  int loop_ub;
  int n;
  unsigned int u;
  Th_data = Th->data;
  xia0_data = xia0->data;
  emxInit_real_T(&Z1, 2);
  /*  % states */
  /*       % spatial steps */
  /*       % position zero */
  /*     % phi zero */
  /*     % intrinsic strain vector */
  /*       % evaluated Theta matrix */
  /*       % state to strain matrix */
  /*      % geometric stiffness */
  /*      % geometric inertia */
  /*     % dampings coefficient */
  /*     % gravitional vector */
  /*        % potential constant */
  /*      % potential energy power (k/((norm(r-r0)-rs)^n)), n must be odd */
  /*  potential constants */
  /*  compute total length */
  i = Z1->size[0] * Z1->size[1];
  Z1->size[0] = 6;
  Z1->size[1] = (int)(2.0 * ((double)x->size[0] - 1.0) + 6.0);
  emxEnsureCapacity_real_T(Z1, i);
  Z1_data = Z1->data;
  loop_ub = 6 * (int)(2.0 * ((double)x->size[0] - 1.0) + 6.0);
  for (i = 0; i < loop_ub; i++) {
    Z1_data[i] = 0.0;
  }
  emxInit_real_T(&Z2, 2);
  i = Z2->size[0] * Z2->size[1];
  Z2->size[0] = x->size[0];
  Z2->size[1] = (int)(3.0 * (double)x->size[0] + 1.0);
  emxEnsureCapacity_real_T(Z2, i);
  Z2_data = Z2->data;
  loop_ub = x->size[0] * (int)(3.0 * (double)x->size[0] + 1.0);
  for (i = 0; i < loop_ub; i++) {
    Z2_data[i] = 0.0;
  }
  for (i = 0; i < 3; i++) {
    Z1_data[6 * i] = Phi0[3 * i];
    Z1_data[6 * i + 1] = Phi0[3 * i + 1];
    Z1_data[6 * i + 2] = Phi0[3 * i + 2];
    Z1_data[i + 18] = p0[i];
  }
  /* NLStiff = false;  */
  i = p->size[0] * p->size[1];
  p->size[0] = 3;
  p->size[1] = (int)((double)Th->size[2] / 2.0);
  emxEnsureCapacity_real_T(p, i);
  p_data = p->data;
  loop_ub = 3 * (int)((double)Th->size[2] / 2.0);
  for (i = 0; i < loop_ub; i++) {
    p_data[i] = 0.0;
  }
  i = Phi->size[0] * Phi->size[1] * Phi->size[2];
  Phi->size[0] = 3;
  Phi->size[1] = 3;
  Phi->size[2] = (int)((double)Th->size[2] / 2.0);
  emxEnsureCapacity_real_T(Phi, i);
  Phi_data = Phi->data;
  loop_ub = 9 * (int)((double)Th->size[2] / 2.0);
  for (i = 0; i < loop_ub; i++) {
    Phi_data[i] = 0.0;
  }
  i = J->size[0] * J->size[1] * J->size[2];
  J->size[0] = 6;
  J->size[1] = x->size[0];
  J->size[2] = (int)((double)Th->size[2] / 2.0);
  emxEnsureCapacity_real_T(J, i);
  J_data = J->data;
  loop_ub = 6 * x->size[0] * (int)((double)Th->size[2] / 2.0);
  for (i = 0; i < loop_ub; i++) {
    J_data[i] = 0.0;
  }
  i = Jt->size[0] * Jt->size[1] * Jt->size[2];
  Jt->size[0] = 6;
  Jt->size[1] = x->size[0];
  Jt->size[2] = (int)((double)Th->size[2] / 2.0);
  emxEnsureCapacity_real_T(Jt, i);
  Jt_data = Jt->data;
  loop_ub = 6 * x->size[0] * (int)((double)Th->size[2] / 2.0);
  for (i = 0; i < loop_ub; i++) {
    Jt_data[i] = 0.0;
  }
  i = (int)((double)Th->size[2] / 2.0);
  if (i - 1 >= 0) {
    a = 0.66666666666666663 * ds;
    a_tmp = 0.25 * ds;
    dv[0] = 0.0;
    dv[4] = 0.0;
    dv[8] = 0.0;
    if (x->size[0] + 4U < 5U) {
      i1 = 0;
      i2 = -1;
    } else {
      i1 = 4;
      i2 = x->size[0] + 3;
    }
    n = i2 - i1;
    s = 2.0 * ((double)x->size[0] - 1.0) + 6.0;
    if (((double)x->size[0] + 6.0) - 1.0 > s) {
      i3 = 0;
      i4 = 0;
    } else {
      i3 = x->size[0] + 4;
      i4 = (int)s;
    }
    b_n = (i4 - i3) - 1;
  }
  emxInit_real_T(&K1Z1, 2);
  emxInit_real_T(&K1Z2, 2);
  emxInit_real_T(&K2Z1, 2);
  emxInit_real_T(&K2Z2, 2);
  emxInit_real_T(&b_Th, 2);
  emxInit_real_T(&b_Z1, 2);
  for (ii = 0; ii < i; ii++) {
    double Rt[9];
    double b_xia0_data[6];
    double d;
    /*  first EL-diff eval */
    b_K1Z2 = (int)((unsigned int)(ii + 1) << 1);
    loop_ub = Th->size[0];
    b_loop_ub = Th->size[1];
    Rt_tmp = b_Th->size[0] * b_Th->size[1];
    b_Th->size[0] = Th->size[0];
    b_Th->size[1] = Th->size[1];
    emxEnsureCapacity_real_T(b_Th, Rt_tmp);
    M_data = b_Th->data;
    for (Rt_tmp = 0; Rt_tmp < b_loop_ub; Rt_tmp++) {
      for (b_k = 0; b_k < loop_ub; b_k++) {
        M_data[b_k + b_Th->size[0] * Rt_tmp] =
            Th_data[(b_k + Th->size[0] * Rt_tmp) +
                    Th->size[0] * Th->size[1] * (b_K1Z2 - 2)];
      }
    }
    for (i5 = 0; i5 < 6; i5++) {
      b_xia0_data[i5] =
          xia0_data[i5 + 6 * ((int)((unsigned int)(ii + 1) << 1) - 2)];
    }
    LagrangianODEX(x, dx, Z1, b_Th, b_xia0_data, Ba, Mtt, Ktt, Gvec, k, npe, r0,
                   rs, K1Z1, K1Z2);
    K1Z2_data = K1Z2->data;
    K1Z1_data = K1Z1->data;
    /*  second EL-diff eval */
    if (Z1->size[1] == K1Z1->size[1]) {
      Rt_tmp = b_Z1->size[0] * b_Z1->size[1];
      b_Z1->size[0] = 6;
      b_Z1->size[1] = Z1->size[1];
      emxEnsureCapacity_real_T(b_Z1, Rt_tmp);
      M_data = b_Z1->data;
      loop_ub = 6 * Z1->size[1];
      for (Rt_tmp = 0; Rt_tmp < loop_ub; Rt_tmp++) {
        M_data[Rt_tmp] = Z1_data[Rt_tmp] + a * K1Z1_data[Rt_tmp];
      }
      loop_ub = Th->size[0];
      b_loop_ub = Th->size[1];
      Rt_tmp = b_Th->size[0] * b_Th->size[1];
      b_Th->size[0] = Th->size[0];
      b_Th->size[1] = Th->size[1];
      emxEnsureCapacity_real_T(b_Th, Rt_tmp);
      M_data = b_Th->data;
      for (Rt_tmp = 0; Rt_tmp < b_loop_ub; Rt_tmp++) {
        for (b_k = 0; b_k < loop_ub; b_k++) {
          M_data[b_k + b_Th->size[0] * Rt_tmp] =
              Th_data[(b_k + Th->size[0] * Rt_tmp) +
                      Th->size[0] * Th->size[1] * (b_K1Z2 - 1)];
        }
      }
      for (i6 = 0; i6 < 6; i6++) {
        b_xia0_data[i6] =
            xia0_data[i6 + 6 * ((int)((unsigned int)(ii + 1) << 1) - 1)];
      }
      LagrangianODEX(x, dx, b_Z1, b_Th, b_xia0_data, Ba, Mtt, Ktt, Gvec, k, npe,
                     r0, rs, K2Z1, K2Z2);
      M_data = K2Z2->data;
      K2Z1_data = K2Z1->data;
    } else {
      c_binary_expand_op(x, dx, Z1, a, K1Z1, Th, ii, xia0, Ba, Mtt, Ktt, Gvec,
                         k, npe, r0, rs, K2Z1, K2Z2);
      M_data = K2Z2->data;
      K2Z1_data = K2Z1->data;
    }
    /*  update integrands */
    if (Z1->size[1] == K1Z1->size[1]) {
      loop_ub = 6 * Z1->size[1];
      b_K1Z2 = Z1->size[0] * Z1->size[1];
      Z1->size[0] = 6;
      emxEnsureCapacity_real_T(Z1, b_K1Z2);
      Z1_data = Z1->data;
      for (b_K1Z2 = 0; b_K1Z2 < loop_ub; b_K1Z2++) {
        Z1_data[b_K1Z2] +=
            a_tmp * (K1Z1_data[b_K1Z2] + 3.0 * K2Z1_data[b_K1Z2]);
      }
    } else {
      b_binary_expand_op(Z1, a_tmp, K1Z1, K2Z1);
      Z1_data = Z1->data;
    }
    if (K1Z2->size[0] == 1) {
      loop_ub = K2Z2->size[0];
    } else {
      loop_ub = K1Z2->size[0];
    }
    if (K1Z2->size[1] == 1) {
      b_K1Z2 = K2Z2->size[1];
    } else {
      b_K1Z2 = K1Z2->size[1];
    }
    if ((K1Z2->size[0] == K2Z2->size[0]) && (K1Z2->size[1] == K2Z2->size[1]) &&
        (Z2->size[0] == loop_ub) && (Z2->size[1] == b_K1Z2)) {
      loop_ub = Z2->size[0] * Z2->size[1];
      for (b_K1Z2 = 0; b_K1Z2 < loop_ub; b_K1Z2++) {
        Z2_data[b_K1Z2] += a_tmp * (K1Z2_data[b_K1Z2] + 3.0 * M_data[b_K1Z2]);
      }
    } else {
      binary_expand_op(Z2, a_tmp, K1Z2, K2Z2);
      Z2_data = Z2->data;
    }
    /*  compute kinematics */
    for (b_K1Z2 = 0; b_K1Z2 < 3; b_K1Z2++) {
      p_data[b_K1Z2 + 3 * ii] = Z1_data[b_K1Z2 + 18];
      Rt_tmp = 3 * b_K1Z2 + 9 * ii;
      Phi_data[Rt_tmp] = Z1_data[6 * b_K1Z2];
      Phi_data[Rt_tmp + 1] = Z1_data[6 * b_K1Z2 + 1];
      Phi_data[Rt_tmp + 2] = Z1_data[6 * b_K1Z2 + 2];
    }
    /* --------------------------------------------------------------------------
     */
    for (b_K1Z2 = 0; b_K1Z2 < 3; b_K1Z2++) {
      Rt_tmp = b_K1Z2 + 9 * ii;
      Rt[3 * b_K1Z2] = Phi_data[Rt_tmp];
      Rt[3 * b_K1Z2 + 1] = Phi_data[Rt_tmp + 3];
      Rt[3 * b_K1Z2 + 2] = Phi_data[Rt_tmp + 6];
    }
    /* --------------------------------------------------------------------------
     */
    memset(&Ai[0], 0, 36U * sizeof(double));
    for (b_K1Z2 = 0; b_K1Z2 < 3; b_K1Z2++) {
      s = Rt[3 * b_K1Z2];
      Ai[6 * b_K1Z2] = s;
      Rt_tmp = 6 * (b_K1Z2 + 3);
      Ai[Rt_tmp + 3] = s;
      s = Rt[3 * b_K1Z2 + 1];
      Ai[6 * b_K1Z2 + 1] = s;
      Ai[Rt_tmp + 4] = s;
      s = Rt[3 * b_K1Z2 + 2];
      Ai[6 * b_K1Z2 + 2] = s;
      Ai[Rt_tmp + 5] = s;
    }
    s = p_data[3 * ii + 2];
    dv[1] = -s;
    d = p_data[3 * ii + 1];
    dv[2] = d;
    dv[3] = s;
    s = p_data[3 * ii];
    dv[5] = -s;
    dv[6] = -d;
    dv[7] = s;
    for (b_K1Z2 = 0; b_K1Z2 < 3; b_K1Z2++) {
      double d1;
      s = Rt[b_K1Z2];
      d = Rt[b_K1Z2 + 3];
      d1 = Rt[b_K1Z2 + 6];
      for (Rt_tmp = 0; Rt_tmp < 3; Rt_tmp++) {
        Ai[(b_K1Z2 + 6 * Rt_tmp) + 3] =
            (s * dv[3 * Rt_tmp] + d * dv[3 * Rt_tmp + 1]) +
            d1 * dv[3 * Rt_tmp + 2];
      }
    }
    b_K1Z2 = K1Z1->size[0] * K1Z1->size[1];
    K1Z1->size[0] = 6;
    K1Z1->size[1] = (i2 - i1) + 1;
    emxEnsureCapacity_real_T(K1Z1, b_K1Z2);
    K1Z1_data = K1Z1->data;
    for (loop_ub = 0; loop_ub <= n; loop_ub++) {
      Rt_tmp = loop_ub * 6;
      for (b_loop_ub = 0; b_loop_ub < 6; b_loop_ub++) {
        s = 0.0;
        for (b_k = 0; b_k < 6; b_k++) {
          b_K1Z2 = Rt_tmp + b_k;
          s += Ai[b_k * 6 + b_loop_ub] *
               Z1_data[b_K1Z2 % 6 + 6 * (i1 + b_K1Z2 / 6)];
        }
        K1Z1_data[Rt_tmp + b_loop_ub] = s;
      }
    }
    loop_ub = K1Z1->size[1];
    for (b_K1Z2 = 0; b_K1Z2 < loop_ub; b_K1Z2++) {
      for (Rt_tmp = 0; Rt_tmp < 6; Rt_tmp++) {
        b_k = Rt_tmp + 6 * b_K1Z2;
        J_data[b_k + 6 * J->size[1] * ii] = K1Z1_data[b_k];
      }
    }
    b_K1Z2 = K1Z1->size[0] * K1Z1->size[1];
    K1Z1->size[0] = 6;
    K1Z1->size[1] = i4 - i3;
    emxEnsureCapacity_real_T(K1Z1, b_K1Z2);
    K1Z1_data = K1Z1->data;
    for (loop_ub = 0; loop_ub <= b_n; loop_ub++) {
      Rt_tmp = loop_ub * 6;
      for (b_loop_ub = 0; b_loop_ub < 6; b_loop_ub++) {
        s = 0.0;
        for (b_k = 0; b_k < 6; b_k++) {
          b_K1Z2 = Rt_tmp + b_k;
          s += Ai[b_k * 6 + b_loop_ub] *
               Z1_data[b_K1Z2 % 6 + 6 * (i3 + b_K1Z2 / 6)];
        }
        K1Z1_data[Rt_tmp + b_loop_ub] = s;
      }
    }
    loop_ub = K1Z1->size[1];
    for (b_K1Z2 = 0; b_K1Z2 < loop_ub; b_K1Z2++) {
      for (Rt_tmp = 0; Rt_tmp < 6; Rt_tmp++) {
        b_k = Rt_tmp + 6 * b_K1Z2;
        Jt_data[b_k + 6 * Jt->size[1] * ii] = K1Z1_data[b_k];
      }
    }
  }
  emxFree_real_T(&b_Z1);
  emxFree_real_T(&b_Th);
  emxFree_real_T(&K2Z2);
  emxFree_real_T(&K2Z1);
  emxFree_real_T(&K1Z2);
  emxFree_real_T(&K1Z1);
  /*  recover the dynamics entities */
  if (x->size[0] < 1) {
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
  M_data = M->data;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      M_data[i1 + M->size[0] * i] = Z2_data[i1 + Z2->size[0] * i];
    }
  }
  if (x->size[0] < 1) {
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
  M_data = C->data;
  for (i1 = 0; i1 < b_loop_ub; i1++) {
    for (i2 = 0; i2 < loop_ub; i2++) {
      M_data[i2 + C->size[0] * i1] = Z2_data[i2 + Z2->size[0] * (i + i1)];
    }
  }
  if (x->size[0] < 1) {
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
  M_data = K->data;
  for (i1 = 0; i1 < b_loop_ub; i1++) {
    for (i2 = 0; i2 < loop_ub; i2++) {
      M_data[i2 + K->size[0] * i1] = Z2_data[i2 + Z2->size[0] * (i + i1)];
    }
  }
  if (x->size[0] < 1) {
    Rt_tmp = 0;
  } else {
    Rt_tmp = x->size[0];
  }
  i1 = (int)(3.0 * (double)x->size[0] + 1.0);
  i2 = G->size[0];
  G->size[0] = Rt_tmp;
  emxEnsureCapacity_real_T(G, i2);
  M_data = G->data;
  for (i2 = 0; i2 < Rt_tmp; i2++) {
    M_data[i2] = Z2_data[i2 + Z2->size[0] * (i1 - 1)];
  }
  *Vg = Z1_data[22];
  *Kin = Z1_data[23];
  i1 = R->size[0] * R->size[1];
  R->size[0] = loop_ub;
  R->size[1] = b_loop_ub;
  emxEnsureCapacity_real_T(R, i1);
  M_data = R->data;
  emxFree_real_T(&Z1);
  for (i1 = 0; i1 < b_loop_ub; i1++) {
    for (i2 = 0; i2 < loop_ub; i2++) {
      M_data[i2 + R->size[0] * i1] =
          Zeta * Z2_data[i2 + Z2->size[0] * (i + i1)];
    }
  }
  emxFree_real_T(&Z2);
}

/* End of code generation (computeLagrangianFast.c) */
