/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * computeForwardKinematicsFast.c
 *
 * Code generation for function 'computeForwardKinematicsFast'
 *
 */

/* Include files */
#include "computeForwardKinematicsFast.h"
#include "computeForwardKinematicsFast_emxutil.h"
#include "computeForwardKinematicsFast_types.h"
#include <string.h>

/* Function Declarations */
static void ForwardODEX(const emxArray_real_T *x, const emxArray_real_T *Z1,
                        const emxArray_real_T *Theta, const double xia0[6],
                        const emxArray_real_T *Ba, emxArray_real_T *dZ1);

static void b_binary_expand_op(emxArray_real_T *in1, const emxArray_real_T *in2,
                               const emxArray_real_T *in3, double in4,
                               const emxArray_real_T *in5,
                               const emxArray_real_T *in6, int in7,
                               const emxArray_real_T *in8,
                               const emxArray_real_T *in9);

static void binary_expand_op(emxArray_real_T *in1, double in2,
                             const emxArray_real_T *in3,
                             const emxArray_real_T *in4);

/* Function Definitions */
static void ForwardODEX(const emxArray_real_T *x, const emxArray_real_T *Z1,
                        const emxArray_real_T *Theta, const double xia0[6],
                        const emxArray_real_T *Ba, emxArray_real_T *dZ1)
{
  emxArray_real_T *BTh;
  emxArray_real_T *C;
  double A[36];
  double Phi_[9];
  double dv[9];
  double XI[6];
  const double *Ba_data;
  const double *Theta_data;
  const double *Z1_data;
  const double *x_data;
  double bkj;
  double *BTh_data;
  double *C_data;
  double *dZ1_data;
  int aoffset;
  int b_i;
  int coffset;
  int i;
  int inner;
  int j;
  int k;
  int nc;
  Ba_data = Ba->data;
  Theta_data = Theta->data;
  Z1_data = Z1->data;
  x_data = x->data;
  /* --------------------------------------------------------------------------
   */
  for (i = 0; i < 3; i++) {
    Phi_[3 * i] = Z1_data[6 * i];
    Phi_[3 * i + 1] = Z1_data[6 * i + 1];
    Phi_[3 * i + 2] = Z1_data[6 * i + 2];
  }
  emxInit_real_T(&BTh, 2);
  inner = Ba->size[1];
  nc = Theta->size[1];
  i = BTh->size[0] * BTh->size[1];
  BTh->size[0] = 6;
  BTh->size[1] = Theta->size[1];
  emxEnsureCapacity_real_T(BTh, i);
  BTh_data = BTh->data;
  for (j = 0; j < nc; j++) {
    int boffset;
    coffset = j * 6;
    boffset = j * Theta->size[0];
    for (b_i = 0; b_i < 6; b_i++) {
      BTh_data[coffset + b_i] = 0.0;
    }
    for (k = 0; k < inner; k++) {
      aoffset = k * 6;
      bkj = Theta_data[boffset + k];
      for (b_i = 0; b_i < 6; b_i++) {
        i = coffset + b_i;
        BTh_data[i] += Ba_data[aoffset + b_i] * bkj;
      }
    }
  }
  inner = BTh->size[1];
  for (b_i = 0; b_i < 6; b_i++) {
    XI[b_i] = 0.0;
  }
  for (k = 0; k < inner; k++) {
    aoffset = k * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      XI[b_i] += BTh_data[aoffset + b_i] * x_data[k];
    }
  }
  for (i = 0; i < 6; i++) {
    XI[i] += xia0[i];
  }
  /*  build forward kin - position */
  /* --------------------------------------------------------------------------
   */
  memset(&A[0], 0, 36U * sizeof(double));
  for (i = 0; i < 3; i++) {
    bkj = Phi_[3 * i];
    A[6 * i] = bkj;
    inner = 6 * (i + 3);
    A[inner + 3] = bkj;
    bkj = Phi_[3 * i + 1];
    A[6 * i + 1] = bkj;
    A[inner + 4] = bkj;
    bkj = Phi_[3 * i + 2];
    A[6 * i + 2] = bkj;
    A[inner + 5] = bkj;
  }
  dv[0] = 0.0;
  dv[3] = -Z1_data[20];
  dv[6] = Z1_data[19];
  dv[1] = Z1_data[20];
  dv[4] = 0.0;
  dv[7] = -Z1_data[18];
  dv[2] = -Z1_data[19];
  dv[5] = Z1_data[18];
  dv[8] = 0.0;
  for (i = 0; i < 3; i++) {
    double d;
    double d1;
    bkj = dv[i];
    d = dv[i + 3];
    d1 = dv[i + 6];
    for (coffset = 0; coffset < 3; coffset++) {
      A[(i + 6 * coffset) + 3] =
          (bkj * Phi_[3 * coffset] + d * Phi_[3 * coffset + 1]) +
          d1 * Phi_[3 * coffset + 2];
    }
  }
  i = dZ1->size[0] * dZ1->size[1];
  dZ1->size[0] = 6;
  dZ1->size[1] = (int)(((double)x->size[0] + 5.0) - 1.0);
  emxEnsureCapacity_real_T(dZ1, i);
  dZ1_data = dZ1->data;
  inner = 6 * (int)(((double)x->size[0] + 5.0) - 1.0);
  for (i = 0; i < inner; i++) {
    dZ1_data[i] = 0.0;
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
    dZ1_data[i + 18] = 0.0;
    for (coffset = 0; coffset < 3; coffset++) {
      dZ1_data[i + 6 * coffset] =
          (Phi_[i] * dv[3 * coffset] + Phi_[i + 3] * dv[3 * coffset + 1]) +
          Phi_[i + 6] * dv[3 * coffset + 2];
      dZ1_data[i + 18] += Phi_[i + 3 * coffset] * XI[coffset + 3];
    }
  }
  if (x->size[0] + 4U < 5U) {
    i = 0;
  } else {
    i = 4;
  }
  emxInit_real_T(&C, 2);
  nc = BTh->size[1];
  coffset = C->size[0] * C->size[1];
  C->size[0] = 6;
  C->size[1] = BTh->size[1];
  emxEnsureCapacity_real_T(C, coffset);
  C_data = C->data;
  for (j = 0; j < nc; j++) {
    inner = j * 6;
    for (b_i = 0; b_i < 6; b_i++) {
      bkj = 0.0;
      for (k = 0; k < 6; k++) {
        bkj += A[k * 6 + b_i] * BTh_data[inner + k];
      }
      C_data[inner + b_i] = bkj;
    }
  }
  emxFree_real_T(&BTh);
  inner = C->size[1];
  for (coffset = 0; coffset < inner; coffset++) {
    for (nc = 0; nc < 6; nc++) {
      dZ1_data[nc + 6 * (i + coffset)] = C_data[nc + 6 * coffset];
    }
  }
  emxFree_real_T(&C);
}

static void b_binary_expand_op(emxArray_real_T *in1, const emxArray_real_T *in2,
                               const emxArray_real_T *in3, double in4,
                               const emxArray_real_T *in5,
                               const emxArray_real_T *in6, int in7,
                               const emxArray_real_T *in8,
                               const emxArray_real_T *in9)
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
  i = (in7 + 1) << 1;
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
  ForwardODEX(in2, b_in3, b_in6, *(double(*)[6]) & in8_data[6 * (i - 1)], in9,
              in1);
  emxFree_real_T(&b_in6);
  emxFree_real_T(&b_in3);
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

void computeForwardKinematicsFast(const emxArray_real_T *x, double ds,
                                  const double p0[3], const double Phi0[9],
                                  const emxArray_real_T *xia0,
                                  const emxArray_real_T *Th,
                                  const emxArray_real_T *Ba,
                                  emxArray_real_T *gtmp, emxArray_real_T *Jtmp)
{
  emxArray_real_T *K1Z1;
  emxArray_real_T *K2Z1;
  emxArray_real_T *Z1;
  emxArray_real_T *b_Th;
  emxArray_real_T *b_Z1;
  double c_a[36];
  const double *Th_data;
  const double *xia0_data;
  double Ns;
  double a;
  double b_a;
  double *Jtmp_data;
  double *K1Z1_data;
  double *K2Z1_data;
  double *Z1_data;
  double *gtmp_data;
  int b_i;
  int b_loop_ub;
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  int i5;
  int ii;
  int k;
  int loop_ub;
  int n;
  Th_data = Th->data;
  xia0_data = xia0->data;
  emxInit_real_T(&Z1, 2);
  /*  % states */
  /*       % spatial steps */
  /*       % position zero */
  /*     % phi zeroclc */
  /*     % intrinsic strain vector */
  /*       % evaluated Theta matrix */
  /*  state to strain matrix         */
  /*  compute total length */
  Ns = (double)Th->size[2] / 2.0;
  i = Z1->size[0] * Z1->size[1];
  Z1->size[0] = 6;
  Z1->size[1] = (int)(((double)x->size[0] + 5.0) - 1.0);
  emxEnsureCapacity_real_T(Z1, i);
  Z1_data = Z1->data;
  loop_ub = 6 * (int)(((double)x->size[0] + 5.0) - 1.0);
  for (i = 0; i < loop_ub; i++) {
    Z1_data[i] = 0.0;
  }
  for (i = 0; i < 3; i++) {
    Z1_data[6 * i] = Phi0[3 * i];
    Z1_data[6 * i + 1] = Phi0[3 * i + 1];
    Z1_data[6 * i + 2] = Phi0[3 * i + 2];
    Z1_data[i + 18] = p0[i];
  }
  i = gtmp->size[0] * gtmp->size[1] * gtmp->size[2];
  gtmp->size[0] = 4;
  gtmp->size[1] = 4;
  i1 = (int)Ns;
  gtmp->size[2] = (int)Ns;
  emxEnsureCapacity_real_T(gtmp, i);
  gtmp_data = gtmp->data;
  i = Jtmp->size[0] * Jtmp->size[1] * Jtmp->size[2];
  Jtmp->size[0] = 6;
  Jtmp->size[1] = x->size[0];
  Jtmp->size[2] = (int)Ns;
  emxEnsureCapacity_real_T(Jtmp, i);
  Jtmp_data = Jtmp->data;
  loop_ub = 6 * x->size[0] * (int)Ns;
  for (i = 0; i < loop_ub; i++) {
    Jtmp_data[i] = 0.0;
  }
  if ((int)Ns - 1 >= 0) {
    a = 0.66666666666666663 * ds;
    b_a = 0.25 * ds;
    if (x->size[0] + 4U < 5U) {
      i2 = 0;
      i3 = -1;
    } else {
      i2 = 4;
      i3 = x->size[0] + 3;
    }
    n = i3 - i2;
  }
  emxInit_real_T(&K1Z1, 2);
  emxInit_real_T(&K2Z1, 2);
  emxInit_real_T(&b_Th, 2);
  emxInit_real_T(&b_Z1, 2);
  for (ii = 0; ii < i1; ii++) {
    double Phi[9];
    double Rt[9];
    double b_xia0_data[6];
    /*  first EL-diff eval */
    i = (ii + 1) << 1;
    loop_ub = Th->size[0];
    b_loop_ub = Th->size[1];
    k = b_Th->size[0] * b_Th->size[1];
    b_Th->size[0] = Th->size[0];
    b_Th->size[1] = Th->size[1];
    emxEnsureCapacity_real_T(b_Th, k);
    K2Z1_data = b_Th->data;
    for (k = 0; k < b_loop_ub; k++) {
      for (b_i = 0; b_i < loop_ub; b_i++) {
        K2Z1_data[b_i + b_Th->size[0] * k] =
            Th_data[(b_i + Th->size[0] * k) +
                    Th->size[0] * Th->size[1] * (i - 2)];
      }
    }
    for (i4 = 0; i4 < 6; i4++) {
      b_xia0_data[i4] = xia0_data[i4 + 6 * (i - 2)];
    }
    ForwardODEX(x, Z1, b_Th, b_xia0_data, Ba, K1Z1);
    K1Z1_data = K1Z1->data;
    /*  second EL-diff eval */
    if (Z1->size[1] == K1Z1->size[1]) {
      k = b_Z1->size[0] * b_Z1->size[1];
      b_Z1->size[0] = 6;
      b_Z1->size[1] = Z1->size[1];
      emxEnsureCapacity_real_T(b_Z1, k);
      K2Z1_data = b_Z1->data;
      loop_ub = 6 * Z1->size[1];
      for (k = 0; k < loop_ub; k++) {
        K2Z1_data[k] = Z1_data[k] + a * K1Z1_data[k];
      }
      loop_ub = Th->size[0];
      b_loop_ub = Th->size[1];
      k = b_Th->size[0] * b_Th->size[1];
      b_Th->size[0] = Th->size[0];
      b_Th->size[1] = Th->size[1];
      emxEnsureCapacity_real_T(b_Th, k);
      K2Z1_data = b_Th->data;
      for (k = 0; k < b_loop_ub; k++) {
        for (b_i = 0; b_i < loop_ub; b_i++) {
          K2Z1_data[b_i + b_Th->size[0] * k] =
              Th_data[(b_i + Th->size[0] * k) +
                      Th->size[0] * Th->size[1] * (i - 1)];
        }
      }
      for (i5 = 0; i5 < 6; i5++) {
        b_xia0_data[i5] = xia0_data[i5 + 6 * (i - 1)];
      }
      ForwardODEX(x, b_Z1, b_Th, b_xia0_data, Ba, K2Z1);
      K2Z1_data = K2Z1->data;
    } else {
      b_binary_expand_op(K2Z1, x, Z1, a, K1Z1, Th, ii, xia0, Ba);
      K2Z1_data = K2Z1->data;
    }
    /*  update integrands */
    if (Z1->size[1] == K1Z1->size[1]) {
      loop_ub = 6 * Z1->size[1];
      i = Z1->size[0] * Z1->size[1];
      Z1->size[0] = 6;
      emxEnsureCapacity_real_T(Z1, i);
      Z1_data = Z1->data;
      for (i = 0; i < loop_ub; i++) {
        Z1_data[i] += b_a * (K1Z1_data[i] + 3.0 * K2Z1_data[i]);
      }
    } else {
      binary_expand_op(Z1, b_a, K1Z1, K2Z1);
      Z1_data = Z1->data;
    }
    /*  recover the kinematics entities */
    for (i = 0; i < 3; i++) {
      Phi[3 * i] = Z1_data[6 * i];
      Phi[3 * i + 1] = Z1_data[6 * i + 1];
      Phi[3 * i + 2] = Z1_data[6 * i + 2];
    }
    /* --------------------------------------------------------------------------
     */
    for (i = 0; i < 4; i++) {
      k = 4 * i + 16 * ii;
      gtmp_data[k] = 0.0;
      gtmp_data[k + 1] = 0.0;
      gtmp_data[k + 2] = 0.0;
      gtmp_data[k + 3] = 0.0;
    }
    gtmp_data[16 * ii + 15] = 1.0;
    /* --------------------------------------------------------------------------
     */
    for (i = 0; i < 3; i++) {
      k = 4 * i + 16 * ii;
      gtmp_data[k] = Phi[3 * i];
      Rt[3 * i] = Phi[i];
      loop_ub = 3 * i + 1;
      gtmp_data[k + 1] = Phi[loop_ub];
      Rt[loop_ub] = Phi[i + 3];
      loop_ub = 3 * i + 2;
      gtmp_data[k + 2] = Phi[loop_ub];
      Rt[loop_ub] = Phi[i + 6];
      gtmp_data[(i + 16 * ii) + 12] = Z1_data[i + 18];
    }
    /* --------------------------------------------------------------------------
     */
    memset(&c_a[0], 0, 36U * sizeof(double));
    for (i = 0; i < 3; i++) {
      Ns = Rt[3 * i];
      c_a[6 * i] = Ns;
      loop_ub = 6 * (i + 3);
      c_a[loop_ub + 3] = Ns;
      Ns = Rt[3 * i + 1];
      c_a[6 * i + 1] = Ns;
      c_a[loop_ub + 4] = Ns;
      Ns = Rt[3 * i + 2];
      c_a[6 * i + 2] = Ns;
      c_a[loop_ub + 5] = Ns;
    }
    Phi[0] = 0.0;
    Phi[1] = -Z1_data[20];
    Phi[2] = Z1_data[19];
    Phi[3] = Z1_data[20];
    Phi[4] = 0.0;
    Phi[5] = -Z1_data[18];
    Phi[6] = -Z1_data[19];
    Phi[7] = Z1_data[18];
    Phi[8] = 0.0;
    for (i = 0; i < 3; i++) {
      double d;
      double d1;
      Ns = Rt[i];
      d = Rt[i + 3];
      d1 = Rt[i + 6];
      for (k = 0; k < 3; k++) {
        c_a[(i + 6 * k) + 3] =
            (Ns * Phi[3 * k] + d * Phi[3 * k + 1]) + d1 * Phi[3 * k + 2];
      }
    }
    i = K1Z1->size[0] * K1Z1->size[1];
    K1Z1->size[0] = 6;
    K1Z1->size[1] = (i3 - i2) + 1;
    emxEnsureCapacity_real_T(K1Z1, i);
    K1Z1_data = K1Z1->data;
    for (b_loop_ub = 0; b_loop_ub <= n; b_loop_ub++) {
      loop_ub = b_loop_ub * 6;
      for (b_i = 0; b_i < 6; b_i++) {
        Ns = 0.0;
        for (k = 0; k < 6; k++) {
          i = loop_ub + k;
          Ns += c_a[k * 6 + b_i] * Z1_data[i % 6 + 6 * (i2 + i / 6)];
        }
        K1Z1_data[loop_ub + b_i] = Ns;
      }
    }
    loop_ub = K1Z1->size[1];
    for (i = 0; i < loop_ub; i++) {
      for (k = 0; k < 6; k++) {
        b_i = k + 6 * i;
        Jtmp_data[b_i + 6 * Jtmp->size[1] * ii] = K1Z1_data[b_i];
      }
    }
  }
  emxFree_real_T(&b_Z1);
  emxFree_real_T(&b_Th);
  emxFree_real_T(&K2Z1);
  emxFree_real_T(&K1Z1);
  emxFree_real_T(&Z1);
}

/* End of code generation (computeForwardKinematicsFast.c) */
