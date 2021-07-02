#include "gjk.h"
#include <immintrin.h>
#include <stdio.h>

// I miss my C++ templates...

/**
 * Dot product of two 3D vectors.
 * */
float dot3D_soa_vectorized(float x0, float x1, float x2, const float *y) {
  return x0 * y[0] + x1 * y[1] + x2 * y[2];
}

/**
 * Cross product of two 3D vectors. Write into res.
 * */
void cross3D_soa_vectorized(const float *x, const float *y, float *res) {
  res[0] = x[1] * y[2] - x[2] * y[1];
  res[1] = x[2] * y[0] - x[0] * y[2];
  res[2] = x[0] * y[1] - x[1] * y[0];
}

/**
 * Among a set of N points in 3D, find the point furthest along dir.
 * */
[[gnu::target("avx512dq")]]
void support3D_soa_v_slow_idx(float** points, const float *dir, float *res, int N) {
  __m512 vertex_x_coord_0, vertex_y_coord_0, vertex_z_coord_0;
  
  // Variables used to compute dot product
  __m512 tmp_0, tmp_1, tmp_2, dp_v;
  
  // Variables used to find argmax
  __m512i tmp_3, tmp_4;
  
  __m512 dpmax_v = _mm512_set1_ps(-1e9);
  __m512i argmax_v = _mm512_set1_epi32(0);
  __m512i curr_idx;

  // Direction vector
  __m512 dir_x = _mm512_set1_ps(dir[0]);
  __m512 dir_y = _mm512_set1_ps(dir[1]);
  __m512 dir_z = _mm512_set1_ps(dir[2]);
  
  // Used to compute dpmax/argmax
  __m512i ones = _mm512_set1_epi32(-1);
  __m512i zeros = _mm512_set1_epi32(0);
  __mmask16 mask_0;
  
  // Used to compute current indices
  __m512i sixteens = _mm512_set1_epi32(16);

  float dp, dpmax = -1e9;
  int argmax = 0;

  int i = 0;
  for (i = 0; i < N - 15; i += 16) {
    // Compute current indices
    curr_idx = _mm512_set_epi32(i + 15, i + 14, i + 13, i + 12, i + 11, i + 10, i + 9, i + 8,
                             i + 7, i + 6, i + 5, i + 4, i + 3, i + 2, i + 1, i);

    // Load vertices
    vertex_x_coord_0 = _mm512_loadu_ps(&points[0][i]);
    vertex_y_coord_0 = _mm512_loadu_ps(&points[1][i]);
    vertex_z_coord_0 = _mm512_loadu_ps(&points[2][i]);
  
    // Compute dot product
    tmp_0 = _mm512_mul_ps(vertex_x_coord_0, dir_x);
    tmp_1 = _mm512_fmadd_ps(dir_y, vertex_y_coord_0, tmp_0);
    dp_v = _mm512_fmadd_ps(dir_z, vertex_z_coord_0, tmp_1);

    // Find local dpmax/argmax
    mask_0 = _mm512_cmp_ps_mask(dp_v, dpmax_v, _CMP_GT_OQ);
    dpmax_v = _mm512_max_ps(dp_v, dpmax_v);
    tmp_3 = _mm512_mask_and_epi32(zeros, mask_0, curr_idx, ones);
    tmp_4 = _mm512_mask_andnot_epi32(argmax_v, mask_0, ones, argmax_v);
    argmax_v = _mm512_or_epi32(tmp_3, tmp_4);
  }

   // Find dpmax/argmax
  if (i != 0) {
    float dps[16];
    int argmaxs[16];

    _mm512_storeu_ps(dps, dpmax_v);
    _mm512_storeu_epi32(argmaxs, argmax_v);
    for (int j = 0; j < 16; j++) {
      if (dps[j] > dpmax) {
        dpmax = dps[j];
        argmax = argmaxs[j];
      }
    }
  }

  for (; i < N; i++) {
    dp = dot3D_soa_vectorized(points[0][i], points[1][i], points[2][i], dir);
    if (dp > dpmax) {
      dpmax = dp;
      argmax = i;
    }
  }

  res[0] = points[0][argmax];
  res[1] = points[1][argmax];
  res[2] = points[2][argmax];
}

[[gnu::target("avx512dq")]]
void support3D_soa_vectorized(float** points, const float *dir, float *res, int N) {
  __m512 vertex_x_coord_0, vertex_y_coord_0, vertex_z_coord_0;
  
  // Variables used to compute dot product
  __m512 tmp_0, tmp_1, tmp_2, dp_v;
  
  // Variables used to compute argmax
  __m512i tmp_3, tmp_4;
  
  __m512 dpmax_v = _mm512_set1_ps(-1e9);
  __m512i argmax_v = _mm512_set1_epi32(0);
  __m512i curr_idx = _mm512_set_epi32(-1, -2, -3, -4, -5, -6, -7, -8,
                                      -9, -10, -11, -12, -13, -14, -15, -16);
   
  // Direction vector
  __m512 dir_x = _mm512_set1_ps(dir[0]);
  __m512 dir_y = _mm512_set1_ps(dir[1]);
  __m512 dir_z = _mm512_set1_ps(dir[2]);
  
  // Variables used to compute argmax
  __m512i ones = _mm512_set1_epi32(-1);
  __m512i zeros = _mm512_set1_epi32(0);
  __mmask16 mask_0;
  
  __m512i sixteens = _mm512_set1_epi32(16);

  float dp, dpmax = -1e9;
  int argmax = 0;

  int i = 0;
  for (i = 0; i < N - 15; i += 16) {
    // Compute current indices
    curr_idx = _mm512_add_epi32(curr_idx, sixteens);

    // Load current vertices
    vertex_x_coord_0 = _mm512_loadu_ps(&points[0][i]);
    vertex_y_coord_0 = _mm512_loadu_ps(&points[1][i]);
    vertex_z_coord_0 = _mm512_loadu_ps(&points[2][i]);

    // Compute dot product
    tmp_0 = _mm512_mul_ps(vertex_x_coord_0, dir_x);
    tmp_1 = _mm512_fmadd_ps(dir_y, vertex_y_coord_0, tmp_0);
    dp_v = _mm512_fmadd_ps(dir_z, vertex_z_coord_0, tmp_1);

    // Find local dpmax/argmax
    mask_0 = _mm512_cmp_ps_mask(dp_v, dpmax_v, _CMP_GT_OQ);
    dpmax_v = _mm512_max_ps(dp_v, dpmax_v);
    tmp_3 = _mm512_mask_and_epi32(zeros, mask_0, curr_idx, ones);
    tmp_4 = _mm512_mask_andnot_epi32(argmax_v, mask_0, ones, argmax_v);
    argmax_v = _mm512_or_epi32(tmp_3, tmp_4);
  }

  // Find dpmax/argmax
  if (i != 0) {
    float dps[16];
    int argmaxs[16];

    _mm512_storeu_ps(dps, dpmax_v);
    _mm512_storeu_epi32(argmaxs, argmax_v);
    for (int j = 0; j < 16; j++) {
      if (dps[j] > dpmax) {
        dpmax = dps[j];
        argmax = argmaxs[j];
      }
    }
  }

  for (; i < N; i++) {
    dp = dot3D_soa_vectorized(points[0][i], points[1][i], points[2][i], dir);
    if (dp > dpmax) {
      dpmax = dp;
      argmax = i;
    }
  }

  res[0] = points[0][argmax];
  res[1] = points[1][argmax];
  res[2] = points[2][argmax];
}

[[gnu::target("avx512dq")]]
void support3D_soa_vectorized_lu2(float** points, const float *dir, float *res, int N) {
  __m512 vertex_x_coord_0, vertex_y_coord_0, vertex_z_coord_0;
  __m512 vertex_x_coord_1, vertex_y_coord_1, vertex_z_coord_1;

  __m512 tmp_0, tmp_1, tmp_2, dp_v;
  __m512i tmp_3, tmp_4;
  __m512 tmp_5, tmp_6, tmp_7, dp_v_1;
  __m512i tmp_8, tmp_9;

  __m512 dpmax_v = _mm512_set1_ps(-1e9);
  __m512i argmax_v = _mm512_set1_epi32(0);
  __m512i curr_idx = _mm512_set_epi32(-17, -18, -19, -20, -21, -22, -23, -24,
                                    -25, -26, -27, -28, -29, -30, -31, -32);
  __m512 dpmax_v_1 = _mm512_set1_ps(-1e9);
  __m512i argmax_v_1 = _mm512_set1_epi32(0);
  __m512i curr_idx_1 = _mm512_set_epi32(-1, -2, -3, -4, -5, -6, -7, -8,
                                      -9, -10, -11, -12, -13, -14, -15, -16);

  __m512 dir_x = _mm512_set1_ps(dir[0]);
  __m512 dir_y = _mm512_set1_ps(dir[1]);
  __m512 dir_z = _mm512_set1_ps(dir[2]);
  __m512i ones = _mm512_set1_epi32(-1);
  __m512i zeros = _mm512_set1_epi32(0);
  __m512i thirtytwos = _mm512_set1_epi32(32);

  __mmask16 mask_0;
  __mmask16 mask_1;

  float dp, dpmax = -1e9;
  int argmax = 0;

  int i = 0;
  for (i = 0; i < N - 31; i += 32) {
    // Compute current indices
    curr_idx = _mm512_add_epi32(curr_idx, thirtytwos);
    curr_idx_1 = _mm512_add_epi32(curr_idx_1, thirtytwos);

    // Load current vertices
    vertex_x_coord_0 = _mm512_loadu_ps(&points[0][i]);
    vertex_y_coord_0 = _mm512_loadu_ps(&points[1][i]);
    vertex_z_coord_0 = _mm512_loadu_ps(&points[2][i]);

    vertex_x_coord_1 = _mm512_loadu_ps(&points[0][i + 16]);
    vertex_y_coord_1 = _mm512_loadu_ps(&points[1][i + 16]);
    vertex_z_coord_1 = _mm512_loadu_ps(&points[2][i + 16]);

    // Compute dot products
    tmp_0 = _mm512_mul_ps(vertex_x_coord_0, dir_x);
    tmp_1 = _mm512_fmadd_ps(dir_y, vertex_y_coord_0, tmp_0);
    dp_v = _mm512_fmadd_ps(dir_z, vertex_z_coord_0, tmp_1);

    tmp_5 = _mm512_mul_ps(vertex_x_coord_1, dir_x);
    tmp_6 = _mm512_fmadd_ps(dir_y, vertex_y_coord_1, tmp_5);
    dp_v_1 = _mm512_fmadd_ps(dir_z, vertex_z_coord_1, tmp_6);

    // Find local dpmax/argmax
    mask_0 = _mm512_cmp_ps_mask(dp_v, dpmax_v, _CMP_GT_OQ);
    dpmax_v = _mm512_max_ps(dp_v, dpmax_v);
    tmp_3 = _mm512_mask_and_epi32(zeros, mask_0, curr_idx, ones);
    tmp_4 = _mm512_mask_andnot_epi32(argmax_v, mask_0, ones, argmax_v);
    argmax_v = _mm512_or_epi32(tmp_3, tmp_4);

    mask_1 = _mm512_cmp_ps_mask(dp_v_1, dpmax_v_1, _CMP_GT_OQ);
    dpmax_v_1 = _mm512_max_ps(dp_v_1, dpmax_v_1);
    tmp_8 = _mm512_mask_and_epi32(zeros, mask_1, curr_idx_1, ones);
    tmp_9 = _mm512_mask_andnot_epi32(argmax_v_1, mask_1, ones, argmax_v_1);
    argmax_v_1 = _mm512_or_epi32(tmp_8, tmp_9);
  }

  // Find dpmax/argmax
  if (i != 0) {
    mask_0 = _mm512_cmp_ps_mask(dpmax_v, dpmax_v_1, _CMP_GT_OQ);
    dpmax_v = _mm512_max_ps(dpmax_v, dpmax_v_1);
    tmp_3 = _mm512_mask_and_epi32(zeros, mask_0, argmax_v, ones);
    tmp_4 = _mm512_mask_andnot_epi32(argmax_v_1, mask_0, ones, argmax_v_1);
    argmax_v = _mm512_or_epi32(tmp_3, tmp_4);

    float dps[16];
    int argmaxs[16];

    _mm512_storeu_ps(dps, dpmax_v);
    _mm512_storeu_epi32(argmaxs, argmax_v);
    for (int j = 0; j < 16; j++) {
      if (dps[j] > dpmax) {
        dpmax = dps[j];
        argmax = argmaxs[j];
      }
    }
  }

  for (; i < N; i++) {
    dp = dot3D_soa_vectorized(points[0][i], points[1][i], points[2][i], dir);
    if (dp > dpmax) {
      dpmax = dp;
      argmax = i;
    }
  }

  res[0] = points[0][argmax];
  res[1] = points[1][argmax];
  res[2] = points[2][argmax];
}

/**
 * "Inverse support"
 * SUBTRACTS THE FOUND POINT FROM res.
 *
 * Among a set of N points in 3D, find the point furthest along NEGATIVE dir.
 * argmax_b (-a)Tb = argmin_b aTb
 * This is so we don't need to invert/negate dir
 * */
[[gnu::target("avx512dq")]]
void invsup3D_soa_v_slow_idx(float** points, const float *dir, float *res, int N) {
  __m512 inv_vertex_x_coord_0, inv_vertex_y_coord_0, inv_vertex_z_coord_0;
  __m512 inv_tmp_0, inv_tmp_1, inv_tmp_2, inv_dp_v;
  __m512i inv_tmp_3, inv_tmp_4;
  __m512 inv_dpmax_v = _mm512_set1_ps(1e9);
  __m512i inv_argmin_v = _mm512_set1_epi32(0);
  __m512i curr_idx;

  __m512 dir_x = _mm512_set1_ps(dir[0]);
  __m512 dir_y = _mm512_set1_ps(dir[1]);
  __m512 dir_z = _mm512_set1_ps(dir[2]);
  __m512i ones = _mm512_set1_epi32(-1);
  __m512i zeros = _mm512_set1_epi32(0);
  __m512i sixteens = _mm512_set1_epi32(16);
  __mmask16 mask_0;

  float inv_dp, inv_dpmax = 1e9;
  int argmin = 0;
  int k;
  for (k = 0; k < N - 15; k += 16) {
    // Compute current indices
    curr_idx = _mm512_set_epi32(k + 15, k + 14, k + 13, k + 12, k + 11, k + 10, k + 9, k + 8,
                                k + 7, k + 6, k + 5, k + 4, k + 3, k + 2, k + 1, k);

    // Load vertices
    inv_vertex_x_coord_0 = _mm512_loadu_ps(&points[0][k]);
    inv_vertex_y_coord_0 = _mm512_loadu_ps(&points[1][k]);
    inv_vertex_z_coord_0 = _mm512_loadu_ps(&points[2][k]);

    // Compute dot products
    inv_tmp_0 = _mm512_mul_ps(inv_vertex_x_coord_0, dir_x);
    inv_tmp_1 = _mm512_fmadd_ps(dir_y, inv_vertex_y_coord_0, inv_tmp_0);
    inv_dp_v = _mm512_fmadd_ps(dir_z, inv_vertex_z_coord_0, inv_tmp_1);

    // Find local argmin/dpmin
    mask_0 = _mm512_cmp_ps_mask(inv_dp_v, inv_dpmax_v, _CMP_LT_OS);
    inv_dpmax_v = _mm512_min_ps(inv_dp_v, inv_dpmax_v);
    inv_tmp_3 = _mm512_mask_and_epi32(zeros, mask_0, curr_idx, ones);
    inv_tmp_4 = _mm512_mask_andnot_epi32(inv_argmin_v, mask_0, ones, inv_argmin_v);
    inv_argmin_v = _mm512_or_epi32(inv_tmp_3, inv_tmp_4);
  }

  // Find dpmin/argmin
  if (k != 0) {
    float inv_dps[16];
    int argmins[16];

    _mm512_storeu_ps(inv_dps, inv_dpmax_v);
    _mm512_storeu_epi32(argmins, inv_argmin_v);

    for (int j = 0; j < 16; j++) {
      if (inv_dps[j] < inv_dpmax) {
        inv_dpmax = inv_dps[j];
        argmin = argmins[j];
      }
    }
  }

  for (; k < N; k++) {
    inv_dp = dot3D_soa_vectorized(points[0][k], points[1][k], points[2][k], dir);
    if (inv_dp < inv_dpmax) {
      inv_dpmax = inv_dp;
      argmin = k;
    }
  }

  res[0] -= points[0][argmin];
  res[1] -= points[1][argmin];
  res[2] -= points[2][argmin];
}

[[gnu::target("avx512dq")]]
void invsup3D_soa_vectorized(float** points, const float *dir, float *res, int N) {
  __m512 inv_vertex_x_coord_0, inv_vertex_y_coord_0, inv_vertex_z_coord_0;
  __m512 inv_tmp_0, inv_tmp_1, inv_tmp_2, inv_dp_v;
  __m512i inv_tmp_3, inv_tmp_4;
  __m512 inv_dpmax_v = _mm512_set1_ps(1e9);
  __m512i inv_argmin_v = _mm512_set1_epi32(0);
  __m512i curr_idx = _mm512_set_epi32(-1, -2, -3, -4, -5, -6, -7, -8,
                                      -9, -10, -11, -12, -13, -14, -15, -16);
  __m512 dir_x = _mm512_set1_ps(dir[0]);
  __m512 dir_y = _mm512_set1_ps(dir[1]);
  __m512 dir_z = _mm512_set1_ps(dir[2]);
  __m512i ones = _mm512_set1_epi32(-1);
  __m512i zeros = _mm512_set1_epi32(0);
  __m512i sixteens = _mm512_set1_epi32(16);
  __mmask16 mask_0;

  float inv_dp, inv_dpmax = 1e9;
  int argmin = 0;
  int k;
  for (k = 0; k < N - 15; k += 16) {
    // Compute current indices
    curr_idx = _mm512_add_epi32(curr_idx, sixteens);

    // Load vertices
    inv_vertex_x_coord_0 = _mm512_loadu_ps(&points[0][k]);
    inv_vertex_y_coord_0 = _mm512_loadu_ps(&points[1][k]);
    inv_vertex_z_coord_0 = _mm512_loadu_ps(&points[2][k]);
  
    // Compute dot product
    inv_tmp_0 = _mm512_mul_ps(inv_vertex_x_coord_0, dir_x);
    inv_tmp_1 = _mm512_fmadd_ps(dir_y, inv_vertex_y_coord_0, inv_tmp_0);
    inv_dp_v = _mm512_fmadd_ps(dir_z, inv_vertex_z_coord_0, inv_tmp_1);

    // Find local dpmin/argmin
    mask_0 = _mm512_cmp_ps_mask(inv_dp_v, inv_dpmax_v, _CMP_LT_OS);
    inv_dpmax_v = _mm512_min_ps(inv_dp_v, inv_dpmax_v);
    inv_tmp_3 = _mm512_mask_and_epi32(zeros, mask_0, curr_idx, ones);
    inv_tmp_4 = _mm512_mask_andnot_epi32(inv_argmin_v, mask_0, ones, inv_argmin_v);
    inv_argmin_v = _mm512_or_epi32(inv_tmp_3, inv_tmp_4);
  }

  // Find dpmin/argmin
  if (k != 0) {
    float inv_dps[16];
    int argmins[16];

    _mm512_storeu_ps(inv_dps, inv_dpmax_v);
    _mm512_storeu_epi32(argmins, inv_argmin_v);

    for (int j = 0; j < 16; j++) {
      if (inv_dps[j] < inv_dpmax) {
        inv_dpmax = inv_dps[j];
        argmin = argmins[j];
      }
    }
  }

  for (; k < N; k++) {
    inv_dp = dot3D_soa_vectorized(points[0][k], points[1][k], points[2][k], dir);
    if (inv_dp < inv_dpmax) {
      inv_dpmax = inv_dp;
      argmin = k;
    }
  }

  res[0] -= points[0][argmin];
  res[1] -= points[1][argmin];
  res[2] -= points[2][argmin];
}

[[gnu::target("avx512dq")]]
void invsup3D_soa_vectorized_lu2(float** points, const float *dir, float *res, int N) {
  __m512 inv_vertex_x_coord_0, inv_vertex_y_coord_0, inv_vertex_z_coord_0;
  __m512 inv_vertex_x_coord_1, inv_vertex_y_coord_1, inv_vertex_z_coord_1;

  __m512 inv_tmp_0, inv_tmp_1, inv_tmp_2, inv_dp_v;
  __m512i inv_tmp_3, inv_tmp_4;
  __m512 inv_tmp_5, inv_tmp_6, inv_tmp_7, inv_dp_v_1;
  __m512i inv_tmp_8, inv_tmp_9;

  __m512 inv_dpmax_v = _mm512_set1_ps(1e9);
  __m512i inv_argmin_v = _mm512_set1_epi32(0);
  __m512i curr_idx = _mm512_set_epi32(-17, -18, -19, -20, -21, -22, -23, -24,
                                      -25, -26, -27, -28, -29, -30, -31, -32);
  __m512 inv_dpmax_v_1 = _mm512_set1_ps(1e9);
  __m512i inv_argmin_v_1 = _mm512_set1_epi32(0);
  __m512i curr_idx_1 = _mm512_set_epi32(-1, -2, -3, -4, -5, -6, -7, -8,
                                        -9, -10, -11, -12, -13, -14, -15, -16);

  __m512 dir_x = _mm512_set1_ps(dir[0]);
  __m512 dir_y = _mm512_set1_ps(dir[1]);
  __m512 dir_z = _mm512_set1_ps(dir[2]);
  __m512i ones = _mm512_set1_epi32(-1);
  __m512i zeros = _mm512_set1_epi32(0);
  __m512i thirtytwos = _mm512_set1_epi32(32);

  __mmask16 mask_0;
  __mmask16 mask_1;

  float inv_dp, inv_dpmax = 1e9;
  int argmin = 0;
  int k;
  for (k = 0; k < N - 31; k += 32) {
    // Compute current indices
    curr_idx = _mm512_add_epi32(curr_idx, thirtytwos);
    curr_idx_1 = _mm512_add_epi32(curr_idx_1, thirtytwos);

    // Load vertices
    inv_vertex_x_coord_0 = _mm512_loadu_ps(&points[0][k]);
    inv_vertex_y_coord_0 = _mm512_loadu_ps(&points[1][k]);
    inv_vertex_z_coord_0 = _mm512_loadu_ps(&points[2][k]);

    inv_vertex_x_coord_1 = _mm512_loadu_ps(&points[0][k + 16]);
    inv_vertex_y_coord_1 = _mm512_loadu_ps(&points[1][k + 16]);
    inv_vertex_z_coord_1 = _mm512_loadu_ps(&points[2][k + 16]);

    // Compute dot products
    inv_tmp_0 = _mm512_mul_ps(inv_vertex_x_coord_0, dir_x);
    inv_tmp_1 = _mm512_fmadd_ps(dir_y, inv_vertex_y_coord_0, inv_tmp_0);
    inv_dp_v = _mm512_fmadd_ps(dir_z, inv_vertex_z_coord_0, inv_tmp_1);

    inv_tmp_5 = _mm512_mul_ps(inv_vertex_x_coord_1, dir_x);
    inv_tmp_6 = _mm512_fmadd_ps(dir_y, inv_vertex_y_coord_1, inv_tmp_5);
    inv_dp_v = _mm512_fmadd_ps(dir_z, inv_vertex_z_coord_1, inv_tmp_6);

    // Find local dpmin/argmin
    mask_0 = _mm512_cmp_ps_mask(inv_dp_v, inv_dpmax_v, _CMP_LT_OQ);
    inv_dpmax_v = _mm512_min_ps(inv_dp_v, inv_dpmax_v);
    inv_tmp_3 = _mm512_mask_and_epi32(zeros, mask_0, curr_idx, ones);
    inv_tmp_4 = _mm512_mask_andnot_epi32(inv_argmin_v, mask_0, ones, inv_argmin_v);
    inv_argmin_v = _mm512_or_epi32(inv_tmp_3, inv_tmp_4);

    mask_1 = _mm512_cmp_ps_mask(inv_dp_v_1, inv_dpmax_v_1, _CMP_LT_OQ);
    inv_dpmax_v_1 = _mm512_min_ps(inv_dp_v_1, inv_dpmax_v_1);
    inv_tmp_8 = _mm512_mask_and_epi32(zeros, mask_1, curr_idx_1, ones);
    inv_tmp_9 = _mm512_mask_andnot_epi32(inv_argmin_v_1, mask_1, ones, inv_argmin_v_1);
    inv_argmin_v_1 = _mm512_or_epi32(inv_tmp_8, inv_tmp_9);
  }

  // Find dpmin/argmin
  if (k != 0) {
    mask_0 = _mm512_cmp_ps_mask(inv_dpmax_v, inv_dpmax_v_1, _CMP_LT_OQ);
    inv_dpmax_v = _mm512_min_ps(inv_dpmax_v, inv_dpmax_v_1);
    inv_tmp_3 = _mm512_mask_and_epi32(zeros, mask_0, inv_argmin_v, ones);
    inv_tmp_4 = _mm512_mask_andnot_epi32(inv_argmin_v_1, mask_0, ones, inv_argmin_v_1);
    inv_argmin_v = _mm512_or_epi32(inv_tmp_3, inv_tmp_4);

    float inv_dps[16];
    int argmins[16];

    _mm512_storeu_ps(inv_dps, inv_dpmax_v);
    _mm512_storeu_epi32(argmins, inv_argmin_v);

    for (int j = 0; j < 16; j++) {
      if (inv_dps[j] < inv_dpmax) {
        inv_dpmax = inv_dps[j];
        argmin = argmins[j];
      }
    }
  }

  for (; k < N; k++) {
    inv_dp = dot3D_soa_vectorized(points[0][k], points[1][k], points[2][k], dir);
    if (inv_dp < inv_dpmax) {
      inv_dpmax = inv_dp;
      argmin = k;
    }
  }

  res[0] -= points[0][argmin];
  res[1] -= points[1][argmin];
  res[2] -= points[2][argmin];
}

int ds_line3D_soa_vectorized(struct Simplex3D* s, float *dir) {
  // WE KNOW A is PAST THE ORIGIN -> new Simplex s is [A, B]
  // Unlike in the video by Casey, I don't think we need to permute A and B.
  // p[0] is B, p[1] is A (the new point)
  const float(* p)[3] = s->p;
  float ab[3] = {p[0][0] - p[1][0], p[0][1] - p[1][1], p[0][2] - p[1][2]};
  float a0[3] = {-p[1][0], -p[1][1], -p[1][2]};

  // Perform the float-cross-product
  float tmp[3];
  cross3D_soa_vectorized(ab, a0, tmp);
  cross3D_soa_vectorized(tmp, ab, dir);

  // We don't (know if we) enclose 0
  return 0;
}

int ds_triangle3D_soa_vectorized(struct Simplex3D* s, float *dir) {
  float(* p)[3] = s->p;
  // p[0] is B, p[1] is C, p[2] is A (the new point)
  float ab[3] = {p[0][0] - p[2][0], p[0][1] - p[2][1], p[0][2] - p[2][2]};
  float ac[3] = {p[1][0] - p[2][0], p[1][1] - p[2][1], p[1][2] - p[2][2]};
  float a0[3] = {-p[2][0], -p[2][1], -p[2][2]};
  float abc[3];
  cross3D_soa_vectorized(ab, ac, abc);

  float tst[3];
  cross3D_soa_vectorized(abc, ac, tst);
  if (dot3D(tst, a0) > 0.0) {
    // [A, C]
    s->num_points = 2;
    p[0][0] = p[2][0];
    p[0][1] = p[2][1];
    p[0][2] = p[2][2];

    cross3D_soa_vectorized(ac, a0, tst);
    cross3D_soa_vectorized(tst, ac, dir);
  } else {
    cross3D_soa_vectorized(ab, abc, tst);
    if (dot3D(tst, a0) > 0.0) {
      // [A, B]
      s->num_points = 2;
      // I don't think we need to permute from [B, A] to [A, B]
      p[1][0] = p[2][0];
      p[1][1] = p[2][1];
      p[1][2] = p[2][2];

      cross3D_soa_vectorized(ab, a0, tst);
      cross3D_soa_vectorized(tst, ab, dir);
    }
    else {
      if (dot3D(abc, a0) > 0.0) {
        // [A,B,C]
        // Permutation [B, C, A] is already correct in this case.
        dir[0] = abc[0];
        dir[1] = abc[1];
        dir[2] = abc[2];
      }
      else {
        // [A,C,B], but we'll use [C, B, A] (swap B and C)
        float temp;  // swap B and C
        SWAP(p[0][0], p[1][0], temp);
        SWAP(p[0][1], p[1][1], temp);
        SWAP(p[0][2], p[1][2], temp);

        dir[0] = -abc[0];
        dir[1] = -abc[1];
        dir[2] = -abc[2];
      }
    }
  }

  // We don't (know if we) enclose 0
  return 0;
}

int ds_tetrahedron3D_soa_vectorized(struct Simplex3D* s, float *dir) {
  float(* p)[3] = s->p;
  // p[0] is B, p[1] is C, p[2] is D, p[3] is A (the new point)
  float ab[3] = {p[0][0] - p[3][0], p[0][1] - p[3][1], p[0][2] - p[3][2]};
  float ac[3] = {p[1][0] - p[3][0], p[1][1] - p[3][1], p[1][2] - p[3][2]};
  float ad[3] = {p[2][0] - p[3][0], p[2][1] - p[3][1], p[2][2] - p[3][2]};
  float a0[3] = {-p[3][0], -p[3][1], -p[3][2]};
  float abc[3], acd[3], adb[3];
  float tst[3];

  // The vectors abc, acd and adb point OUTSIDE the tetrahedron.
  cross3D_soa_vectorized(ab, ac, abc);
  if (dot3D(abc, a0) > 0.0) {
    cross3D_soa_vectorized(ac, ad, acd);
    if (dot3D(acd, a0) > 0.0) {
      s->num_points = 2;
      // [A, C]
      p[0][0] = p[3][0];
      p[0][1] = p[3][1];
      p[0][2] = p[3][2];

      // Like in the line case (I think ?)
      cross3D_soa_vectorized(ac, a0, tst);
      cross3D_soa_vectorized(tst, ac, dir);
    }
    else {
      cross3D_soa_vectorized(ad, ab, adb);
      if (dot3D(adb, a0) > 0.0) {
        s->num_points = 2;
        // [A, B]  (tho we return [B, A] as line perm doesn't matter).
        p[1][0] = p[3][0];
        p[1][1] = p[3][1];
        p[1][2] = p[3][2];

        // Like in the line case (I think ?)
        cross3D_soa_vectorized(ab, a0, tst);
        cross3D_soa_vectorized(tst, ab, dir);
      }
      else {
        s->num_points = 3;
        // [A, B, C]  (we return [B, C, A], an equivalent perm).
        p[2][0] = p[3][0];
        p[2][1] = p[3][1];
        p[2][2] = p[3][2];

        dir[0] = abc[0];
        dir[1] = abc[1];
        dir[2] = abc[2];
      }
    }
  }
  else {
    cross3D_soa_vectorized(ac, ad, acd);
    cross3D_soa_vectorized(ad, ab, adb);
    if (dot3D(acd, a0) > 0.0) {
      if (dot3D(adb, a0) > 0.0) {
        s->num_points = 2;
        // [A, D]  ... For once we actually need to move two vectors ;-;
        p[0][0] = p[3][0];
        p[0][1] = p[3][1];
        p[0][2] = p[3][2];

        p[1][0] = p[2][0];
        p[1][1] = p[2][1];
        p[1][2] = p[2][2];

        // Like in the line case
        cross3D_soa_vectorized(ad, a0, tst);
        cross3D_soa_vectorized(tst, ad, dir);
      }
      else {
        s->num_points = 3;
        // [A, C, D]
        p[0][0] = p[3][0];
        p[0][1] = p[3][1];
        p[0][2] = p[3][2];

        dir[0] = acd[0];
        dir[1] = acd[1];
        dir[2] = acd[2];
      }
    }
    else {
      if (dot3D(adb, a0) > 0.0) {
        s->num_points = 3;
        // [A, D, B]  (we return [B, A, D], an equivalent perm).
        p[1][0] = p[3][0];
        p[1][1] = p[3][1];
        p[1][2] = p[3][2];

        dir[0] = adb[0];
        dir[1] = adb[1];
        dir[2] = adb[2];
      }
      else {  // WE HAVE FOUND A COLLISION
        return 1;
      }
    }
  }
  return 0;  // Didn't prove a collision, if we reach this.
}

// if AB ''means'' if AB x AO > 0
int do_simplex3D_soa_vectorized(struct Simplex3D* s, float *dir) {
  switch (s->num_points) {
    case 2:return ds_line3D_soa_vectorized(s, dir);
    case 3:return ds_triangle3D_soa_vectorized(s, dir);
    case 4:return ds_tetrahedron3D_soa_vectorized(s, dir);
    default:;
  }
  return -1;  // unreachable unless s points to an invalid Simplex
}

int do_intersect3D_soa_v_slow_idx(const struct CHObject_soa* obj1, const struct CHObject_soa* obj2) {
  // The search direction
  float d[3] = {1.f, 1.f, 1.f};  // could also be random
  // The point we examine in each iteration
  float a[3];

  // Our simplex that do_simplex modifies.
  struct Simplex3D s;
  s.num_points = 1;

  // S = Support(A-B) = Support(A) - Invsup(B)
  support3D_soa_v_slow_idx(obj1->points, d, s.p[0], obj1->num_points);
  invsup3D_soa_v_slow_idx(obj2->points, d, s.p[0], obj2->num_points);

  // d = -S
  d[0] = -s.p[0][0];
  d[1] = -s.p[0][1];
  d[2] = -s.p[0][2];

  int max_iter = 100;  // most likely we need less
  while (max_iter--) {
    // a = Support(A-B) = Support(A) - Invsup(B)
    support3D_soa_v_slow_idx(obj1->points, d, a, obj1->num_points);
    invsup3D_soa_v_slow_idx(obj2->points, d, a, obj2->num_points);

    // Are we past the origin? If no, we will never enclose it.
    if (dot3D(a, d) < 0.) return 0;

    // Add newly found point to Simplex
    // AT THIS POINT WE KNOW a IS PAST THE ORIGIN
    s.p[s.num_points][0] = a[0];
    s.p[s.num_points][1] = a[1];
    s.p[s.num_points][2] = a[2];
    s.num_points++;

    if (do_simplex3D_soa_vectorized(&s, d)) return 1;
  }

  return 1; // Most likely not reachable. Say we intersect by default.
}

int do_intersect3D_soa_vectorized(const struct CHObject_soa* obj1, const struct CHObject_soa* obj2) {
  // The search direction
  float d[3] = {1.f, 1.f, 1.f};  // could also be random
  // The point we examine in each iteration
  float a[3];

  // Our simplex that do_simplex modifies.
  struct Simplex3D s;
  s.num_points = 1;

  // S = Support(A-B) = Support(A) - Invsup(B)
  support3D_soa_vectorized(obj1->points, d, s.p[0], obj1->num_points);
  invsup3D_soa_vectorized(obj2->points, d, s.p[0], obj2->num_points);

  // d = -S
  d[0] = -s.p[0][0];
  d[1] = -s.p[0][1];
  d[2] = -s.p[0][2];

  int max_iter = 100;  // most likely we need less
  while (max_iter--) {
    // a = Support(A-B) = Support(A) - Invsup(B)
    support3D_soa_vectorized(obj1->points, d, a, obj1->num_points);
    invsup3D_soa_vectorized(obj2->points, d, a, obj2->num_points);

    // Are we past the origin? If no, we will never enclose it.
    if (dot3D(a, d) < 0.) return 0;

    // Add newly found point to Simplex
    // AT THIS POINT WE KNOW a IS PAST THE ORIGIN
    s.p[s.num_points][0] = a[0];
    s.p[s.num_points][1] = a[1];
    s.p[s.num_points][2] = a[2];
    s.num_points++;

    if (do_simplex3D_soa_vectorized(&s, d)) return 1;
  }

  return 1; // Most likely not reachable. Say we intersect by default.
}

int do_intersect3D_soa_vectorized_lu2(const struct CHObject_soa* obj1, const struct CHObject_soa* obj2) {
  // The search direction
  float d[3] = {1.f, 1.f, 1.f};  // could also be random
  // The point we examine in each iteration
  float a[3];

  // Our simplex that do_simplex modifies.
  struct Simplex3D s;
  s.num_points = 1;

  // S = Support(A-B) = Support(A) - Invsup(B)
  support3D_soa_vectorized_lu2(obj1->points, d, s.p[0], obj1->num_points);
  invsup3D_soa_vectorized_lu2(obj2->points, d, s.p[0], obj2->num_points);

  // d = -S
  d[0] = -s.p[0][0];
  d[1] = -s.p[0][1];
  d[2] = -s.p[0][2];

  int max_iter = 100;  // most likely we need less
  while (max_iter--) {
    // a = Support(A-B) = Support(A) - Invsup(B)
    support3D_soa_vectorized_lu2(obj1->points, d, a, obj1->num_points);
    invsup3D_soa_vectorized_lu2(obj2->points, d, a, obj2->num_points);

    // Are we past the origin? If no, we will never enclose it.
    if (dot3D(a, d) < 0.) return 0;

    // Add newly found point to Simplex
    // AT THIS POINT WE KNOW a IS PAST THE ORIGIN
    s.p[s.num_points][0] = a[0];
    s.p[s.num_points][1] = a[1];
    s.p[s.num_points][2] = a[2];
    s.num_points++;

    if (do_simplex3D_soa_vectorized(&s, d)) return 1;
  }

  return 1; // Most likely not reachable. Say we intersect by default.
}

[[gnu::target("avx512dq")]]
int do_intersect3D_soa_v_slow_idx_inlined(const struct CHObject_soa* obj1, const struct CHObject_soa* obj2) {
  // The search direction
  float d[3] = {1.f, 1.f, 1.f};  // could also be random
  // The point we examine in each iteration
  float a[3];

  // Our simplex that do_simplex modifies.
  struct Simplex3D s;
  s.num_points = 1;

  int no_points_min = obj1->num_points > obj2->num_points ? obj2->num_points : obj1->num_points;
  float dp, dpmax = -1e9;
  int argmax = 0;
  float inv_dp, inv_dpmax = 1e9;
  int argmin = 0;

  __m512 vertex_x_coord_0, vertex_y_coord_0, vertex_z_coord_0;
  __m512 inv_vertex_x_coord_0, inv_vertex_y_coord_0, inv_vertex_z_coord_0;

  __m512 tmp_0, tmp_1, tmp_2, dp_v;
  __m512i tmp_3, tmp_4;

  __m512 inv_tmp_0, inv_tmp_1, inv_tmp_2, inv_dp_v;
  __m512i inv_tmp_3, inv_tmp_4;

  __m512 dpmax_v = _mm512_set1_ps(-1e9);
  __m512 inv_dpmax_v = _mm512_set1_ps(1e9);

  __m512i argmax_v = _mm512_set1_epi32(0);
  __m512i inv_argmin_v = _mm512_set1_epi32(0);

  __m512i curr_idx;
  __mmask16 mask_0, mask_1;

  __m512 dir_x = _mm512_set1_ps(d[0]);
  __m512 dir_y = _mm512_set1_ps(d[1]);
  __m512 dir_z = _mm512_set1_ps(d[2]);
  __m512i ones = _mm512_set1_epi32(-1);
  __m512i zeros = _mm512_set1_epi32(0);
  __m512i sixteens = _mm512_set1_epi32(16);

  int i;
  // Joint loop
  for (i = 0; i < no_points_min - 15; i += 16) {
    // Compute current indices
    curr_idx = _mm512_set_epi32(i + 15, i + 14, i + 13, i + 12, i + 11, i + 10, i + 9, i + 8,
                                i + 7, i + 6, i + 5, i + 4, i + 3, i + 2, i + 1, i);

    // Load vertices
    vertex_x_coord_0 = _mm512_loadu_ps(&obj1->points[0][i]);
    vertex_y_coord_0 = _mm512_loadu_ps(&obj1->points[1][i]);
    vertex_z_coord_0 = _mm512_loadu_ps(&obj1->points[2][i]);

    inv_vertex_x_coord_0 = _mm512_loadu_ps(&obj2->points[0][i]);
    inv_vertex_y_coord_0 = _mm512_loadu_ps(&obj2->points[1][i]);
    inv_vertex_z_coord_0 = _mm512_loadu_ps(&obj2->points[2][i]);

    // Compute dot products
    tmp_0 = _mm512_mul_ps(vertex_x_coord_0, dir_x);
    tmp_1 = _mm512_fmadd_ps(dir_y, vertex_y_coord_0, tmp_0);
    dp_v = _mm512_fmadd_ps(dir_z, vertex_z_coord_0, tmp_1);

    inv_tmp_0 = _mm512_mul_ps(inv_vertex_x_coord_0, dir_x);
    inv_tmp_1 = _mm512_fmadd_ps(dir_y, inv_vertex_y_coord_0, inv_tmp_0);
    inv_dp_v = _mm512_fmadd_ps(dir_z, inv_vertex_z_coord_0, inv_tmp_1);

    // Find local dpmax/argmax
    mask_0 = _mm512_cmp_ps_mask(dp_v, dpmax_v, _CMP_GT_OQ);
    dpmax_v = _mm512_max_ps(dp_v, dpmax_v);
    tmp_3 = _mm512_mask_and_epi32(zeros, mask_0, curr_idx, ones);
    tmp_4 = _mm512_mask_andnot_epi32(argmax_v, mask_0, ones, argmax_v);
    argmax_v = _mm512_or_epi32(tmp_3, tmp_4);

    // Find local dpmin/argmin
    mask_1 = _mm512_cmp_ps_mask(inv_dp_v, inv_dpmax_v, _CMP_LT_OS);
    inv_dpmax_v = _mm512_min_ps(inv_dp_v, inv_dpmax_v);
    inv_tmp_3 = _mm512_mask_and_epi32(zeros, mask_1, curr_idx, ones);
    inv_tmp_4 = _mm512_mask_andnot_epi32(inv_argmin_v, mask_1, ones, inv_argmin_v);
    inv_argmin_v = _mm512_or_epi32(inv_tmp_3, inv_tmp_4);
  }

  int k = i;
  // Iterate over the remaining points of the first object
  for (; i < obj1->num_points - 15; i += 16) {
    // Compute current indices
    curr_idx = _mm512_set_epi32(i + 15, i + 14, i + 13, i + 12, i + 11, i + 10, i + 9, i + 8,
                                i + 7, i + 6, i + 5, i + 4, i + 3, i + 2, i + 1, i);

    // Load vertices
    vertex_x_coord_0 = _mm512_loadu_ps(&obj1->points[0][i]);
    vertex_y_coord_0 = _mm512_loadu_ps(&obj1->points[1][i]);
    vertex_z_coord_0 = _mm512_loadu_ps(&obj1->points[2][i]);

    // Compute dot products
    tmp_0 = _mm512_mul_ps(vertex_x_coord_0, dir_x);
    tmp_1 = _mm512_fmadd_ps(dir_y, vertex_y_coord_0, tmp_0);
    dp_v = _mm512_fmadd_ps(dir_z, vertex_z_coord_0, tmp_1);

    // Find local dpmax/argmax
    mask_0 = _mm512_cmp_ps_mask(dp_v, dpmax_v, _CMP_GT_OQ);
    dpmax_v = _mm512_max_ps(dp_v, dpmax_v);
    tmp_3 = _mm512_mask_and_epi32(zeros, mask_0, curr_idx, ones);
    tmp_4 = _mm512_mask_andnot_epi32(argmax_v, mask_0, ones, argmax_v);
    argmax_v = _mm512_or_epi32(tmp_3, tmp_4);
  }

  // Find dpmax/argmax
  if (i != 0) {
    float dps[16];
    int argmaxs[16];

    _mm512_storeu_ps(dps, dpmax_v);
    _mm512_storeu_epi32(argmaxs, argmax_v);
    for (int j = 0; j < 16; j++) {
      if (dps[j] > dpmax) {
        dpmax = dps[j];
        argmax = argmaxs[j];
      }
    }
  }

  // Handle remaining cases
  for (; i < obj1->num_points; i++) {
    dp = dot3D_soa_vectorized(obj1->points[0][i], obj1->points[1][i], obj1->points[2][i], d);
    if (dp > dpmax) {
      dpmax = dp;
      argmax = i;
    }
  }

  // Iterate over the remaining points of the second object
  for (; k < obj2->num_points - 15; k += 16) {
    // Compute current indices
    curr_idx = _mm512_set_epi32(k + 15, k + 14, k + 13, k + 12, k + 11, k + 10, k + 9, k + 8,
                                k + 7, k + 6, k + 5, k + 4, k + 3, k + 2, k + 1, k);

    // Load vertices
    inv_vertex_x_coord_0 = _mm512_loadu_ps(&obj2->points[0][k]);
    inv_vertex_y_coord_0 = _mm512_loadu_ps(&obj2->points[1][k]);
    inv_vertex_z_coord_0 = _mm512_loadu_ps(&obj2->points[2][k]);

    // Compute dot product
    inv_tmp_0 = _mm512_mul_ps(inv_vertex_x_coord_0, dir_x);
    inv_tmp_1 = _mm512_fmadd_ps(dir_y, inv_vertex_y_coord_0, inv_tmp_0);
    inv_dp_v = _mm512_fmadd_ps(dir_z, inv_vertex_z_coord_0, inv_tmp_1);

    // Find local dpmin/argmin
    mask_1 = _mm512_cmp_ps_mask(inv_dp_v, inv_dpmax_v, _CMP_LT_OS);
    inv_dpmax_v = _mm512_min_ps(inv_dp_v, inv_dpmax_v);
    inv_tmp_3 = _mm512_mask_and_epi32(zeros, mask_1, curr_idx, ones);
    inv_tmp_4 = _mm512_mask_andnot_epi32(inv_argmin_v, mask_1, ones, inv_argmin_v);
    inv_argmin_v = _mm512_or_epi32(inv_tmp_3, inv_tmp_4);
  }

  // Find dpmin/argmin
  if (k != 0) {
    float inv_dps[16];
    int argmins[16];

    _mm512_storeu_ps(inv_dps, inv_dpmax_v);
    _mm512_storeu_epi32(argmins, inv_argmin_v);

    for (int j = 0; j < 16; j++) {
      if (inv_dps[j] < inv_dpmax) {
        inv_dpmax = inv_dps[j];
        argmin = argmins[j];
      }
    }
  }

  // Handle remaining cases
  for (; k < obj2->num_points; k++) {
    inv_dp = dot3D_soa_vectorized(obj2->points[0][k], obj2->points[1][k], obj2->points[2][k], d);
    if (inv_dp < inv_dpmax) {
      inv_dpmax = inv_dp;
      argmin = k;
    }
  }

  // S = Support(A-B) = Support(A) - Invsup(B)
  s.p[0][0] = obj1->points[0][argmax] - obj2->points[0][argmin];
  s.p[0][1] = obj1->points[1][argmax] - obj2->points[1][argmin];
  s.p[0][2] = obj1->points[2][argmax] - obj2->points[2][argmin];

  // d = -S
  d[0] = -s.p[0][0];
  d[1] = -s.p[0][1];
  d[2] = -s.p[0][2];

  int max_iter = 100;  // most likely we need less
  while (max_iter--) {
    // a = Support(A-B) = Support(A) - Invsup(B)
    dpmax = -1e9;
    argmax = 0;
    inv_dpmax = 1e9;
    argmin = 0;

    dpmax_v = _mm512_set1_ps(-1e9);
    inv_dpmax_v = _mm512_set1_ps(1e9);

    argmax_v = _mm512_set1_epi32(0);
    inv_argmin_v = _mm512_set1_epi32(0);

    dir_x = _mm512_set1_ps(d[0]);
    dir_y = _mm512_set1_ps(d[1]);
    dir_z = _mm512_set1_ps(d[2]);

    // Joint loop
    for (i = 0; i < no_points_min - 15; i += 16) {
      // Compute current indices
      curr_idx = _mm512_set_epi32(i + 15, i + 14, i + 13, i + 12, i + 11, i + 10, i + 9, i + 8,
                                  i + 7, i + 6, i + 5, i + 4, i + 3, i + 2, i + 1, i);

      // Load vertices
      vertex_x_coord_0 = _mm512_loadu_ps(&obj1->points[0][i]);
      vertex_y_coord_0 = _mm512_loadu_ps(&obj1->points[1][i]);
      vertex_z_coord_0 = _mm512_loadu_ps(&obj1->points[2][i]);

      inv_vertex_x_coord_0 = _mm512_loadu_ps(&obj2->points[0][i]);
      inv_vertex_y_coord_0 = _mm512_loadu_ps(&obj2->points[1][i]);
      inv_vertex_z_coord_0 = _mm512_loadu_ps(&obj2->points[2][i]);

      // Compute dot products
      tmp_0 = _mm512_mul_ps(vertex_x_coord_0, dir_x);
      tmp_1 = _mm512_fmadd_ps(dir_y, vertex_y_coord_0, tmp_0);
      dp_v = _mm512_fmadd_ps(dir_z, vertex_z_coord_0, tmp_1);

      inv_tmp_0 = _mm512_mul_ps(inv_vertex_x_coord_0, dir_x);
      inv_tmp_1 = _mm512_fmadd_ps(dir_y, inv_vertex_y_coord_0, inv_tmp_0);
      inv_dp_v = _mm512_fmadd_ps(dir_z, inv_vertex_z_coord_0, inv_tmp_1);

      // Find local dpmax/argmax
      mask_0 = _mm512_cmp_ps_mask(dp_v, dpmax_v, _CMP_GT_OQ);
      dpmax_v = _mm512_max_ps(dp_v, dpmax_v);
      tmp_3 = _mm512_mask_and_epi32(zeros, mask_0, curr_idx, ones);
      tmp_4 = _mm512_mask_andnot_epi32(argmax_v, mask_0, ones, argmax_v);
      argmax_v = _mm512_or_epi32(tmp_3, tmp_4);

      // Find local dpmin/argmin
      mask_1 = _mm512_cmp_ps_mask(inv_dp_v, inv_dpmax_v, _CMP_LT_OS);
      inv_dpmax_v = _mm512_min_ps(inv_dp_v, inv_dpmax_v);
      inv_tmp_3 = _mm512_mask_and_epi32(zeros, mask_1, curr_idx, ones);
      inv_tmp_4 = _mm512_mask_andnot_epi32(inv_argmin_v, mask_1, ones, inv_argmin_v);
      inv_argmin_v = _mm512_or_epi32(inv_tmp_3, inv_tmp_4);
    }

    k = i;
    // Iterate over the remaining points of the first object
    for (; i < obj1->num_points - 15; i += 16) {
      // Compute current indices
      curr_idx = _mm512_set_epi32(i + 15, i + 14, i + 13, i + 12, i + 11, i + 10, i + 9, i + 8,
                                  i + 7, i + 6, i + 5, i + 4, i + 3, i + 2, i + 1, i);

      // Load vertices
      vertex_x_coord_0 = _mm512_loadu_ps(&obj1->points[0][i]);
      vertex_y_coord_0 = _mm512_loadu_ps(&obj1->points[1][i]);
      vertex_z_coord_0 = _mm512_loadu_ps(&obj1->points[2][i]);

      // Compute dot products
      tmp_0 = _mm512_mul_ps(vertex_x_coord_0, dir_x);
      tmp_1 = _mm512_fmadd_ps(dir_y, vertex_y_coord_0, tmp_0);
      dp_v = _mm512_fmadd_ps(dir_z, vertex_z_coord_0, tmp_1);

      // Find local dpmax/argmax
      mask_0 = _mm512_cmp_ps_mask(dp_v, dpmax_v, _CMP_GT_OQ);
      dpmax_v = _mm512_max_ps(dp_v, dpmax_v);
      tmp_3 = _mm512_mask_and_epi32(zeros, mask_0, curr_idx, ones);
      tmp_4 = _mm512_mask_andnot_epi32(argmax_v, mask_0, ones, argmax_v);
      argmax_v = _mm512_or_epi32(tmp_3, tmp_4);
    }

    // Find dpmax/argmax
    if (i != 0) {
      float dps[16];
      int argmaxs[16];

      _mm512_storeu_ps(dps, dpmax_v);
      _mm512_storeu_epi32(argmaxs, argmax_v);
      for (int j = 0; j < 16; j++) {
        if (dps[j] > dpmax) {
          dpmax = dps[j];
          argmax = argmaxs[j];
        }
      }
    }

    // Handle remaining cases
    for (; i < obj1->num_points; i++) {
      dp = dot3D_soa_vectorized(obj1->points[0][i], obj1->points[1][i], obj1->points[2][i], d);
      if (dp > dpmax) {
        dpmax = dp;
        argmax = i;
      }
    }

    // Iterate over the remaining points of the second object
    for (; k < obj2->num_points - 15; k += 16) {
      // Compute the current indices
      curr_idx = _mm512_set_epi32(k + 15, k + 14, k + 13, k + 12, k + 11, k + 10, k + 9, k + 8,
                                  k + 7, k + 6, k + 5, k + 4, k + 3, k + 2, k + 1, k);

      // Load the vertices
      inv_vertex_x_coord_0 = _mm512_loadu_ps(&obj2->points[0][k]);
      inv_vertex_y_coord_0 = _mm512_loadu_ps(&obj2->points[1][k]);
      inv_vertex_z_coord_0 = _mm512_loadu_ps(&obj2->points[2][k]);

      // Compute the dot product
      inv_tmp_0 = _mm512_mul_ps(inv_vertex_x_coord_0, dir_x);
      inv_tmp_1 = _mm512_fmadd_ps(dir_y, inv_vertex_y_coord_0, inv_tmp_0);
      inv_dp_v = _mm512_fmadd_ps(dir_z, inv_vertex_z_coord_0, inv_tmp_1);

      // Find local dpmin/argmin
      mask_1 = _mm512_cmp_ps_mask(inv_dp_v, inv_dpmax_v, _CMP_LT_OS);
      inv_dpmax_v = _mm512_min_ps(inv_dp_v, inv_dpmax_v);
      inv_tmp_3 = _mm512_mask_and_epi32(zeros, mask_1, curr_idx, ones);
      inv_tmp_4 = _mm512_mask_andnot_epi32(inv_argmin_v, mask_1, ones, inv_argmin_v);
      inv_argmin_v = _mm512_or_epi32(inv_tmp_3, inv_tmp_4);
    }

    // Find dpmin/argmin 
    if (k != 0) {
      float inv_dps[16];
      int argmins[16];

      _mm512_storeu_ps(inv_dps, inv_dpmax_v);
      _mm512_storeu_epi32(argmins, inv_argmin_v);

      for (int j = 0; j < 16; j++) {
        if (inv_dps[j] < inv_dpmax) {
          inv_dpmax = inv_dps[j];
          argmin = argmins[j];
        }
      }
    }

    // Handle remaining cases
    for (; k < obj2->num_points; k++) {
      inv_dp = dot3D_soa_vectorized(obj2->points[0][k], obj2->points[1][k], obj2->points[2][k], d);
      if (inv_dp < inv_dpmax) {
        inv_dpmax = inv_dp;
        argmin = k;
      }
    }
    
    // a = Support(A-B) = Support(A) - Invsup(B)
    a[0] = obj1->points[0][argmax] - obj2->points[0][argmin];
    a[1] = obj1->points[1][argmax] - obj2->points[1][argmin];
    a[2] = obj1->points[2][argmax] - obj2->points[2][argmin];

    // Are we past the origin? If no, we will never enclose it.
    if (dot3D(a, d) < 0.) return 0;

    // Add newly found point to Simplex
    // AT THIS POINT WE KNOW a IS PAST THE ORIGIN
    s.p[s.num_points][0] = a[0];
    s.p[s.num_points][1] = a[1];
    s.p[s.num_points][2] = a[2];
    s.num_points++;

    if (do_simplex3D_soa_vectorized(&s, d)) return 1;
  }

  return 1; // Most likely not reachable. Say we intersect by default.
}

[[gnu::target("avx512dq")]]
int do_intersect3D_soa_vectorized_inlined(const struct CHObject_soa* obj1, const struct CHObject_soa* obj2) {
  // The search direction
  float d[3] = {1.f, 1.f, 1.f};  // could also be random
  // The point we examine in each iteration
  float a[3];

  // Our simplex that do_simplex modifies.
  struct Simplex3D s;
  s.num_points = 1;

  int no_points_min = obj1->num_points > obj2->num_points ? obj2->num_points : obj1->num_points;
  float dp, dpmax = -1e9;
  int argmax = 0;
  float inv_dp, inv_dpmax = 1e9;
  int argmin = 0;

  __m512 vertex_x_coord_0, vertex_y_coord_0, vertex_z_coord_0;
  __m512 inv_vertex_x_coord_0, inv_vertex_y_coord_0, inv_vertex_z_coord_0;

  __m512 tmp_0, tmp_1, tmp_2, dp_v;
  __m512i tmp_3, tmp_4;

  __m512 inv_tmp_0, inv_tmp_1, inv_tmp_2, inv_dp_v;
  __m512i inv_tmp_3, inv_tmp_4;

  __m512 dpmax_v = _mm512_set1_ps(-1e9);
  __m512 inv_dpmax_v = _mm512_set1_ps(1e9);

  __m512i argmax_v = _mm512_set1_epi32(0);
  __m512i inv_argmin_v = _mm512_set1_epi32(0);

  __m512i curr_idx = _mm512_set_epi32(-1, -2, -3, -4, -5, -6, -7, -8,
                                      -9, -10, -11, -12, -13, -14, -15, -16);
  __mmask16 mask_0, mask_1;

  __m512 dir_x = _mm512_set1_ps(d[0]);
  __m512 dir_y = _mm512_set1_ps(d[1]);
  __m512 dir_z = _mm512_set1_ps(d[2]);
  __m512i ones = _mm512_set1_epi32(-1);
  __m512i zeros = _mm512_set1_epi32(0);
  __m512i sixteens = _mm512_set1_epi32(16);

  int i;
  // Joint loop
  for (i = 0; i < no_points_min - 15; i += 16) {
    // Compute the current indices
    curr_idx = _mm512_add_epi32(curr_idx, sixteens);

    // Load the vertices
    vertex_x_coord_0 = _mm512_loadu_ps(&obj1->points[0][i]);
    vertex_y_coord_0 = _mm512_loadu_ps(&obj1->points[1][i]);
    vertex_z_coord_0 = _mm512_loadu_ps(&obj1->points[2][i]);

    inv_vertex_x_coord_0 = _mm512_loadu_ps(&obj2->points[0][i]);
    inv_vertex_y_coord_0 = _mm512_loadu_ps(&obj2->points[1][i]);
    inv_vertex_z_coord_0 = _mm512_loadu_ps(&obj2->points[2][i]);

    // Compute dot products
    tmp_0 = _mm512_mul_ps(vertex_x_coord_0, dir_x);
    tmp_1 = _mm512_fmadd_ps(dir_y, vertex_y_coord_0, tmp_0);
    dp_v = _mm512_fmadd_ps(dir_z, vertex_z_coord_0, tmp_1);

    inv_tmp_0 = _mm512_mul_ps(inv_vertex_x_coord_0, dir_x);
    inv_tmp_1 = _mm512_fmadd_ps(dir_y, inv_vertex_y_coord_0, inv_tmp_0);
    inv_dp_v = _mm512_fmadd_ps(dir_z, inv_vertex_z_coord_0, inv_tmp_1);

    // Find local dpmax/argmax
    mask_0 = _mm512_cmp_ps_mask(dp_v, dpmax_v, _CMP_GT_OQ);
    dpmax_v = _mm512_max_ps(dp_v, dpmax_v);
    tmp_3 = _mm512_mask_and_epi32(zeros, mask_0, curr_idx, ones);
    tmp_4 = _mm512_mask_andnot_epi32(argmax_v, mask_0, ones, argmax_v);
    argmax_v = _mm512_or_epi32(tmp_3, tmp_4);

    // Find local dpmin/argmin
    mask_1 = _mm512_cmp_ps_mask(inv_dp_v, inv_dpmax_v, _CMP_LT_OS);
    inv_dpmax_v = _mm512_min_ps(inv_dp_v, inv_dpmax_v);
    inv_tmp_3 = _mm512_mask_and_epi32(zeros, mask_1, curr_idx, ones);
    inv_tmp_4 = _mm512_mask_andnot_epi32(inv_argmin_v, mask_1, ones, inv_argmin_v);
    inv_argmin_v = _mm512_or_epi32(inv_tmp_3, inv_tmp_4);
  }

  int k = i;
  // Iterate over the remaining points of the first object
  for (; i < obj1->num_points - 15; i += 16) {
    // Compute the current indices
    curr_idx = _mm512_add_epi32(curr_idx, sixteens);

    // Load the vertices
    vertex_x_coord_0 = _mm512_loadu_ps(&obj1->points[0][i]);
    vertex_y_coord_0 = _mm512_loadu_ps(&obj1->points[1][i]);
    vertex_z_coord_0 = _mm512_loadu_ps(&obj1->points[2][i]);

    // Compute the dot products
    tmp_0 = _mm512_mul_ps(vertex_x_coord_0, dir_x);
    tmp_1 = _mm512_fmadd_ps(dir_y, vertex_y_coord_0, tmp_0);
    dp_v = _mm512_fmadd_ps(dir_z, vertex_z_coord_0, tmp_1);

    // Find the local dpmax/argmax
    mask_0 = _mm512_cmp_ps_mask(dp_v, dpmax_v, _CMP_GT_OQ);
    dpmax_v = _mm512_max_ps(dp_v, dpmax_v);
    tmp_3 = _mm512_mask_and_epi32(zeros, mask_0, curr_idx, ones);
    tmp_4 = _mm512_mask_andnot_epi32(argmax_v, mask_0, ones, argmax_v);
    argmax_v = _mm512_or_epi32(tmp_3, tmp_4);
  }

  // Find the dpmax/argmax
  if (i != 0) {
    float dps[16];
    int argmaxs[16];

    _mm512_storeu_ps(dps, dpmax_v);
    _mm512_storeu_epi32(argmaxs, argmax_v);
    for (int j = 0; j < 16; j++) {
      if (dps[j] > dpmax) {
        dpmax = dps[j];
        argmax = argmaxs[j];
      }
    }
  }

  // Handle remaining cases
  for (; i < obj1->num_points; i++) {
    dp = dot3D_soa_vectorized(obj1->points[0][i], obj1->points[1][i], obj1->points[2][i], d);
    if (dp > dpmax) {
      dpmax = dp;
      argmax = i;
    }
  }

  // Iterate over the remaining points of the second object
  for (; k < obj2->num_points - 15; k += 16) {
    // Compute the current indices
    curr_idx = _mm512_add_epi32(curr_idx, sixteens);

    // Load the vertices
    inv_vertex_x_coord_0 = _mm512_loadu_ps(&obj2->points[0][k]);
    inv_vertex_y_coord_0 = _mm512_loadu_ps(&obj2->points[1][k]);
    inv_vertex_z_coord_0 = _mm512_loadu_ps(&obj2->points[2][k]);

    // Compute the dot products
    inv_tmp_0 = _mm512_mul_ps(inv_vertex_x_coord_0, dir_x);
    inv_tmp_1 = _mm512_fmadd_ps(dir_y, inv_vertex_y_coord_0, inv_tmp_0);
    inv_dp_v = _mm512_fmadd_ps(dir_z, inv_vertex_z_coord_0, inv_tmp_1);

    // Find the local dpmin/argmin
    mask_1 = _mm512_cmp_ps_mask(inv_dp_v, inv_dpmax_v, _CMP_LT_OS);
    inv_dpmax_v = _mm512_min_ps(inv_dp_v, inv_dpmax_v);
    inv_tmp_3 = _mm512_mask_and_epi32(zeros, mask_1, curr_idx, ones);
    inv_tmp_4 = _mm512_mask_andnot_epi32(inv_argmin_v, mask_1, ones, inv_argmin_v);
    inv_argmin_v = _mm512_or_epi32(inv_tmp_3, inv_tmp_4);
  }

  // Find dpmin/argmin
  if (k != 0) {
    float inv_dps[16];
    int argmins[16];

    _mm512_storeu_ps(inv_dps, inv_dpmax_v);
    _mm512_storeu_epi32(argmins, inv_argmin_v);

    for (int j = 0; j < 16; j++) {
      if (inv_dps[j] < inv_dpmax) {
        inv_dpmax = inv_dps[j];
        argmin = argmins[j];
      }
    }
  }

  // Handle remaining cases
  for (; k < obj2->num_points; k++) {
    inv_dp = dot3D_soa_vectorized(obj2->points[0][k], obj2->points[1][k], obj2->points[2][k], d);
    if (inv_dp < inv_dpmax) {
      inv_dpmax = inv_dp;
      argmin = k;
    }
  }

  // S = Support(A-B) = Support(A) - Invsup(B)
  s.p[0][0] = obj1->points[0][argmax] - obj2->points[0][argmin];
  s.p[0][1] = obj1->points[1][argmax] - obj2->points[1][argmin];
  s.p[0][2] = obj1->points[2][argmax] - obj2->points[2][argmin];

  // d = -S
  d[0] = -s.p[0][0];
  d[1] = -s.p[0][1];
  d[2] = -s.p[0][2];

  int max_iter = 100;  // most likely we need less
  while (max_iter--) {
    dpmax = -1e9;
    argmax = 0;
    inv_dpmax = 1e9;
    argmin = 0;

    dpmax_v = _mm512_set1_ps(-1e9);
    inv_dpmax_v = _mm512_set1_ps(1e9);

    argmax_v = _mm512_set1_epi32(0);
    inv_argmin_v = _mm512_set1_epi32(0);

    dir_x = _mm512_set1_ps(d[0]);
    dir_y = _mm512_set1_ps(d[1]);
    dir_z = _mm512_set1_ps(d[2]);

    curr_idx = _mm512_set_epi32(-1, -2, -3, -4, -5, -6, -7, -8,
                                -9, -10, -11, -12, -13, -14, -15, -16);

    // Joint loop
    for (i = 0; i < no_points_min - 15; i += 16) {
      // Compute current indices
      curr_idx = _mm512_add_epi32(curr_idx, sixteens);
  
      // Load the vertices
      vertex_x_coord_0 = _mm512_loadu_ps(&obj1->points[0][i]);
      vertex_y_coord_0 = _mm512_loadu_ps(&obj1->points[1][i]);
      vertex_z_coord_0 = _mm512_loadu_ps(&obj1->points[2][i]);

      inv_vertex_x_coord_0 = _mm512_loadu_ps(&obj2->points[0][i]);
      inv_vertex_y_coord_0 = _mm512_loadu_ps(&obj2->points[1][i]);
      inv_vertex_z_coord_0 = _mm512_loadu_ps(&obj2->points[2][i]);

      // Compute the dot products
      tmp_0 = _mm512_mul_ps(vertex_x_coord_0, dir_x);
      tmp_1 = _mm512_fmadd_ps(dir_y, vertex_y_coord_0, tmp_0);
      dp_v = _mm512_fmadd_ps(dir_z, vertex_z_coord_0, tmp_1);

      inv_tmp_0 = _mm512_mul_ps(inv_vertex_x_coord_0, dir_x);
      inv_tmp_1 = _mm512_fmadd_ps(dir_y, inv_vertex_y_coord_0, inv_tmp_0);
      inv_dp_v = _mm512_fmadd_ps(dir_z, inv_vertex_z_coord_0, inv_tmp_1);

      // Find the local dpmax/argmax
      mask_0 = _mm512_cmp_ps_mask(dp_v, dpmax_v, _CMP_GT_OQ);
      dpmax_v = _mm512_max_ps(dp_v, dpmax_v);
      tmp_3 = _mm512_mask_and_epi32(zeros, mask_0, curr_idx, ones);
      tmp_4 = _mm512_mask_andnot_epi32(argmax_v, mask_0, ones, argmax_v);
      argmax_v = _mm512_or_epi32(tmp_3, tmp_4);

      // Find the local dpmin/argmin
      mask_1 = _mm512_cmp_ps_mask(inv_dp_v, inv_dpmax_v, _CMP_LT_OS);
      inv_dpmax_v = _mm512_min_ps(inv_dp_v, inv_dpmax_v);
      inv_tmp_3 = _mm512_mask_and_epi32(zeros, mask_1, curr_idx, ones);
      inv_tmp_4 = _mm512_mask_andnot_epi32(inv_argmin_v, mask_1, ones, inv_argmin_v);
      inv_argmin_v = _mm512_or_epi32(inv_tmp_3, inv_tmp_4);
    }

    k = i;
    // Iterate over the remaining points of the first object
    for (; i < obj1->num_points - 15; i += 16) {
      // Compute the current indices
      curr_idx = _mm512_add_epi32(curr_idx, sixteens);

      // Load the vertices
      vertex_x_coord_0 = _mm512_loadu_ps(&obj1->points[0][i]);
      vertex_y_coord_0 = _mm512_loadu_ps(&obj1->points[1][i]);
      vertex_z_coord_0 = _mm512_loadu_ps(&obj1->points[2][i]);

      // Compute the dot products
      tmp_0 = _mm512_mul_ps(vertex_x_coord_0, dir_x);
      tmp_1 = _mm512_fmadd_ps(dir_y, vertex_y_coord_0, tmp_0);
      dp_v = _mm512_fmadd_ps(dir_z, vertex_z_coord_0, tmp_1);

      // Find local dpmax/argmax
      mask_0 = _mm512_cmp_ps_mask(dp_v, dpmax_v, _CMP_GT_OQ);
      dpmax_v = _mm512_max_ps(dp_v, dpmax_v);
      tmp_3 = _mm512_mask_and_epi32(zeros, mask_0, curr_idx, ones);
      tmp_4 = _mm512_mask_andnot_epi32(argmax_v, mask_0, ones, argmax_v);
      argmax_v = _mm512_or_epi32(tmp_3, tmp_4);
    }

    // Find dpmax/argmax
    if (i != 0) {
      float dps[16];
      int argmaxs[16];

      _mm512_storeu_ps(dps, dpmax_v);
      _mm512_storeu_epi32(argmaxs, argmax_v);
      for (int j = 0; j < 16; j++) {
        if (dps[j] > dpmax) {
          dpmax = dps[j];
          argmax = argmaxs[j];
        }
      }
    }

    // Handle remaining cases
    for (; i < obj1->num_points; i++) {
      dp = dot3D_soa_vectorized(obj1->points[0][i], obj1->points[1][i], obj1->points[2][i], d);
      if (dp > dpmax) {
        dpmax = dp;
        argmax = i;
      }
    }

    // Iterate over the remaining points of the second object
    for (; k < obj2->num_points - 15; k += 16) {
      // Compute the current indices
      curr_idx = _mm512_add_epi32(curr_idx, sixteens);

      // Load the vertices
      inv_vertex_x_coord_0 = _mm512_loadu_ps(&obj2->points[0][k]);
      inv_vertex_y_coord_0 = _mm512_loadu_ps(&obj2->points[1][k]);
      inv_vertex_z_coord_0 = _mm512_loadu_ps(&obj2->points[2][k]);

      // Compute the dot products
      inv_tmp_0 = _mm512_mul_ps(inv_vertex_x_coord_0, dir_x);
      inv_tmp_1 = _mm512_fmadd_ps(dir_y, inv_vertex_y_coord_0, inv_tmp_0);
      inv_dp_v = _mm512_fmadd_ps(dir_z, inv_vertex_z_coord_0, inv_tmp_1);

      // Find local dpmin/argmin
      mask_1 = _mm512_cmp_ps_mask(inv_dp_v, inv_dpmax_v, _CMP_LT_OS);
      inv_dpmax_v = _mm512_min_ps(inv_dp_v, inv_dpmax_v);
      inv_tmp_3 = _mm512_mask_and_epi32(zeros, mask_1, curr_idx, ones);
      inv_tmp_4 = _mm512_mask_andnot_epi32(inv_argmin_v, mask_1, ones, inv_argmin_v);
      inv_argmin_v = _mm512_or_epi32(inv_tmp_3, inv_tmp_4);
    }

    // Find dpmin/argmin
    if (k != 0) {
      float inv_dps[16];
      int argmins[16];

      _mm512_storeu_ps(inv_dps, inv_dpmax_v);
      _mm512_storeu_epi32(argmins, inv_argmin_v);

      for (int j = 0; j < 16; j++) {
        if (inv_dps[j] < inv_dpmax) {
          inv_dpmax = inv_dps[j];
          argmin = argmins[j];
        }
      }
    }

    // Handle remaining cases
    for (; k < obj2->num_points; k++) {
      inv_dp = dot3D_soa_vectorized(obj2->points[0][k], obj2->points[1][k], obj2->points[2][k], d);
      if (inv_dp < inv_dpmax) {
        inv_dpmax = inv_dp;
        argmin = k;
      }
    }

    // a = Support(A-B) = Support(A) - Invsup(B)
    a[0] = obj1->points[0][argmax] - obj2->points[0][argmin];
    a[1] = obj1->points[1][argmax] - obj2->points[1][argmin];
    a[2] = obj1->points[2][argmax] - obj2->points[2][argmin];

    // Are we past the origin? If no, we will never enclose it.
    if (dot3D(a, d) < 0.) return 0;

    // Add newly found point to Simplex
    // AT THIS POINT WE KNOW a IS PAST THE ORIGIN
    s.p[s.num_points][0] = a[0];
    s.p[s.num_points][1] = a[1];
    s.p[s.num_points][2] = a[2];
    s.num_points++;

    if (do_simplex3D_soa_vectorized(&s, d)) return 1;
  }

  return 1; // Most likely not reachable. Say we intersect by default.
}
