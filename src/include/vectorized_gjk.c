#include "gjk.h"
#include <immintrin.h>
#include <stdio.h>

/**
 * Dot product of two 3D vectors.
 * */
float dot3D_vectorized(const float *x, const float *y) {
  return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

/**
 * Cross product of two 3D vectors. Write into res.
 * */
void cross3D_vectorized(const float *x, const float *y, float *res) {
  res[0] = x[1] * y[2] - x[2] * y[1];
  res[1] = x[2] * y[0] - x[0] * y[2];
  res[2] = x[0] * y[1] - x[1] * y[0];
}

/**
 * Among a set of N points in 3D, find the point furthest along dir.
 * TODO can we do this more efficiently than brute force?
 * */
[[gnu::target("avx512vl")]]
void support3D_vectorized(const float (* points)[3], const float *dir, float *res, int N) {
  __m128 vertex_0, vertex_1, vertex_2, vertex_3, vertex_4, vertex_5, vertex_6, vertex_7;
  __m128 vertex_8, vertex_9, vertex_10, vertex_11, vertex_12, vertex_13, vertex_14, vertex_15;

  __m256 vertex_0_1, vertex_2_3, vertex_4_5, vertex_6_7;
  __m256 vertex_8_9, vertex_10_11, vertex_12_13, vertex_14_15;

  __m256 dp_01, dp_23, dp_45, dp_67, dp0_7;
  __m256 dp_89, dp_1011, dp_1213, dp_1415, dp8_15;

  __m256 dpmax_v = _mm256_set1_ps(-1e9);
  __m256 dpmax_v_1 = _mm256_set1_ps(-1e9);

  __m256i argmax_v = _mm256_set1_epi32(0);
  __m256i argmax_v_1 = _mm256_set1_epi32(0);

  __m256i curr_ind_v = _mm256_set_epi32(-9, -11, -13, -15, -10, -12, -14, -16);
  __m256i curr_ind_v_1 = _mm256_set_epi32(-1, -3, -5, -7, -2, -4, -6, -8);
  __m256 dir_v = _mm256_set_ps(0.0f, dir[2], dir[1], dir[0], 0.0f, dir[2], dir[1], dir[0]);

  __m256 tmp1, tmp2, mask;
  __m256 tmp5, tmp6, mask1;

  __m256i tmp3, tmp4;
  __m256i tmp7, tmp8;

  __m256i sixteens = _mm256_set1_epi32(16);

  float dp, dpmax = -1e9;
  int argmax = 0;
  int i = 0;

  for (i = 0; i < N - 15; i += 16) {
    curr_ind_v = _mm256_add_epi32(curr_ind_v, sixteens);
    curr_ind_v_1 = _mm256_add_epi32(curr_ind_v_1, sixteens);

    vertex_0 = _mm_loadu_ps(points[i]);
    vertex_1 = _mm_loadu_ps(points[i + 1]);
    vertex_2 = _mm_loadu_ps(points[i + 2]);
    vertex_3 = _mm_loadu_ps(points[i + 3]);
    vertex_4 = _mm_loadu_ps(points[i + 4]);
    vertex_5 = _mm_loadu_ps(points[i + 5]);
    vertex_6 = _mm_loadu_ps(points[i + 6]);
    vertex_7 = _mm_loadu_ps(points[i + 7]);

    vertex_8 = _mm_loadu_ps(points[i + 8]);
    vertex_9 = _mm_loadu_ps(points[i + 9]);
    vertex_10 = _mm_loadu_ps(points[i + 10]);
    vertex_11 = _mm_loadu_ps(points[i + 11]);
    vertex_12 = _mm_loadu_ps(points[i + 12]);
    vertex_13 = _mm_loadu_ps(points[i + 13]);
    vertex_14 = _mm_loadu_ps(points[i + 14]);
    vertex_15 = _mm_set_ps(0.0f, points[i + 15][2], points[i + 15][1], points[i + 15][0]);

    vertex_0_1 = _mm256_set_m128(vertex_1, vertex_0);
    vertex_2_3 = _mm256_set_m128(vertex_3, vertex_2);
    vertex_4_5 = _mm256_set_m128(vertex_5, vertex_4);
    vertex_6_7 = _mm256_set_m128(vertex_7, vertex_6);

    vertex_8_9 = _mm256_set_m128(vertex_9, vertex_8);
    vertex_10_11 = _mm256_set_m128(vertex_11, vertex_10);
    vertex_12_13 = _mm256_set_m128(vertex_13, vertex_12);
    vertex_14_15 = _mm256_set_m128(vertex_15, vertex_14);

    dp_01 = _mm256_dp_ps(vertex_0_1, dir_v, 113);
    dp_23 = _mm256_dp_ps(vertex_2_3, dir_v, 114);
    dp_45 = _mm256_dp_ps(vertex_4_5, dir_v, 116);
    dp_67 = _mm256_dp_ps(vertex_6_7, dir_v, 120);

    dp_89 = _mm256_dp_ps(vertex_8_9, dir_v, 113);
    dp_1011 = _mm256_dp_ps(vertex_10_11, dir_v, 114);
    dp_1213 = _mm256_dp_ps(vertex_12_13, dir_v, 116);
    dp_1415 = _mm256_dp_ps(vertex_14_15, dir_v, 120);

    tmp1 = _mm256_or_ps(dp_01, dp_23);
    tmp2 = _mm256_or_ps(dp_45, dp_67);
    dp0_7 = _mm256_or_ps(tmp1, tmp2);

    tmp5 = _mm256_or_ps(dp_89, dp_1011);
    tmp6 = _mm256_or_ps(dp_1213, dp_1415);
    dp8_15 = _mm256_or_ps(tmp5, tmp6);

    mask = _mm256_cmp_ps(dp0_7, dpmax_v, _CMP_GT_OS);
    dpmax_v = _mm256_max_ps(dp0_7, dpmax_v);
    tmp3 = _mm256_castps_si256(_mm256_and_ps(mask, _mm256_castsi256_ps(curr_ind_v)));
    tmp4 = _mm256_castps_si256(_mm256_andnot_ps(mask, _mm256_castsi256_ps(argmax_v)));
    argmax_v = _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(tmp3), _mm256_castsi256_ps(tmp4)));

    mask1 = _mm256_cmp_ps(dp8_15, dpmax_v_1, _CMP_GT_OS);
    dpmax_v_1 = _mm256_max_ps(dp8_15, dpmax_v_1);
    tmp7 = _mm256_castps_si256(_mm256_and_ps(mask1, _mm256_castsi256_ps(curr_ind_v_1)));
    tmp8 = _mm256_castps_si256(_mm256_andnot_ps(mask1, _mm256_castsi256_ps(argmax_v_1)));
    argmax_v_1 = _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(tmp7), _mm256_castsi256_ps(tmp8)));
  }

  mask = _mm256_cmp_ps(dpmax_v, dpmax_v_1, _CMP_GT_OS);
  dpmax_v = _mm256_max_ps(dpmax_v, dpmax_v_1);
  tmp3 = _mm256_castps_si256(_mm256_and_ps(mask, _mm256_castsi256_ps(argmax_v)));
  tmp4 = _mm256_castps_si256(_mm256_andnot_ps(mask, _mm256_castsi256_ps(argmax_v_1)));
  argmax_v = _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(tmp3), _mm256_castsi256_ps(tmp4)));

  if (i != 0) {
    float dps[8];
    int argmaxs[8];
    _mm256_storeu_ps(dps, dpmax_v);
    _mm256_storeu_epi32(argmaxs, argmax_v);

    for (int j = 0; j < 8; j++) {
      if (dps[j] > dpmax) {
        dpmax = dps[j];
        argmax = argmaxs[j];
      }
    }
  }

  // Handle the remaining cases
  for (; i < N; i++) {
    dp = dot3D_vectorized(points[i], dir);
    if (dp > dpmax) {
      dpmax = dp;
      argmax = i;
    }
  }

  res[0] = points[argmax][0];
  res[1] = points[argmax][1];
  res[2] = points[argmax][2];
}

/**
 * "Inverse support"
 * SUBTRACTS THE FOUND POINT FROM res.
 *
 * Among a set of N points in 3D, find the point furthest along NEGATIVE dir.
 * argmax_b (-a)Tb = argmin_b aTb
 * This is so we don't need to invert/negate dir
 * */
[[gnu::target("avx512vl")]]
void invsup3D_vectorized(const float (* points)[3], const float *dir, float *res, int N) {
  __m128 vertex_0, vertex_1, vertex_2, vertex_3, vertex_4, vertex_5, vertex_6, vertex_7;
  __m256 vertex_0_1, vertex_2_3, vertex_4_5, vertex_6_7;
  __m256 dp_01, dp_23, dp_45, dp_67, dp0_7;

  __m256 dpmax_v = _mm256_set1_ps(1e9);
  __m256i argmin_v = _mm256_set1_epi32(0);
  __m256i curr_ind_v = _mm256_set_epi32(-1, -3, -5, -7, -2, -4, -6, -8);

  __m256 dir_v = _mm256_set_ps(0.0f, dir[2], dir[1], dir[0], 0.0f, dir[2], dir[1], dir[0]);

  __m256 tmp1, tmp2, mask;
  __m256 tmp5, tmp6, mask1;
  __m256i tmp3, tmp4;
  __m256i tmp7, tmp8;

  __m256i eights = _mm256_set1_epi32(8);

  float dp, dpmax = 1e9;
  int argmin = 0;
  int i;
  for (i = 0; i < N - 7; i += 8) {
    curr_ind_v = _mm256_add_epi32(curr_ind_v, eights);

    vertex_0 = _mm_loadu_ps(points[i]);
    vertex_1 = _mm_loadu_ps(points[i + 1]);
    vertex_2 = _mm_loadu_ps(points[i + 2]);
    vertex_3 = _mm_loadu_ps(points[i + 3]);
    vertex_4 = _mm_loadu_ps(points[i + 4]);
    vertex_5 = _mm_loadu_ps(points[i + 5]);
    vertex_6 = _mm_loadu_ps(points[i + 6]);
    vertex_7 = _mm_set_ps(0.0f, points[i + 7][2], points[i + 7][1], points[i + 7][0]);

    vertex_0_1 = _mm256_set_m128(vertex_1, vertex_0);
    vertex_2_3 = _mm256_set_m128(vertex_3, vertex_2);
    vertex_4_5 = _mm256_set_m128(vertex_5, vertex_4);
    vertex_6_7 = _mm256_set_m128(vertex_7, vertex_6);

    dp_01 = _mm256_dp_ps(vertex_0_1, dir_v, 113);
    dp_23 = _mm256_dp_ps(vertex_2_3, dir_v, 114);
    dp_45 = _mm256_dp_ps(vertex_4_5, dir_v, 116);
    dp_67 = _mm256_dp_ps(vertex_6_7, dir_v, 120);

    tmp1 = _mm256_or_ps(dp_01, dp_23);
    tmp2 = _mm256_or_ps(dp_45, dp_67);
    dp0_7 = _mm256_or_ps(tmp1, tmp2);

    mask = _mm256_cmp_ps(dp0_7, dpmax_v, _CMP_LT_OS);
    dpmax_v = _mm256_min_ps(dp0_7, dpmax_v);
    tmp3 = _mm256_castps_si256(_mm256_and_ps(mask, _mm256_castsi256_ps(curr_ind_v)));
    tmp4 = _mm256_castps_si256(_mm256_andnot_ps(mask, _mm256_castsi256_ps(argmin_v)));
    argmin_v = _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(tmp3), _mm256_castsi256_ps(tmp4)));
  }

  if (i != 0) {
    float dps[8];
    int argmins[8];
    _mm256_storeu_ps(dps, dpmax_v);
    _mm256_storeu_epi32(argmins, argmin_v);

    for (int j = 0; j < 8; j++) {
      if (dps[j] < dpmax) {
        dpmax = dps[j];
        argmin = argmins[j];
      }
    }
  }

  // Handle the remaining cases
  for (; i < N; i++) {
    dp = dot3D_vectorized(points[i], dir);
    if (dp < dpmax) {
      dpmax = dp;
      argmin = i;
    }
  }

  res[0] -= points[argmin][0];
  res[1] -= points[argmin][1];
  res[2] -= points[argmin][2];
}

int ds_line3D_vectorized(struct Simplex3D* s, float *dir) {
  // WE KNOW A is PAST THE ORIGIN -> new Simplex s is [A, B]
  // Unlike in the video by Casey, I don't think we need to permute A and B.
  // p[0] is B, p[1] is A (the new point)
  const float(* p)[3] = s->p;
  float ab[3] = {p[0][0] - p[1][0], p[0][1] - p[1][1], p[0][2] - p[1][2]};
  float a0[3] = {-p[1][0], -p[1][1], -p[1][2]};

  // Perform the float-cross-product
  float tmp[3];
  cross3D_vectorized(ab, a0, tmp);
  cross3D_vectorized(tmp, ab, dir);

  // We don't (know if we) enclose 0
  return 0;
}

int ds_triangle3D_vectorized(struct Simplex3D* s, float *dir) {
  float(* p)[3] = s->p;
  // p[0] is B, p[1] is C, p[2] is A (the new point)
  float ab[3] = {p[0][0] - p[2][0], p[0][1] - p[2][1], p[0][2] - p[2][2]};
  float ac[3] = {p[1][0] - p[2][0], p[1][1] - p[2][1], p[1][2] - p[2][2]};
  float a0[3] = {-p[2][0], -p[2][1], -p[2][2]};
  float abc[3];
  cross3D_vectorized(ab, ac, abc);

  float tst[3];
  cross3D_vectorized(abc, ac, tst);
  if (dot3D_vectorized(tst, a0) > 0.0) {
    // [A, C]
    s->num_points = 2;
    p[0][0] = p[2][0];
    p[0][1] = p[2][1];
    p[0][2] = p[2][2];

    cross3D_vectorized(ac, a0, tst);
    cross3D_vectorized(tst, ac, dir);
  } else {
    cross3D_vectorized(ab, abc, tst);
    if (dot3D_vectorized(tst, a0) > 0.0) {
      // [A, B]
      s->num_points = 2;
      // I don't think we need to permute from [B, A] to [A, B]
      p[1][0] = p[2][0];
      p[1][1] = p[2][1];
      p[1][2] = p[2][2];

      cross3D_vectorized(ab, a0, tst);
      cross3D_vectorized(tst, ab, dir);
    }
    else {
      if (dot3D_vectorized(abc, a0) > 0.0) {
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

int ds_tetrahedron3D_vectorized(struct Simplex3D* s, float *dir) {
  float(* p)[3] = s->p;
  // p[0] is B, p[1] is C, p[2] is D, p[3] is A (the new point)
  float ab[3] = {p[0][0] - p[3][0], p[0][1] - p[3][1], p[0][2] - p[3][2]};
  float ac[3] = {p[1][0] - p[3][0], p[1][1] - p[3][1], p[1][2] - p[3][2]};
  float ad[3] = {p[2][0] - p[3][0], p[2][1] - p[3][1], p[2][2] - p[3][2]};
  float a0[3] = {-p[3][0], -p[3][1], -p[3][2]};
  float abc[3], acd[3], adb[3];
  float tst[3];

  // The vectors abc, acd and adb point OUTSIDE the tetrahedron.
  cross3D_vectorized(ab, ac, abc);
  if (dot3D_vectorized(abc, a0) > 0.0) {
    cross3D_vectorized(ac, ad, acd);
    if (dot3D_vectorized(acd, a0) > 0.0) {
      s->num_points = 2;
      // [A, C]
      p[0][0] = p[3][0];
      p[0][1] = p[3][1];
      p[0][2] = p[3][2];

      // Like in the line case (I think ?)
      cross3D_vectorized(ac, a0, tst);
      cross3D_vectorized(tst, ac, dir);
    }
    else {
      cross3D_vectorized(ad, ab, adb);
      if (dot3D_vectorized(adb, a0) > 0.0) {
        s->num_points = 2;
        // [A, B]  (tho we return [B, A] as line perm doesn't matter).
        p[1][0] = p[3][0];
        p[1][1] = p[3][1];
        p[1][2] = p[3][2];

        // Like in the line case (I think ?)
        cross3D_vectorized(ab, a0, tst);
        cross3D_vectorized(tst, ab, dir);
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
    cross3D_vectorized(ac, ad, acd);
    cross3D_vectorized(ad, ab, adb);
    if (dot3D_vectorized(acd, a0) > 0.0) {
      if (dot3D_vectorized(adb, a0) > 0.0) {
        s->num_points = 2;
        // [A, D]  ... For once we actually need to move two vectors ;-;
        p[0][0] = p[3][0];
        p[0][1] = p[3][1];
        p[0][2] = p[3][2];

        p[1][0] = p[2][0];
        p[1][1] = p[2][1];
        p[1][2] = p[2][2];

        // Like in the line case
        cross3D_vectorized(ad, a0, tst);
        cross3D_vectorized(tst, ad, dir);
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
      if (dot3D_vectorized(adb, a0) > 0.0) {
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
int do_simplex3D_vectorized(struct Simplex3D* s, float *dir) {
  switch (s->num_points) {
    case 2:return ds_line3D_vectorized(s, dir);
    case 3:return ds_triangle3D_vectorized(s, dir);
    case 4:return ds_tetrahedron3D_vectorized(s, dir);
    default:;
  }
  return -1;  // unreachable unless s points to an invalid Simplex
}

int do_intersect3D_vectorized(const struct CHObject* obj1, const struct CHObject* obj2) {
  // The search direction
  float d[3] = {1.f, 1.f, 1.f};  // could also be random
  // The point we examine in each iteration
  float a[3];

  // Our simplex that do_simplex modifies.
  struct Simplex3D s;
  s.num_points = 1;

  // S = Support(A-B) = Support(A) - Invsup(B)
  support3D_vectorized(obj1->points, d, s.p[0], obj1->num_points);
  invsup3D_vectorized(obj2->points, d, s.p[0], obj2->num_points);

  // d = -S
  d[0] = -s.p[0][0];
  d[1] = -s.p[0][1];
  d[2] = -s.p[0][2];

  int max_iter = 100;  // most likely we need less
  while (max_iter--) {
    // a = Support(A-B) = Support(A) - Invsup(B)
    support3D_vectorized(obj1->points, d, a, obj1->num_points);
    invsup3D_vectorized(obj2->points, d, a, obj2->num_points);

    // Are we past the origin? If no, we will never enclose it.
    if (dot3D_vectorized(a, d) < 0.) return 0;

    // Add newly found point to Simplex
    // AT THIS POINT WE KNOW a IS PAST THE ORIGIN
    s.p[s.num_points][0] = a[0];
    s.p[s.num_points][1] = a[1];
    s.p[s.num_points][2] = a[2];
    s.num_points++;

    if (do_simplex3D_vectorized(&s, d)) return 1;
  }

  return 1; // Most likely not reachable. Say we intersect by default.
}

[[gnu::target("avx512vl")]]
int do_intersect3D_vectorized_inlined(const struct CHObject* obj1, const struct CHObject* obj2) {
  // The search direction
  float d[3] = {1.f, 1.f, 1.f};  // could also be random
  // The point we examine in each iteration
  float a[3];

  // Our simplex that do_simplex modifies.
  struct Simplex3D s;
  s.num_points = 1;

  // Method inlining for support function
  // Object with minimum no of points determine the no. of iterations for the joint loop.
  int no_points_min = obj1->num_points > obj2->num_points ? obj2->num_points : obj1->num_points;

  // Vertices of object 1 accessed during one iteration of the loop.
  __m128 vertex_0, vertex_1, vertex_2, vertex_3, vertex_4, vertex_5, vertex_6, vertex_7;
  __m128 vertex_8, vertex_9, vertex_10, vertex_11, vertex_12, vertex_13, vertex_14, vertex_15;

  // Vertices of object 1 will be combined into a 256 bit vector.
  __m256 vertex_0_1, vertex_2_3, vertex_4_5, vertex_6_7;
  __m256 vertex_8_9, vertex_10_11, vertex_12_13, vertex_14_15;

  // Vertices of object 2 accessed during one iteration of the loop.
  __m128 inv_vertex_0, inv_vertex_1, inv_vertex_2, inv_vertex_3, inv_vertex_4, inv_vertex_5, inv_vertex_6, inv_vertex_7;

  // Vertices of object 2 will be combined into a 256 bit vector.
  __m256 inv_vertex_0_1, inv_vertex_2_3, inv_vertex_4_5, inv_vertex_6_7;

  // Variables for holding the local and maximum dot products of first object.
  __m256 dp_01, dp_23, dp_45, dp_67, dp0_7;
  __m256 dp_89, dp_1011, dp_1213, dp_1415, dp8_15;

  __m256 dpmax_v = _mm256_set1_ps(-1e9);
  __m256 dpmax_v_1 = _mm256_set1_ps(-1e9);

  __m256i argmax_v = _mm256_set1_epi32(0);
  __m256i argmax_v_1 = _mm256_set1_epi32(0);

  __m256i curr_ind_v = _mm256_set_epi32(-1, -3, -5, -7, -2, -4, -6, -8);
  __m256i curr_ind_v_1;

  // Variables for holding the local and maximum dot products of second object.
  __m256 inv_dp_01, inv_dp_23, inv_dp_45, inv_dp_67, inv_dp0_7;

  __m256 inv_dpmax_v = _mm256_set1_ps(1e9);
  __m256i inv_argmin_v = _mm256_set1_epi32(0);

  __m256 dir_v = _mm256_set_ps(0.0f, d[2], d[1], d[0], 0.0f, d[2], d[1], d[0]);

  __m256 tmp1, tmp2, mask;
  __m256 tmp5, tmp6, mask1;
  __m256i tmp3, tmp4;
  __m256i tmp7, tmp8;

  __m256i eights = _mm256_set1_epi32(8);
  __m256i sixteens = _mm256_set1_epi32(16);

  float dp, dpmax = -1e9;
  int argmax = 0;
  int i;

  for (i = 0; i < no_points_min - 7; i += 8) {
    curr_ind_v = _mm256_add_epi32(curr_ind_v, eights);

    vertex_0 = _mm_load_ps(obj1->points[i]);
    vertex_1 = _mm_loadu_ps(obj1->points[i + 1]);
    vertex_2 = _mm_loadu_ps(obj1->points[i + 2]);
    vertex_3 = _mm_loadu_ps(obj1->points[i + 3]);
    vertex_4 = _mm_loadu_ps(obj1->points[i + 4]);
    vertex_5 = _mm_loadu_ps(obj1->points[i + 5]);
    vertex_6 = _mm_loadu_ps(obj1->points[i + 6]);
    vertex_7 = _mm_set_ps(0.0f, obj1->points[i + 7][2], obj1->points[i + 7][1], obj1->points[i + 7][0]);

    inv_vertex_0 = _mm_load_ps(obj2->points[i]);
    inv_vertex_1 = _mm_loadu_ps(obj2->points[i + 1]);
    inv_vertex_2 = _mm_loadu_ps(obj2->points[i + 2]);
    inv_vertex_3 = _mm_loadu_ps(obj2->points[i + 3]);
    inv_vertex_4 = _mm_loadu_ps(obj2->points[i + 4]);
    inv_vertex_5 = _mm_loadu_ps(obj2->points[i + 5]);
    inv_vertex_6 = _mm_loadu_ps(obj2->points[i + 6]);
    inv_vertex_7 = _mm_set_ps(0.0f, obj2->points[i + 7][2], obj2->points[i + 7][1], obj2->points[i + 7][0]);

    vertex_0_1 = _mm256_set_m128(vertex_1, vertex_0);
    vertex_2_3 = _mm256_set_m128(vertex_3, vertex_2);
    vertex_4_5 = _mm256_set_m128(vertex_5, vertex_4);
    vertex_6_7 = _mm256_set_m128(vertex_7, vertex_6);

    inv_vertex_0_1 = _mm256_set_m128(inv_vertex_1, inv_vertex_0);
    inv_vertex_2_3 = _mm256_set_m128(inv_vertex_3, inv_vertex_2);
    inv_vertex_4_5 = _mm256_set_m128(inv_vertex_5, inv_vertex_4);
    inv_vertex_6_7 = _mm256_set_m128(inv_vertex_7, inv_vertex_6);

    dp_01 = _mm256_dp_ps(vertex_0_1, dir_v, 113);
    dp_23 = _mm256_dp_ps(vertex_2_3, dir_v, 114);
    dp_45 = _mm256_dp_ps(vertex_4_5, dir_v, 116);
    dp_67 = _mm256_dp_ps(vertex_6_7, dir_v, 120);

    inv_dp_01 = _mm256_dp_ps(inv_vertex_0_1, dir_v, 113);
    inv_dp_23 = _mm256_dp_ps(inv_vertex_2_3, dir_v, 114);
    inv_dp_45 = _mm256_dp_ps(inv_vertex_4_5, dir_v, 116);
    inv_dp_67 = _mm256_dp_ps(inv_vertex_6_7, dir_v, 120);

    tmp1 = _mm256_or_ps(dp_01, dp_23);
    tmp2 = _mm256_or_ps(dp_45, dp_67);
    dp0_7 = _mm256_or_ps(tmp1, tmp2);

    tmp5 = _mm256_or_ps(inv_dp_01, inv_dp_23);
    tmp6 = _mm256_or_ps(inv_dp_45, inv_dp_67);
    inv_dp0_7 = _mm256_or_ps(tmp5, tmp6);

    mask = _mm256_cmp_ps(dp0_7, dpmax_v, _CMP_GT_OS);
    dpmax_v = _mm256_max_ps(dp0_7, dpmax_v);
    tmp3 = _mm256_castps_si256(_mm256_and_ps(mask, _mm256_castsi256_ps(curr_ind_v)));
    tmp4 = _mm256_castps_si256(_mm256_andnot_ps(mask, _mm256_castsi256_ps(argmax_v)));
    argmax_v = _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(tmp3), _mm256_castsi256_ps(tmp4)));

    mask1 = _mm256_cmp_ps(inv_dp0_7, inv_dpmax_v, _CMP_LT_OS);
    inv_dpmax_v = _mm256_min_ps(inv_dp0_7, inv_dpmax_v);
    tmp7 = _mm256_castps_si256(_mm256_and_ps(mask1, _mm256_castsi256_ps(curr_ind_v)));
    tmp8 = _mm256_castps_si256(_mm256_andnot_ps(mask1, _mm256_castsi256_ps(inv_argmin_v)));
    inv_argmin_v = _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(tmp7), _mm256_castsi256_ps(tmp8)));
  }

  int k = i;
  curr_ind_v = _mm256_set_epi32( i - 9, i - 11, i - 13, i - 15, i - 10, i - 12, i - 14, i - 16);
  curr_ind_v_1 = _mm256_set_epi32( i - 1, i - 3, i - 5, i - 7, i - 2, i - 4, i - 6, i - 8);
  for (; i < obj1->num_points - 15; i += 16) {
    curr_ind_v = _mm256_add_epi32(curr_ind_v, sixteens);
    curr_ind_v_1 = _mm256_add_epi32(curr_ind_v_1, sixteens);

    vertex_0 = _mm_load_ps(obj1->points[i]);
    vertex_1 = _mm_loadu_ps(obj1->points[i + 1]);
    vertex_2 = _mm_loadu_ps(obj1->points[i + 2]);
    vertex_3 = _mm_loadu_ps(obj1->points[i + 3]);
    vertex_4 = _mm_loadu_ps(obj1->points[i + 4]);
    vertex_5 = _mm_loadu_ps(obj1->points[i + 5]);
    vertex_6 = _mm_loadu_ps(obj1->points[i + 6]);
    vertex_7 = _mm_loadu_ps(obj1->points[i + 7]);

    vertex_8 = _mm_load_ps(obj1->points[i + 8]);
    vertex_9 = _mm_loadu_ps(obj1->points[i + 9]);
    vertex_10 = _mm_loadu_ps(obj1->points[i + 10]);
    vertex_11 = _mm_loadu_ps(obj1->points[i + 11]);
    vertex_12 = _mm_loadu_ps(obj1->points[i + 12]);
    vertex_13 = _mm_loadu_ps(obj1->points[i + 13]);
    vertex_14 = _mm_loadu_ps(obj1->points[i + 14]);
    vertex_15 = _mm_set_ps(0.0f, obj1->points[i + 15][2], obj1->points[i + 15][1], obj1->points[i + 15][0]);

    vertex_0_1 = _mm256_set_m128(vertex_1, vertex_0);
    vertex_2_3 = _mm256_set_m128(vertex_3, vertex_2);
    vertex_4_5 = _mm256_set_m128(vertex_5, vertex_4);
    vertex_6_7 = _mm256_set_m128(vertex_7, vertex_6);

    vertex_8_9 = _mm256_set_m128(vertex_9, vertex_8);
    vertex_10_11 = _mm256_set_m128(vertex_11, vertex_10);
    vertex_12_13 = _mm256_set_m128(vertex_13, vertex_12);
    vertex_14_15 = _mm256_set_m128(vertex_15, vertex_14);

    dp_01 = _mm256_dp_ps(vertex_0_1, dir_v, 113);
    dp_23 = _mm256_dp_ps(vertex_2_3, dir_v, 114);
    dp_45 = _mm256_dp_ps(vertex_4_5, dir_v, 116);
    dp_67 = _mm256_dp_ps(vertex_6_7, dir_v, 120);

    dp_89 = _mm256_dp_ps(vertex_8_9, dir_v, 113);
    dp_1011 = _mm256_dp_ps(vertex_10_11, dir_v, 114);
    dp_1213 = _mm256_dp_ps(vertex_12_13, dir_v, 116);
    dp_1415 = _mm256_dp_ps(vertex_14_15, dir_v, 120);

    tmp1 = _mm256_or_ps(dp_01, dp_23);
    tmp2 = _mm256_or_ps(dp_45, dp_67);
    dp0_7 = _mm256_or_ps(tmp1, tmp2);

    tmp5 = _mm256_or_ps(dp_89, dp_1011);
    tmp6 = _mm256_or_ps(dp_1213, dp_1415);
    dp8_15 = _mm256_or_ps(tmp5, tmp6);

    mask = _mm256_cmp_ps(dp0_7, dpmax_v, _CMP_GT_OS);
    dpmax_v = _mm256_max_ps(dp0_7, dpmax_v);
    tmp3 = _mm256_castps_si256(_mm256_and_ps(mask, _mm256_castsi256_ps(curr_ind_v)));
    tmp4 = _mm256_castps_si256(_mm256_andnot_ps(mask, _mm256_castsi256_ps(argmax_v)));
    argmax_v = _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(tmp3), _mm256_castsi256_ps(tmp4)));

    mask1 = _mm256_cmp_ps(dp8_15, dpmax_v_1, _CMP_GT_OS);
    dpmax_v_1 = _mm256_max_ps(dp8_15, dpmax_v_1);
    tmp7 = _mm256_castps_si256(_mm256_and_ps(mask1, _mm256_castsi256_ps(curr_ind_v_1)));
    tmp8 = _mm256_castps_si256(_mm256_andnot_ps(mask1, _mm256_castsi256_ps(argmax_v_1)));
    argmax_v_1 = _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(tmp7), _mm256_castsi256_ps(tmp8)));
  }

  mask = _mm256_cmp_ps(dpmax_v, dpmax_v_1, _CMP_GT_OS);
  dpmax_v = _mm256_max_ps(dpmax_v, dpmax_v_1);
  tmp3 = _mm256_castps_si256(_mm256_and_ps(mask, _mm256_castsi256_ps(argmax_v)));
  tmp4 = _mm256_castps_si256(_mm256_andnot_ps(mask, _mm256_castsi256_ps(argmax_v_1)));
  argmax_v = _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(tmp3), _mm256_castsi256_ps(tmp4)));

  float dps[8];
  int argmaxs[8];
  if (i != 0) {
    _mm256_store_ps(dps, dpmax_v);
    _mm256_store_epi32(argmaxs, argmax_v);

    for (int j = 0; j < 8; j++) {
      if (dps[j] > dpmax) {
        dpmax = dps[j];
        argmax = argmaxs[j];
      }
    }
  }

  // Handle the remaining cases
  for (; i < obj1->num_points; i++) {
    dp = dot3D_vectorized(obj1->points[i], d);
    if (dp > dpmax) {
      dpmax = dp;
      argmax = i;
    }
  }

  float inv_dp, inv_dpmax = 1e9;
  int argmin = 0;
  curr_ind_v = _mm256_set_epi32( k - 1, k - 3, k - 5, k - 7, k - 2, k - 4, k - 6, k - 8);
  for (; k < obj2->num_points - 7; k += 8) {
    curr_ind_v = _mm256_add_epi32(curr_ind_v, eights);

    inv_vertex_0 = _mm_load_ps(obj2->points[k]);
    inv_vertex_1 = _mm_loadu_ps(obj2->points[k + 1]);
    inv_vertex_2 = _mm_loadu_ps(obj2->points[k + 2]);
    inv_vertex_3 = _mm_loadu_ps(obj2->points[k + 3]);
    inv_vertex_4 = _mm_loadu_ps(obj2->points[k + 4]);
    inv_vertex_5 = _mm_loadu_ps(obj2->points[k + 5]);
    inv_vertex_6 = _mm_loadu_ps(obj2->points[k + 6]);
    inv_vertex_7 = _mm_set_ps(0.0f, obj2->points[k + 7][2], obj2->points[k + 7][1], obj2->points[k + 7][0]);

    inv_vertex_0_1 = _mm256_set_m128(inv_vertex_1, inv_vertex_0);
    inv_vertex_2_3 = _mm256_set_m128(inv_vertex_3, inv_vertex_2);
    inv_vertex_4_5 = _mm256_set_m128(inv_vertex_5, inv_vertex_4);
    inv_vertex_6_7 = _mm256_set_m128(inv_vertex_7, inv_vertex_6);

    inv_dp_01 = _mm256_dp_ps(inv_vertex_0_1, dir_v, 113);
    inv_dp_23 = _mm256_dp_ps(inv_vertex_2_3, dir_v, 114);
    inv_dp_45 = _mm256_dp_ps(inv_vertex_4_5, dir_v, 116);
    inv_dp_67 = _mm256_dp_ps(inv_vertex_6_7, dir_v, 120);

    tmp1 = _mm256_or_ps(inv_dp_01, inv_dp_23);
    tmp2 = _mm256_or_ps(inv_dp_45, inv_dp_67);
    inv_dp0_7 = _mm256_or_ps(tmp1, tmp2);

    mask = _mm256_cmp_ps(inv_dp0_7, inv_dpmax_v, _CMP_LT_OS);
    inv_dpmax_v = _mm256_min_ps(inv_dp0_7, inv_dpmax_v);
    tmp3 = _mm256_castps_si256(_mm256_and_ps(mask, _mm256_castsi256_ps(curr_ind_v)));
    tmp4 = _mm256_castps_si256(_mm256_andnot_ps(mask, _mm256_castsi256_ps(inv_argmin_v)));
    inv_argmin_v = _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(tmp3), _mm256_castsi256_ps(tmp4)));
  }

  int argmins[8];
  if (k != 0) {
    _mm256_store_ps(dps, inv_dpmax_v);
    _mm256_store_epi32(argmins, inv_argmin_v);

    for (int j = 0; j < 8; j++) {
      if (dps[j] < inv_dpmax) {
        inv_dpmax = dps[j];
        argmin = argmins[j];
      }
    }
  }

  // Handle the remaining cases
  for (; k < obj2->num_points; k++) {
    inv_dp = dot3D_vectorized(obj2->points[k], d);
    if (inv_dp < inv_dpmax) {
      inv_dpmax = inv_dp;
      argmin = k;
    }
  }

  // S = Support(A-B) = Support(A) - Invsup(B)
  s.p[0][0] = obj1->points[argmax][0] - obj2->points[argmin][0];
  s.p[0][1] = obj1->points[argmax][1] - obj2->points[argmin][1];
  s.p[0][2] = obj1->points[argmax][2] - obj2->points[argmin][2];

  // d = -S
  d[0] = -s.p[0][0];
  d[1] = -s.p[0][1];
  d[2] = -s.p[0][2];

  int max_iter = 100;  // most likely we need less
  while (max_iter--) {
    // a = Support(A-B) = Support(A) - Invsup(B)

    // Reset variables
    dpmax_v = _mm256_set1_ps(-1e9);
    dpmax_v_1 = _mm256_set1_ps(-1e9);
    dpmax = -1e9;

    argmax_v = _mm256_set1_epi32(0);
    argmax_v_1 = _mm256_set1_epi32(0);
    argmax = 0;

    inv_dpmax_v = _mm256_set1_ps(1e9);
    inv_dpmax = 1e9;
    inv_argmin_v = _mm256_set1_epi32(0);
    argmin = 0;

    dir_v = _mm256_set_ps(0.0f, d[2], d[1], d[0], 0.0f, d[2], d[1], d[0]);
    curr_ind_v = _mm256_set_epi32(-1, -3, -5, -7, -2, -4, -6, -8);

    // Joint loop for
    for (i = 0; i < no_points_min - 7; i += 8) {
      curr_ind_v = _mm256_add_epi32(curr_ind_v, eights);

      vertex_0 = _mm_load_ps(obj1->points[i]);
      vertex_1 = _mm_loadu_ps(obj1->points[i + 1]);
      vertex_2 = _mm_loadu_ps(obj1->points[i + 2]);
      vertex_3 = _mm_loadu_ps(obj1->points[i + 3]);
      vertex_4 = _mm_loadu_ps(obj1->points[i + 4]);
      vertex_5 = _mm_loadu_ps(obj1->points[i + 5]);
      vertex_6 = _mm_loadu_ps(obj1->points[i + 6]);
      vertex_7 = _mm_set_ps(0.0f, obj1->points[i + 7][2], obj1->points[i + 7][1], obj1->points[i + 7][0]);//_mm_loadu_ps(obj1->points[i + 7]);

      inv_vertex_0 = _mm_load_ps(obj2->points[i]);
      inv_vertex_1 = _mm_loadu_ps(obj2->points[i + 1]);
      inv_vertex_2 = _mm_loadu_ps(obj2->points[i + 2]);
      inv_vertex_3 = _mm_loadu_ps(obj2->points[i + 3]);
      inv_vertex_4 = _mm_loadu_ps(obj2->points[i + 4]);
      inv_vertex_5 = _mm_loadu_ps(obj2->points[i + 5]);
      inv_vertex_6 = _mm_loadu_ps(obj2->points[i + 6]);
      inv_vertex_7 = _mm_set_ps(0.0f, obj2->points[i + 7][2], obj2->points[i + 7][1], obj2->points[i + 7][0]);

      vertex_0_1 = _mm256_set_m128(vertex_1, vertex_0);
      vertex_2_3 = _mm256_set_m128(vertex_3, vertex_2);
      vertex_4_5 = _mm256_set_m128(vertex_5, vertex_4);
      vertex_6_7 = _mm256_set_m128(vertex_7, vertex_6);

      inv_vertex_0_1 = _mm256_set_m128(inv_vertex_1, inv_vertex_0);
      inv_vertex_2_3 = _mm256_set_m128(inv_vertex_3, inv_vertex_2);
      inv_vertex_4_5 = _mm256_set_m128(inv_vertex_5, inv_vertex_4);
      inv_vertex_6_7 = _mm256_set_m128(inv_vertex_7, inv_vertex_6);

      dp_01 = _mm256_dp_ps(vertex_0_1, dir_v, 113);
      dp_23 = _mm256_dp_ps(vertex_2_3, dir_v, 114);
      dp_45 = _mm256_dp_ps(vertex_4_5, dir_v, 116);
      dp_67 = _mm256_dp_ps(vertex_6_7, dir_v, 120);

      inv_dp_01 = _mm256_dp_ps(inv_vertex_0_1, dir_v, 113);
      inv_dp_23 = _mm256_dp_ps(inv_vertex_2_3, dir_v, 114);
      inv_dp_45 = _mm256_dp_ps(inv_vertex_4_5, dir_v, 116);
      inv_dp_67 = _mm256_dp_ps(inv_vertex_6_7, dir_v, 120);

      tmp1 = _mm256_or_ps(dp_01, dp_23);
      tmp2 = _mm256_or_ps(dp_45, dp_67);
      dp0_7 = _mm256_or_ps(tmp1, tmp2);

      tmp5 = _mm256_or_ps(inv_dp_01, inv_dp_23);
      tmp6 = _mm256_or_ps(inv_dp_45, inv_dp_67);
      inv_dp0_7 = _mm256_or_ps(tmp5, tmp6);

      mask = _mm256_cmp_ps(dp0_7, dpmax_v, _CMP_GT_OS);
      dpmax_v = _mm256_max_ps(dp0_7, dpmax_v);
      tmp3 = _mm256_castps_si256(_mm256_and_ps(mask, _mm256_castsi256_ps(curr_ind_v)));
      tmp4 = _mm256_castps_si256(_mm256_andnot_ps(mask, _mm256_castsi256_ps(argmax_v)));
      argmax_v = _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(tmp3), _mm256_castsi256_ps(tmp4)));

      mask1 = _mm256_cmp_ps(inv_dp0_7, inv_dpmax_v, _CMP_LT_OS);
      inv_dpmax_v = _mm256_min_ps(inv_dp0_7, inv_dpmax_v);
      tmp7 = _mm256_castps_si256(_mm256_and_ps(mask1, _mm256_castsi256_ps(curr_ind_v)));
      tmp8 = _mm256_castps_si256(_mm256_andnot_ps(mask1, _mm256_castsi256_ps(inv_argmin_v)));
      inv_argmin_v = _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(tmp7), _mm256_castsi256_ps(tmp8)));
    }

    k = i;
    curr_ind_v = _mm256_set_epi32( i - 9, i - 11, i - 13, i - 15, i - 10, i - 12, i - 14, i - 16);
    curr_ind_v_1 = _mm256_set_epi32( i - 1, i - 3, i - 5, i - 7, i - 2, i - 4, i - 6, i - 8);
    for (; i < obj1->num_points - 15; i += 16) {
      curr_ind_v = _mm256_add_epi32(curr_ind_v, sixteens);
      curr_ind_v_1 = _mm256_add_epi32(curr_ind_v_1, sixteens);

      vertex_0 = _mm_load_ps(obj1->points[i]);
      vertex_1 = _mm_loadu_ps(obj1->points[i + 1]);
      vertex_2 = _mm_loadu_ps(obj1->points[i + 2]);
      vertex_3 = _mm_loadu_ps(obj1->points[i + 3]);
      vertex_4 = _mm_loadu_ps(obj1->points[i + 4]);
      vertex_5 = _mm_loadu_ps(obj1->points[i + 5]);
      vertex_6 = _mm_loadu_ps(obj1->points[i + 6]);
      vertex_7 = _mm_loadu_ps(obj1->points[i + 7]);

      vertex_8 = _mm_load_ps(obj1->points[i + 8]);
      vertex_9 = _mm_loadu_ps(obj1->points[i + 9]);
      vertex_10 = _mm_loadu_ps(obj1->points[i + 10]);
      vertex_11 = _mm_loadu_ps(obj1->points[i + 11]);
      vertex_12 = _mm_loadu_ps(obj1->points[i + 12]);
      vertex_13 = _mm_loadu_ps(obj1->points[i + 13]);
      vertex_14 = _mm_loadu_ps(obj1->points[i + 14]);
      vertex_15 = _mm_set_ps(0.0f, obj1->points[i + 15][2], obj1->points[i + 15][1], obj1->points[i + 15][0]);

      vertex_0_1 = _mm256_set_m128(vertex_1, vertex_0);
      vertex_2_3 = _mm256_set_m128(vertex_3, vertex_2);
      vertex_4_5 = _mm256_set_m128(vertex_5, vertex_4);
      vertex_6_7 = _mm256_set_m128(vertex_7, vertex_6);

      vertex_8_9 = _mm256_set_m128(vertex_9, vertex_8);
      vertex_10_11 = _mm256_set_m128(vertex_11, vertex_10);
      vertex_12_13 = _mm256_set_m128(vertex_13, vertex_12);
      vertex_14_15 = _mm256_set_m128(vertex_15, vertex_14);

      dp_01 = _mm256_dp_ps(vertex_0_1, dir_v, 113);
      dp_23 = _mm256_dp_ps(vertex_2_3, dir_v, 114);
      dp_45 = _mm256_dp_ps(vertex_4_5, dir_v, 116);
      dp_67 = _mm256_dp_ps(vertex_6_7, dir_v, 120);

      dp_89 = _mm256_dp_ps(vertex_8_9, dir_v, 113);
      dp_1011 = _mm256_dp_ps(vertex_10_11, dir_v, 114);
      dp_1213 = _mm256_dp_ps(vertex_12_13, dir_v, 116);
      dp_1415 = _mm256_dp_ps(vertex_14_15, dir_v, 120);

      tmp1 = _mm256_or_ps(dp_01, dp_23);
      tmp2 = _mm256_or_ps(dp_45, dp_67);
      dp0_7 = _mm256_or_ps(tmp1, tmp2);

      tmp5 = _mm256_or_ps(dp_89, dp_1011);
      tmp6 = _mm256_or_ps(dp_1213, dp_1415);
      dp8_15 = _mm256_or_ps(tmp5, tmp6);

      mask = _mm256_cmp_ps(dp0_7, dpmax_v, _CMP_GT_OS);
      dpmax_v = _mm256_max_ps(dp0_7, dpmax_v);
      tmp3 = _mm256_castps_si256(_mm256_and_ps(mask, _mm256_castsi256_ps(curr_ind_v)));
      tmp4 = _mm256_castps_si256(_mm256_andnot_ps(mask, _mm256_castsi256_ps(argmax_v)));
      argmax_v = _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(tmp3), _mm256_castsi256_ps(tmp4)));

      mask1 = _mm256_cmp_ps(dp8_15, dpmax_v_1, _CMP_GT_OS);
      dpmax_v_1 = _mm256_max_ps(dp8_15, dpmax_v_1);
      tmp7 = _mm256_castps_si256(_mm256_and_ps(mask1, _mm256_castsi256_ps(curr_ind_v_1)));
      tmp8 = _mm256_castps_si256(_mm256_andnot_ps(mask1, _mm256_castsi256_ps(argmax_v_1)));
      argmax_v_1 = _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(tmp7), _mm256_castsi256_ps(tmp8)));
    }

    if (i != 0) {
      mask = _mm256_cmp_ps(dpmax_v, dpmax_v_1, _CMP_GT_OS);
      dpmax_v = _mm256_max_ps(dpmax_v, dpmax_v_1);
      tmp3 = _mm256_castps_si256(_mm256_and_ps(mask, _mm256_castsi256_ps(argmax_v)));
      tmp4 = _mm256_castps_si256(_mm256_andnot_ps(mask, _mm256_castsi256_ps(argmax_v_1)));
      argmax_v = _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(tmp3), _mm256_castsi256_ps(tmp4)));

      _mm256_store_ps(dps, dpmax_v);
      _mm256_store_epi32(argmaxs, argmax_v);

      for (int j = 0; j < 8; j++) {
        if (dps[j] > dpmax) {
          dpmax = dps[j];
          argmax = argmaxs[j];
        }
      }
    }

    // Handle the remaining cases
    for (; i < obj1->num_points; i++) {
      dp = dot3D_vectorized(obj1->points[i], d);
      if (dp > dpmax) {
        dpmax = dp;
        argmax = i;
      }
    }

    curr_ind_v = _mm256_set_epi32( k - 1, k - 3, k - 5, k - 7, k - 2, k - 4, k - 6, k - 8);
    for (; k < obj2->num_points - 7; k += 8) {
      curr_ind_v = _mm256_add_epi32(curr_ind_v, eights);

      inv_vertex_0 = _mm_load_ps(obj2->points[k]);
      inv_vertex_1 = _mm_loadu_ps(obj2->points[k + 1]);
      inv_vertex_2 = _mm_loadu_ps(obj2->points[k + 2]);
      inv_vertex_3 = _mm_loadu_ps(obj2->points[k + 3]);
      inv_vertex_4 = _mm_loadu_ps(obj2->points[k + 4]);
      inv_vertex_5 = _mm_loadu_ps(obj2->points[k + 5]);
      inv_vertex_6 = _mm_loadu_ps(obj2->points[k + 6]);
      inv_vertex_7 = _mm_set_ps(0.0f, obj2->points[i + 7][2], obj2->points[i + 7][1], obj2->points[i + 7][0]);

      inv_vertex_0_1 = _mm256_set_m128(inv_vertex_1, inv_vertex_0);
      inv_vertex_2_3 = _mm256_set_m128(inv_vertex_3, inv_vertex_2);
      inv_vertex_4_5 = _mm256_set_m128(inv_vertex_5, inv_vertex_4);
      inv_vertex_6_7 = _mm256_set_m128(inv_vertex_7, inv_vertex_6);

      inv_dp_01 = _mm256_dp_ps(inv_vertex_0_1, dir_v, 113);
      inv_dp_23 = _mm256_dp_ps(inv_vertex_2_3, dir_v, 114);
      inv_dp_45 = _mm256_dp_ps(inv_vertex_4_5, dir_v, 116);
      inv_dp_67 = _mm256_dp_ps(inv_vertex_6_7, dir_v, 120);

      tmp1 = _mm256_or_ps(inv_dp_01, inv_dp_23);
      tmp2 = _mm256_or_ps(inv_dp_45, inv_dp_67);
      inv_dp0_7 = _mm256_or_ps(tmp1, tmp2);

      mask = _mm256_cmp_ps(inv_dp0_7, inv_dpmax_v, _CMP_LT_OS);
      inv_dpmax_v = _mm256_min_ps(inv_dp0_7, inv_dpmax_v);
      tmp3 = _mm256_castps_si256(_mm256_and_ps(mask, _mm256_castsi256_ps(curr_ind_v)));
      tmp4 = _mm256_castps_si256(_mm256_andnot_ps(mask, _mm256_castsi256_ps(inv_argmin_v)));
      inv_argmin_v = _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(tmp3), _mm256_castsi256_ps(tmp4)));
    }

    if (k != 0) {
      int argmins[8];
      _mm256_store_ps(dps, inv_dpmax_v);
      _mm256_store_epi32(argmins, inv_argmin_v);

      for (int j = 0; j < 8; j++) {
        if (dps[j] < inv_dpmax) {
          inv_dpmax = dps[j];
          argmin = argmins[j];
        }
      }
    }

    // Handle the remaining cases
    for (; k < obj2->num_points; k++) {
      inv_dp = dot3D_vectorized(obj2->points[k], d);
      if (inv_dp < inv_dpmax) {
        inv_dpmax = inv_dp;
        argmin = k;
      }
    }

    a[0] = obj1->points[argmax][0] - obj2->points[argmin][0];
    a[1] = obj1->points[argmax][1] - obj2->points[argmin][1];
    a[2] = obj1->points[argmax][2] - obj2->points[argmin][2];

    // Are we past the origin? If no, we will never enclose it.
    if (dot3D_vectorized(a, d) < 0.) return 0;

    // Add newly found point to Simplex
    // AT THIS POINT WE KNOW a IS PAST THE ORIGIN
    s.p[s.num_points][0] = a[0];
    s.p[s.num_points][1] = a[1];
    s.p[s.num_points][2] = a[2];
    s.num_points++;

    if (do_simplex3D_vectorized(&s, d)) return 1;
  }

  return 1; // Most likely not reachable. Say we intersect by default.
}