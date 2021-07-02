#include "gjk.h"
#include <stdio.h>

/**
 * Dot product of two 3D vectors.
 * */
float dot3D_optimized(const float x[3], const float y[3]) {
  return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

/**
 * Cross product of two 3D vectors. Write into res.
 * */
void cross3D_optimized(const float x[3], const float y[3], float res[3]) {
  // Scalar replacement - Decrease memory accesses
  // - probably won't affect performance, values are already in the cache
  float x0, x1, x2, y0, y1, y2;

  x0 = x[0];
  x1 = x[1];
  x2 = x[2];
  y0 = y[0];
  y1 = y[1];
  y2 = y[2];

  res[0] = x1 * y2 - x2 * y1;
  res[1] = x2 * y0 - x0 * y2;
  res[2] = x0 * y1 - x1 * y0;
}

int ds_line3D_optimized(struct Simplex3D* s, float dir[3]) {
  // WE KNOW A is PAST THE ORIGIN -> new Simplex s is [A, B]
  // Unlike in the video by Casey, I don't think we need to permute A and B.
  // p[0] is B, p[1] is A (the new point)
  const float(* p)[3] = s->p;
  float ab[3] = {p[0][0] - p[1][0], p[0][1] - p[1][1], p[0][2] - p[1][2]};
  float a0[3] = {-p[1][0], -p[1][1], -p[1][2]};

  // Perform the float-cross-product
  float tmp[3];
  cross3D_optimized(ab, a0, tmp);
  cross3D_optimized(tmp, ab, dir);

  // We don't (know if we) enclose 0
  return 0;
}

int ds_triangle3D_optimized(struct Simplex3D* s, float dir[3]) {
  float(* p)[3] = s->p;
  // p[0] is B, p[1] is C, p[2] is A (the new point)
  float ab[3] = {p[0][0] - p[2][0], p[0][1] - p[2][1], p[0][2] - p[2][2]};
  float ac[3] = {p[1][0] - p[2][0], p[1][1] - p[2][1], p[1][2] - p[2][2]};
  float a0[3] = {-p[2][0], -p[2][1], -p[2][2]};
  float abc[3];
  cross3D_optimized(ab, ac, abc);

  float tst[3];
  cross3D_optimized(abc, ac, tst);
  if (dot3D_optimized(tst, a0) > 0.0) {
    // [A, C]
    s->num_points = 2;
    p[0][0] = p[2][0];
    p[0][1] = p[2][1];
    p[0][2] = p[2][2];

    cross3D_optimized(ac, a0, tst);
    cross3D_optimized(tst, ac, dir);
  } else {
    cross3D_optimized(ab, abc, tst);
    if (dot3D_optimized(tst, a0) > 0.0) {
      // [A, B]
      s->num_points = 2;
      // I don't think we need to permute from [B, A] to [A, B]
      p[1][0] = p[2][0];
      p[1][1] = p[2][1];
      p[1][2] = p[2][2];

      cross3D_optimized(ab, a0, tst);
      cross3D_optimized(tst, ab, dir);
    }
    else {
      if (dot3D_optimized(abc, a0) > 0.0) {
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

int ds_tetrahedron3D_optimized(struct Simplex3D* s, float dir[3]) {
  float(* p)[3] = s->p;
  // p[0] is B, p[1] is C, p[2] is D, p[3] is A (the new point)
  float ab[3] = {p[0][0] - p[3][0], p[0][1] - p[3][1], p[0][2] - p[3][2]};
  float ac[3] = {p[1][0] - p[3][0], p[1][1] - p[3][1], p[1][2] - p[3][2]};
  float ad[3] = {p[2][0] - p[3][0], p[2][1] - p[3][1], p[2][2] - p[3][2]};
  float a0[3] = {-p[3][0], -p[3][1], -p[3][2]};
  float abc[3], acd[3], adb[3];
  float tst[3];

  // The vectors abc, acd and adb point OUTSIDE the tetrahedron.
  cross3D_optimized(ab, ac, abc);
  if (dot3D_optimized(abc, a0) > 0.0) {
    cross3D_optimized(ac, ad, acd);
    if (dot3D_optimized(acd, a0) > 0.0) {
      s->num_points = 2;
      // [A, C]
      p[0][0] = p[3][0];
      p[0][1] = p[3][1];
      p[0][2] = p[3][2];

      // Like in the line case (I think ?)
      cross3D_optimized(ac, a0, tst);
      cross3D_optimized(tst, ac, dir);
    }
    else {
      cross3D_optimized(ad, ab, adb);
      if (dot3D_optimized(adb, a0) > 0.0) {
        s->num_points = 2;
        // [A, B]  (tho we return [B, A] as line perm doesn't matter).
        p[1][0] = p[3][0];
        p[1][1] = p[3][1];
        p[1][2] = p[3][2];

        // Like in the line case (I think ?)
        cross3D_optimized(ab, a0, tst);
        cross3D_optimized(tst, ab, dir);
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
    cross3D_optimized(ac, ad, acd);
    cross3D_optimized(ad, ab, adb);
    if (dot3D_optimized(acd, a0) > 0.0) {
      if (dot3D_optimized(adb, a0) > 0.0) {
        s->num_points = 2;
        // [A, D]  ... For once we actually need to move two vectors ;-;
        p[0][0] = p[3][0];
        p[0][1] = p[3][1];
        p[0][2] = p[3][2];

        p[1][0] = p[2][0];
        p[1][1] = p[2][1];
        p[1][2] = p[2][2];

        // Like in the line case
        cross3D_optimized(ad, a0, tst);
        cross3D_optimized(tst, ad, dir);
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
      if (dot3D_optimized(adb, a0) > 0.0) {
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
int do_simplex3D_optimized(struct Simplex3D* s, float dir[3]) {
  switch (s->num_points) {
    case 2:return ds_line3D_optimized(s, dir);
    case 3:return ds_triangle3D_optimized(s, dir);
    case 4:return ds_tetrahedron3D_optimized(s, dir);
    default:;
  }
  return -1;  // unreachable unless s points to an invalid Simplex
}

//-----------------------------OPTIMIZED IMPLEMENTATION HERE-------------------------------
/**
 * Optimized Implementation of support3D_vectorized function.
 * */
void support3D_optimized(const float (* points)[3], const float dir[3], float res[3], int N) {
  float dp, dp_1, dp_2, dp_3, dp_4, dp_5, dp_6, dp_7;
  float dpmax = -1e9, dpmax_1 = -1e9, dpmax_2 = -1e9,
    dpmax_3 = -1e9, dpmax_4 = -1e9, dpmax_5 = -1e9, dpmax_6 = -1e9, dpmax_7 = -1e9;
  int argmax = 0, argmax_1 = 0, argmax_2 = 0,
    argmax_3 = 0, argmax_4 = 0, argmax_5 = 0, argmax_6 = 0, argmax_7 = 0;

  // Scalar replacement variables for input points
  float pts_i0, pts_i1, pts_i2;
  float pts_i0_1, pts_i1_1, pts_i2_1;
  float pts_i0_2, pts_i1_2, pts_i2_2;
  float pts_i0_3, pts_i1_3, pts_i2_3;
  float pts_i0_4, pts_i1_4, pts_i2_4;
  float pts_i0_5, pts_i1_5, pts_i2_5;
  float pts_i0_6, pts_i1_6, pts_i2_6;
  float pts_i0_7, pts_i1_7, pts_i2_7;

  // Scalar replacement variables for products and summations in dot product
  float prd_0, prd_1, prd_2, sum_0;
  float prd_0_1, prd_1_1, prd_2_1, sum_0_1;
  float prd_0_2, prd_1_2, prd_2_2, sum_0_2;
  float prd_0_3, prd_1_3, prd_2_3, sum_0_3;
  float prd_0_4, prd_1_4, prd_2_4, sum_0_4;
  float prd_0_5, prd_1_5, prd_2_5, sum_0_5;
  float prd_0_6, prd_1_6, prd_2_6, sum_0_6;
  float prd_0_7, prd_1_7, prd_2_7, sum_0_7;

  // Reduce memory access by factoring out locations accessed in every iteration.
  float dir_0 = dir[0];
  float dir_1 = dir[1];
  float dir_2 = dir[2];

  int i;
  // Loop unroll 8 times
  for (i = 0; i < N - 7; i += 8) {
    // Load
    pts_i0 = points[i][0];
    pts_i1 = points[i][1];
    pts_i2 = points[i][2];

    pts_i0_1 = points[i + 1][0];
    pts_i1_1 = points[i + 1][1];
    pts_i2_1 = points[i + 1][2];

    pts_i0_2 = points[i + 2][0];
    pts_i1_2 = points[i + 2][1];
    pts_i2_2 = points[i + 2][2];

    pts_i0_3 = points[i + 3][0];
    pts_i1_3 = points[i + 3][1];
    pts_i2_3 = points[i + 3][2];

    pts_i0_4 = points[i + 4][0];
    pts_i1_4 = points[i + 4][1];
    pts_i2_4 = points[i + 4][2];

    pts_i0_5 = points[i + 5][0];
    pts_i1_5 = points[i + 5][1];
    pts_i2_5 = points[i + 5][2];

    pts_i0_6 = points[i + 6][0];
    pts_i1_6 = points[i + 6][1];
    pts_i2_6 = points[i + 6][2];

    pts_i0_7 = points[i + 7][0];
    pts_i1_7 = points[i + 7][1];
    pts_i2_7 = points[i + 7][2];

    // Compute
    prd_0 = pts_i0 * dir_0;
    prd_1 = pts_i1 * dir_1;
    sum_0 = prd_0 + prd_1;
    prd_2 = pts_i2 * dir_2;

    prd_0_1 = pts_i0_1 * dir_0;
    prd_1_1 = pts_i1_1 * dir_1;
    sum_0_1 = prd_0_1 + prd_1_1;
    prd_2_1 = pts_i2_1 * dir_2;

    prd_0_2 = pts_i0_2 * dir_0;
    prd_1_2 = pts_i1_2 * dir_1;
    sum_0_2 = prd_0_2 + prd_1_2;
    prd_2_2 = pts_i2_2 * dir_2;

    prd_0_3 = pts_i0_3 * dir_0;
    prd_1_3 = pts_i1_3 * dir_1;
    sum_0_3 = prd_0_3 + prd_1_3;
    prd_2_3 = pts_i2_3 * dir_2;

    prd_0_4 = pts_i0_4 * dir_0;
    prd_1_4 = pts_i1_4 * dir_1;
    sum_0_4 = prd_0_4 + prd_1_4;
    prd_2_4 = pts_i2_4 * dir_2;

    prd_0_5 = pts_i0_5 * dir_0;
    prd_1_5 = pts_i1_5 * dir_1;
    sum_0_5 = prd_0_5 + prd_1_5;
    prd_2_5 = pts_i2_5 * dir_2;

    prd_0_6 = pts_i0_6 * dir_0;
    prd_1_6 = pts_i1_6 * dir_1;
    sum_0_6 = prd_0_6 + prd_1_6;
    prd_2_6 = pts_i2_6 * dir_2;

    prd_0_7 = pts_i0_7 * dir_0;
    prd_1_7 = pts_i1_7 * dir_1;
    sum_0_7 = prd_0_7 + prd_1_7;
    prd_2_7 = pts_i2_7 * dir_2;

    // Store
    dp = sum_0 + prd_2;
    dp_1 = sum_0_1 + prd_2_1;
    dp_2 = sum_0_2 + prd_2_2;

    dp_3 = sum_0_3 + prd_2_3;
    dp_4 = sum_0_4 + prd_2_4;
    dp_5 = sum_0_5 + prd_2_5;
    dp_6 = sum_0_6 + prd_2_6;
    dp_7 = sum_0_7 + prd_2_7;

    // Find the local dpmax and argmax
    if (dp > dpmax) {
      dpmax = dp;
      argmax = i;
    }

    if (dp_1 > dpmax_1) {
      dpmax_1 = dp_1;
      argmax_1 = i + 1;
    }

    if (dp_2 > dpmax_2) {
      dpmax_2 = dp_2;
      argmax_2 = i + 2;
    }

    if (dp_3 > dpmax_3) {
      dpmax_3 = dp_3;
      argmax_3 = i + 3;
    }

    if (dp_4 > dpmax_4) {
      dpmax_4 = dp_4;
      argmax_4 = i + 4;
    }

    if (dp_5 > dpmax_5) {
      dpmax_5 = dp_5;
      argmax_5 = i + 5;
    }

    if (dp_6 > dpmax_6) {
      dpmax_6 = dp_6;
      argmax_6 = i + 6;
    }

    if (dp_7 > dpmax_7) {
      dpmax_7 = dp_7;
      argmax_7 = i + 7;
    }
  }

  // Find the global dpmax and argmax
  if (dpmax < dpmax_1) {
    dpmax = dpmax_1;
    argmax = argmax_1;
  }

  if (dpmax < dpmax_2) {
    dpmax = dpmax_2;
    argmax = argmax_2;
  }

  if (dpmax < dpmax_3) {
    dpmax = dpmax_3;
    argmax = argmax_3;
  }

  if (dpmax < dpmax_4) {
    dpmax = dpmax_4;
    argmax = argmax_4;
  }

  if (dpmax < dpmax_5) {
    dpmax = dpmax_5;
    argmax = argmax_5;
  }

  if (dpmax < dpmax_6) {
    dpmax = dpmax_6;
    argmax = argmax_6;
  }

  if (dpmax < dpmax_7) {
    dpmax = dpmax_7;
    argmax = argmax_7;
  }

  // Handle the remaining cases
  for (; i < N; i++) {
    // Load
    pts_i0 = points[i][0];
    pts_i1 = points[i][1];
    pts_i2 = points[i][2];

    // Calculate
    prd_0 = pts_i0 * dir_0;
    prd_1 = pts_i1 * dir_1;
    prd_2 = pts_i2 * dir_2;
    sum_0 = prd_0 + prd_1;

    // Store
    dp = sum_0 + prd_2;

    if (dp > dpmax) {
      dpmax = dp;
      argmax = i;
    }
  }

  res[0] = points[argmax][0];
  res[1] = points[argmax][1];
  res[2] = points[argmax][2];
}

void invsup3D_optimized(const float (* points)[3], const float dir[3], float res[3], int N) {
  float dp, dp_1, dp_2, dp_3, dp_4, dp_5, dp_6, dp_7;
  float dpmax = 1e9, dpmax_1 = 1e9, dpmax_2 = 1e9,
    dpmax_3 = 1e9, dpmax_4 = 1e9, dpmax_5 = 1e9, dpmax_6 = 1e9, dpmax_7 = 1e9;
  int argmin = 0, argmin_1 = 0, argmin_2 = 0,
    argmin_3 = 0, argmin_4 = 0, argmin_5 = 0, argmin_6 = 0, argmin_7 = 0;

  // Scalar replacement variables for input points
  float pts_i0, pts_i1, pts_i2;
  float pts_i0_1, pts_i1_1, pts_i2_1;
  float pts_i0_2, pts_i1_2, pts_i2_2;
  float pts_i0_3, pts_i1_3, pts_i2_3;
  float pts_i0_4, pts_i1_4, pts_i2_4;
  float pts_i0_5, pts_i1_5, pts_i2_5;
  float pts_i0_6, pts_i1_6, pts_i2_6;
  float pts_i0_7, pts_i1_7, pts_i2_7;

  // Scalar replacement variables for products and summations in dot products
  float prd_0, prd_1, prd_2, sum_0;
  float prd_0_1, prd_1_1, prd_2_1, sum_0_1;
  float prd_0_2, prd_1_2, prd_2_2, sum_0_2;
  float prd_0_3, prd_1_3, prd_2_3, sum_0_3;
  float prd_0_4, prd_1_4, prd_2_4, sum_0_4;
  float prd_0_5, prd_1_5, prd_2_5, sum_0_5;
  float prd_0_6, prd_1_6, prd_2_6, sum_0_6;
  float prd_0_7, prd_1_7, prd_2_7, sum_0_7;

  // Reduce memory access by factoring out locations accessed in every iteration.
  float dir_0 = dir[0];
  float dir_1 = dir[1];
  float dir_2 = dir[2];

  int i;
  // Loop unroll 4 times 
  for (i = 0; i < N - 3; i += 4) {
    // Load
    pts_i0 = points[i][0];
    pts_i1 = points[i][1];
    pts_i2 = points[i][2];

    pts_i0_1 = points[i + 1][0];
    pts_i1_1 = points[i + 1][1];
    pts_i2_1 = points[i + 1][2];

    pts_i0_2 = points[i + 2][0];
    pts_i1_2 = points[i + 2][1];
    pts_i2_2 = points[i + 2][2];

    pts_i0_3 = points[i + 3][0];
    pts_i1_3 = points[i + 3][1];
    pts_i2_3 = points[i + 3][2];

    // Compute
    prd_0 = pts_i0 * dir_0;
    prd_1 = pts_i1 * dir_1;
    sum_0 = prd_0 + prd_1;
    prd_2 = pts_i2 * dir_2;

    prd_0_1 = pts_i0_1 * dir_0;
    prd_1_1 = pts_i1_1 * dir_1;
    sum_0_1 = prd_0_1 + prd_1_1;
    prd_2_1 = pts_i2_1 * dir_2;

    prd_0_2 = pts_i0_2 * dir_0;
    prd_1_2 = pts_i1_2 * dir_1;
    sum_0_2 = prd_0_2 + prd_1_2;
    prd_2_2 = pts_i2_2 * dir_2;

    prd_0_3 = pts_i0_3 * dir_0;
    prd_1_3 = pts_i1_3 * dir_1;
    sum_0_3 = prd_0_3 + prd_1_3;
    prd_2_3 = pts_i2_3 * dir_2;

    // Store
    dp = sum_0 + prd_2;
    dp_1 = sum_0_1 + prd_2_1;
    dp_2 = sum_0_2 + prd_2_2;
    dp_3 = sum_0_3 + prd_2_3;

    // Compute local dpmax and argmax
    if (dp < dpmax) {
      dpmax = dp;
      argmin = i;
    }

    if (dp_1 < dpmax_1) {
      dpmax_1 = dp_1;
      argmin_1 = i + 1;
    }

    if (dp_2 < dpmax_2) {
      dpmax_2 = dp_2;
      argmin_2 = i + 2;
    }

    if (dp_3 < dpmax_3) {
      dpmax_3 = dp_3;
      argmin_3 = i + 3;
    }
  }

  // Find the dpmax and argmax
  if (dpmax > dpmax_1) {
    dpmax = dpmax_1;
    argmin = argmin_1;
  }

  if (dpmax > dpmax_2) {
    dpmax = dpmax_2;
    argmin = argmin_2;
  }

  if (dpmax > dpmax_3) {
    dpmax = dpmax_3;
    argmin = argmin_3;
  }

  // Handle the remaining cases
  for (; i < N; i++) {
    // Load
    pts_i0 = points[i][0];
    pts_i1 = points[i][1];
    pts_i2 = points[i][2];

    // Calculate
    prd_0 = pts_i0 * dir_0;
    prd_1 = pts_i1 * dir_1;
    prd_2 = pts_i2 * dir_2;
    sum_0 = prd_0 + prd_1;

    // Store
    dp = sum_0 + prd_2;

    if (dp < dpmax) {
      dpmax = dp;
      argmin = i;
    }
  }

  res[0] -= points[argmin][0];
  res[1] -= points[argmin][1];
  res[2] -= points[argmin][2];
}

int do_intersect3D_optimized(const struct CHObject* obj1, const struct CHObject* obj2) {
  // The search direction
  float d[3] = {1.f, 1.f, 1.f};  // could also be random
  // The point we examine in each iteration
  float a[3];

  // Our simplex that do_simplex modifies.
  struct Simplex3D s;
  s.num_points = 1;

  // S = Support(A-B) = Support(A) - Invsup(B)
  support3D_optimized(obj1->points, d, s.p[0], obj1->num_points);
  invsup3D_optimized(obj2->points, d, s.p[0], obj2->num_points);

  // d = -S
  d[0] = -s.p[0][0];
  d[1] = -s.p[0][1];
  d[2] = -s.p[0][2];

  int max_iter = 100;  // most likely we need less
  while (max_iter--) {
    // a = Support(A-B) = Support(A) - Invsup(B)
    support3D_optimized(obj1->points, d, a, obj1->num_points);
    invsup3D_optimized(obj2->points, d, a, obj2->num_points);

    // Are we past the origin? If no, we will never enclose it.
    if (dot3D_optimized(a, d) < 0.) return 0;

    // Add newly found point to Simplex
    // AT THIS POINT WE KNOW a IS PAST THE ORIGIN
    s.p[s.num_points][0] = a[0];
    s.p[s.num_points][1] = a[1];
    s.p[s.num_points][2] = a[2];
    s.num_points++;

    if (do_simplex3D_optimized(&s, d)) return 1;
  }

  return 1; // Most likely not reachable. Say we intersect by default.
}

int do_intersect3D_optimized_inlined(const struct CHObject* obj1, const struct CHObject* obj2) {
  // The search direction
  float d[3] = {1.f, 1.f, 1.f};  // could also be random
  // The point we examine in each iteration
  float a[3];

  // Our simplex that do_simplex modifies.
  struct Simplex3D s;
  s.num_points = 1;

  // Method inlining for support function
  int no_points_min = obj1->num_points > obj2->num_points ? obj2->num_points : obj1->num_points;

  float dp, dp_1, dp_2, dpmax = -1e9;
  int i, argmax = 0;

  float invsup_dp, invsup_dp_1, invsup_dp_2, invsup_dpmax = 1e9;
  int argmin = 0;

  float dir_0;
  float dir_1;
  float dir_2;

  float pts_i0, pts_i1, pts_i2;
  float pts_i0_1, pts_i1_1, pts_i2_1;
  float pts_i0_2, pts_i1_2, pts_i2_2;

  float invsup_pts_i0, invsup_pts_i1, invsup_pts_i2;
  float invsup_pts_i0_1, invsup_pts_i1_1, invsup_pts_i2_1;
  float invsup_pts_i0_2, invsup_pts_i1_2, invsup_pts_i2_2;

  float prd_0, prd_1, prd_2, sum_0;
  float prd_0_1, prd_1_1, prd_2_1, sum_0_1;
  float prd_0_2, prd_1_2, prd_2_2, sum_0_2;

  float invsup_prd_0, invsup_prd_1, invsup_prd_2, invsup_sum_0;
  float invsup_prd_0_1, invsup_prd_1_1, invsup_prd_2_1, invsup_sum_0_1;
  float invsup_prd_0_2, invsup_prd_1_2, invsup_prd_2_2, invsup_sum_0_2;

  dir_0 = d[0];
  dir_1 = d[1];
  dir_2 = d[2];

  // Joint loop that calculates the dot product and argmax of both points at the same time.
  for (i = 0; i < no_points_min - 2; i += 3) {
    // Load
    pts_i0 = obj1->points[i][0];
    pts_i1 = obj1->points[i][1];
    pts_i2 = obj1->points[i][2];

    pts_i0_1 = obj1->points[i + 1][0];
    pts_i1_1 = obj1->points[i + 1][1];
    pts_i2_1 = obj1->points[i + 1][2];

    pts_i0_2 = obj1->points[i + 2][0];
    pts_i1_2 = obj1->points[i + 2][1];
    pts_i2_2 = obj1->points[i + 2][2];

    invsup_pts_i0 = obj2->points[i][0];
    invsup_pts_i1 = obj2->points[i][1];
    invsup_pts_i2 = obj2->points[i][2];

    invsup_pts_i0_1 = obj2->points[i + 1][0];
    invsup_pts_i1_1 = obj2->points[i + 1][1];
    invsup_pts_i2_1 = obj2->points[i + 1][2];

    invsup_pts_i0_2 = obj2->points[i + 2][0];
    invsup_pts_i1_2 = obj2->points[i + 2][1];
    invsup_pts_i2_2 = obj2->points[i + 2][2];

    // Calculate
    prd_0 = pts_i0 * dir_0;
    prd_1 = pts_i1 * dir_1;
    sum_0 = prd_0 + prd_1;
    prd_2 = pts_i2 * dir_2;

    prd_0_1 = pts_i0_1 * dir_0;
    prd_1_1 = pts_i1_1 * dir_1;
    sum_0_1 = prd_0_1 + prd_1_1;
    prd_2_1 = pts_i2_1 * dir_2;

    prd_0_2 = pts_i0_2 * dir_0;
    prd_1_2 = pts_i1_2 * dir_1;
    sum_0_2 = prd_0_2 + prd_1_2;
    prd_2_2 = pts_i2_2 * dir_2;

    invsup_prd_0 = invsup_pts_i0 * dir_0;
    invsup_prd_1 = invsup_pts_i1 * dir_1;
    invsup_sum_0 = invsup_prd_0 + invsup_prd_1;
    invsup_prd_2 = invsup_pts_i2 * dir_2;

    invsup_prd_0_1 = invsup_pts_i0_1 * dir_0;
    invsup_prd_1_1 = invsup_pts_i1_1 * dir_1;
    invsup_sum_0_1 = invsup_prd_0_1 + invsup_prd_1_1;
    invsup_prd_2_1 = invsup_pts_i2_1 * dir_2;

    invsup_prd_0_2 = invsup_pts_i0_2 * dir_0;
    invsup_prd_1_2 = invsup_pts_i1_2 * dir_1;
    invsup_sum_0_2 = invsup_prd_0_2 + invsup_prd_1_2;
    invsup_prd_2_2 = invsup_pts_i2_2 * dir_2;

    // Store
    dp = sum_0 + prd_2;
    dp_1 = sum_0_1 + prd_2_1;
    dp_2 = sum_0_2 + prd_2_2;

    invsup_dp = invsup_sum_0 + invsup_prd_2;
    invsup_dp_1 = invsup_sum_0_1 + invsup_prd_2_1;
    invsup_dp_2 = invsup_sum_0_2 + invsup_prd_2_2;

    // Find the dpmax and argmax
    if (dp > dpmax) {
      dpmax = dp;
      argmax = i;
    }

    if (dp_1 > dpmax) {
      dpmax = dp_1;
      argmax = i + 1;
    }

    if (dp_2 > dpmax) {
      dpmax = dp_2;
      argmax = i + 2;
    }

    // Find dpmin and argmin
    if (invsup_dp < invsup_dpmax) {
      invsup_dpmax = invsup_dp;
      argmin = i;
    }

    if (invsup_dp_1 < invsup_dpmax) {
      invsup_dpmax = invsup_dp_1;
      argmin = i + 1;
    }

    if (invsup_dp_2 < invsup_dpmax) {
      invsup_dpmax = invsup_dp_2;
      argmin = i + 2;
    }
  }

  int j = i;
  // Iterate over the remaining points of first object.
  for (; i < obj1->num_points - 2; i += 3) {
    // Load
    pts_i0 = obj1->points[i][0];
    pts_i1 = obj1->points[i][1];
    pts_i2 = obj1->points[i][2];

    pts_i0_1 = obj1->points[i + 1][0];
    pts_i1_1 = obj1->points[i + 1][1];
    pts_i2_1 = obj1->points[i + 1][2];

    pts_i0_2 = obj1->points[i + 2][0];
    pts_i1_2 = obj1->points[i + 2][1];
    pts_i2_2 = obj1->points[i + 2][2];

    // Compute
    prd_0 = pts_i0 * dir_0;
    prd_1 = pts_i1 * dir_1;
    sum_0 = prd_0 + prd_1;
    prd_2 = pts_i2 * dir_2;

    prd_0_1 = pts_i0_1 * dir_0;
    prd_1_1 = pts_i1_1 * dir_1;
    sum_0_1 = prd_0_1 + prd_1_1;
    prd_2_1 = pts_i2_1 * dir_2;

    prd_0_2 = pts_i0_2 * dir_0;
    prd_1_2 = pts_i1_2 * dir_1;
    sum_0_2 = prd_0_2 + prd_1_2;
    prd_2_2 = pts_i2_2 * dir_2;

    // Store
    dp = sum_0 + prd_2;
    dp_1 = sum_0_1 + prd_2_1;
    dp_2 = sum_0_2 + prd_2_2;

    // Find the dpmax and argmax
    if (dp > dpmax) {
      dpmax = dp;
      argmax = i;
    }

    if (dp_1 > dpmax) {
      dpmax = dp_1;
      argmax = i + 1;
    }

    if (dp_2 > dpmax) {
      dpmax = dp_2;
      argmax = i + 2;
    }
  }

  // Handle the remaining cases
  for (; i < obj1->num_points; i++) {
    // Load
    pts_i0 = obj1->points[i][0];
    pts_i1 = obj1->points[i][1];
    pts_i2 = obj1->points[i][2];

    // Calculate
    prd_0 = pts_i0 * dir_0;
    prd_1 = pts_i1 * dir_1;
    prd_2 = pts_i2 * dir_2;
    sum_0 = prd_0 + prd_1;

    // Store
    dp = sum_0 + prd_2;

    if (dp > dpmax) {
      dpmax = dp;
      argmax = i;
    }
  }

  // Iterate over the remaining points of second object.
  for (; j < obj2->num_points - 2; j += 3) {
    // Load
    invsup_pts_i0 = obj2->points[j][0];
    invsup_pts_i1 = obj2->points[j][1];
    invsup_pts_i2 = obj2->points[j][2];

    invsup_pts_i0_1 = obj2->points[j + 1][0];
    invsup_pts_i1_1 = obj2->points[j + 1][1];
    invsup_pts_i2_1 = obj2->points[j + 1][2];

    invsup_pts_i0_2 = obj2->points[j + 2][0];
    invsup_pts_i1_2 = obj2->points[j + 2][1];
    invsup_pts_i2_2 = obj2->points[j + 2][2];

    // Compute
    invsup_prd_0 = invsup_pts_i0 * dir_0;
    invsup_prd_1 = invsup_pts_i1 * dir_1;
    invsup_sum_0 = invsup_prd_0 + invsup_prd_1;
    invsup_prd_2 = invsup_pts_i2 * dir_2;

    invsup_prd_0_1 = invsup_pts_i0_1 * dir_0;
    invsup_prd_1_1 = invsup_pts_i1_1 * dir_1;
    invsup_sum_0_1 = invsup_prd_0_1 + invsup_prd_1_1;
    invsup_prd_2_1 = invsup_pts_i2_1 * dir_2;

    invsup_prd_0_2 = invsup_pts_i0_2 * dir_0;
    invsup_prd_1_2 = invsup_pts_i1_2 * dir_1;
    invsup_sum_0_2 = invsup_prd_0_2 + invsup_prd_1_2;
    invsup_prd_2_2 = invsup_pts_i2_2 * dir_2;

    // Store
    invsup_dp = invsup_sum_0 + invsup_prd_2;
    invsup_dp_1 = invsup_sum_0_1 + invsup_prd_2_1;
    invsup_dp_2 = invsup_sum_0_2 + invsup_prd_2_2;

    // Find the dpmin and argmin
    if (invsup_dp < invsup_dpmax) {
      invsup_dpmax = invsup_dp;
      argmin = j;
    }

    if (invsup_dp_1 < invsup_dpmax) {
      invsup_dpmax = invsup_dp_1;
      argmin = j + 1;
    }

    if (invsup_dp_2 < invsup_dpmax) {
      invsup_dpmax = invsup_dp_2;
      argmin = j + 2;
    }
  }

  // Handle the remaining cases
  for (; j < obj2->num_points; j++) {
    // Load
    invsup_pts_i0 = obj2->points[j][0];
    invsup_pts_i1 = obj2->points[j][1];
    invsup_pts_i2 = obj2->points[j][2];

    // Calculate
    invsup_prd_0 = invsup_pts_i0 * dir_0;
    invsup_prd_1 = invsup_pts_i1 * dir_1;
    invsup_prd_2 = invsup_pts_i2 * dir_2;
    invsup_sum_0 = invsup_prd_0 + invsup_prd_1;

    // Store
    invsup_dp = invsup_sum_0 + invsup_prd_2;

    if (invsup_dp < invsup_dpmax) {
      invsup_dpmax = invsup_dp;
      argmin = j;
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

    // Reset values
    dpmax = -1e9;
    invsup_dpmax = 1e9;
    argmax = 0;
    argmin = 0;

    // Set new search direction
    dir_0 = d[0];
    dir_1 = d[1];
    dir_2 = d[2];
    
    // Joint loop that iterates over both objects.
    for (i = 0; i < no_points_min - 2; i += 3) {
      // Load
      pts_i0 = obj1->points[i][0];
      pts_i1 = obj1->points[i][1];
      pts_i2 = obj1->points[i][2];

      pts_i0_1 = obj1->points[i + 1][0];
      pts_i1_1 = obj1->points[i + 1][1];
      pts_i2_1 = obj1->points[i + 1][2];

      pts_i0_2 = obj1->points[i + 2][0];
      pts_i1_2 = obj1->points[i + 2][1];
      pts_i2_2 = obj1->points[i + 2][2];

      invsup_pts_i0 = obj2->points[i][0];
      invsup_pts_i1 = obj2->points[i][1];
      invsup_pts_i2 = obj2->points[i][2];

      invsup_pts_i0_1 = obj2->points[i + 1][0];
      invsup_pts_i1_1 = obj2->points[i + 1][1];
      invsup_pts_i2_1 = obj2->points[i + 1][2];

      invsup_pts_i0_2 = obj2->points[i + 2][0];
      invsup_pts_i1_2 = obj2->points[i + 2][1];
      invsup_pts_i2_2 = obj2->points[i + 2][2];

      // Compute
      prd_0 = pts_i0 * dir_0;
      prd_1 = pts_i1 * dir_1;
      sum_0 = prd_0 + prd_1;
      prd_2 = pts_i2 * dir_2;

      prd_0_1 = pts_i0_1 * dir_0;
      prd_1_1 = pts_i1_1 * dir_1;
      sum_0_1 = prd_0_1 + prd_1_1;
      prd_2_1 = pts_i2_1 * dir_2;

      prd_0_2 = pts_i0_2 * dir_0;
      prd_1_2 = pts_i1_2 * dir_1;
      sum_0_2 = prd_0_2 + prd_1_2;
      prd_2_2 = pts_i2_2 * dir_2;

      invsup_prd_0 = invsup_pts_i0 * dir_0;
      invsup_prd_1 = invsup_pts_i1 * dir_1;
      invsup_sum_0 = invsup_prd_0 + invsup_prd_1;
      invsup_prd_2 = invsup_pts_i2 * dir_2;

      invsup_prd_0_1 = invsup_pts_i0_1 * dir_0;
      invsup_prd_1_1 = invsup_pts_i1_1 * dir_1;
      invsup_sum_0_1 = invsup_prd_0_1 + invsup_prd_1_1;
      invsup_prd_2_1 = invsup_pts_i2_1 * dir_2;

      invsup_prd_0_2 = invsup_pts_i0_2 * dir_0;
      invsup_prd_1_2 = invsup_pts_i1_2 * dir_1;
      invsup_sum_0_2 = invsup_prd_0_2 + invsup_prd_1_2;
      invsup_prd_2_2 = invsup_pts_i2_2 * dir_2;

      // Store
      dp = sum_0 + prd_2;
      dp_1 = sum_0_1 + prd_2_1;
      dp_2 = sum_0_2 + prd_2_2;

      invsup_dp = invsup_sum_0 + invsup_prd_2;
      invsup_dp_1 = invsup_sum_0_1 + invsup_prd_2_1;
      invsup_dp_2 = invsup_sum_0_2 + invsup_prd_2_2;

      // Find local dpmax and argmax
      if (dp > dpmax) {
        dpmax = dp;
        argmax = i;
      }

      if (dp_1 > dpmax) {
        dpmax = dp_1;
        argmax = i + 1;
      }

      if (dp_2 > dpmax) {
        dpmax = dp_2;
        argmax = i + 2;
      }

      // Find dpmin and argmin
      if (invsup_dp < invsup_dpmax) {
        invsup_dpmax = invsup_dp;
        argmin = i;
      }

      if (invsup_dp_1 < invsup_dpmax) {
        invsup_dpmax = invsup_dp_1;
        argmin = i + 1;
      }

      if (invsup_dp_2 < invsup_dpmax) {
        invsup_dpmax = invsup_dp_2;
        argmin = i + 2;
      }
    }

    int j = i;
    // Finish remaining points of the first object
    for (; i < obj1->num_points - 2; i += 3) {
      // Load
      pts_i0 = obj1->points[i][0];
      pts_i1 = obj1->points[i][1];
      pts_i2 = obj1->points[i][2];

      pts_i0_1 = obj1->points[i + 1][0];
      pts_i1_1 = obj1->points[i + 1][1];
      pts_i2_1 = obj1->points[i + 1][2];

      pts_i0_2 = obj1->points[i + 2][0];
      pts_i1_2 = obj1->points[i + 2][1];
      pts_i2_2 = obj1->points[i + 2][2];

      // Calculate
      prd_0 = pts_i0 * dir_0;
      prd_1 = pts_i1 * dir_1;
      sum_0 = prd_0 + prd_1;
      prd_2 = pts_i2 * dir_2;

      prd_0_1 = pts_i0_1 * dir_0;
      prd_1_1 = pts_i1_1 * dir_1;
      sum_0_1 = prd_0_1 + prd_1_1;
      prd_2_1 = pts_i2_1 * dir_2;

      prd_0_2 = pts_i0_2 * dir_0;
      prd_1_2 = pts_i1_2 * dir_1;
      sum_0_2 = prd_0_2 + prd_1_2;
      prd_2_2 = pts_i2_2 * dir_2;

      // Store
      dp = sum_0 + prd_2;
      dp_1 = sum_0_1 + prd_2_1;
      dp_2 = sum_0_2 + prd_2_2;

      // Find dpmax and argmax
      if (dp > dpmax) {
        dpmax = dp;
        argmax = i;
      }

      if (dp_1 > dpmax) {
        dpmax = dp_1;
        argmax = i + 1;
      }

      if (dp_2 > dpmax) {
        dpmax = dp_2;
        argmax = i + 2;
      }
    }

    // Handle the remaining cases
    for (; i < obj1->num_points; i++) {
      // Load
      pts_i0 = obj1->points[i][0];
      pts_i1 = obj1->points[i][1];
      pts_i2 = obj1->points[i][2];

      // Calculate
      prd_0 = pts_i0 * dir_0;
      prd_1 = pts_i1 * dir_1;
      prd_2 = pts_i2 * dir_2;
      sum_0 = prd_0 + prd_1;

      // Store
      dp = sum_0 + prd_2;

      if (dp > dpmax) {
        dpmax = dp;
        argmax = i;
      }
    }

    // Finish remaining points of the first object
    for (; j < obj2->num_points - 2; j += 3) {
      // Load
      invsup_pts_i0 = obj2->points[j][0];
      invsup_pts_i1 = obj2->points[j][1];
      invsup_pts_i2 = obj2->points[j][2];

      invsup_pts_i0_1 = obj2->points[j + 1][0];
      invsup_pts_i1_1 = obj2->points[j + 1][1];
      invsup_pts_i2_1 = obj2->points[j + 1][2];

      invsup_pts_i0_2 = obj2->points[j + 2][0];
      invsup_pts_i1_2 = obj2->points[j + 2][1];
      invsup_pts_i2_2 = obj2->points[j + 2][2];

      // Compute
      invsup_prd_0 = invsup_pts_i0 * dir_0;
      invsup_prd_1 = invsup_pts_i1 * dir_1;
      invsup_sum_0 = invsup_prd_0 + invsup_prd_1;
      invsup_prd_2 = invsup_pts_i2 * dir_2;

      invsup_prd_0_1 = invsup_pts_i0_1 * dir_0;
      invsup_prd_1_1 = invsup_pts_i1_1 * dir_1;
      invsup_sum_0_1 = invsup_prd_0_1 + invsup_prd_1_1;
      invsup_prd_2_1 = invsup_pts_i2_1 * dir_2;

      invsup_prd_0_2 = invsup_pts_i0_2 * dir_0;
      invsup_prd_1_2 = invsup_pts_i1_2 * dir_1;
      invsup_sum_0_2 = invsup_prd_0_2 + invsup_prd_1_2;
      invsup_prd_2_2 = invsup_pts_i2_2 * dir_2;

      // Store
      invsup_dp = invsup_sum_0 + invsup_prd_2;
      invsup_dp_1 = invsup_sum_0_1 + invsup_prd_2_1;
      invsup_dp_2 = invsup_sum_0_2 + invsup_prd_2_2;

      // Find local dpmin and argmin
      if (invsup_dp < invsup_dpmax) {
        invsup_dpmax = invsup_dp;
        argmin = j;
      }

      if (invsup_dp_1 < invsup_dpmax) {
        invsup_dpmax = invsup_dp_1;
        argmin = j + 1;
      }

      if (invsup_dp_2 < invsup_dpmax) {
        invsup_dpmax = invsup_dp_2;
        argmin = j + 2;
      }
    }

    // Handle the remaining cases
    for (; j < obj2->num_points; j++) {
      // Load
      invsup_pts_i0 = obj2->points[j][0];
      invsup_pts_i1 = obj2->points[j][1];
      invsup_pts_i2 = obj2->points[j][2];

      // Calculate
      invsup_prd_0 = invsup_pts_i0 * dir_0;
      invsup_prd_1 = invsup_pts_i1 * dir_1;
      invsup_prd_2 = invsup_pts_i2 * dir_2;
      invsup_sum_0 = invsup_prd_0 + invsup_prd_1;

      // Store
      invsup_dp = invsup_sum_0 + invsup_prd_2;

      if (invsup_dp < invsup_dpmax) {
        invsup_dpmax = invsup_dp;
        argmin = j;
      }
    }

    // a = Support(A-B) = Support(A) - Invsup(B)
    a[0] = obj1->points[argmax][0] - obj2->points[argmin][0];
    a[1] = obj1->points[argmax][1] - obj2->points[argmin][1];
    a[2] = obj1->points[argmax][2] - obj2->points[argmin][2];

    // Are we past the origin? If no, we will never enclose it.
    if (dot3D_optimized(a, d) < 0.) return 0;

    // Add newly found point to Simplex
    // AT THIS POINT WE KNOW a IS PAST THE ORIGIN
    s.p[s.num_points][0] = a[0];
    s.p[s.num_points][1] = a[1];
    s.p[s.num_points][2] = a[2];
    s.num_points++;

    if (do_simplex3D_optimized(&s, d)) return 1;
  }

  return 1; // Most likely not reachable. Say we intersect by default.
}

// Loop unroll 8 times.
int do_intersect3D_o_i_lu8(const struct CHObject* obj1, const struct CHObject* obj2) {
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

  // Variables for holding the current and maximum dot products of first object.
  float dp, dp_1, dp_2, dp_3, dp_4, dp_5, dp_6, dp_7;
  float dpmax = -1e9, dpmax_1 = -1e9, dpmax_2 = -1e9,
    dpmax_3 = -1e9, dpmax_4 = -1e9, dpmax_5 = -1e9, dpmax_6 = -1e9, dpmax_7 = -1e9;
  int argmax = 0, argmax_1 = 0, argmax_2 = 0,
    argmax_3 = 0, argmax_4 = 0, argmax_5 = 0, argmax_6 = 0, argmax_7 = 0;

  // Variables for holding the current and maximum dot products of second object.
  float invsup_dp, invsup_dp_1, invsup_dp_2,
    invsup_dp_3, invsup_dp_4, invsup_dp_5, invsup_dp_6, invsup_dp_7;
  float invsup_dpmax = 1e9, invsup_dpmax_1 = 1e9, invsup_dpmax_2 = 1e9,
    invsup_dpmax_3 = 1e9, invsup_dpmax_4 = 1e9, invsup_dpmax_5 = 1e9,
    invsup_dpmax_6 = 1e9, invsup_dpmax_7 = 1e9;
  int argmin = 0, argmin_1 = 0, argmin_2 = 0,
    argmin_3 = 0, argmin_4 = 0, argmin_5 = 0, argmin_6 = 0, argmin_7 = 0;

  // Vertices of object 1 accessed during one iteration of the loop.
  float pts_i0, pts_i1, pts_i2;
  float pts_i0_1, pts_i1_1, pts_i2_1;
  float pts_i0_2, pts_i1_2, pts_i2_2;
  float pts_i0_3, pts_i1_3, pts_i2_3;

  // Not used in the joint loop
  float pts_i0_4, pts_i1_4, pts_i2_4;
  float pts_i0_5, pts_i1_5, pts_i2_5;
  float pts_i0_6, pts_i1_6, pts_i2_6;
  float pts_i0_7, pts_i1_7, pts_i2_7;

  // Vertices of object 2 accessed during one iteration of the loop.
  float invsup_pts_i0, invsup_pts_i1, invsup_pts_i2;
  float invsup_pts_i0_1, invsup_pts_i1_1, invsup_pts_i2_1;
  float invsup_pts_i0_2, invsup_pts_i1_2, invsup_pts_i2_2;
  float invsup_pts_i0_3, invsup_pts_i1_3, invsup_pts_i2_3;

  // Not used in the joint loop
  float invsup_pts_i0_4, invsup_pts_i1_4, invsup_pts_i2_4;
  float invsup_pts_i0_5, invsup_pts_i1_5, invsup_pts_i2_5;
  float invsup_pts_i0_6, invsup_pts_i1_6, invsup_pts_i2_6;
  float invsup_pts_i0_7, invsup_pts_i1_7, invsup_pts_i2_7;

  // Variables that hold the values of individual operations for dot product (object 1).
  float prd_0, prd_1, prd_2, sum_0;
  float prd_0_1, prd_1_1, prd_2_1, sum_0_1;
  float prd_0_2, prd_1_2, prd_2_2, sum_0_2;
  float prd_0_3, prd_1_3, prd_2_3, sum_0_3;

  // Not used during joint loop.
  float prd_0_4, prd_1_4, prd_2_4, sum_0_4;
  float prd_0_5, prd_1_5, prd_2_5, sum_0_5;
  float prd_0_6, prd_1_6, prd_2_6, sum_0_6;
  float prd_0_7, prd_1_7, prd_2_7, sum_0_7;

  // Variables that hold the values of individual operations for dot product (object 2).
  float invsup_prd_0, invsup_prd_1, invsup_prd_2, invsup_sum_0;
  float invsup_prd_0_1, invsup_prd_1_1, invsup_prd_2_1, invsup_sum_0_1;
  float invsup_prd_0_2, invsup_prd_1_2, invsup_prd_2_2, invsup_sum_0_2;
  float invsup_prd_0_3, invsup_prd_1_3, invsup_prd_2_3, invsup_sum_0_3;

  // Not used during joint loop.
  float invsup_prd_0_4, invsup_prd_1_4, invsup_prd_2_4, invsup_sum_0_4;
  float invsup_prd_0_5, invsup_prd_1_5, invsup_prd_2_5, invsup_sum_0_5;
  float invsup_prd_0_6, invsup_prd_1_6, invsup_prd_2_6, invsup_sum_0_6;
  float invsup_prd_0_7, invsup_prd_1_7, invsup_prd_2_7, invsup_sum_0_7;

  // Scalar replacement variables for the direction array to reduce memory access.
  float dir_0 = d[0];
  float dir_1 = d[1];
  float dir_2 = d[2];

  int i;
  // Joint loop that iterates both object 1 and object 2 at the same time to reduce total no. of iterations.
  for (i = 0; i < no_points_min - 3; i += 4) {
    // Load first object
    pts_i0 = obj1->points[i][0];
    pts_i1 = obj1->points[i][1];
    pts_i2 = obj1->points[i][2];

    pts_i0_1 = obj1->points[i + 1][0];
    pts_i1_1 = obj1->points[i + 1][1];
    pts_i2_1 = obj1->points[i + 1][2];

    pts_i0_2 = obj1->points[i + 2][0];
    pts_i1_2 = obj1->points[i + 2][1];
    pts_i2_2 = obj1->points[i + 2][2];

    pts_i0_3 = obj1->points[i + 3][0];
    pts_i1_3 = obj1->points[i + 3][1];
    pts_i2_3 = obj1->points[i + 3][2];

    // Load second object
    invsup_pts_i0 = obj2->points[i][0];
    invsup_pts_i1 = obj2->points[i][1];
    invsup_pts_i2 = obj2->points[i][2];

    invsup_pts_i0_1 = obj2->points[i + 1][0];
    invsup_pts_i1_1 = obj2->points[i + 1][1];
    invsup_pts_i2_1 = obj2->points[i + 1][2];

    invsup_pts_i0_2 = obj2->points[i + 2][0];
    invsup_pts_i1_2 = obj2->points[i + 2][1];
    invsup_pts_i2_2 = obj2->points[i + 2][2];

    invsup_pts_i0_3 = obj2->points[i + 3][0];
    invsup_pts_i1_3 = obj2->points[i + 3][1];
    invsup_pts_i2_3 = obj2->points[i + 3][2];

    // Calculate support
    prd_0 = pts_i0 * dir_0;
    prd_1 = pts_i1 * dir_1;
    sum_0 = prd_0 + prd_1;
    prd_2 = pts_i2 * dir_2;

    prd_0_1 = pts_i0_1 * dir_0;
    prd_1_1 = pts_i1_1 * dir_1;
    sum_0_1 = prd_0_1 + prd_1_1;
    prd_2_1 = pts_i2_1 * dir_2;

    prd_0_2 = pts_i0_2 * dir_0;
    prd_1_2 = pts_i1_2 * dir_1;
    sum_0_2 = prd_0_2 + prd_1_2;
    prd_2_2 = pts_i2_2 * dir_2;

    prd_0_3 = pts_i0_3 * dir_0;
    prd_1_3 = pts_i1_3 * dir_1;
    sum_0_3 = prd_0_3 + prd_1_3;
    prd_2_3 = pts_i2_3 * dir_2;

    // Calculate the inverse support
    invsup_prd_0 = invsup_pts_i0 * dir_0;
    invsup_prd_1 = invsup_pts_i1 * dir_1;
    invsup_sum_0 = invsup_prd_0 + invsup_prd_1;
    invsup_prd_2 = invsup_pts_i2 * dir_2;

    invsup_prd_0_1 = invsup_pts_i0_1 * dir_0;
    invsup_prd_1_1 = invsup_pts_i1_1 * dir_1;
    invsup_sum_0_1 = invsup_prd_0_1 + invsup_prd_1_1;
    invsup_prd_2_1 = invsup_pts_i2_1 * dir_2;

    invsup_prd_0_2 = invsup_pts_i0_2 * dir_0;
    invsup_prd_1_2 = invsup_pts_i1_2 * dir_1;
    invsup_sum_0_2 = invsup_prd_0_2 + invsup_prd_1_2;
    invsup_prd_2_2 = invsup_pts_i2_2 * dir_2;

    invsup_prd_0_3 = invsup_pts_i0_3 * dir_0;
    invsup_prd_1_3 = invsup_pts_i1_3 * dir_1;
    invsup_sum_0_3 = invsup_prd_0_3 + invsup_prd_1_3;
    invsup_prd_2_3 = invsup_pts_i2_3 * dir_2;

    // Store first object's dot products
    dp = sum_0 + prd_2;
    dp_1 = sum_0_1 + prd_2_1;
    dp_2 = sum_0_2 + prd_2_2;
    dp_3 = sum_0_3 + prd_2_3;

    // Store second object's dot products
    invsup_dp = invsup_sum_0 + invsup_prd_2;
    invsup_dp_1 = invsup_sum_0_1 + invsup_prd_2_1;
    invsup_dp_2 = invsup_sum_0_2 + invsup_prd_2_2;
    invsup_dp_3 = invsup_sum_0_3 + invsup_prd_2_3;

    // Find local dpmax and argmax
    if (dp > dpmax) {
      dpmax = dp;
      argmax = i;
    }

    if (dp_1 > dpmax_1) {
      dpmax_1 = dp_1;
      argmax_1 = i + 1;
    }

    if (dp_2 > dpmax_2) {
      dpmax_2 = dp_2;
      argmax_2 = i + 2;
    }

    if (dp_3 > dpmax_3) {
      dpmax_3 = dp_3;
      argmax_3 = i + 3;
    }

    // Find local dpmin and argmin
    if (invsup_dp < invsup_dpmax) {
      invsup_dpmax = invsup_dp;
      argmin = i;
    }

    if (invsup_dp_1 < invsup_dpmax_1) {
      invsup_dpmax_1 = invsup_dp_1;
      argmin_1 = i + 1;
    }

    if (invsup_dp_2 < invsup_dpmax_2) {
      invsup_dpmax_2 = invsup_dp_2;
      argmin_2 = i + 2;
    }

    if (invsup_dp_3 < invsup_dpmax_3) {
      invsup_dpmax_3 = invsup_dp_3;
      argmin_3 = i + 3;
    }
  }

  // Joint loop ended - Complete the search for support of first object if it wasn't completed
  int j = i;
  for (; i < obj1->num_points - 7; i += 8) {
    // Load
    pts_i0 = obj1->points[i][0];
    pts_i1 = obj1->points[i][1];
    pts_i2 = obj1->points[i][2];

    pts_i0_1 = obj1->points[i + 1][0];
    pts_i1_1 = obj1->points[i + 1][1];
    pts_i2_1 = obj1->points[i + 1][2];

    pts_i0_2 = obj1->points[i + 2][0];
    pts_i1_2 = obj1->points[i + 2][1];
    pts_i2_2 = obj1->points[i + 2][2];

    pts_i0_3 = obj1->points[i + 3][0];
    pts_i1_3 = obj1->points[i + 3][1];
    pts_i2_3 = obj1->points[i + 3][2];

    pts_i0_3 = obj1->points[i + 3][0];
    pts_i1_3 = obj1->points[i + 3][1];
    pts_i2_3 = obj1->points[i + 3][2];

    pts_i0_4 = obj1->points[i + 4][0];
    pts_i1_4 = obj1->points[i + 4][1];
    pts_i2_4 = obj1->points[i + 4][2];

    pts_i0_5 = obj1->points[i + 5][0];
    pts_i1_5 = obj1->points[i + 5][1];
    pts_i2_5 = obj1->points[i + 5][2];

    pts_i0_6 = obj1->points[i + 6][0];
    pts_i1_6 = obj1->points[i + 6][1];
    pts_i2_6 = obj1->points[i + 6][2];

    pts_i0_7 = obj1->points[i + 7][0];
    pts_i1_7 = obj1->points[i + 7][1];
    pts_i2_7 = obj1->points[i + 7][2];

    // Calculate
    prd_0 = pts_i0 * dir_0;
    prd_1 = pts_i1 * dir_1;
    sum_0 = prd_0 + prd_1;
    prd_2 = pts_i2 * dir_2;

    prd_0_1 = pts_i0_1 * dir_0;
    prd_1_1 = pts_i1_1 * dir_1;
    sum_0_1 = prd_0_1 + prd_1_1;
    prd_2_1 = pts_i2_1 * dir_2;

    prd_0_2 = pts_i0_2 * dir_0;
    prd_1_2 = pts_i1_2 * dir_1;
    sum_0_2 = prd_0_2 + prd_1_2;
    prd_2_2 = pts_i2_2 * dir_2;

    prd_0_3 = pts_i0_3 * dir_0;
    prd_1_3 = pts_i1_3 * dir_1;
    sum_0_3 = prd_0_3 + prd_1_3;
    prd_2_3 = pts_i2_3 * dir_2;

    prd_0_4 = pts_i0_4 * dir_0;
    prd_1_4 = pts_i1_4 * dir_1;
    sum_0_4 = prd_0_4 + prd_1_4;
    prd_2_4 = pts_i2_4 * dir_2;

    prd_0_5 = pts_i0_5 * dir_0;
    prd_1_5 = pts_i1_5 * dir_1;
    sum_0_5 = prd_0_5 + prd_1_5;
    prd_2_5 = pts_i2_5 * dir_2;

    prd_0_6 = pts_i0_6 * dir_0;
    prd_1_6 = pts_i1_6 * dir_1;
    sum_0_6 = prd_0_6 + prd_1_6;
    prd_2_6 = pts_i2_6 * dir_2;

    prd_0_7 = pts_i0_7 * dir_0;
    prd_1_7 = pts_i1_7 * dir_1;
    sum_0_7 = prd_0_7 + prd_1_7;
    prd_2_7 = pts_i2_7 * dir_2;

    // Store
    dp = sum_0 + prd_2;
    dp_1 = sum_0_1 + prd_2_1;
    dp_2 = sum_0_2 + prd_2_2;
    dp_3 = sum_0_3 + prd_2_3;
    dp_4 = sum_0_4 + prd_2_4;
    dp_5 = sum_0_5 + prd_2_5;
    dp_6 = sum_0_6 + prd_2_6;
    dp_7 = sum_0_7 + prd_2_7;

    // Find local dpmax and argmax
    if (dp > dpmax) {
      dpmax = dp;
      argmax = i;
    }

    if (dp_1 > dpmax_1) {
      dpmax_1 = dp_1;
      argmax_1 = i + 1;
    }

    if (dp_2 > dpmax_2) {
      dpmax_2 = dp_2;
      argmax_2 = i + 2;
    }

    if (dp_3 > dpmax_3) {
      dpmax_3 = dp_3;
      argmax_3 = i + 3;
    }

    if (dp_4 > dpmax_4) {
      dpmax_4 = dp_4;
      argmax_4 = i + 4;
    }

    if (dp_5 > dpmax_5) {
      dpmax_5 = dp_5;
      argmax_5 = i + 5;
    }

    if (dp_6 > dpmax_6) {
      dpmax_6 = dp_6;
      argmax_6 = i + 6;
    }

    if (dp_7 > dpmax_7) {
      dpmax_7 = dp_7;
      argmax_7 = i + 7;
    }
  }

  // Find dpmax and argmax
  if (dpmax < dpmax_1) {
    dpmax = dpmax_1;
    argmax = argmax_1;
  }

  if (dpmax < dpmax_2) {
    dpmax = dpmax_2;
    argmax = argmax_2;
  }

  if (dpmax < dpmax_3) {
    dpmax = dpmax_3;
    argmax = argmax_3;
  }

  if (dpmax < dpmax_4) {
    dpmax = dpmax_4;
    argmax = argmax_4;
  }

  if (dpmax < dpmax_5) {
    dpmax = dpmax_5;
    argmax = argmax_5;
  }

  if (dpmax < dpmax_6) {
    dpmax = dpmax_6;
    argmax = argmax_6;
  }

  if (dpmax < dpmax_7) {
    dpmax = dpmax_7;
    argmax = argmax_7;
  }

  // Handle the remaining cases
  for (; i < obj1->num_points; i++) {
    // Load
    pts_i0 = obj1->points[i][0];
    pts_i1 = obj1->points[i][1];
    pts_i2 = obj1->points[i][2];

    // Calculate
    prd_0 = pts_i0 * dir_0;
    prd_1 = pts_i1 * dir_1;
    prd_2 = pts_i2 * dir_2;
    sum_0 = prd_0 + prd_1;

    // Store
    dp = sum_0 + prd_2;

    if (dp > dpmax) {
      dpmax = dp;
      argmax = i;
    }
  }

  // Joint loop ended - Complete the search for support of second object if it wasn't completed
  for (; j < obj2->num_points - 7; j += 8) {
    // Load
    invsup_pts_i0 = obj2->points[j][0];
    invsup_pts_i1 = obj2->points[j][1];
    invsup_pts_i2 = obj2->points[j][2];

    invsup_pts_i0_1 = obj2->points[j + 1][0];
    invsup_pts_i1_1 = obj2->points[j + 1][1];
    invsup_pts_i2_1 = obj2->points[j + 1][2];

    invsup_pts_i0_2 = obj2->points[j + 2][0];
    invsup_pts_i1_2 = obj2->points[j + 2][1];
    invsup_pts_i2_2 = obj2->points[j + 2][2];

    invsup_pts_i0_3 = obj2->points[j + 3][0];
    invsup_pts_i1_3 = obj2->points[j + 3][1];
    invsup_pts_i2_3 = obj2->points[j + 3][2];

    invsup_pts_i0_4 = obj2->points[j + 4][0];
    invsup_pts_i1_4 = obj2->points[j + 4][1];
    invsup_pts_i2_4 = obj2->points[j + 4][2];

    invsup_pts_i0_5 = obj2->points[j + 5][0];
    invsup_pts_i1_5 = obj2->points[j + 5][1];
    invsup_pts_i2_5 = obj2->points[j + 5][2];

    invsup_pts_i0_6 = obj2->points[j + 6][0];
    invsup_pts_i1_6 = obj2->points[j + 6][1];
    invsup_pts_i2_6 = obj2->points[j + 6][2];

    invsup_pts_i0_7 = obj2->points[j + 7][0];
    invsup_pts_i1_7 = obj2->points[j + 7][1];
    invsup_pts_i2_7 = obj2->points[j + 7][2];

    // Calculate
    invsup_prd_0 = invsup_pts_i0 * dir_0;
    invsup_prd_1 = invsup_pts_i1 * dir_1;
    invsup_sum_0 = invsup_prd_0 + invsup_prd_1;
    invsup_prd_2 = invsup_pts_i2 * dir_2;

    invsup_prd_0_1 = invsup_pts_i0_1 * dir_0;
    invsup_prd_1_1 = invsup_pts_i1_1 * dir_1;
    invsup_sum_0_1 = invsup_prd_0_1 + invsup_prd_1_1;
    invsup_prd_2_1 = invsup_pts_i2_1 * dir_2;

    invsup_prd_0_2 = invsup_pts_i0_2 * dir_0;
    invsup_prd_1_2 = invsup_pts_i1_2 * dir_1;
    invsup_sum_0_2 = invsup_prd_0_2 + invsup_prd_1_2;
    invsup_prd_2_2 = invsup_pts_i2_2 * dir_2;

    invsup_prd_0_3 = invsup_pts_i0_3 * dir_0;
    invsup_prd_1_3 = invsup_pts_i1_3 * dir_1;
    invsup_sum_0_3 = invsup_prd_0_3 + invsup_prd_1_3;
    invsup_prd_2_3 = invsup_pts_i2_3 * dir_2;

    invsup_prd_0_4 = invsup_pts_i0_4 * dir_0;
    invsup_prd_1_4 = invsup_pts_i1_4 * dir_1;
    invsup_sum_0_4 = invsup_prd_0_4 + invsup_prd_1_4;
    invsup_prd_2_4 = invsup_pts_i2_4 * dir_2;

    invsup_prd_0_5 = invsup_pts_i0_5 * dir_0;
    invsup_prd_1_5 = invsup_pts_i1_5 * dir_1;
    invsup_sum_0_5 = invsup_prd_0_5 + invsup_prd_1_5;
    invsup_prd_2_5 = invsup_pts_i2_5 * dir_2;

    invsup_prd_0_6 = invsup_pts_i0_6 * dir_0;
    invsup_prd_1_6 = invsup_pts_i1_6 * dir_1;
    invsup_sum_0_6 = invsup_prd_0_6 + invsup_prd_1_6;
    invsup_prd_2_6 = invsup_pts_i2_6 * dir_2;

    invsup_prd_0_7 = invsup_pts_i0_7 * dir_0;
    invsup_prd_1_7 = invsup_pts_i1_7 * dir_1;
    invsup_sum_0_7 = invsup_prd_0_7 + invsup_prd_1_7;
    invsup_prd_2_7 = invsup_pts_i2_7 * dir_2;

    // Store
    invsup_dp = invsup_sum_0 + invsup_prd_2;
    invsup_dp_1 = invsup_sum_0_1 + invsup_prd_2_1;
    invsup_dp_2 = invsup_sum_0_2 + invsup_prd_2_2;
    invsup_dp_3 = invsup_sum_0_3 + invsup_prd_2_3;
    invsup_dp_4 = invsup_sum_0_4 + invsup_prd_2_4;
    invsup_dp_5 = invsup_sum_0_5 + invsup_prd_2_5;
    invsup_dp_6 = invsup_sum_0_6 + invsup_prd_2_6;
    invsup_dp_7 = invsup_sum_0_7 + invsup_prd_2_7;

    // Find local dpmin and argmin
    if (invsup_dp < invsup_dpmax) {
      invsup_dpmax = invsup_dp;
      argmin = j;
    }

    if (invsup_dp_1 < invsup_dpmax_1) {
      invsup_dpmax_1 = invsup_dp_1;
      argmin_1 = j + 1;
    }

    if (invsup_dp_2 < invsup_dpmax_2) {
      invsup_dpmax_2 = invsup_dp_2;
      argmin_2 = j + 2;
    }

    if (invsup_dp_3 < invsup_dpmax_3) {
      invsup_dpmax_3 = invsup_dp_3;
      argmin_3 = j + 3;
    }

    if (invsup_dp_4 < invsup_dpmax_4) {
      invsup_dpmax_4 = invsup_dp_4;
      argmin_4 = j + 4;
    }

    if (invsup_dp_5 < invsup_dpmax_5) {
      invsup_dpmax_5 = invsup_dp_5;
      argmin_5 = j + 5;
    }

    if (invsup_dp_6 < invsup_dpmax_6) {
      invsup_dpmax_6 = invsup_dp_6;
      argmin_6 = j + 6;
    }

    if (invsup_dp_7 < invsup_dpmax_7) {
      invsup_dpmax_7 = invsup_dp_7;
      argmin_7 = j + 7;
    }
  }

  // Find dpmin and argmin
  if (invsup_dpmax > invsup_dpmax_1) {
    invsup_dpmax = invsup_dpmax_1;
    argmin = argmin_1;
  }

  if (invsup_dpmax > invsup_dpmax_2) {
    invsup_dpmax = invsup_dpmax_2;
    argmin = argmin_2;
  }

  if (invsup_dpmax > invsup_dpmax_3) {
    invsup_dpmax = invsup_dpmax_3;
    argmin = argmin_3;
  }

  if (invsup_dpmax > invsup_dpmax_4) {
    invsup_dpmax = invsup_dpmax_4;
    argmin = argmin_4;
  }

  if (invsup_dpmax > invsup_dpmax_5) {
    invsup_dpmax = invsup_dpmax_5;
    argmin = argmin_5;
  }

  if (invsup_dpmax > invsup_dpmax_6) {
    invsup_dpmax = invsup_dpmax_6;
    argmin = argmin_6;
  }

  if (invsup_dpmax > invsup_dpmax_7) {
    invsup_dpmax = invsup_dpmax_7;
    argmin = argmin_7;
  }

  // Handle the remaining cases
  for (; j < obj2->num_points; j++) {
    // Load
    invsup_pts_i0 = obj2->points[j][0];
    invsup_pts_i1 = obj2->points[j][1];
    invsup_pts_i2 = obj2->points[j][2];

    // Calculate
    invsup_prd_0 = invsup_pts_i0 * dir_0;
    invsup_prd_1 = invsup_pts_i1 * dir_1;
    invsup_prd_2 = invsup_pts_i2 * dir_2;
    invsup_sum_0 = invsup_prd_0 + invsup_prd_1;

    // Store
    invsup_dp = invsup_sum_0 + invsup_prd_2;

    if (invsup_dp < invsup_dpmax) {
      invsup_dpmax = invsup_dp;
      argmin = j;
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
    dpmax = -1e9, dpmax_1 = -1e9, dpmax_2 = -1e9,
      dpmax_3 = -1e9, dpmax_4 = -1e9, dpmax_5 = -1e9, dpmax_6 = -1e9, dpmax_7 = -1e9;
    argmax = 0, argmax_1 = 0, argmax_2 = 0,
      argmax_3 = 0, argmax_4 = 0, argmax_5 = 0, argmax_6 = 0, argmax_7 = 0;

    invsup_dpmax = 1e9, invsup_dpmax_1 = 1e9, invsup_dpmax_2 = 1e9,
          invsup_dpmax_3 = 1e9, invsup_dpmax_4 = 1e9, invsup_dpmax_5 = 1e9,
          invsup_dpmax_6 = 1e9, invsup_dpmax_7 = 1e9;
    argmin = 0, argmin_1 = 0, argmin_2 = 0,
      argmin_3 = 0, argmin_4 = 0, argmin_5 = 0, argmin_6 = 0, argmin_7 = 0;

    dir_0 = d[0];
    dir_1 = d[1];
    dir_2 = d[2];

    // Joint loop that iterates over both object to calculate dp and argmin/argmax.
    for (i = 0; i < no_points_min - 3; i += 4) {
      // Load first object
      pts_i0 = obj1->points[i][0];
      pts_i1 = obj1->points[i][1];
      pts_i2 = obj1->points[i][2];

      pts_i0_1 = obj1->points[i + 1][0];
      pts_i1_1 = obj1->points[i + 1][1];
      pts_i2_1 = obj1->points[i + 1][2];

      pts_i0_2 = obj1->points[i + 2][0];
      pts_i1_2 = obj1->points[i + 2][1];
      pts_i2_2 = obj1->points[i + 2][2];

      pts_i0_3 = obj1->points[i + 3][0];
      pts_i1_3 = obj1->points[i + 3][1];
      pts_i2_3 = obj1->points[i + 3][2];

      // Load second object
      invsup_pts_i0 = obj2->points[i][0];
      invsup_pts_i1 = obj2->points[i][1];
      invsup_pts_i2 = obj2->points[i][2];

      invsup_pts_i0_1 = obj2->points[i + 1][0];
      invsup_pts_i1_1 = obj2->points[i + 1][1];
      invsup_pts_i2_1 = obj2->points[i + 1][2];

      invsup_pts_i0_2 = obj2->points[i + 2][0];
      invsup_pts_i1_2 = obj2->points[i + 2][1];
      invsup_pts_i2_2 = obj2->points[i + 2][2];

      invsup_pts_i0_3 = obj2->points[i + 3][0];
      invsup_pts_i1_3 = obj2->points[i + 3][1];
      invsup_pts_i2_3 = obj2->points[i + 3][2];

      // Calculate support
      prd_0 = pts_i0 * dir_0;
      prd_1 = pts_i1 * dir_1;
      sum_0 = prd_0 + prd_1;
      prd_2 = pts_i2 * dir_2;

      prd_0_1 = pts_i0_1 * dir_0;
      prd_1_1 = pts_i1_1 * dir_1;
      sum_0_1 = prd_0_1 + prd_1_1;
      prd_2_1 = pts_i2_1 * dir_2;

      prd_0_2 = pts_i0_2 * dir_0;
      prd_1_2 = pts_i1_2 * dir_1;
      sum_0_2 = prd_0_2 + prd_1_2;
      prd_2_2 = pts_i2_2 * dir_2;

      prd_0_3 = pts_i0_3 * dir_0;
      prd_1_3 = pts_i1_3 * dir_1;
      sum_0_3 = prd_0_3 + prd_1_3;
      prd_2_3 = pts_i2_3 * dir_2;

      // Calculate the inverse support
      invsup_prd_0 = invsup_pts_i0 * dir_0;
      invsup_prd_1 = invsup_pts_i1 * dir_1;
      invsup_sum_0 = invsup_prd_0 + invsup_prd_1;
      invsup_prd_2 = invsup_pts_i2 * dir_2;

      invsup_prd_0_1 = invsup_pts_i0_1 * dir_0;
      invsup_prd_1_1 = invsup_pts_i1_1 * dir_1;
      invsup_sum_0_1 = invsup_prd_0_1 + invsup_prd_1_1;
      invsup_prd_2_1 = invsup_pts_i2_1 * dir_2;

      invsup_prd_0_2 = invsup_pts_i0_2 * dir_0;
      invsup_prd_1_2 = invsup_pts_i1_2 * dir_1;
      invsup_sum_0_2 = invsup_prd_0_2 + invsup_prd_1_2;
      invsup_prd_2_2 = invsup_pts_i2_2 * dir_2;

      invsup_prd_0_3 = invsup_pts_i0_3 * dir_0;
      invsup_prd_1_3 = invsup_pts_i1_3 * dir_1;
      invsup_sum_0_3 = invsup_prd_0_3 + invsup_prd_1_3;
      invsup_prd_2_3 = invsup_pts_i2_3 * dir_2;

      // Store first object's dot products
      dp = sum_0 + prd_2;
      dp_1 = sum_0_1 + prd_2_1;
      dp_2 = sum_0_2 + prd_2_2;
      dp_3 = sum_0_3 + prd_2_3;

      // Store second object's dot products
      invsup_dp = invsup_sum_0 + invsup_prd_2;
      invsup_dp_1 = invsup_sum_0_1 + invsup_prd_2_1;
      invsup_dp_2 = invsup_sum_0_2 + invsup_prd_2_2;
      invsup_dp_3 = invsup_sum_0_3 + invsup_prd_2_3;

      // Find local dpmax and argmax
      if (dp > dpmax) {
        dpmax = dp;
        argmax = i;
      }

      if (dp_1 > dpmax_1) {
        dpmax_1 = dp_1;
        argmax_1 = i + 1;
      }

      if (dp_2 > dpmax_2) {
        dpmax_2 = dp_2;
        argmax_2 = i + 2;
      }

      if (dp_3 > dpmax_3) {
        dpmax_3 = dp_3;
        argmax_3 = i + 3;
      }

      // Find local dpmin and argmin
      if (invsup_dp < invsup_dpmax) {
        invsup_dpmax = invsup_dp;
        argmin = i;
      }

      if (invsup_dp_1 < invsup_dpmax_1) {
        invsup_dpmax_1 = invsup_dp_1;
        argmin_1 = i + 1;
      }

      if (invsup_dp_2 < invsup_dpmax_2) {
        invsup_dpmax_2 = invsup_dp_2;
        argmin_2 = i + 2;
      }

      if (invsup_dp_3 < invsup_dpmax_3) {
        invsup_dpmax_3 = invsup_dp_3;
        argmin_3 = i + 3;
      }
    }

    // Joint loop ended - Complete the search for support of first object if it wasn't completed
    j = i;
    for (; i < obj1->num_points - 7; i += 8) {
      // Load
      pts_i0 = obj1->points[i][0];
      pts_i1 = obj1->points[i][1];
      pts_i2 = obj1->points[i][2];

      pts_i0_1 = obj1->points[i + 1][0];
      pts_i1_1 = obj1->points[i + 1][1];
      pts_i2_1 = obj1->points[i + 1][2];

      pts_i0_2 = obj1->points[i + 2][0];
      pts_i1_2 = obj1->points[i + 2][1];
      pts_i2_2 = obj1->points[i + 2][2];

      pts_i0_3 = obj1->points[i + 3][0];
      pts_i1_3 = obj1->points[i + 3][1];
      pts_i2_3 = obj1->points[i + 3][2];

      pts_i0_3 = obj1->points[i + 3][0];
      pts_i1_3 = obj1->points[i + 3][1];
      pts_i2_3 = obj1->points[i + 3][2];

      pts_i0_4 = obj1->points[i + 4][0];
      pts_i1_4 = obj1->points[i + 4][1];
      pts_i2_4 = obj1->points[i + 4][2];

      pts_i0_5 = obj1->points[i + 5][0];
      pts_i1_5 = obj1->points[i + 5][1];
      pts_i2_5 = obj1->points[i + 5][2];

      pts_i0_6 = obj1->points[i + 6][0];
      pts_i1_6 = obj1->points[i + 6][1];
      pts_i2_6 = obj1->points[i + 6][2];

      pts_i0_7 = obj1->points[i + 7][0];
      pts_i1_7 = obj1->points[i + 7][1];
      pts_i2_7 = obj1->points[i + 7][2];

      // Calculate
      prd_0 = pts_i0 * dir_0;
      prd_1 = pts_i1 * dir_1;
      sum_0 = prd_0 + prd_1;
      prd_2 = pts_i2 * dir_2;

      prd_0_1 = pts_i0_1 * dir_0;
      prd_1_1 = pts_i1_1 * dir_1;
      sum_0_1 = prd_0_1 + prd_1_1;
      prd_2_1 = pts_i2_1 * dir_2;

      prd_0_2 = pts_i0_2 * dir_0;
      prd_1_2 = pts_i1_2 * dir_1;
      sum_0_2 = prd_0_2 + prd_1_2;
      prd_2_2 = pts_i2_2 * dir_2;

      prd_0_3 = pts_i0_3 * dir_0;
      prd_1_3 = pts_i1_3 * dir_1;
      sum_0_3 = prd_0_3 + prd_1_3;
      prd_2_3 = pts_i2_3 * dir_2;

      prd_0_4 = pts_i0_4 * dir_0;
      prd_1_4 = pts_i1_4 * dir_1;
      sum_0_4 = prd_0_4 + prd_1_4;
      prd_2_4 = pts_i2_4 * dir_2;

      prd_0_5 = pts_i0_5 * dir_0;
      prd_1_5 = pts_i1_5 * dir_1;
      sum_0_5 = prd_0_5 + prd_1_5;
      prd_2_5 = pts_i2_5 * dir_2;

      prd_0_6 = pts_i0_6 * dir_0;
      prd_1_6 = pts_i1_6 * dir_1;
      sum_0_6 = prd_0_6 + prd_1_6;
      prd_2_6 = pts_i2_6 * dir_2;

      prd_0_7 = pts_i0_7 * dir_0;
      prd_1_7 = pts_i1_7 * dir_1;
      sum_0_7 = prd_0_7 + prd_1_7;
      prd_2_7 = pts_i2_7 * dir_2;

      // Store
      dp = sum_0 + prd_2;
      dp_1 = sum_0_1 + prd_2_1;
      dp_2 = sum_0_2 + prd_2_2;
      dp_3 = sum_0_3 + prd_2_3;
      dp_4 = sum_0_4 + prd_2_4;
      dp_5 = sum_0_5 + prd_2_5;
      dp_6 = sum_0_6 + prd_2_6;
      dp_7 = sum_0_7 + prd_2_7;

      // Find local dpmax and argmax
      if (dp > dpmax) {
        dpmax = dp;
        argmax = i;
      }

      if (dp_1 > dpmax_1) {
        dpmax_1 = dp_1;
        argmax_1 = i + 1;
      }

      if (dp_2 > dpmax_2) {
        dpmax_2 = dp_2;
        argmax_2 = i + 2;
      }

      if (dp_3 > dpmax_3) {
        dpmax_3 = dp_3;
        argmax_3 = i + 3;
      }

      if (dp_4 > dpmax_4) {
        dpmax_4 = dp_4;
        argmax_4 = i + 4;
      }

      if (dp_5 > dpmax_5) {
        dpmax_5 = dp_5;
        argmax_5 = i + 5;
      }

      if (dp_6 > dpmax_6) {
        dpmax_6 = dp_6;
        argmax_6 = i + 6;
      }

      if (dp_7 > dpmax_7) {
        dpmax_7 = dp_7;
        argmax_7 = i + 7;
      }
    }

    // Find dpmax and argmax
    if (dpmax < dpmax_1) {
      dpmax = dpmax_1;
      argmax = argmax_1;
    }

    if (dpmax < dpmax_2) {
      dpmax = dpmax_2;
      argmax = argmax_2;
    }

    if (dpmax < dpmax_3) {
      dpmax = dpmax_3;
      argmax = argmax_3;
    }

    if (dpmax < dpmax_4) {
      dpmax = dpmax_4;
      argmax = argmax_4;
    }

    if (dpmax < dpmax_5) {
      dpmax = dpmax_5;
      argmax = argmax_5;
    }

    if (dpmax < dpmax_6) {
      dpmax = dpmax_6;
      argmax = argmax_6;
    }

    if (dpmax < dpmax_7) {
      dpmax = dpmax_7;
      argmax = argmax_7;
    }

    // Handle the remaining cases
    for (; i < obj1->num_points; i++) {
      // Load
      pts_i0 = obj1->points[i][0];
      pts_i1 = obj1->points[i][1];
      pts_i2 = obj1->points[i][2];

      // Calculate
      prd_0 = pts_i0 * dir_0;
      prd_1 = pts_i1 * dir_1;
      prd_2 = pts_i2 * dir_2;
      sum_0 = prd_0 + prd_1;

      // Store
      dp = sum_0 + prd_2;

      if (dp > dpmax) {
        dpmax = dp;
        argmax = i;
      }
    }

    // Joint loop ended - Complete the search for support of second object if it wasn't completed
    for (; j < obj2->num_points - 7; j += 8) {
      // Load
      invsup_pts_i0 = obj2->points[j][0];
      invsup_pts_i1 = obj2->points[j][1];
      invsup_pts_i2 = obj2->points[j][2];

      invsup_pts_i0_1 = obj2->points[j + 1][0];
      invsup_pts_i1_1 = obj2->points[j + 1][1];
      invsup_pts_i2_1 = obj2->points[j + 1][2];

      invsup_pts_i0_2 = obj2->points[j + 2][0];
      invsup_pts_i1_2 = obj2->points[j + 2][1];
      invsup_pts_i2_2 = obj2->points[j + 2][2];

      invsup_pts_i0_3 = obj2->points[j + 3][0];
      invsup_pts_i1_3 = obj2->points[j + 3][1];
      invsup_pts_i2_3 = obj2->points[j + 3][2];

      invsup_pts_i0_4 = obj2->points[j + 4][0];
      invsup_pts_i1_4 = obj2->points[j + 4][1];
      invsup_pts_i2_4 = obj2->points[j + 4][2];

      invsup_pts_i0_5 = obj2->points[j + 5][0];
      invsup_pts_i1_5 = obj2->points[j + 5][1];
      invsup_pts_i2_5 = obj2->points[j + 5][2];

      invsup_pts_i0_6 = obj2->points[j + 6][0];
      invsup_pts_i1_6 = obj2->points[j + 6][1];
      invsup_pts_i2_6 = obj2->points[j + 6][2];

      invsup_pts_i0_7 = obj2->points[j + 7][0];
      invsup_pts_i1_7 = obj2->points[j + 7][1];
      invsup_pts_i2_7 = obj2->points[j + 7][2];

      // Calculate
      invsup_prd_0 = invsup_pts_i0 * dir_0;
      invsup_prd_1 = invsup_pts_i1 * dir_1;
      invsup_sum_0 = invsup_prd_0 + invsup_prd_1;
      invsup_prd_2 = invsup_pts_i2 * dir_2;

      invsup_prd_0_1 = invsup_pts_i0_1 * dir_0;
      invsup_prd_1_1 = invsup_pts_i1_1 * dir_1;
      invsup_sum_0_1 = invsup_prd_0_1 + invsup_prd_1_1;
      invsup_prd_2_1 = invsup_pts_i2_1 * dir_2;

      invsup_prd_0_2 = invsup_pts_i0_2 * dir_0;
      invsup_prd_1_2 = invsup_pts_i1_2 * dir_1;
      invsup_sum_0_2 = invsup_prd_0_2 + invsup_prd_1_2;
      invsup_prd_2_2 = invsup_pts_i2_2 * dir_2;

      invsup_prd_0_3 = invsup_pts_i0_3 * dir_0;
      invsup_prd_1_3 = invsup_pts_i1_3 * dir_1;
      invsup_sum_0_3 = invsup_prd_0_3 + invsup_prd_1_3;
      invsup_prd_2_3 = invsup_pts_i2_3 * dir_2;

      invsup_prd_0_4 = invsup_pts_i0_4 * dir_0;
      invsup_prd_1_4 = invsup_pts_i1_4 * dir_1;
      invsup_sum_0_4 = invsup_prd_0_4 + invsup_prd_1_4;
      invsup_prd_2_4 = invsup_pts_i2_4 * dir_2;

      invsup_prd_0_5 = invsup_pts_i0_5 * dir_0;
      invsup_prd_1_5 = invsup_pts_i1_5 * dir_1;
      invsup_sum_0_5 = invsup_prd_0_5 + invsup_prd_1_5;
      invsup_prd_2_5 = invsup_pts_i2_5 * dir_2;

      invsup_prd_0_6 = invsup_pts_i0_6 * dir_0;
      invsup_prd_1_6 = invsup_pts_i1_6 * dir_1;
      invsup_sum_0_6 = invsup_prd_0_6 + invsup_prd_1_6;
      invsup_prd_2_6 = invsup_pts_i2_6 * dir_2;

      invsup_prd_0_7 = invsup_pts_i0_7 * dir_0;
      invsup_prd_1_7 = invsup_pts_i1_7 * dir_1;
      invsup_sum_0_7 = invsup_prd_0_7 + invsup_prd_1_7;
      invsup_prd_2_7 = invsup_pts_i2_7 * dir_2;

      // Store
      invsup_dp = invsup_sum_0 + invsup_prd_2;
      invsup_dp_1 = invsup_sum_0_1 + invsup_prd_2_1;
      invsup_dp_2 = invsup_sum_0_2 + invsup_prd_2_2;
      invsup_dp_3 = invsup_sum_0_3 + invsup_prd_2_3;
      invsup_dp_4 = invsup_sum_0_4 + invsup_prd_2_4;
      invsup_dp_5 = invsup_sum_0_5 + invsup_prd_2_5;
      invsup_dp_6 = invsup_sum_0_6 + invsup_prd_2_6;
      invsup_dp_7 = invsup_sum_0_7 + invsup_prd_2_7;

      // Find local dpmin and argmin
      if (invsup_dp < invsup_dpmax) {
        invsup_dpmax = invsup_dp;
        argmin = j;
      }

      if (invsup_dp_1 < invsup_dpmax_1) {
        invsup_dpmax_1 = invsup_dp_1;
        argmin_1 = j + 1;
      }

      if (invsup_dp_2 < invsup_dpmax_2) {
        invsup_dpmax_2 = invsup_dp_2;
        argmin_2 = j + 2;
      }

      if (invsup_dp_3 < invsup_dpmax_3) {
        invsup_dpmax_3 = invsup_dp_3;
        argmin_3 = j + 3;
      }

      if (invsup_dp_4 < invsup_dpmax_4) {
        invsup_dpmax_4 = invsup_dp_4;
        argmin_4 = j + 4;
      }

      if (invsup_dp_5 < invsup_dpmax_5) {
        invsup_dpmax_5 = invsup_dp_5;
        argmin_5 = j + 5;
      }

      if (invsup_dp_6 < invsup_dpmax_6) {
        invsup_dpmax_6 = invsup_dp_6;
        argmin_6 = j + 6;
      }

      if (invsup_dp_7 < invsup_dpmax_7) {
        invsup_dpmax_7 = invsup_dp_7;
        argmin_7 = j + 7;
      }
    }

    // Find dpmin and argmin
    if (invsup_dpmax > invsup_dpmax_1) {
      invsup_dpmax = invsup_dpmax_1;
      argmin = argmin_1;
    }

    if (invsup_dpmax > invsup_dpmax_2) {
      invsup_dpmax = invsup_dpmax_2;
      argmin = argmin_2;
    }

    if (invsup_dpmax > invsup_dpmax_3) {
      invsup_dpmax = invsup_dpmax_3;
      argmin = argmin_3;
    }

    if (invsup_dpmax > invsup_dpmax_4) {
      invsup_dpmax = invsup_dpmax_4;
      argmin = argmin_4;
    }

    if (invsup_dpmax > invsup_dpmax_5) {
      invsup_dpmax = invsup_dpmax_5;
      argmin = argmin_5;
    }

    if (invsup_dpmax > invsup_dpmax_6) {
      invsup_dpmax = invsup_dpmax_6;
      argmin = argmin_6;
    }

    if (invsup_dpmax > invsup_dpmax_7) {
      invsup_dpmax = invsup_dpmax_7;
      argmin = argmin_7;
    }

    // Handle the remaining cases
    for (; j < obj2->num_points; j++) {
      // Load
      invsup_pts_i0 = obj2->points[j][0];
      invsup_pts_i1 = obj2->points[j][1];
      invsup_pts_i2 = obj2->points[j][2];

      // Calculate
      invsup_prd_0 = invsup_pts_i0 * dir_0;
      invsup_prd_1 = invsup_pts_i1 * dir_1;
      invsup_prd_2 = invsup_pts_i2 * dir_2;
      invsup_sum_0 = invsup_prd_0 + invsup_prd_1;

      // Store
      invsup_dp = invsup_sum_0 + invsup_prd_2;

      if (invsup_dp < invsup_dpmax) {
        invsup_dpmax = invsup_dp;
        argmin = j;
      }
    }

    // a = Support(A-B) = Support(A) - Invsup(B)
    a[0] = obj1->points[argmax][0] - obj2->points[argmin][0];
    a[1] = obj1->points[argmax][1] - obj2->points[argmin][1];
    a[2] = obj1->points[argmax][2] - obj2->points[argmin][2];

    // Are we past the origin? If no, we will never enclose it.
    if (dot3D_optimized(a, d) < 0.) return 0;

    // Add newly found point to Simplex
    // AT THIS POINT WE KNOW a IS PAST THE ORIGIN
    s.p[s.num_points][0] = a[0];
    s.p[s.num_points][1] = a[1];
    s.p[s.num_points][2] = a[2];
    s.num_points++;

    if (do_simplex3D_optimized(&s, d)) return 1;
  }

  return 1; // Most likely not reachable. Say we intersect by default.
}

// Seperate argmin and dot product operations.
int do_intersect3D_o_i_branch_later(const struct CHObject* obj1, const struct CHObject* obj2) {
  // The search direction
  float d[3] = {1.f, 1.f, 1.f};  // could also be random
  // The point we examine in each iteration
  float a[3];

  // Our simplex that do_simplex modifies.
  struct Simplex3D s;
  s.num_points = 1;

  // Method inlining for support function
  int no_points_min = obj1->num_points > obj2->num_points ? obj2->num_points : obj1->num_points;
  float* mvm_result_1 = (float*) malloc(obj1->num_points * sizeof(float));
  float* mvm_result_2 = (float*) malloc(obj2->num_points * sizeof(float));

  float dp, dp_1, dp_2, dpmax = -1e9;
  int i, argmax = 0;

  float invsup_dp, invsup_dp_1, invsup_dp_2, invsup_dpmax = 1e9;
  int argmin = 0;

  float dir_0;
  float dir_1;
  float dir_2;

  float pts_i0, pts_i1, pts_i2;
  float pts_i0_1, pts_i1_1, pts_i2_1;
  float pts_i0_2, pts_i1_2, pts_i2_2;

  float invsup_pts_i0, invsup_pts_i1, invsup_pts_i2;
  float invsup_pts_i0_1, invsup_pts_i1_1, invsup_pts_i2_1;
  float invsup_pts_i0_2, invsup_pts_i1_2, invsup_pts_i2_2;

  float prd_0, prd_1, prd_2, sum_0;
  float prd_0_1, prd_1_1, prd_2_1, sum_0_1;
  float prd_0_2, prd_1_2, prd_2_2, sum_0_2;

  float invsup_prd_0, invsup_prd_1, invsup_prd_2, invsup_sum_0;
  float invsup_prd_0_1, invsup_prd_1_1, invsup_prd_2_1, invsup_sum_0_1;
  float invsup_prd_0_2, invsup_prd_1_2, invsup_prd_2_2, invsup_sum_0_2;

  dir_0 = d[0];
  dir_1 = d[1];
  dir_2 = d[2];

  for (i = 0; i < no_points_min - 2; i += 3) {
    // Load
    pts_i0 = obj1->points[i][0];
    pts_i1 = obj1->points[i][1];
    pts_i2 = obj1->points[i][2];

    pts_i0_1 = obj1->points[i + 1][0];
    pts_i1_1 = obj1->points[i + 1][1];
    pts_i2_1 = obj1->points[i + 1][2];

    pts_i0_2 = obj1->points[i + 2][0];
    pts_i1_2 = obj1->points[i + 2][1];
    pts_i2_2 = obj1->points[i + 2][2];

    invsup_pts_i0 = obj2->points[i][0];
    invsup_pts_i1 = obj2->points[i][1];
    invsup_pts_i2 = obj2->points[i][2];

    invsup_pts_i0_1 = obj2->points[i + 1][0];
    invsup_pts_i1_1 = obj2->points[i + 1][1];
    invsup_pts_i2_1 = obj2->points[i + 1][2];

    invsup_pts_i0_2 = obj2->points[i + 2][0];
    invsup_pts_i1_2 = obj2->points[i + 2][1];
    invsup_pts_i2_2 = obj2->points[i + 2][2];

    // Calculate
    prd_0 = pts_i0 * dir_0;
    prd_1 = pts_i1 * dir_1;
    sum_0 = prd_0 + prd_1;
    prd_2 = pts_i2 * dir_2;

    prd_0_1 = pts_i0_1 * dir_0;
    prd_1_1 = pts_i1_1 * dir_1;
    sum_0_1 = prd_0_1 + prd_1_1;
    prd_2_1 = pts_i2_1 * dir_2;

    prd_0_2 = pts_i0_2 * dir_0;
    prd_1_2 = pts_i1_2 * dir_1;
    sum_0_2 = prd_0_2 + prd_1_2;
    prd_2_2 = pts_i2_2 * dir_2;

    invsup_prd_0 = invsup_pts_i0 * dir_0;
    invsup_prd_1 = invsup_pts_i1 * dir_1;
    invsup_sum_0 = invsup_prd_0 + invsup_prd_1;
    invsup_prd_2 = invsup_pts_i2 * dir_2;

    invsup_prd_0_1 = invsup_pts_i0_1 * dir_0;
    invsup_prd_1_1 = invsup_pts_i1_1 * dir_1;
    invsup_sum_0_1 = invsup_prd_0_1 + invsup_prd_1_1;
    invsup_prd_2_1 = invsup_pts_i2_1 * dir_2;

    invsup_prd_0_2 = invsup_pts_i0_2 * dir_0;
    invsup_prd_1_2 = invsup_pts_i1_2 * dir_1;
    invsup_sum_0_2 = invsup_prd_0_2 + invsup_prd_1_2;
    invsup_prd_2_2 = invsup_pts_i2_2 * dir_2;

    // Store
    mvm_result_1[i] = sum_0 + prd_2;
    mvm_result_1[i + 1] = sum_0_1 + prd_2_1;
    mvm_result_1[i + 2] = sum_0_2 + prd_2_2;

    mvm_result_2[i] = invsup_sum_0 + invsup_prd_2;
    mvm_result_2[i + 1] = invsup_sum_0_1 + invsup_prd_2_1;
    mvm_result_2[i + 2] = invsup_sum_0_2 + invsup_prd_2_2;
  }

  int j = i;
  for (; i < obj1->num_points - 2; i += 3) {
    // Load
    pts_i0 = obj1->points[i][0];
    pts_i1 = obj1->points[i][1];
    pts_i2 = obj1->points[i][2];

    pts_i0_1 = obj1->points[i + 1][0];
    pts_i1_1 = obj1->points[i + 1][1];
    pts_i2_1 = obj1->points[i + 1][2];

    pts_i0_2 = obj1->points[i + 2][0];
    pts_i1_2 = obj1->points[i + 2][1];
    pts_i2_2 = obj1->points[i + 2][2];

    // Calculate
    prd_0 = pts_i0 * dir_0;
    prd_1 = pts_i1 * dir_1;
    sum_0 = prd_0 + prd_1;
    prd_2 = pts_i2 * dir_2;

    prd_0_1 = pts_i0_1 * dir_0;
    prd_1_1 = pts_i1_1 * dir_1;
    sum_0_1 = prd_0_1 + prd_1_1;
    prd_2_1 = pts_i2_1 * dir_2;

    prd_0_2 = pts_i0_2 * dir_0;
    prd_1_2 = pts_i1_2 * dir_1;
    sum_0_2 = prd_0_2 + prd_1_2;
    prd_2_2 = pts_i2_2 * dir_2;

    // Store
    mvm_result_1[i] = sum_0 + prd_2;
    mvm_result_1[i + 1] = sum_0_1 + prd_2_1;
    mvm_result_1[i + 2] = sum_0_2 + prd_2_2;
  }

  // Handle the remaining cases
  for (; i < obj1->num_points; i++) {
    // Load
    pts_i0 = obj1->points[i][0];
    pts_i1 = obj1->points[i][1];
    pts_i2 = obj1->points[i][2];

    // Calculate
    prd_0 = pts_i0 * dir_0;
    prd_1 = pts_i1 * dir_1;
    prd_2 = pts_i2 * dir_2;
    sum_0 = prd_0 + prd_1;

    // Store
    mvm_result_1[i] = sum_0 + prd_2;
  }

  for (; j < obj2->num_points - 2; j += 3) {
    /**
    * 1) Function Inlining
    * 2) Scalar Replacement
    * 3) Loop unrolling
    * */

    // Load
    invsup_pts_i0 = obj2->points[j][0];
    invsup_pts_i1 = obj2->points[j][1];
    invsup_pts_i2 = obj2->points[j][2];

    invsup_pts_i0_1 = obj2->points[j + 1][0];
    invsup_pts_i1_1 = obj2->points[j + 1][1];
    invsup_pts_i2_1 = obj2->points[j + 1][2];

    invsup_pts_i0_2 = obj2->points[j + 2][0];
    invsup_pts_i1_2 = obj2->points[j + 2][1];
    invsup_pts_i2_2 = obj2->points[j + 2][2];

    // Calculate
    invsup_prd_0 = invsup_pts_i0 * dir_0;
    invsup_prd_1 = invsup_pts_i1 * dir_1;
    invsup_sum_0 = invsup_prd_0 + invsup_prd_1;
    invsup_prd_2 = invsup_pts_i2 * dir_2;

    invsup_prd_0_1 = invsup_pts_i0_1 * dir_0;
    invsup_prd_1_1 = invsup_pts_i1_1 * dir_1;
    invsup_sum_0_1 = invsup_prd_0_1 + invsup_prd_1_1;
    invsup_prd_2_1 = invsup_pts_i2_1 * dir_2;

    invsup_prd_0_2 = invsup_pts_i0_2 * dir_0;
    invsup_prd_1_2 = invsup_pts_i1_2 * dir_1;
    invsup_sum_0_2 = invsup_prd_0_2 + invsup_prd_1_2;
    invsup_prd_2_2 = invsup_pts_i2_2 * dir_2;

    // Store
    mvm_result_2[j] = invsup_sum_0 + invsup_prd_2;
    mvm_result_2[j + 1] = invsup_sum_0_1 + invsup_prd_2_1;
    mvm_result_2[j + 2] = invsup_sum_0_2 + invsup_prd_2_2;
  }

  // Handle the remaining cases
  for (; j < obj2->num_points; j++) {
    // Load
    invsup_pts_i0 = obj2->points[j][0];
    invsup_pts_i1 = obj2->points[j][1];
    invsup_pts_i2 = obj2->points[j][2];

    // Calculate
    invsup_prd_0 = invsup_pts_i0 * dir_0;
    invsup_prd_1 = invsup_pts_i1 * dir_1;
    invsup_prd_2 = invsup_pts_i2 * dir_2;
    invsup_sum_0 = invsup_prd_0 + invsup_prd_1;

    // Store
    mvm_result_2[j] = invsup_sum_0 + invsup_prd_2;
  }

  for (i = 0; i < no_points_min - 2; i += 3) {
    dp = mvm_result_1[i];
    dp_1 = mvm_result_1[i + 1];
    dp_2 = mvm_result_1[i + 2];

    invsup_dp = mvm_result_2[i];
    invsup_dp_1 = mvm_result_2[i + 1];
    invsup_dp_2 = mvm_result_2[i + 2];

    if (dp > dpmax) {
      dpmax = dp;
      argmax = i;
    }

    if (dp_1 > dpmax) {
      dpmax = dp_1;
      argmax = i + 1;
    }

    if (dp_2 > dpmax) {
      dpmax = dp_2;
      argmax = i + 2;
    }

    if (invsup_dp < invsup_dpmax) {
      invsup_dpmax = invsup_dp;
      argmin = i;
    }

    if (invsup_dp_1 < invsup_dpmax) {
      invsup_dpmax = invsup_dp_1;
      argmin = i + 1;
    }

    if (invsup_dp_2 < invsup_dpmax) {
      invsup_dpmax = invsup_dp_2;
      argmin = i + 2;
    }
  }

  j = i;
  for (; i < obj1->num_points - 2; i += 3) {
    dp = mvm_result_1[i];
    dp_1 = mvm_result_1[i + 1];
    dp_2 = mvm_result_1[i + 2];

    if (dp > dpmax) {
      dpmax = dp;
      argmax = i;
    }

    if (dp_1 > dpmax) {
      dpmax = dp_1;
      argmax = i + 1;
    }

    if (dp_2 > dpmax) {
      dpmax = dp_2;
      argmax = i + 2;
    }
  }

  for (; i < obj1->num_points; i++) {
    dp = mvm_result_1[i];

    if (dp > dpmax) {
      dpmax = dp;
      argmax = i;
    }
  }

  for (; j < obj2->num_points - 2; j += 3) {
    invsup_dp = mvm_result_2[j];
    invsup_dp_1 = mvm_result_2[j + 1];
    invsup_dp_2 = mvm_result_2[j + 2];

    if (invsup_dp < invsup_dpmax) {
      invsup_dpmax = invsup_dp;
      argmin = j;
    }

    if (invsup_dp_1 < invsup_dpmax) {
      invsup_dpmax = invsup_dp_1;
      argmin = j + 1;
    }

    if (invsup_dp_2 < invsup_dpmax) {
      invsup_dpmax = invsup_dp_2;
      argmin = j + 2;
    }
  }

  for (; j < obj2->num_points; j++) {
    invsup_dp = mvm_result_2[j];

    if (invsup_dp < invsup_dpmax) {
      invsup_dpmax = invsup_dp;
      argmin = j;
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
    dpmax = -1e9;
    invsup_dpmax = 1e9;
    argmax = 0;
    argmin = 0;

    dir_0 = d[0];
    dir_1 = d[1];
    dir_2 = d[2];
    for (i = 0; i < no_points_min - 2; i += 3) {
      /**
      * 1) Function Inlining
      * 2) Scalar Replacement
      * 3) Loop unrolling
      * */

      // Load
      pts_i0 = obj1->points[i][0];
      pts_i1 = obj1->points[i][1];
      pts_i2 = obj1->points[i][2];

      pts_i0_1 = obj1->points[i + 1][0];
      pts_i1_1 = obj1->points[i + 1][1];
      pts_i2_1 = obj1->points[i + 1][2];

      pts_i0_2 = obj1->points[i + 2][0];
      pts_i1_2 = obj1->points[i + 2][1];
      pts_i2_2 = obj1->points[i + 2][2];

      invsup_pts_i0 = obj2->points[i][0];
      invsup_pts_i1 = obj2->points[i][1];
      invsup_pts_i2 = obj2->points[i][2];

      invsup_pts_i0_1 = obj2->points[i + 1][0];
      invsup_pts_i1_1 = obj2->points[i + 1][1];
      invsup_pts_i2_1 = obj2->points[i + 1][2];

      invsup_pts_i0_2 = obj2->points[i + 2][0];
      invsup_pts_i1_2 = obj2->points[i + 2][1];
      invsup_pts_i2_2 = obj2->points[i + 2][2];

      // Calculate
      prd_0 = pts_i0 * dir_0;
      prd_1 = pts_i1 * dir_1;
      sum_0 = prd_0 + prd_1;
      prd_2 = pts_i2 * dir_2;

      prd_0_1 = pts_i0_1 * dir_0;
      prd_1_1 = pts_i1_1 * dir_1;
      sum_0_1 = prd_0_1 + prd_1_1;
      prd_2_1 = pts_i2_1 * dir_2;

      prd_0_2 = pts_i0_2 * dir_0;
      prd_1_2 = pts_i1_2 * dir_1;
      sum_0_2 = prd_0_2 + prd_1_2;
      prd_2_2 = pts_i2_2 * dir_2;

      invsup_prd_0 = invsup_pts_i0 * dir_0;
      invsup_prd_1 = invsup_pts_i1 * dir_1;
      invsup_sum_0 = invsup_prd_0 + invsup_prd_1;
      invsup_prd_2 = invsup_pts_i2 * dir_2;

      invsup_prd_0_1 = invsup_pts_i0_1 * dir_0;
      invsup_prd_1_1 = invsup_pts_i1_1 * dir_1;
      invsup_sum_0_1 = invsup_prd_0_1 + invsup_prd_1_1;
      invsup_prd_2_1 = invsup_pts_i2_1 * dir_2;

      invsup_prd_0_2 = invsup_pts_i0_2 * dir_0;
      invsup_prd_1_2 = invsup_pts_i1_2 * dir_1;
      invsup_sum_0_2 = invsup_prd_0_2 + invsup_prd_1_2;
      invsup_prd_2_2 = invsup_pts_i2_2 * dir_2;

      // Store
      mvm_result_1[i] = sum_0 + prd_2;
      mvm_result_1[i + 1] = sum_0_1 + prd_2_1;
      mvm_result_1[i + 2] = sum_0_2 + prd_2_2;

      mvm_result_2[i] = invsup_sum_0 + invsup_prd_2;
      mvm_result_2[i + 1] = invsup_sum_0_1 + invsup_prd_2_1;
      mvm_result_2[i + 2] = invsup_sum_0_2 + invsup_prd_2_2;
    }

    int j = i;
    for (; i < obj1->num_points - 2; i += 3) {
      /**
      * 1) Function Inlining
      * 2) Scalar Replacement
      * 3) Loop unrolling
      * */

      // Load
      pts_i0 = obj1->points[i][0];
      pts_i1 = obj1->points[i][1];
      pts_i2 = obj1->points[i][2];

      pts_i0_1 = obj1->points[i + 1][0];
      pts_i1_1 = obj1->points[i + 1][1];
      pts_i2_1 = obj1->points[i + 1][2];

      pts_i0_2 = obj1->points[i + 2][0];
      pts_i1_2 = obj1->points[i + 2][1];
      pts_i2_2 = obj1->points[i + 2][2];

      // Calculate
      prd_0 = pts_i0 * dir_0;
      prd_1 = pts_i1 * dir_1;
      sum_0 = prd_0 + prd_1;
      prd_2 = pts_i2 * dir_2;

      prd_0_1 = pts_i0_1 * dir_0;
      prd_1_1 = pts_i1_1 * dir_1;
      sum_0_1 = prd_0_1 + prd_1_1;
      prd_2_1 = pts_i2_1 * dir_2;

      prd_0_2 = pts_i0_2 * dir_0;
      prd_1_2 = pts_i1_2 * dir_1;
      sum_0_2 = prd_0_2 + prd_1_2;
      prd_2_2 = pts_i2_2 * dir_2;

      // Store
      mvm_result_1[i] = sum_0 + prd_2;
      mvm_result_1[i + 1] = sum_0_1 + prd_2_1;
      mvm_result_1[i + 2] = sum_0_2 + prd_2_2;
    }

    // Handle the remaining cases
    for (; i < obj1->num_points; i++) {
      // Load
      pts_i0 = obj1->points[i][0];
      pts_i1 = obj1->points[i][1];
      pts_i2 = obj1->points[i][2];

      // Calculate
      prd_0 = pts_i0 * dir_0;
      prd_1 = pts_i1 * dir_1;
      prd_2 = pts_i2 * dir_2;
      sum_0 = prd_0 + prd_1;

      // Store
      mvm_result_1[i] = sum_0 + prd_2;
    }

    for (; j < obj2->num_points - 2; j += 3) {
      /**
      * 1) Function Inlining
      * 2) Scalar Replacement
      * 3) Loop unrolling
      * */

      // Load
      invsup_pts_i0 = obj2->points[j][0];
      invsup_pts_i1 = obj2->points[j][1];
      invsup_pts_i2 = obj2->points[j][2];

      invsup_pts_i0_1 = obj2->points[j + 1][0];
      invsup_pts_i1_1 = obj2->points[j + 1][1];
      invsup_pts_i2_1 = obj2->points[j + 1][2];

      invsup_pts_i0_2 = obj2->points[j + 2][0];
      invsup_pts_i1_2 = obj2->points[j + 2][1];
      invsup_pts_i2_2 = obj2->points[j + 2][2];

      // Calculate
      invsup_prd_0 = invsup_pts_i0 * dir_0;
      invsup_prd_1 = invsup_pts_i1 * dir_1;
      invsup_sum_0 = invsup_prd_0 + invsup_prd_1;
      invsup_prd_2 = invsup_pts_i2 * dir_2;

      invsup_prd_0_1 = invsup_pts_i0_1 * dir_0;
      invsup_prd_1_1 = invsup_pts_i1_1 * dir_1;
      invsup_sum_0_1 = invsup_prd_0_1 + invsup_prd_1_1;
      invsup_prd_2_1 = invsup_pts_i2_1 * dir_2;

      invsup_prd_0_2 = invsup_pts_i0_2 * dir_0;
      invsup_prd_1_2 = invsup_pts_i1_2 * dir_1;
      invsup_sum_0_2 = invsup_prd_0_2 + invsup_prd_1_2;
      invsup_prd_2_2 = invsup_pts_i2_2 * dir_2;

      // Store
      mvm_result_2[j] = invsup_sum_0 + invsup_prd_2;
      mvm_result_2[j + 1] = invsup_sum_0_1 + invsup_prd_2_1;
      mvm_result_2[j + 2] = invsup_sum_0_2 + invsup_prd_2_2;
    }

    // Handle the remaining cases
    for (; j < obj2->num_points; j++) {
      // Load
      invsup_pts_i0 = obj2->points[j][0];
      invsup_pts_i1 = obj2->points[j][1];
      invsup_pts_i2 = obj2->points[j][2];

      // Calculate
      invsup_prd_0 = invsup_pts_i0 * dir_0;
      invsup_prd_1 = invsup_pts_i1 * dir_1;
      invsup_prd_2 = invsup_pts_i2 * dir_2;
      invsup_sum_0 = invsup_prd_0 + invsup_prd_1;

      // Store
      mvm_result_2[j] = invsup_sum_0 + invsup_prd_2;
    }

    for (i = 0; i < no_points_min - 2; i += 3) {
      dp = mvm_result_1[i];
      dp_1 = mvm_result_1[i + 1];
      dp_2 = mvm_result_1[i + 2];

      invsup_dp = mvm_result_2[i];
      invsup_dp_1 = mvm_result_2[i + 1];
      invsup_dp_2 = mvm_result_2[i + 2];

      if (dp > dpmax) {
        dpmax = dp;
        argmax = i;
      }

      if (dp_1 > dpmax) {
        dpmax = dp_1;
        argmax = i + 1;
      }

      if (dp_2 > dpmax) {
        dpmax = dp_2;
        argmax = i + 2;
      }

      if (invsup_dp < invsup_dpmax) {
        invsup_dpmax = invsup_dp;
        argmin = i;
      }

      if (invsup_dp_1 < invsup_dpmax) {
        invsup_dpmax = invsup_dp_1;
        argmin = i + 1;
      }

      if (invsup_dp_2 < invsup_dpmax) {
        invsup_dpmax = invsup_dp_2;
        argmin = i + 2;
      }
    }

    j = i;
    for (; i < obj1->num_points - 2; i += 3) {
      dp = mvm_result_1[i];
      dp_1 = mvm_result_1[i + 1];
      dp_2 = mvm_result_1[i + 2];

      if (dp > dpmax) {
        dpmax = dp;
        argmax = i;
      }

      if (dp_1 > dpmax) {
        dpmax = dp_1;
        argmax = i + 1;
      }

      if (dp_2 > dpmax) {
        dpmax = dp_2;
        argmax = i + 2;
      }
    }

    for (; i < obj1->num_points; i++) {
      dp = mvm_result_1[i];

      if (dp > dpmax) {
        dpmax = dp;
        argmax = i;
      }
    }

    for (; j < obj2->num_points - 2; j += 3) {
      invsup_dp = mvm_result_2[j];
      invsup_dp_1 = mvm_result_2[j + 1];
      invsup_dp_2 = mvm_result_2[j + 2];

      if (invsup_dp < invsup_dpmax) {
        invsup_dpmax = invsup_dp;
        argmin = j;
      }

      if (invsup_dp_1 < invsup_dpmax) {
        invsup_dpmax = invsup_dp_1;
        argmin = j + 1;
      }

      if (invsup_dp_2 < invsup_dpmax) {
        invsup_dpmax = invsup_dp_2;
        argmin = j + 2;
      }
    }

    for (; j < obj2->num_points; j++) {
      invsup_dp = mvm_result_2[j];

      if (invsup_dp < invsup_dpmax) {
        invsup_dpmax = invsup_dp;
        argmin = j;
      }
    }

    a[0] = obj1->points[argmax][0] - obj2->points[argmin][0];
    a[1] = obj1->points[argmax][1] - obj2->points[argmin][1];
    a[2] = obj1->points[argmax][2] - obj2->points[argmin][2];

    // Are we past the origin? If no, we will never enclose it.
    if (dot3D_optimized(a, d) < 0.) {
      free(mvm_result_1);
      free(mvm_result_2);
      return 0;
    }

    // Add newly found point to Simplex
    // AT THIS POINT WE KNOW a IS PAST THE ORIGIN
    s.p[s.num_points][0] = a[0];
    s.p[s.num_points][1] = a[1];
    s.p[s.num_points][2] = a[2];
    s.num_points++;

    if (do_simplex3D_optimized(&s, d)) {
      free(mvm_result_1);
      free(mvm_result_2);
      return 1;
    }
  }

  return 1; // Most likely not reachable. Say we intersect by default.
}
