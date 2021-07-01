#include "gjk.h"
#include <stdio.h>

// I miss my C++ templates...
//#define SWAP(a,b,tmp) tmp = b; b = a; a = tmp

/**
 * Dot product of two 3D vectors.
 * */
float dot3D_soa(float x0, float x1, float x2, const float *y) {
  return x0 * y[0] + x1 * y[1] + x2 * y[2];
}

/**
 * Cross product of two 3D vectors. Write into res.
 * */
void cross3D_soa(const float *x, const float *y, float *res) {
  res[0] = x[1] * y[2] - x[2] * y[1];
  res[1] = x[2] * y[0] - x[0] * y[2];
  res[2] = x[0] * y[1] - x[1] * y[0];
}

/**
 * Among a set of N points in 3D, find the point furthest along dir.
 * TODO can we do this more efficiently than brute force?
 * */
void support3D_soa(float** points, const float *dir, float *res, int N) {
  float dp, dpmax = -1e9;
  int argmax = 0;
  for (int i = 0; i < N; ++i) {
    dp = dot3D_soa(points[0][i], points[1][i], points[2][i], dir);
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
void invsup3D_soa(float** points, const float *dir, float *res, int N) {
  float dp, dpmax = 1e9;
  int argmin = 0;
  for (int i = 0; i < N; ++i) {
    dp = dot3D_soa(points[0][i], points[1][i], points[2][i], dir);
    if (dp < dpmax) {
      dpmax = dp;
      argmin = i;
    }
  }
  res[0] -= points[0][argmin];
  res[1] -= points[1][argmin];
  res[2] -= points[2][argmin];
}

int ds_line3D_soa(struct Simplex3D* s, float *dir) {
  // WE KNOW A is PAST THE ORIGIN -> new Simplex s is [A, B]
  // Unlike in the video by Casey, I don't think we need to permute A and B.
  // p[0] is B, p[1] is A (the new point)
  const float(* p)[3] = s->p;
  float ab[3] = {p[0][0] - p[1][0], p[0][1] - p[1][1], p[0][2] - p[1][2]};
  float a0[3] = {-p[1][0], -p[1][1], -p[1][2]};

  // Perform the float-cross-product
  float tmp[3];
  cross3D_soa(ab, a0, tmp);
  cross3D_soa(tmp, ab, dir);

  // We don't (know if we) enclose 0
  return 0;
}

int ds_triangle3D_soa(struct Simplex3D* s, float *dir) {
  float(* p)[3] = s->p;
  // p[0] is B, p[1] is C, p[2] is A (the new point)
  float ab[3] = {p[0][0] - p[2][0], p[0][1] - p[2][1], p[0][2] - p[2][2]};
  float ac[3] = {p[1][0] - p[2][0], p[1][1] - p[2][1], p[1][2] - p[2][2]};
  float a0[3] = {-p[2][0], -p[2][1], -p[2][2]};
  float abc[3];
  cross3D_soa(ab, ac, abc);

  float tst[3];
  cross3D_soa(abc, ac, tst);
  if (dot3D(tst, a0) > 0.0) {
    // [A, C]
    s->num_points = 2;
    p[0][0] = p[2][0];
    p[0][1] = p[2][1];
    p[0][2] = p[2][2];

    cross3D_soa(ac, a0, tst);
    cross3D_soa(tst, ac, dir);
  } else {
    cross3D_soa(ab, abc, tst);
    if (dot3D(tst, a0) > 0.0) {
      // [A, B]
      s->num_points = 2;
      // I don't think we need to permute from [B, A] to [A, B]
      p[1][0] = p[2][0];
      p[1][1] = p[2][1];
      p[1][2] = p[2][2];

      cross3D_soa(ab, a0, tst);
      cross3D_soa(tst, ab, dir);
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

int ds_tetrahedron3D_soa(struct Simplex3D* s, float *dir) {
  float(* p)[3] = s->p;
  // p[0] is B, p[1] is C, p[2] is D, p[3] is A (the new point)
  float ab[3] = {p[0][0] - p[3][0], p[0][1] - p[3][1], p[0][2] - p[3][2]};
  float ac[3] = {p[1][0] - p[3][0], p[1][1] - p[3][1], p[1][2] - p[3][2]};
  float ad[3] = {p[2][0] - p[3][0], p[2][1] - p[3][1], p[2][2] - p[3][2]};
  float a0[3] = {-p[3][0], -p[3][1], -p[3][2]};
  float abc[3], acd[3], adb[3];
  float tst[3];

  // The vectors abc, acd and adb point OUTSIDE the tetrahedron.
  cross3D_soa(ab, ac, abc);
  if (dot3D(abc, a0) > 0.0) {
    cross3D_soa(ac, ad, acd);
    if (dot3D(acd, a0) > 0.0) {
      s->num_points = 2;
      // [A, C]
      p[0][0] = p[3][0];
      p[0][1] = p[3][1];
      p[0][2] = p[3][2];

      // Like in the line case (I think ?)
      cross3D_soa(ac, a0, tst);
      cross3D_soa(tst, ac, dir);
    }
    else {
      cross3D_soa(ad, ab, adb);
      if (dot3D(adb, a0) > 0.0) {
        s->num_points = 2;
        // [A, B]  (tho we return [B, A] as line perm doesn't matter).
        p[1][0] = p[3][0];
        p[1][1] = p[3][1];
        p[1][2] = p[3][2];

        // Like in the line case (I think ?)
        cross3D_soa(ab, a0, tst);
        cross3D_soa(tst, ab, dir);
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
    cross3D_soa(ac, ad, acd);
    cross3D_soa(ad, ab, adb);
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
        cross3D_soa(ad, a0, tst);
        cross3D_soa(tst, ad, dir);
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
int do_simplex3D_soa(struct Simplex3D* s, float *dir) {
  switch (s->num_points) {
    case 2:return ds_line3D_soa(s, dir);
    case 3:return ds_triangle3D_soa(s, dir);
    case 4:return ds_tetrahedron3D_soa(s, dir);
    default:;
  }
  return -1;  // unreachable unless s points to an invalid Simplex
}

int do_intersect3D_soa(const struct CHObject_soa* obj1, const struct CHObject_soa* obj2) {
  // The search direction
  float d[3] = {1.f, 1.f, 1.f};  // could also be random
  // The point we examine in each iteration
  float a[3];

  // Our simplex that do_simplex modifies.
  struct Simplex3D s;
  s.num_points = 1;

  // S = Support(A-B) = Support(A) - Invsup(B)
  support3D_soa(obj1->points, d, s.p[0], obj1->num_points);
  invsup3D_soa(obj2->points, d, s.p[0], obj2->num_points);

  // d = -S
  d[0] = -s.p[0][0];
  d[1] = -s.p[0][1];
  d[2] = -s.p[0][2];

  int max_iter = 100;  // most likely we need less
  while (max_iter--) {
    // a = Support(A-B) = Support(A) - Invsup(B)
    support3D_soa(obj1->points, d, a, obj1->num_points);
    invsup3D_soa(obj2->points, d, a, obj2->num_points);

    // Are we past the origin? If no, we will never enclose it.
    if (dot3D(a, d) < 0.) return 0;

    // Add newly found point to Simplex
    // AT THIS POINT WE KNOW a IS PAST THE ORIGIN
    s.p[s.num_points][0] = a[0];
    s.p[s.num_points][1] = a[1];
    s.p[s.num_points][2] = a[2];
    s.num_points++;

    if (do_simplex3D_soa(&s, d)) return 1;
  }

  return 1; // Most likely not reachable. Say we intersect by default.
}

