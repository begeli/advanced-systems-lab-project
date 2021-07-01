#pragma once

#include <float.h>
#include <stdint.h>
#include <math.h>

///**
// * Dot product of two 3D vectors.
// * */
//__always_inline float dot3D(const float x[3], const float y[3]);
//
///**
// * Dot product of two 2D vectors.
// * */
//__always_inline float dot2D(const float x[2], const float y[2]);
//
///**
// * Cross product of two 3D vectors. Write into res.
// * */
//__always_inline void cross3D(const float x[3], const float y[3], float res[3]);
//
///**
// * Maximum element of a vector x of n float elements
// * */
//__always_inline float vmax(const float*, int);
//
///**
// * Index of maximum element of a vector x of n float elements
// * */
//__always_inline int vargmax(const float*, int);
//
///**
// * Matrix-Vector Product of an (n x 3)-dim. matrix and a 3D vector
// * n is passed as a parameter
// * */
//__always_inline void matvec_nx3(const float (* mat)[3], const float vec[3], float res[3], int n);
//

struct Vec3f {
  float x, y, z;
} __attribute__((packed));

struct Vec3f_SoA {
  float* x;
  float* y;
  float* z;
};

/**
 * Dot product of two 3D vectors.
 * */
__always_inline float vdot3D(const struct Vec3f a, const struct Vec3f b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

/**
 * Dot product of two 2D vectors.
 * */
__always_inline float vdot2D(const float x[2], const float y[2]) {
  return x[0] * y[0] + x[1] * y[1];
}

/**
 * Cross product of two 3D vectors. Write into res.
 * */
__always_inline struct Vec3f vcross3D(const struct Vec3f a, const struct Vec3f b) {
  struct Vec3f res;
  res.x = a.y * b.z - a.z * b.y;
  res.y = a.z * b.x - a.x * b.z;
  res.z = a.x * b.y - a.y * b.x;
  return res;
}

/**
 * Absolute value of a 3D float vector
 * */
__always_inline struct Vec3f vabs3D(const struct Vec3f v) {
  struct Vec3f res;
  res.x = fabsf(v.x);
  res.y = fabsf(v.y);
  res.z = fabsf(v.z);
  return res;
}

__always_inline struct Vec3f vadd(const struct Vec3f a, const struct Vec3f b) {
  struct Vec3f res = {a.x + b.x, a.y + b.y, a.z + b.z};
  return res;
}

__always_inline struct Vec3f vsub(const struct Vec3f a, const struct Vec3f b) {
  struct Vec3f res = {a.x - b.x, a.y - b.y, a.z - b.z};
  return res;
}

__always_inline struct Vec3f vmul(const struct Vec3f a, const struct Vec3f b) {
  struct Vec3f res = {a.x * b.x, a.y * b.y, a.z * b.z};
  return res;
}

__always_inline struct Vec3f dviv(const struct Vec3f a, const struct Vec3f b) {
  struct Vec3f res = {a.x / b.x, a.y / b.y, a.z / b.z};
  return res;
}

__always_inline struct Vec3f sadd(const struct Vec3f a, float f) {
  struct Vec3f res = {a.x + f, a.y + f, a.z + f};
  return res;
}

__always_inline struct Vec3f ssub(const struct Vec3f a, float f) {
  struct Vec3f res = {a.x - f, a.y - f, a.z - f};
  return res;
}

__always_inline struct Vec3f smul(const struct Vec3f a, float f) {
  struct Vec3f res = {a.x * f, a.y * f, a.z * f};
  return res;
}

__always_inline struct Vec3f sdiv(const struct Vec3f a, float f) {
  struct Vec3f res = {a.x / f, a.y / f, a.z / f};
  return res;
}

/**
 * Distance of two vectors
 * */
__always_inline float l2_dist3D(const struct Vec3f a, const struct Vec3f b) {
  struct Vec3f r = vsub(a, b);
  return vdot3D(r, r);
}

/**
 * Maximum element of a vector x of n float elements
 * */
__always_inline float vmax(const float* x, int n) {
  float m = FLT_MIN;
  for (int i = 0; i < n; ++i) {
    m = fmaxf(m, x[i]);
  }
  return m;
}

/**
 * Index of maximum element of a vector x of n float elements
 * */
__always_inline int vargmax(const float* x, int n) {
  float m = FLT_MIN;
  int argmax = 0;
  for (int i = 0; i < n; ++i) {
    if (x[i] > m) {
      m = x[i];
      argmax = i;
    }
  }
  return argmax;
}

/**
 * Skinny Matrix Vector multiplication.
 * */
__always_inline void matvec_nx3(const struct Vec3f* mat, const struct Vec3f vec, float* res, int n) {
  for (int i = 0; i < n; ++i) {
    res[i] = vdot3D(mat[i], vec);
  }
}
