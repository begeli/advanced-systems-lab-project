#include "LinAlg.h"
#include <float.h>
#include <math.h>

///**
// * Dot product of two 3D vectors.
// * */
//__always_inline float dot3D(const float x[3], const float y[3]) {
//  return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
//}
//
///**
// * Dot product of two 2D vectors.
// * */
//__always_inline float dot2D(const float x[2], const float y[2]) {
//  return x[0] * y[0] + x[1] * y[1];
//}
//
///**
// * Cross product of two 3D vectors. Write into res.
// * */
//__always_inline void cross3D(const float x[3], const float y[3], float res[3]) {
//  res[0] = x[1] * y[2] - x[2] * y[1];
//  res[1] = x[2] * y[0] - x[0] * y[2];
//  res[2] = x[0] * y[1] - x[1] * y[0];
//}
//
///**
// * Maximum element of a vector x of n float elements
// * */
//__always_inline float vmax(const float* x, int n) {
//  float m = FLT_MIN;
//  for (int i = 0; i < n; ++i) {
//    m = fmaxf(m, x[i]);
//  }
//  return m;
//}
//
///**
// * Index of maximum element of a vector x of n float elements
// * */
//__always_inline int vargmax(const float* x, int n) {
//  float m = FLT_MIN;
//  int argmax = 0;
//  for (int i = 0; i < n; ++i) {
//    if (x[i] > m) {
//      m = x[i];
//      argmax = i;
//    }
//  }
//  return argmax;
//}
//
///**
// * Skinny Matrix Vector multiplication.
// * */
//__always_inline void matvec_nx3(const float (* mat)[3], const float vec[3], float res[3], int n) {
//  for (int i = 0; i < n; ++i) {
//    res[i] = dot3D(mat[i], vec);
//  }
//}
