#pragma once

#include "src/geometry/CHObject.h"
/// This is the main entry point into our GJK library.
/// Big sad that we have to use C. It hurts my soul not to use modern C++

/// TODO array layout: (dim x n) or (n x dim)?
/// TODO understand the importance of permuting points in the Simplex

#define SWAP(a,b,tmp) tmp = b; b = a; a = tmp

/**
 * Represents our search points in 3D.
 * If 0 is enclosed by this shape we have an intersection.
 * */
struct Simplex3D {
    /// num_points will be:
    /// 1: never
    /// 2: if it's a line
    /// 3: if it's a triangle
    /// 4: if it's a tetrahedron
    int num_points;
    float p[4][3];
};

/**
 * Represents our search points in 2D.
 * If 0 is enclosed by this shape we have an intersection.
 * */
struct Simplex2D {
    /// num_points will be:
    /// 1: never
    /// 2: if it's a line
    /// 3: if it's a triangle
    int num_points;
    float p[3][2];
};

/**
 * Same struct as Simplex3D but in SoA format instead AoS.
 * If 0 is enclosed by this shape we have an intersection.
 * Might be very unnecessary.
 * */
struct Simplex3D_soa {
    /// num_points will be:
    /// 1: never
    /// 2: if it's a line
    /// 3: if it's a triangle
    /// 4: if it's a tetrahedron
    int num_points;
    float p[3][4];
};

struct CHObject_soa{
    /// The number of vertices comprising the convex hull of this object.
    int num_points;
    float** points;
};

/**
 * Do two three-dimensional Objects intersect?
 * This will be our entry point to the 3D boolean GJK Algorithm.
 * @returns 1 if the objects intersect, 0 otherwise.
 * */
int do_intersect3D(const struct CHObject* obj1, const struct CHObject* obj2);

struct Object2D{
  /// The number of vertices comprising the convex hull of this object.
  int num_points;
  const float (* points)[2];
};

/**
 * Do two two-dimensional Objects intersect?
 * This will be our entry point to the 2D boolean GJK Algorithm.
 * @returns 1 if the objects intersect, 0 otherwise.
 * */
int do_intersect2D(const struct Object2D* obj1, const struct Object2D* obj2);

/**
 * TODO scoping: should we pursue gjk distance?
 * This function, if implemented, shall return the distance of two 3D objects.
 * Uses the GJK algorithm, with the improved distance subroutine from:
 * "Improving the GJK algorithm for faster and more reliable
 * distance queries between convex objects" by M. Montanari and E.Barbieri
 * */
float gjk_distance3D(const struct CHObject* obj1, const struct CHObject* obj2);