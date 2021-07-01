#pragma once

#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include "bvhstats.h"
#include "../geometry/LinAlg.h"
#include "../geometry/TriMesh.h"

/// Here lie dragons.
/// The interface to our BVH implementation
/// Can collide two or more TriMesh objects
/// Once again, C++ is dearly missed... dearly

#define BVH_kDOP 32
#define BVH_DEGREE 4

/**
 * Stores a tagged pointer to aligned BVHNode, BVHNode_f16c or leaf data.
 * Pointed-to element MUST be at least 16-byte aligned.
 * */
struct BVHPointer {
  uintptr_t ptr;
};

// 4*8 + 4*32*4 = 544 bytes = 136 bytes per child
struct BVHNode {
  float kdop[BVH_kDOP][BVH_DEGREE];
  struct BVHPointer children[BVH_DEGREE];
};

// 4*8 + 4*16*4 = 288 = 72 bytes per child
struct BVHNode_f16c {
  // kDOP intervals stored using half-precision F16C.
  unsigned short kdop[BVH_kDOP][BVH_DEGREE];
  struct BVHPointer children[BVH_DEGREE];
};

// Object-Oriented C sucks

/// ====================== BVHPointer related methods ======================
/// The minimum byte alignment of either BVHNode or any primitive
const size_t minAlignment = 16;
/// The resulting mask for all bits guaranteed to be 0 by minAlignment
const size_t alignMask = minAlignment - 1;
/// Negative of this mask
const size_t alignNegativeMask = ~alignMask;
/// Leaf Tag
const size_t leafTag = 8;
/// Info Mask
const size_t numberMask = leafTag - 1;

__always_inline struct BVHNode* get_bvh(const struct BVHPointer b) {
  return (struct BVHNode*) (b.ptr & alignNegativeMask);
}

__always_inline struct BVHNode_f16c* get_bvh_f16c(const struct BVHPointer b) {
  return (struct BVHNode_f16c*) (b.ptr & alignNegativeMask);
}

__always_inline void* get_leaf(const struct BVHPointer b) {
  return (void*) (b.ptr & alignNegativeMask);
}

__always_inline size_t get_n(const struct BVHPointer b) {
  return (b.ptr) & numberMask;
}

__always_inline void clear_info(struct BVHPointer* b) {
  b->ptr &= alignNegativeMask;
}

__always_inline size_t is_leaf(struct BVHPointer b) {
  return b.ptr & leafTag;
}

__always_inline struct BVHPointer ref_leaf(void* prim, size_t num_primitives) {
  assert(!((uintptr_t) prim & alignMask));
  struct BVHPointer b = {((uintptr_t) prim) | (leafTag + num_primitives)};
  return b;
}

__always_inline struct BVHPointer ref_node(struct BVHNode* target) {
  assert(!((uintptr_t) target & alignMask));
  struct BVHPointer b = {(uintptr_t) target};
  return b;
}

__always_inline struct BVHPointer ref_node_f16c(struct BVHNode_f16c* target) {
  assert(!((uintptr_t) target & alignMask));
  struct BVHPointer b = {(uintptr_t) target};
  return b;
}

/// =============== BVHNode related methods (including f16c) ===============
/**
 * The main entry point: constructs a kDOP-BVH given a TriMesh
 * Cannot be inline, as it references the static function _mm_malloc()
 *
 * Guarantees to never return a leaf, as required by do_intersect_bvh().
 * */
struct BVHPointer construct(struct TriMesh*);
struct BVHPointer construct_f16c(struct TriMesh*);

/**
 * Free up a tree's memory.
 *
 * NOTE: CANNOT HAVE A LEAF NODE AS INPUT!
 * [The construct() method guarantees to never return a leaf.]
 *
 * Calling this method will be necessary, as C lacks C++'s destructors.
 * Will recursively destroy all children of the pointed-to node and then the node itself
 * */
void destruct(struct BVHPointer);
void destruct_f16c(struct BVHPointer);

/**
 * Intersection test two BVH Trees
 *
 * NOTE: CANNOT HAVE LEAF NODES AS INPUT!
 * [It's fine to call this method with any BVHPointers returned by construct()]
 *
 * - Independent of DOP-directions, but
 * - Implicitly assumes the same DOP-directions for left and right.
 * */
int do_intersect_bvh (struct BVHPointer left, struct BVHPointer right);
int do_intersect_bvh_vec (struct BVHPointer left, struct BVHPointer right);

/**
 * Intersection test two BVH Trees with half-precision kDOPs
 *
 * NOTE: CANNOT HAVE FULL-PRECISION BVH TREES or LEAF NODES AS INPUT!
 * [This method can be called with any BVHPointers returned by construct_f16c()]
 *
 * - Independent of DOP-directions, but
 * - Implicitly assumes the same DOP-directions for left and right.
 * */
int do_intersect_bvh_f16c (struct BVHPointer left, struct BVHPointer right);


/// ========== Relevant only to the internals of bvh.c and bvh_f16.c ==========
// Should not be in the header ideally,
// but everything below is common to both .c files.

#define GOLDEN_RATIO 1.618033988749895f
#define ICOSMALL 0.88705799822f
#define SGR 1.4352899911124037f
#define ISGR 0.5482319928924037f
#define T_MAX 100.f
#define ROUND_MAX 100
#define BVH_EPS 1e-5f
// New for vectorized code: BVHNode has to be 64-byte aligned
#define BVH_NODE_ALIGNMENT 64

#ifndef NDEBUG
#define BVH_VERBOSE
#endif

// https://stackoverflow.com/questions/40447195/can-i-hint-the-optimizer-by-giving-the-range-of-an-integer
#define assume(cond) do { if (!(cond)) __builtin_unreachable(); } while (0)

// https://stackoverflow.com/questions/35011579/fill-constant-floats-in-avx-intrinsics-vec
#ifdef __GNUC__
#define ALIGN(x) x __attribute__((aligned(64)))
#elif defined(_MSC_VER)
#define ALIGN(x) __declspec(align(32))
#endif

/// The vertices of the Pentakis Dodecahedron will form the 32 DOP directions.
/// We have pentakis[i][j] == -pentakis[(i + 16) % 32][j], as needed for kDOPs.
/// The lowest energy configuration of 32 electrons
/// according to Thomson's "Plum Pudding" atomic model.
const struct Vec3f pentakis[] = {
    {ICOSMALL, ICOSMALL, ICOSMALL},
    {ICOSMALL, ICOSMALL, -ICOSMALL},
    {-ICOSMALL, ICOSMALL, ICOSMALL},
    {ICOSMALL, -ICOSMALL, ICOSMALL},
    {0.f, 1.f, GOLDEN_RATIO},
    {GOLDEN_RATIO, 0.f, 1.f},
    {1.f, GOLDEN_RATIO, 0.f},
    {0.f, 1.f, -GOLDEN_RATIO},
    {-GOLDEN_RATIO, 0.f, 1.f},
    {1.f, -GOLDEN_RATIO, 0.f},
    {SGR, ISGR, 0.f},
    {0.f, SGR, ISGR},
    {ISGR, 0.f, SGR},
    {SGR, -ISGR, 0.f},
    {0.f, SGR, -ISGR},
    {-ISGR, 0.f, SGR},
    {-ICOSMALL, -ICOSMALL, -ICOSMALL},
    {-ICOSMALL, -ICOSMALL, ICOSMALL},
    {ICOSMALL, -ICOSMALL, -ICOSMALL},
    {-ICOSMALL, ICOSMALL, -ICOSMALL},
    {0.f, -1.f, -GOLDEN_RATIO},
    {-GOLDEN_RATIO, 0.f, -1.f},
    {-1.f, -GOLDEN_RATIO, 0.f},
    {0.f, -1.f, GOLDEN_RATIO},
    {GOLDEN_RATIO, 0.f, -1.f},
    {-1.f, GOLDEN_RATIO, 0.f},
    {-SGR, -ISGR, 0.f},
    {0.f, -SGR, -ISGR},
    {-ISGR, 0.f, -SGR},
    {-SGR, ISGR, 0.f},
    {0.f, -SGR, ISGR},
    {ISGR, 0.f, -SGR}
};

/// Stores center point and index of one triangle face
struct BVHPrimitive {
  struct Vec3f center;
  int triangle_index;
};
