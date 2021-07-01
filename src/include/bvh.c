#include <assert.h>
#include <immintrin.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/param.h>
#include "bvh.h"
#include "src/geometry/TriangleTests.h"
#include <stdio.h>


/************************ kDOP related [Start] ************************/

__always_inline void compute_kDOP_SoA(struct BVHPrimitive* prims[BVH_DEGREE],
                                      struct Vec3f (* const faces)[3],
                                      float result[BVH_kDOP][BVH_DEGREE],
                                      const int n_prims[BVH_DEGREE]) {
  float support, inv_sup, res1, res2, res3, ms1, m23, is1, i23;
  const int half_kDOP = BVH_kDOP / 2;
  for (int i = 0; i < half_kDOP; ++i) {
    for (int j = 0; j < BVH_DEGREE; ++j) {
      support = FLT_MIN;
      inv_sup = FLT_MAX;
      result[i][j] = support;
      result[half_kDOP + i][j] = inv_sup;
      if (n_prims[j] <= 0) {
        printf("KDOP CLUSTER FUCKED UP!!!!!!\n");
      }
      for (int k = 0; k < n_prims[j]; ++k) {
        int ik = prims[j][k].triangle_index;
        res1 = vdot3D(pentakis[i], faces[ik][0]);
        res2 = vdot3D(pentakis[i], faces[ik][1]);
        res3 = vdot3D(pentakis[i], faces[ik][2]);

        // Compute the max value along this direction
        ms1 = fmaxf(support, res1);
        m23 = fmaxf(res2, res3);
        support = fmaxf(ms1, m23);

        // Compute the min value along this same direction
        // (the negative of the max. value along the negative direction)
        is1 = fminf(inv_sup, res1);
        i23 = fminf(res2, res3);
        inv_sup = fminf(is1, i23);
      }
      result[i][j] = support;
      result[half_kDOP + i][j] = inv_sup;
    }
  }
}

/************************ kDOP related [End] ************************/

/************************ Neural Gas Algorithm [Start] ************************/

__always_inline float h_lambda(float rank, float lambda) {
  return fmaxf(expf(-rank / lambda), 1e-6f);
}

__always_inline void ng_step(struct BVHPrimitive* prims,
                             struct Vec3f prototypes[BVH_DEGREE],
                             float (* res)[BVH_DEGREE],
                             float (* ranks)[BVH_DEGREE],
                             float (* hl)[BVH_DEGREE],
                             int n_prims,
                             float t) {
  // Compute prototype <--> datapoint distances
  for (int i = 0; i < BVH_DEGREE; ++i) {
    for (int j = 0; j < n_prims; ++j) {
      res[j][i] = l2_dist3D(prototypes[i], prims[j].center);
    }
  }

  // Find the ranks of the prototypes w.r.t. each data point
  for (int i = 0; i < n_prims; ++i) {
    ranks[i][0] = 0.f;
    ranks[i][1] = 0.f;
    ranks[i][2] = 0.f;
    ranks[i][3] = 0.f;
    for (int j = 0; j < BVH_DEGREE; ++j) {
      for (int k = j + 1; k < BVH_DEGREE; ++k) {
        if (res[i][j] > res[i][k]) ranks[i][j] = ranks[i][j] + 1.f;
        else ranks[i][k] = ranks[i][k] + 1.f;
      }
    }
  }

  // Compute current lambda
  float n = (float) n_prims;
  float lambda_t = 0.5f * n * powf(0.02f / n, t / T_MAX);
  lambda_t = fmaxf(lambda_t, 1e-6f);

  // Compute the h_lambda values
  for (int i = 0; i < n_prims; ++i) {
    for (int j = 0; j < BVH_DEGREE; ++j) {
      hl[i][j] = h_lambda(ranks[i][j], lambda_t);
    }
  }

  // Compute and write the new prototype positions
  float normalizer[BVH_DEGREE] = {0.f, 0.f, 0.f, 0.f};
  for (int i = 0; i < BVH_DEGREE; ++i) {
    prototypes[i].x = 0.f;
    prototypes[i].y = 0.f;
    prototypes[i].z = 0.f;
  }
  for (int i = 0; i < n_prims; ++i) {
    for (int j = 0; j < BVH_DEGREE; ++j) {
      prototypes[j] = vadd(prototypes[j], smul(prims[i].center, hl[i][j]));
      normalizer[j] += hl[i][j];
    }
  }
  for (int j = 0; j < BVH_DEGREE; ++j) {
    prototypes[j] = sdiv(prototypes[j], normalizer[j]);
  }
}

__always_inline void neural_gas(struct BVHPrimitive* prims,
                                struct BVHPrimitive* spare,
                                struct Vec3f prototypes[BVH_DEGREE],
                                float (* res)[BVH_DEGREE],
                                float (* ranks)[BVH_DEGREE],
                                float (* hl)[BVH_DEGREE],
                                int cluster_sizes[BVH_DEGREE],
                                int n_prims) {
  struct Vec3f previous_protos[BVH_DEGREE];

  // init the prototypes: first point, last point, second point, second-to-last point
  prototypes[0] = prims[0].center;
  prototypes[1] = prims[n_prims / 3].center;
  prototypes[2] = prims[2 * n_prims / 3].center;
  prototypes[3] = prims[n_prims - 1].center;

  for (int round = 0; round < ROUND_MAX; ++round) {
    memcpy(previous_protos, prototypes, BVH_DEGREE * sizeof(struct Vec3f));
    ng_step(prims, prototypes, res, ranks, hl, n_prims, (float) round);

    float movep1 = l2_dist3D(previous_protos[0], prototypes[0]);
    float movep2 = l2_dist3D(previous_protos[1], prototypes[1]);
    float movep3 = l2_dist3D(previous_protos[2], prototypes[2]);
    float movep4 = l2_dist3D(previous_protos[3], prototypes[3]);
    float movement = movep1 + movep2 + movep3 + movep4;
    // TODO in the paper eps = 10^-5 * BoundingBoxSize, our value is static BVH_EPS.
    if (movement < BVH_EPS) {

#ifdef BVH_VERBOSE
      printf("\nn=%d Converged in %d round(s).\n", n_prims, 1 + round);
#endif

      goto done;
    }
  }
  done:

  // Compute prototype <--> datapoint distances
  for (int i = 0; i < BVH_DEGREE; ++i) {
    for (int j = 0; j < n_prims; ++j) {
      res[j][i] = l2_dist3D(prototypes[i], prims[j].center);
    }
  }

  // Find the ranks of the prototypes w.r.t. each data point
  for (int i = 0; i < n_prims; ++i) {
    ranks[i][0] = 0.f;
    ranks[i][1] = 0.f;
    ranks[i][2] = 0.f;
    ranks[i][3] = 0.f;
    for (int j = 0; j < BVH_DEGREE; ++j) {
      for (int k = j + 1; k < BVH_DEGREE; ++k) {
        if (res[i][j] > res[i][k]) ranks[i][j] = ranks[i][j] + 1.f;
        else ranks[i][k] = ranks[i][k] + 1.f;
      }
    }
  }

  // Sort the datapoints according to their clusters (closest prototypes)
  int left = 0, right = n_prims - 1;
  int c1, c2, c3, c4, old_r; // count no. of elements in each cluster
  for (int i = 0; i < n_prims; ++i) {
    if (ranks[i][0] == 0.f) {
      spare[left++] = prims[i];
    } else if (ranks[i][3] == 0.f) {
      spare[right--] = prims[i];
    }
  }
  c1 = left;
  old_r = right;
  for (int i = 0; i < n_prims; ++i) {
    if (ranks[i][1] == 0.f) {
      spare[left++] = prims[i];
    } else if (ranks[i][2] == 0.f) {
      spare[right--] = prims[i];
    }
  }
  c2 = left - c1;
  c3 = old_r - right;
  c4 = n_prims - 1 - old_r;

#ifdef BVH_VERBOSE
  printf("c1p=%.2f (%d)  ", ((double)c1) / ((double)n_prims), c1);
  printf("c2p=%.2f (%d)  ", ((double)c2) / ((double)n_prims), c2);
  printf("c3p=%.2f (%d)  ", ((double)c3) / ((double)n_prims), c3);
  printf("c4p=%.2f (%d)\n", ((double)c4) / ((double)n_prims), c4);
#endif

  // Edge-Case handling. Goodluck vectorizing this crap :-(
  while (c1 <= 3) {
    if (c2) --c2;
    else if (c3) --c3;
    else --c4;
    ++c1;
  }
  // c1 now at least 4
  while (c2 <= 2) {
    --c1;
    ++c2;
  }
  // c2 now at least 3
  while (c3 <= 1) {
    --c2;
    ++c3;
  }
  // c3 now at least 2
  if (c4 == 0) {
    --c3;
    ++c4;
  }
  // c4 now at least 1

  // All of them now guaranteed to be >= 1. Even if Neural Gas doesn't converge.
  assert(c1 > 0);
  assert(c2 > 0);
  assert(c3 > 0);
  assert(c4 > 0);
  assert(c1 + c2 + c3 + c4 == n_prims);
  cluster_sizes[0] = c1;
  cluster_sizes[1] = c2;
  cluster_sizes[2] = c3;
  cluster_sizes[3] = c4;

#ifdef BVH_VERBOSE
  printf("c1p=%.2f (%d)  ", ((double)c1) / ((double)n_prims), c1);
  printf("c2p=%.2f (%d)  ", ((double)c2) / ((double)n_prims), c2);
  printf("c3p=%.2f (%d)  ", ((double)c3) / ((double)n_prims), c3);
  printf("c4p=%.2f (%d)\n", ((double)c4) / ((double)n_prims), c4);
#endif

  // Copy the newly sorted primitives back to the original array
  memcpy(prims, spare, n_prims * sizeof(struct BVHPrimitive));
}

/************************ Neural Gas Algorithm [End] ************************/

/************************ BVH construction [Start] ************************/

// TODO optimization: see if recursion or explicit stack is faster
struct BVHPointer construct_rec(struct BVHPrimitive* prims,
                                struct BVHPrimitive* spare,
                                struct Vec3f (* const faces)[3],
                                struct Vec3f prototypes[BVH_DEGREE],
                                float (* res)[BVH_DEGREE],
                                float (* ranks)[BVH_DEGREE],
                                float (* hl)[BVH_DEGREE],
                                int n_prims) {

  // Base case: We are a leaf if we construct over n <= 4 primitives.
  if (n_prims <= BVH_DEGREE) {
    // Set up the leaf region of the idx array
    struct Vec3f (* leafmem)[3] = (struct Vec3f (*)[3]) _mm_malloc(n_prims * sizeof *leafmem, 16);
    for (int i = 0; i < n_prims; ++i) {
      leafmem[i][0] = faces[prims[i].triangle_index][0];
      leafmem[i][1] = faces[prims[i].triangle_index][1];
      leafmem[i][2] = faces[prims[i].triangle_index][2];
    }

    // Separate contiguous memory for each leaf
    return ref_leaf(leafmem, n_prims);
  }

  // Recursive case: build our BVH node, with exactly BVH_DEGREE children
  // Step 1: Run Neural Gas clustering to find child splits
  int cluster_sizes[BVH_DEGREE];
  neural_gas(prims, spare, prototypes, res, ranks, hl, cluster_sizes, n_prims);

  // Step 2: Allocate space for a new BVHNode: the parent of all four children
  struct BVHNode* bnode = (struct BVHNode*) _mm_malloc(sizeof(struct BVHNode), BVH_NODE_ALIGNMENT);

  // Step 3: Split the children. Now the primitive region is contiguous for each child.
  struct BVHPrimitive* subregions[BVH_DEGREE];
  int offset = 0;
  for (int i = 0; i < BVH_DEGREE; ++i) {
    subregions[i] = prims + offset;
    bnode->children[i] = construct_rec(subregions[i],
                                       spare,
                                       faces,
                                       prototypes,
                                       res,
                                       ranks,
                                       hl,
                                       cluster_sizes[i]);
    // Compute new offset
    offset = offset + cluster_sizes[i];
  }

  // Step 4: Compute the bounding kDOP for all children
  compute_kDOP_SoA(subregions, faces, bnode->kdop, cluster_sizes);

  // Return the ready node
  return ref_node(bnode);
}

struct BVHPointer construct(struct TriMesh* mesh) {
  int num_faces = MAX(mesh->num_faces, 1 + BVH_DEGREE);
  struct BVHPrimitive* prims = (struct BVHPrimitive*) _mm_malloc(num_faces * sizeof(struct BVHPrimitive), 16);
  struct BVHPrimitive* spare = (struct BVHPrimitive*) _mm_malloc(num_faces * sizeof(struct BVHPrimitive), 16);
  float (* res)[BVH_DEGREE] = (float (*)[BVH_DEGREE]) _mm_malloc(num_faces * sizeof *res, 16);
  float (* ranks)[BVH_DEGREE] = (float (*)[BVH_DEGREE]) _mm_malloc(num_faces * sizeof *ranks, 16);
  float (* hl)[BVH_DEGREE] = (float (*)[BVH_DEGREE]) _mm_malloc(num_faces * sizeof *hl, 16);

  // Initialize the primitives with the indices and the centers of the triangles
  int i;
  for (i = 0; i < mesh->num_faces; ++i) {
    // center = (face[i][0] + face[i][1] + face[i][2]) / 3.f
    prims[i].center = smul(vadd(mesh->faces[i][0], vadd(mesh->faces[i][1], mesh->faces[i][2])), (1.f / 3.f));
    prims[i].triangle_index = i;
  }

  // Hack: copy the last triangle again (i is fixed as the last valid idx)
  // This means we call construct_rec with at least 5 primitives,
  // which guarantees we never return a leaf from this method.
  for (int j = i; j < num_faces; ++j) {
    // center = (face[i][0] + face[i][1] + face[i][2]) / 3.f
    prims[j].center = smul(vadd(mesh->faces[i][0], vadd(mesh->faces[i][1], mesh->faces[i][2])), (1.f / 3.f));
    prims[j].triangle_index = i;
  }

  // Construct the tree using neural gas clustering as splitting criterion
  struct Vec3f prototypes[BVH_DEGREE];
  struct BVHPointer root = construct_rec(prims,
                                         spare,
                                         mesh->faces,
                                         prototypes,
                                         res,
                                         ranks,
                                         hl,
                                         num_faces);

  // Free the memory we don't need anymore
  _mm_free(prims);
  _mm_free(spare);
  _mm_free(res);
  _mm_free(ranks);
  _mm_free(hl);
  return root;
}

/**
 * Free up a tree's memory.
 *
 * NOTE: CANNOT HAVE A LEAF NODE AS INPUT!
 * [The construct() method guarantees to never return a leaf.]
 *
 * Calling this method will be necessary, as C lacks C++'s destructors.
 * Will recursively destroy all children of the pointed-to node and then the node itself
 * */
void destruct(struct BVHPointer b) {
  assert(!is_leaf(b));
  struct BVHNode* node = get_bvh(b);
  for (int i = 0; i < BVH_DEGREE; ++i) {
    struct BVHPointer child = node->children[i];
    if (!is_leaf(child)) destruct(child);
    else _mm_free(get_leaf(child)); // free the leaf memory
  }
  _mm_free(node);
}

/************************ BVH construction [End] ************************/

/************************ Intersection testing [Start] ************************/
//int forcezero = 0;
__always_inline int intersect_triangle_triangle(struct Vec3f l[3], struct Vec3f r[3]) {
#ifdef BVH_TRICOUNT
  ++bvh_trianglecount;
#endif
#ifdef BVH_COUNT_BYTES
  bvh_movedbytecount += 2 * 3 * sizeof(struct Vec3f);
#endif
// flops counted inside NoDivTriTriIsect()
  //return forcezero;
  // TODO own implementation, maybe use faster algo (for non-coplanar case) from
  //  [Olivier Devillers, Philippe Guigue. Faster Triangle-Triangle Intersection Tests. RR-4488, INRIA. 2002.]
  // TODO vectorize (but for this we need to refactor the whole method)
  return NoDivTriTriIsect(l[0], l[1], l[2], r[0], r[1], r[2]);
}

__always_inline void intersect_leaf_inner(float leafparent_kdop[BVH_kDOP][BVH_DEGREE],
                                          float right[BVH_kDOP][BVH_DEGREE],
                                          int res[BVH_DEGREE],
                                          size_t leaf_index) {
#ifdef BVH_COUNT_BYTES
  bvh_movedbytecount += 4 * BVH_kDOP * BVH_DEGREE;
#endif
  // TODO maybe stream (cache bypassing) loads for leafparent kdop
  const int half_kDOP = BVH_kDOP / 2;
  for (int i = 0; i < BVH_DEGREE; ++i) {
    for (int k = 0; k < half_kDOP; ++k) {
#ifdef BVH_COUNT_FLOPS
      bvh_flopcount += 2;
#endif
      // Max. of l < Min. of right OR Max. of right < Min. of l
      if (leafparent_kdop[k][leaf_index] < right[half_kDOP + k][i]
          || right[k][i] < leafparent_kdop[half_kDOP + k][leaf_index]) {
        res[i] = 0;
        goto cont_inner;
      }
    }
    res[i] = 1; // overlap
    cont_inner:;
  }
}

__always_inline void intersect_kDOP_kDOP(float left[BVH_kDOP][BVH_DEGREE],
                                         float right[BVH_kDOP][BVH_DEGREE],
                                         int res[BVH_DEGREE][BVH_DEGREE]) {
#ifdef BVH_COUNT_BYTES
  bvh_movedbytecount += 4 * 2 * BVH_kDOP * BVH_DEGREE;
#endif

  const int half_kDOP = BVH_kDOP / 2;
  for (int i = 0; i < BVH_DEGREE; ++i) {
    for (int j = 0; j < BVH_DEGREE; ++j) {
      for (int k = 0; k < half_kDOP; ++k) {
#ifdef BVH_COUNT_FLOPS
        bvh_flopcount += 2;
#endif
        // Max. of left < Min. of right OR Max. of right < Min. of left
        if (left[k][i] < right[half_kDOP + k][j] || right[k][j] < left[half_kDOP + k][i]) {
          res[i][j] = 0;
          goto cont_inner;
        }
      }
      res[i][j] = 1; // overlap
      cont_inner:;
    }
  }
}

/**
 * Intersection test two BVH Trees
 * This is the recursive (for now?) helper method.
 * - Independent of DOP-directions, but
 * - Implicitly assumes the same DOP-directions for left and right.
 * */
int do_intersect_rec(struct BVHPointer left,
                     struct BVHPointer right,
                     struct BVHPointer left_parent,
                     struct BVHPointer right_parent,
                     size_t left_index,
                     size_t right_index) {
  if (is_leaf(left) && is_leaf(right)) {
    // triangle - triangle
    size_t left_n = get_n(left), right_n = get_n(right);
    struct Vec3f (* ltri)[3] = (struct Vec3f (*)[3]) get_leaf(left);
    struct Vec3f (* rtri)[3] = (struct Vec3f (*)[3]) get_leaf(right);

    for (size_t i = 0; i < left_n; ++i) {
      for (size_t j = 0; j < right_n; ++j) {
        if (intersect_triangle_triangle(ltri[i], rtri[j]))
          return 1;
      }
    }
    return 0;
  } // return within if block
  if (is_leaf(left)) {
    int intersection_vector[BVH_DEGREE];
    struct BVHNode* lparent = get_bvh(left_parent);
    struct BVHNode* rnode = get_bvh(right);
    intersect_leaf_inner(lparent->kdop, rnode->kdop, intersection_vector, left_index);
    for (int i = 0; i < BVH_DEGREE; ++i) {
      if (intersection_vector[i]) {
        if (do_intersect_rec(left, rnode->children[i], left_parent, right, left_index, i))
          return 1;
      }
    }
    return 0;
  } // return within if block
  if (is_leaf(right)) {
    int intersection_vector[BVH_DEGREE];
    struct BVHNode* rparent = get_bvh(right_parent);
    struct BVHNode* lnode = get_bvh(left);
    intersect_leaf_inner(rparent->kdop, lnode->kdop, intersection_vector, right_index);
    for (int i = 0; i < BVH_DEGREE; ++i) {
      if (intersection_vector[i]) {
        if (do_intersect_rec(right, lnode->children[i], right_parent, left, right_index, i))
          return 1;
      }
    }
    return 0;
  } // return within if block
  // kDOP - kDOP
  struct BVHNode* lnode = get_bvh(left);
  struct BVHNode* rnode = get_bvh(right);
  int intersection_matrix[BVH_DEGREE][BVH_DEGREE];
  intersect_kDOP_kDOP(lnode->kdop, rnode->kdop, intersection_matrix);
  for (int j = 0; j < BVH_DEGREE; ++j) {
    for (int i = 0; i < BVH_DEGREE; ++i) {
      if (intersection_matrix[i][j]) {
        if (do_intersect_rec(lnode->children[i], rnode->children[j], left, right, i, j))
          return 1;
      }
    }
  }
  return 0;
}

/**
 * Intersection test two BVH Trees
 *
 * NOTE: CANNOT HAVE LEAF NODES AS INPUT!
 * [The construct() method guarantees to never return a leaf because of this.]
 *
 * - Independent of DOP-directions, but
 * - Implicitly assumes the same DOP-directions for left and right.
 * */
int do_intersect_bvh(struct BVHPointer left, struct BVHPointer right) {
#ifdef BVH_TRICOUNT
  bvh_trianglecount = 0;
#endif
#ifdef BVH_COUNT_BYTES
  bvh_movedbytecount = 0;
#endif
#ifdef BVH_COUNT_FLOPS
  bvh_flopcount = 0;
#endif

  assert(!is_leaf(left));
  assert(!is_leaf(right));
  // The top recursion level thus never accesses left_parent or right_parent
  int res = do_intersect_rec(left, right, left, right, 0, 0);
#ifdef BVH_TRICOUNT
  printf("Needed %ld triangle-triangle tests.\n", bvh_trianglecount);
#endif
#ifdef BVH_COUNT_BYTES
  printf("Moved %ld bytes Main memory -> cache.\n", bvh_movedbytecount);
#endif
#ifdef BVH_COUNT_FLOPS
  printf("Performed %ld float ops.\n", bvh_flopcount);
#endif
  return res;
}

/************************ Intersection testing [End] ************************/

/************************ Vectorized testing [Start] ************************/
#ifdef NDEBUG

/** Madness begins: as we don't have templates, we need macro. */

#define LEAF_INNER_LOOP(x) for (int i = 0; i < BVH_kDOP / 2; i += 4) { \
/* 8 cycles, in parallel */ \
__m512 left_cl_hi = _mm512_load_ps(leafparent_kdop[i]); \
__m512 left_cl_lo = _mm512_load_ps(leafparent_kdop[(BVH_kDOP / 2) + i]); \
/* 8 cycles + 1 waiting for ports, in parallel */ \
__m512 r_dops_hi = _mm512_load_ps(right[i]); \
__m512 r_dops_lo = _mm512_load_ps(right[(BVH_kDOP / 2) + i]); \
/* 0000,1111,2222,3333, 1 cycle each */ \
__m512 l_dops_hi = _mm512_permute_ps(left_cl_hi, (x)); \
__m512 l_dops_lo = _mm512_permute_ps(left_cl_lo, (x)); \
/* Compute overlaps, 3 cycles each */ \
__mmask16 overlap_lhrl = _mm512_cmp_ps_mask(l_dops_hi, r_dops_lo, _CMP_GE_OS); \
__mmask16 overlap_llrh = _mm512_cmp_ps_mask(l_dops_lo, r_dops_hi, _CMP_LE_OS); \
/* Which left and right children overlap along this kdop direction? */ \
__mmask16 overlap_current_axis = overlap_lhrl & overlap_llrh; \
result = result & overlap_current_axis; \
}

/** Madness ends. */

__always_inline void intersect_leaf_inner_vec(float leafparent_kdop[BVH_kDOP][BVH_DEGREE],
                                              float right[BVH_kDOP][BVH_DEGREE],
                                              int res[BVH_DEGREE],
                                              int leaf_index) {
#ifdef BVH_COUNT_BYTES
  bvh_movedbytecount += 4 * BVH_kDOP * BVH_DEGREE;
#endif
#ifdef BVH_COUNT_FLOPS
  bvh_flopcount += 2 * 4 * 16;
#endif

  // Help the compiler optimize our function better.
  assume(0 <= leaf_index && leaf_index <= 3);

  __mmask16 result = 0xFFFF;
  const __mmask16 b1 = 0x1111;
  const __mmask16 b2 = 0x2222;
  const __mmask16 b3 = 0x4444;
  const __mmask16 b4 = 0x8888;

  // We need this madness as only immediates are accepted for permute_ps
  switch (leaf_index) {
    case 0: LEAF_INNER_LOOP(0);
      break;
    case 1: LEAF_INNER_LOOP(0x55);
      break;
    case 2: LEAF_INNER_LOOP(0xAA);
      break;
    case 3: LEAF_INNER_LOOP(0xFF);
      break;
    default:; // don't merge case 3 and default, that adds 10% runtime.
  }

  // mask has 16 bit. every FOURTH bit has to be exactly 0xF for the res to be 1
  res[0] = _ktestc_mask16_u8(result, b1);
  res[1] = _ktestc_mask16_u8(result, b2);
  res[2] = _ktestc_mask16_u8(result, b3);
  res[3] = _ktestc_mask16_u8(result, b4);
}

/**
 * The vectorized version of the kDOP-kDOP intersection test.
 * The main approach taken is the vectorized 4-vs-4 approach from the paper:
 * "SIMDop: SIMD Optimized Bounding Volume Hierarchies for Collision Detection" by T.Tan, R. Weller and G. Zachmann
 * */
__always_inline void intersect_kDOP_kDOP_vec(float left[BVH_kDOP][BVH_DEGREE],
                                             float right[BVH_kDOP][BVH_DEGREE],
                                             int res[BVH_DEGREE][BVH_DEGREE]) {
#ifdef BVH_COUNT_BYTES
  bvh_movedbytecount += 4 * 2 * BVH_kDOP * BVH_DEGREE;
#endif
#ifdef BVH_COUNT_FLOPS
  bvh_flopcount += 4 * 4 * 2 * 16;
#endif

  __m512i r_mask = _mm512_set_epi32(3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0);
  __m512i l_mask = _mm512_set_epi32(3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0);
  const __m512i offset = _mm512_set_epi32(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4);
  __mmask16 result = 0xFFFF;

  for (int i = 0; i < BVH_kDOP / 2; i += 4) {
    // 8 cycles, in parallel
    __m512 l_hi_group = _mm512_load_ps(left[i]);
    __m512 l_lo_group = _mm512_load_ps(left[(BVH_kDOP / 2) + i]);

    // 8 cycles + 1 waiting for ports, in parallel
    __m512 r_hi_group = _mm512_load_ps(right[i]);
    __m512 r_lo_group = _mm512_load_ps(right[(BVH_kDOP / 2) + i]);

    for (int j = 0; j < 4; ++j) {
      // DOP value permutation like in Fig. 4 of the SIMDop paper.
      __m512 l_hi = _mm512_permutexvar_ps(l_mask, l_hi_group);
      __m512 r_lo = _mm512_permutexvar_ps(r_mask, r_lo_group);

      // DOP value permutation like in Fig. 4 of the SIMDop paper.
      __m512 l_lo = _mm512_permutexvar_ps(l_mask, l_lo_group);
      __m512 r_hi = _mm512_permutexvar_ps(r_mask, r_hi_group);

      // Compute overlaps
      __mmask16 overlap_lhrl = _mm512_cmp_ps_mask(l_hi, r_lo, _CMP_GE_OS);
      __mmask16 overlap_llrh = _mm512_cmp_ps_mask(l_lo, r_hi, _CMP_LE_OS);

      // Which left and right children overlap along this kdop direction?
      __mmask16 overlap_current_axis = overlap_lhrl & overlap_llrh;
      result = result & overlap_current_axis;

      // Get the new masks for permutexvar
      l_mask = _mm512_add_epi32(l_mask, offset);
      r_mask = _mm512_add_epi32(r_mask, offset);
    }
  }
  __m512i res_vec = _mm512_movm_epi32(result);
  _mm512_store_epi32(res, res_vec);
}

/**** Specialized recursion Begin ****/

/**
 * Handle the case where both nodes are leaf nodes.
 * */
__always_inline int case_leaf_leaf(struct BVHPointer left,
                                   struct BVHPointer right) {
  size_t left_n = get_n(left), right_n = get_n(right);
  struct Vec3f (* ltri)[3] = (struct Vec3f (*)[3]) get_leaf(left);
  struct Vec3f (* rtri)[3] = (struct Vec3f (*)[3]) get_leaf(right);

  for (size_t i = 0; i < left_n; ++i) {
    for (size_t j = 0; j < right_n; ++j) {
      if (intersect_triangle_triangle(ltri[i], rtri[j]))
        return 1;
    }
  }
  return 0;
}

// pre-declare for access in case_leaf_inner()
int do_intersect_rec_left_leaf_vec(struct BVHPointer, struct BVHPointer, int);

/**
 * Handle the case where one node is a leaf and the other is an inner node.
 * */
__always_inline int case_leaf_inner(struct BVHPointer leaf_parent,
                                    struct BVHPointer right_node,
                                    int leaf_index) {
#ifdef BVH_COUNT_BYTES
  bvh_movedbytecount += 4 * BVH_kDOP * BVH_DEGREE;
#endif
#ifdef BVH_COUNT_FLOPS
  bvh_flopcount += 4 * 2 * 16;
#endif

  struct BVHNode* lparent = get_bvh(leaf_parent);
  struct BVHNode* rnode = get_bvh(right_node);
  float (* leafparent_kdop)[BVH_DEGREE] = lparent->kdop;
  float (* right)[BVH_DEGREE] = rnode->kdop;

  // Help the compiler optimize our function better.
  assume(0 <= leaf_index && leaf_index <= 3);

  __mmask16 result = 0xFFFF;
  const __mmask16 b1 = 0x1111;
  const __mmask16 b2 = 0x2222;
  const __mmask16 b3 = 0x4444;
  const __mmask16 b4 = 0x8888;

  // We need this madness as only immediates are accepted for permute_ps
  switch (leaf_index) {
    case 0: LEAF_INNER_LOOP(0);
      break;
    case 1: LEAF_INNER_LOOP(0x55);
      break;
    case 2: LEAF_INNER_LOOP(0xAA);
      break;
    case 3: LEAF_INNER_LOOP(0xFF);
      break;
    default:; // don't merge case 3 and default, that adds 10% runtime.
  }

  // mask has 16 bit. every FOURTH bit has to be exactly 0xF for the res to be 1
  if (_ktestc_mask16_u8(result, b1) && do_intersect_rec_left_leaf_vec(leaf_parent, rnode->children[0], leaf_index))
    return 1;
  if (_ktestc_mask16_u8(result, b2) && do_intersect_rec_left_leaf_vec(leaf_parent, rnode->children[1], leaf_index))
    return 1;
  if (_ktestc_mask16_u8(result, b3) && do_intersect_rec_left_leaf_vec(leaf_parent, rnode->children[2], leaf_index))
    return 1;
  if (_ktestc_mask16_u8(result, b4) && do_intersect_rec_left_leaf_vec(leaf_parent, rnode->children[3], leaf_index))
    return 1;
  return 0;
}

/**
 * Specialized recursive traversal:
 * At this point, we know we will never encounter the inner-inner case again.
 * */
int do_intersect_rec_left_leaf_vec(struct BVHPointer left_parent,
                                   struct BVHPointer right,
                                   int left_index) {
  // leaf - leaf
  if (is_leaf(right))
    return case_leaf_leaf(get_bvh(left_parent)->children[left_index], right);
  // leaf - inner
  return case_leaf_inner(left_parent, right, left_index);
}

// pre-declare for access in case_inner_inner()
int do_intersect_rec_vec(struct BVHPointer left,
                         struct BVHPointer right,
                         struct BVHPointer left_parent,
                         struct BVHPointer right_parent,
                         int left_index,
                         int right_index);

/**
 * Handle the case where both nodes left and right are inner nodes.
 * The main approach taken is the vectorized 4-vs-4 approach from the paper:
 * "SIMDop: SIMD Optimized Bounding Volume Hierarchies for Collision Detection" by T.Tan, R. Weller and G. Zachmann
 * */
__always_inline int case_inner_inner(struct BVHPointer left_node,
                                     struct BVHPointer right_node) {
#ifdef BVH_COUNT_BYTES
  bvh_movedbytecount += 4 * 2 * BVH_kDOP * BVH_DEGREE;
#endif
#ifdef BVH_COUNT_FLOPS
  bvh_flopcount += 4 * 4 * 2 * 16;
#endif

  struct BVHNode* lnode = get_bvh(left_node);
  struct BVHNode* rnode = get_bvh(right_node);
  float (* left)[BVH_DEGREE] = lnode->kdop;
  float (* right)[BVH_DEGREE] = rnode->kdop;

  __m512i l_mask = _mm512_set_epi32(3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0);
  __m512i r_mask = _mm512_set_epi32(3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0);
  const __m512i offset = _mm512_set_epi32(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4);
  __mmask16 result = 0xFFFF;

  for (int i = 0; i < BVH_kDOP / 2; i += 4) {
    // 8 cycles, in parallel
    __m512 l_hi_group = _mm512_load_ps(left[i]);
    __m512 l_lo_group = _mm512_load_ps(left[(BVH_kDOP / 2) + i]);

    // 8 cycles + 1 waiting for ports, in parallel
    __m512 r_hi_group = _mm512_load_ps(right[i]);
    __m512 r_lo_group = _mm512_load_ps(right[(BVH_kDOP / 2) + i]);

    for (int j = 0; j < 4; ++j) {
      // DOP value permutation like in Fig. 4 of the SIMDop paper.
      __m512 l_hi = _mm512_permutexvar_ps(l_mask, l_hi_group);
      __m512 r_lo = _mm512_permutexvar_ps(r_mask, r_lo_group);

      // DOP value permutation like in Fig. 4 of the SIMDop paper.
      __m512 l_lo = _mm512_permutexvar_ps(l_mask, l_lo_group);
      __m512 r_hi = _mm512_permutexvar_ps(r_mask, r_hi_group);

      // Compute overlaps
      __mmask16 overlap_lhrl = _mm512_cmp_ps_mask(l_hi, r_lo, _CMP_GE_OS);
      __mmask16 overlap_llrh = _mm512_cmp_ps_mask(l_lo, r_hi, _CMP_LE_OS);

      // Which left and right children overlap along this kdop direction?
      __mmask16 overlap_current_axis = overlap_lhrl & overlap_llrh;
      result = result & overlap_current_axis;

      // Get the new masks for permutexvar
      l_mask = _mm512_add_epi32(l_mask, offset);
      r_mask = _mm512_add_epi32(r_mask, offset);
    }
  }

  for (int i = 0; i < BVH_DEGREE; ++i) {
    for (int j = 0; j < BVH_DEGREE; ++j) {
      if (result & 0x1) {
        if (do_intersect_rec_vec(lnode->children[j], rnode->children[i], left_node, right_node, j, i))
          return 1;
      }
      result >>= 0x1;
    }
  }
  return 0;
}

/**
 * Intersection test two BVH Trees
 * This is the recursive (for now?) helper method.
 * - Independent of DOP-directions, but
 * - Implicitly assumes the same DOP-directions for left and right.
 * */
int do_intersect_rec_vec(struct BVHPointer left,
                         struct BVHPointer right,
                         struct BVHPointer left_parent,
                         struct BVHPointer right_parent,
                         int left_index,
                         int right_index) {
  // Could be just left leaf, or both
  if (is_leaf(left))
    return do_intersect_rec_left_leaf_vec(left_parent, right, left_index);
  // right is the only leaf
  if (is_leaf(right))
    return case_leaf_inner(right_parent, left, right_index);
  // both are inner nodes
  return case_inner_inner(left, right);
}

/**** Specialized recursion End ****/

/**
 * Slow version: does not call specialized methods, only itself.
 * Intersection test two BVH Trees
 * This is the recursive (for now?) helper method.
 * - Independent of DOP-directions, but
 * - Implicitly assumes the same DOP-directions for left and right.
 * */
int do_intersect_rec_vec_slow(struct BVHPointer left,
                              struct BVHPointer right,
                              struct BVHPointer left_parent,
                              struct BVHPointer right_parent,
                              int left_index,
                              int right_index) {
  if (is_leaf(left) && is_leaf(right)) {
    // triangle - triangle
    size_t left_n = get_n(left), right_n = get_n(right);
    struct Vec3f (* ltri)[3] = (struct Vec3f (*)[3]) get_leaf(left);
    struct Vec3f (* rtri)[3] = (struct Vec3f (*)[3]) get_leaf(right);

    for (size_t i = 0; i < left_n; ++i) {
      for (size_t j = 0; j < right_n; ++j) {
        if (intersect_triangle_triangle(ltri[i], rtri[j]))
          return 1;
      }
    }
    return 0;
  } // return within if block
  if (is_leaf(left)) {
    int intersection_vector[BVH_DEGREE];
    struct BVHNode* lparent = get_bvh(left_parent);
    struct BVHNode* rnode = get_bvh(right);
    intersect_leaf_inner_vec(lparent->kdop, rnode->kdop, intersection_vector, left_index);
    for (int i = 0; i < BVH_DEGREE; ++i) {
      if (intersection_vector[i]) {
        if (do_intersect_rec_vec_slow(left, rnode->children[i], left_parent, right, left_index, i))
          return 1;
      }
    }
    return 0;
  } // return within if block
  if (is_leaf(right)) {
    int intersection_vector[BVH_DEGREE];
    struct BVHNode* rparent = get_bvh(right_parent);
    struct BVHNode* lnode = get_bvh(left);
    intersect_leaf_inner_vec(rparent->kdop, lnode->kdop, intersection_vector, right_index);
    for (int i = 0; i < BVH_DEGREE; ++i) {
      if (intersection_vector[i]) {
        if (do_intersect_rec_vec_slow(right, lnode->children[i], right_parent, left, right_index, i))
          return 1;
      }
    }
    return 0;
  } // return within if block
  // kDOP - kDOP
  struct BVHNode* lnode = get_bvh(left);
  struct BVHNode* rnode = get_bvh(right);
  ALIGN(int intersection_matrix[BVH_DEGREE][BVH_DEGREE]);
  intersect_kDOP_kDOP_vec(lnode->kdop, rnode->kdop, intersection_matrix);
  for (int j = 0; j < BVH_DEGREE; ++j) {
    for (int i = 0; i < BVH_DEGREE; ++i) {
      if (intersection_matrix[i][j]) {
        if (do_intersect_rec_vec_slow(lnode->children[i], rnode->children[j], left, right, i, j))
          return 1;
      }
    }
  }
  return 0;
}

/**
 * Intersection test two BVH Trees
 *
 * NOTE: CANNOT HAVE LEAF NODES AS INPUT!
 * [The construct() method guarantees to never return a leaf because of this.]
 *
 * - Independent of DOP-directions, but
 * - Implicitly assumes the same DOP-directions for left and right.
 * */
int do_intersect_bvh_vec(struct BVHPointer left, struct BVHPointer right) {
#ifdef BVH_TRICOUNT
  bvh_trianglecount = 0;
#endif
#ifdef BVH_COUNT_BYTES
  bvh_movedbytecount = 0;
#endif
#ifdef BVH_COUNT_FLOPS
  bvh_flopcount = 0;
#endif

  assert(!is_leaf(left));
  assert(!is_leaf(right));
  // The top recursion level thus never accesses left_parent or right_parent
  int res = do_intersect_rec_vec(left, right, left, right, 0, 0);

#ifdef BVH_TRICOUNT
  printf("Needed %ld triangle-triangle tests.\n", bvh_trianglecount);
#endif
#ifdef BVH_COUNT_BYTES
  printf("Moved %ld bytes Main memory -> cache.\n", bvh_movedbytecount);
#endif
#ifdef BVH_COUNT_FLOPS
  printf("Performed %ld float ops.\n", bvh_flopcount);
#endif
  return res;
}

/************************ Vectorized testing [End] ************************/
#endif
