#pragma once

/**
 * Enable or disable several statistics on the BVH algorithm.
 *
 * All statistics are transparent when disabled,
 * but add additional runtime overhead when enabled.
 * */
//#define BVH_TRACK_STATS

/// Count how many triangle tests are performed
#ifdef BVH_TRACK_STATS
#define BVH_TRICOUNT
#endif
long bvh_trianglecount = 0;

/// Measure float memory movement (in bytes) RAM <-> cache
#ifdef BVH_TRACK_STATS
#define BVH_COUNT_BYTES
#endif
long bvh_movedbytecount = 0;

/// Count how many float operations are performed
#ifdef BVH_TRACK_STATS
#define BVH_COUNT_FLOPS
#endif
long bvh_flopcount = 0;
