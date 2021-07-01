#include "src/include/gjk.c"
#include "src/include/gjk_soa.c"
#include "src/include/optimized_gjk.c"
#include "src/include/vectorized_gjk.c"
#include "src/include/vectorized_gjk_soa.c"
#include "benchmarking/common.h"

/**
 * Called by the driver (main.cpp) to register your functions
 * Use add_function(func, description) to add your own functions
 * This function assumes that the optimized gjk implementations are in
 * different files, so, as you implement newer versions, include the
 * appropriate files above.
 * */
void register_functions() {
  add_gjk_function(&do_intersect3D, "do_intersect3D", 1);
  add_gjk_function(&do_intersect3D_optimized, "do_intersect3D_optimized", 1);
  add_gjk_function(&do_intersect3D_optimized_inlined, "do_intersect3D_optimized_inlined", 1);
  add_gjk_function(&do_intersect3D_o_i_lu8, "do_intersect3D_o_i_lu8", 1);
  add_gjk_function(&do_intersect3D_o_i_branch_later, "do_intersect3D_o_i_branch_later", 1);
  add_gjk_function(&do_intersect3D_vectorized, "do_intersect3D_vectorized", 1);
  add_gjk_function(&do_intersect3D_vectorized_inlined, "do_intersect3D_vectorized_inlined", 1);
  add_gjk_soa_function(&do_intersect3D_soa, "do_intersect3D_soa", 1);
  add_gjk_soa_function(&do_intersect3D_soa_v_slow_idx, "do_intersect3D_soa_v_slow_idx", 1);
  add_gjk_soa_function(&do_intersect3D_soa_v_slow_idx_inlined, "do_intersect3D_soa_v_slow_idx_inlined", 1);
  add_gjk_soa_function(&do_intersect3D_soa_vectorized, "do_intersect3D_soa_vectorized", 1);
  add_gjk_soa_function(&do_intersect3D_soa_vectorized_lu2, "do_intersect3D_soa_vectorized_lu2", 1);
  add_gjk_soa_function(&do_intersect3D_soa_vectorized_inlined, "do_intersect3D_soa_vectorized_inlined", 1);
}

