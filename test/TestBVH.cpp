#include <string>
#include <cmath>

#include "test/utils/TestUtils.hpp"
#include "test/utils/TriMesh_Loader.hpp"
#include "src/include/bvh.c"
#include "src/include/bvh_f16c.c"

//#define MEASURE

#define bench(left, right) \
auto cycles = test::rdtsc_bvh( &do_intersect_bvh, (left), (right)); \
std::cout << cycles << " cycles" << std::endl

#define bench_vec(left, right) \
auto cycles = test::rdtsc_bvh( &do_intersect_bvh_vec, (left), (right)); \
std::cout << cycles << " cycles" << std::endl

#define bench_f16c(left, right) \
auto cycles = test::rdtsc_bvh( &do_intersect_bvh_f16c, (left), (right)); \
std::cout << cycles << " cycles" << std::endl

#define BENCH_BVH

#ifdef NDEBUG
#define BVH_TEST_VEC
#define BVH_TEST_F16C
#define BVH_TEST_THICC
#endif
//#define CALIBRATE

//#define BVH_TEST_VERBOSE

/**
 * Test BVH implementation
 *
 * Also includes cycle count measurements
 *
 * Additionally, this group of tests shall help to verify that our
 * conversion from .obj to TriMesh works as expected and as intended.
 *
 * The tests in this file created by us by exporting objects from Blender 3D
 */

namespace test {

namespace {
int result = -1; // global so that function calls are not optimized away.

#include <cstdint>
/* ==================== GNU C and possibly other UNIX compilers ===================== */
#if !defined(WIN32) || defined(__GNUC__)

#if defined(__GNUC__) || defined(__linux__)
#define VOLATILE __volatile__
#define ASM __asm__
#else
/* if we're neither compiling with gcc or under linux, we can hope
         * the following lines work, they probably won't */
#define ASM asm
#define VOLATILE
#endif

#define myInt64 unsigned long long
#define INT32 unsigned int

/* ======================== WIN32 ======================= */
#else

#define myInt64 signed __int64
#define INT32 unsigned __int32

#endif

/* This is the RDTSC timer.
 * RDTSC is an instruction on several Intel and compatible CPUs that Reads the
 * Time Stamp Counter. The Intel manuals contain more information.
 */


#define COUNTER_LO(a) ((a).int32.lo)
#define COUNTER_HI(a) ((a).int32.hi)
#define COUNTER_VAL(a) ((a).int64)

#define COUNTER(a) \
    ((unsigned long long)COUNTER_VAL(a))

#define COUNTER_DIFF(a, b) \
    (COUNTER(a)-COUNTER(b))

/* ==================== GNU C and possibly other UNIX compilers ===================== */
#if !defined(WIN32) || defined(__GNUC__)

typedef union {
  myInt64 int64;
  struct { INT32 lo, hi; } int32;
} tsc_counter;

#define RDTSC(cpu_c) \
      ASM VOLATILE ("rdtsc" : "=a" ((cpu_c).int32.lo), "=d"((cpu_c).int32.hi))
#define CPUID() \
        ASM VOLATILE ("cpuid" : : "a" (0) : "bx", "cx", "dx" )

/* ======================== WIN32 ======================= */
#else

typedef union
    {       myInt64 int64;
            struct {INT32 lo, hi;} int32;
    } tsc_counter;

#define RDTSC(cpu_c)   \
    {       __asm rdtsc    \
            __asm mov (cpu_c).int32.lo,eax  \
            __asm mov (cpu_c).int32.hi,edx  \
    }

#define CPUID() \
    { \
        __asm mov eax, 0 \
        __asm cpuid \
    }

#endif

void init_tsc() {
  ; // no need to initialize anything for x86
}

myInt64 start_tsc(void) {
  tsc_counter start;
  CPUID();
  RDTSC(start);
  return COUNTER_VAL(start);
}

myInt64 stop_tsc(myInt64 start) {
  tsc_counter end;
  RDTSC(end);
  CPUID();
  return COUNTER_VAL(end) - start;
}
} // namespace

using bvh_func = int (*)(struct BVHPointer, struct BVHPointer);

/**
 * Measures the runtime of gjk function
 * */
myInt64 rdtsc_bvh(bvh_func func, struct BVHPointer left, struct BVHPointer right) {
  int i, num_runs;
  myInt64 cycles;
  myInt64 start;
  num_runs = 1;

  /*
   * The CPUID instruction serializes the pipeline.
   * Using it, we can create execution barriers around the code we want to time.
   * The calibrate section is used to make the computation large enough so as to
   * avoid measurements bias due to the timing overhead.
   */
#ifdef CALIBRATE
  while (num_runs < (1 << 20)) {
    start = start_tsc();
    for (i = 0; i < num_runs; ++i) {
      result = func(left, right);
    }
    cycles = stop_tsc(start);

    if (cycles >= 10000000) break;

    num_runs *= 2;
  }
#endif

  start = start_tsc();
  for (i = 0; i < num_runs; ++i) {
    result = func(left, right);
  }

  cycles = stop_tsc(start) / num_runs;
  return cycles;
}

auto runBVH(TriMesh* m1, TriMesh* m2) -> int {  // NOLINT
  struct BVHPointer b1 = construct(m1);
#ifdef BVH_TEST_VERBOSE
  std::cout << "runBVH finished construction 1" << std::endl;
#endif
  struct BVHPointer b2 = construct(m2);
#ifdef BVH_TEST_VERBOSE
  std::cout << "runBVH finished construction 2" << std::endl;
#endif

#ifdef BENCH_BVH
  bench(b1, b2);
  auto res = result;
#else
  auto res = do_intersect_bvh(b1, b2);
#endif

  // Free the memory again
  destruct(b1);
  destruct(b2);

  return res;
}

#ifdef NDEBUG
auto runBVH_vec(TriMesh* m1, TriMesh* m2) -> int {  // NOLINT
  struct BVHPointer b1 = construct(m1);
#ifdef BVH_TEST_VERBOSE
  std::cout << "runBVH_vec finished construction 1" << std::endl;
#endif
  struct BVHPointer b2 = construct(m2);
#ifdef BVH_TEST_VERBOSE
  std::cout << "runBVH_vec finished construction 2" << std::endl;
#endif

#ifdef BENCH_BVH
  bench_vec(b1, b2);
  auto res = result;
#else
  auto res = do_intersect_bvh_vec(b1, b2);
#endif

  // Free the memory again
  destruct(b1);
  destruct(b2);

  return res;
}
#endif

#ifdef NDEBUG
auto runBVH_f16c(TriMesh* m1, TriMesh* m2) -> int {  // NOLINT
  struct BVHPointer b1 = construct_f16c(m1);
#ifdef BVH_TEST_VERBOSE
  std::cout << "runBVH_f16c finished construction 1" << std::endl;
#endif
  struct BVHPointer b2 = construct_f16c(m2);
#ifdef BVH_TEST_VERBOSE
  std::cout << "runBVH_f16c finished construction 2" << std::endl;
#endif

#ifdef BENCH_BVH
  bench_f16c(b1, b2);
  auto res = result;
#else
  auto res = do_intersect_bvh_f16c(b1, b2);
#endif

  // Free the memory again
  destruct_f16c(b1);
  destruct_f16c(b2);

  return res;
}
#endif

TEST(TestBVH, OverlappingCubes) {
  std::cout << "OverlappingCubes" << std::endl;
  std::vector<TriMesh>
      tm = tri_m_direct::load_meshes_to_TriMesh("../../Data/objects/cubes/overlapping/overlapping_cubes_triangulated.obj");

  TriMesh* m1 = &tm[0];
  TriMesh* m2 = &tm[1];

  auto res = runBVH(m1, m2);
//    std::cout << "BVH terminated with output " << res << std::endl;
  EXPECT_EQ(res, true); // the cubes intersect.

  free(m1->faces);
  free(m2->faces);
}

TEST(TestBVH, NonOverlappingIcoInUvSphere) {
  std::cout << "NonOverlappingIcoInUvSphere" << std::endl;
  std::string file = "../../Data/objects/spheres/non-overlapping_ico_in_concave_uv_sphere.obj";
  std::vector<TriMesh> tm = tri_m_direct::load_meshes_to_TriMesh(file);
  TriMesh* uvSphere = &tm[0];
  TriMesh* icoSphere = &tm[1];

  auto res = runBVH(uvSphere, icoSphere);
  EXPECT_EQ(res, false); // the spheres do NOT intersect.

  free(uvSphere->faces);
  free(icoSphere->faces);
}

TEST(TestBVH, TouchingIcoInUvSphere) {
  std::cout << "TouchingIcoInUvSphere" << std::endl;
  std::vector<TriMesh>
      tm =
      tri_m_direct::load_meshes_to_TriMesh(
          "../../Data/objects/spheres/touching_1_vertex_ico_in_concave_uv_sphere_scaled_x100.obj");

  TriMesh* uvSphere = &tm[0];
  TriMesh* icoSphere = &tm[1];

  auto res = runBVH(uvSphere, icoSphere);
//    std::cout << "BVH terminated with output " << res << std::endl;
  EXPECT_EQ(res, true); // the spheres touch in one vertex, so they intersect.

  free(uvSphere->faces);
  free(icoSphere->faces);
}

TEST(TestBVH, OverlappingBarelyIcoInUvSphere) {
  std::cout << "OverlappingBarelyIcoInUvSphere" << std::endl;
  std::vector<TriMesh>
      tm =
      tri_m_direct::load_meshes_to_TriMesh(
          "../../Data/objects/spheres/barely_overlapping_ico_in_concave_uv_sphere_scaled_x100.obj");

  TriMesh* uvSphere = &tm[0];
  TriMesh* icoSphere = &tm[1];

  auto res = runBVH(uvSphere, icoSphere);
//    std::cout << "BVH terminated with output " << res << std::endl;
  EXPECT_EQ(res, true); // the spheres touch in one vertex, so they intersect.

  free(uvSphere->faces);
  free(icoSphere->faces);
}

#ifdef BVH_TEST_VEC

TEST(TestBVH, VecOverlappingCubes) {
  std::cout << "VecOverlappingCubes" << std::endl;
  std::vector<TriMesh>
      tm = tri_m_direct::load_meshes_to_TriMesh("../../Data/objects/cubes/overlapping/overlapping_cubes_triangulated.obj");

  TriMesh* m1 = &tm[0];
  TriMesh* m2 = &tm[1];

  auto res = runBVH_vec(m1, m2);
//    std::cout << "BVH terminated with output " << res << std::endl;
  EXPECT_EQ(res, true); // the cubes intersect.

  free(m1->faces);
  free(m2->faces);
}

TEST(TestBVH, VecNonOverlappingIcoInUvSphere) {
  std::cout << "VecNonOverlappingIcoInUvSphere" << std::endl;
  std::string file = "../../Data/objects/spheres/non-overlapping_ico_in_concave_uv_sphere.obj";
  std::vector<TriMesh> tm = tri_m_direct::load_meshes_to_TriMesh(file);
  TriMesh* uvSphere = &tm[0];
  TriMesh* icoSphere = &tm[1];

  auto res = runBVH_vec(uvSphere, icoSphere);
//    std::cout << "BVH terminated with output " << res << std::endl;
  EXPECT_EQ(res, false); // the spheres do NOT intersect.

  free(uvSphere->faces);
  free(icoSphere->faces);
}

TEST(TestBVH, VecTouchingIcoInUvSphere) {
  std::cout << "VecTouchingIcoInUvSphere" << std::endl;
  std::vector<TriMesh>
      tm =
      tri_m_direct::load_meshes_to_TriMesh("../../Data/objects/spheres/touching_1_vertex_ico_in_concave_uv_sphere.obj");

  TriMesh* uvSphere = &tm[0];
  TriMesh* icoSphere = &tm[1];

  auto res = runBVH_vec(uvSphere, icoSphere);
//    std::cout << "BVH terminated with output " << res << std::endl;
  EXPECT_EQ(res, true); // the spheres touch in one vertex, so they intersect.

  free(uvSphere->faces);
  free(icoSphere->faces);
}

TEST(TestBVH, VecOverlappingBarelyIcoInUvSphere) {
  std::cout << "VecOverlappingBarelyIcoInUvSphere" << std::endl;
  std::vector<TriMesh>
      tm =
      tri_m_direct::load_meshes_to_TriMesh("../../Data/objects/spheres/barely_overlapping_ico_in_concave_uv_sphere.obj");

  TriMesh* uvSphere = &tm[0];
  TriMesh* icoSphere = &tm[1];

  auto res = runBVH_vec(uvSphere, icoSphere);
//    std::cout << "BVH terminated with output " << res << std::endl;
  EXPECT_EQ(res, true); // the spheres touch in one vertex, so they intersect.

  free(uvSphere->faces);
  free(icoSphere->faces);
}

#endif

TEST(TestBVH, PlaneBodyTouchingWings) {
  std::cout << "PlaneBodyTouchingWings" << std::endl;
  std::vector<TriMesh>
      tm =
      tri_m_direct::load_meshes_to_TriMesh("../../Data/objects/complex/airplane_v2_triangularized.obj");

  TriMesh* planebody = &tm[0];
  TriMesh* planewings = &tm[1];

  auto res = runBVH(planebody, planewings);
//    std::cout << "BVH terminated with output " << res << std::endl;
  EXPECT_EQ(res, true);

  free(planebody->faces);
  free(planewings->faces);
}

TEST(TestBVH, PlaneBodyWingsExploded) {
  std::cout << "PlaneBodyWingsExploded" << std::endl;
  std::vector<TriMesh>
      tm =
      tri_m_direct::load_meshes_to_TriMesh("../../Data/objects/complex/airplane_slightly_exploded.obj");

  TriMesh* planebody = &tm[0];
  TriMesh* planewings = &tm[1];

  auto res = runBVH(planebody, planewings);
//    std::cout << "BVH terminated with output " << res << std::endl;
  EXPECT_EQ(res, false);

  free(planebody->faces);
  free(planewings->faces);
}

#ifdef BVH_TEST_VEC

TEST(TestBVH, VecPlaneBodyTouchingWings) {
  std::cout << "VecPlaneBodyTouchingWings" << std::endl;
  std::vector<TriMesh>
      tm =
      tri_m_direct::load_meshes_to_TriMesh("../../Data/objects/complex/airplane_v2_triangularized.obj");

  TriMesh* planebody = &tm[0];
  TriMesh* planewings = &tm[1];

  auto res = runBVH_vec(planebody, planewings);
//    std::cout << "BVH terminated with output " << res << std::endl;
  EXPECT_EQ(res, true);

  free(planebody->faces);
  free(planewings->faces);
}

TEST(TestBVH, VecPlaneBodyWingsExploded) {
  std::cout << "VecPlaneBodyWingsExploded" << std::endl;
  std::vector<TriMesh>
      tm =
      tri_m_direct::load_meshes_to_TriMesh("../../Data/objects/complex/airplane_slightly_exploded.obj");

  TriMesh* planebody = &tm[0];
  TriMesh* planewings = &tm[1];

  auto res = runBVH_vec(planebody, planewings);
//    std::cout << "BVH terminated with output " << res << std::endl;
  EXPECT_EQ(res, false);

  free(planebody->faces);
  free(planewings->faces);
}

#endif

#ifdef BVH_TEST_F16C

TEST(TestBVH, F16COverlappingCubes) {
  std::cout << "F16COverlappingCubes" << std::endl;
  std::vector<TriMesh>
      tm = tri_m_direct::load_meshes_to_TriMesh("../../Data/objects/cubes/overlapping/overlapping_cubes_triangulated.obj");

  TriMesh* m1 = &tm[0];
  TriMesh* m2 = &tm[1];

  auto res = runBVH_f16c(m1, m2);
//    std::cout << "BVH terminated with output " << res << std::endl;
  EXPECT_EQ(res, true); // the cubes intersect.

  free(m1->faces);
  free(m2->faces);
}

TEST(TestBVH, F16CNonOverlappingIcoInUvSphere) {
  std::cout << "F16CNonOverlappingIcoInUvSphere" << std::endl;
  std::string file = "../../Data/objects/spheres/non-overlapping_ico_in_concave_uv_sphere.obj";
  std::vector<TriMesh> tm = tri_m_direct::load_meshes_to_TriMesh(file);
  TriMesh* uvSphere = &tm[0];
  TriMesh* icoSphere = &tm[1];

  auto res = runBVH_f16c(uvSphere, icoSphere);
//    std::cout << "BVH terminated with output " << res << std::endl;
  EXPECT_EQ(res, false); // the spheres do NOT intersect.

  free(uvSphere->faces);
  free(icoSphere->faces);
}

TEST(TestBVH, F16CTouchingIcoInUvSphere) {
  std::cout << "F16CTouchingIcoInUvSphere" << std::endl;
  std::vector<TriMesh>
      tm =
      tri_m_direct::load_meshes_to_TriMesh("../../Data/objects/spheres/touching_1_vertex_ico_in_concave_uv_sphere.obj");

  TriMesh* uvSphere = &tm[0];
  TriMesh* icoSphere = &tm[1];

  auto res = runBVH_f16c(uvSphere, icoSphere);
//    std::cout << "BVH terminated with output " << res << std::endl;
  EXPECT_EQ(res, true); // the spheres touch in one vertex, so they intersect.

  free(uvSphere->faces);
  free(icoSphere->faces);
}

TEST(TestBVH, F16COverlappingBarelyIcoInUvSphere) {
  std::cout << "F16COverlappingBarelyIcoInUvSphere" << std::endl;
  std::vector<TriMesh>
      tm =
      tri_m_direct::load_meshes_to_TriMesh("../../Data/objects/spheres/barely_overlapping_ico_in_concave_uv_sphere.obj");

  TriMesh* uvSphere = &tm[0];
  TriMesh* icoSphere = &tm[1];

  auto res = runBVH_f16c(uvSphere, icoSphere);
//    std::cout << "BVH terminated with output " << res << std::endl;.
  EXPECT_EQ(res, true); // the spheres touch in one vertex, so they intersect.

  free(uvSphere->faces);
  free(icoSphere->faces);
}

TEST(TestBVH, F16CPlaneBodyTouchingWings) {
  std::cout << "F16CPlaneBodyTouchingWings" << std::endl;
  std::vector<TriMesh>
      tm =
      tri_m_direct::load_meshes_to_TriMesh("../../Data/objects/complex/airplane_v2_triangularized.obj");

  TriMesh* planebody = &tm[0];
  TriMesh* planewings = &tm[1];

  auto res = runBVH_f16c(planebody, planewings);
//    std::cout << "BVH terminated with output " << res << std::endl;
  EXPECT_EQ(res, true);

  free(planebody->faces);
  free(planewings->faces);
}

TEST(TestBVH, F16CPlaneBodyWingsExploded) {
  std::cout << "F16CPlaneBodyWingsExploded" << std::endl;
  std::vector<TriMesh>
      tm =
      tri_m_direct::load_meshes_to_TriMesh("../../Data/objects/complex/airplane_slightly_exploded.obj");

  TriMesh* planebody = &tm[0];
  TriMesh* planewings = &tm[1];

  auto res = runBVH_f16c(planebody, planewings);
//    std::cout << "BVH terminated with output " << res << std::endl;
  EXPECT_EQ(res, false);

  free(planebody->faces);
  free(planewings->faces);
}

#endif

#ifdef BVH_TEST_THICC

TEST(TestBVH, ToyotaHoveringAboveBracelet) {
  std::cout << "ToyotaHoveringAboveBracelet" << std::endl;
  std::vector<TriMesh>
      tm =
      tri_m_direct::load_meshes_to_TriMesh("../../Data/thicc/toyota.obj");

  TriMesh* bracelet = &tm[0];
  TriMesh* toyota = &tm[1];

  // 5.4M triangles total
  std::cout << "Baseline" << std::endl;
  auto res_base = runBVH(bracelet, toyota);
  std::cout << "Vectorized" << std::endl;
  auto res_vec = runBVH_vec(bracelet, toyota);
  std::cout << "Vectorized + Quantized" << std::endl;
  auto res_f16c = runBVH_f16c(bracelet, toyota);
  EXPECT_EQ(res_base, false);
  EXPECT_EQ(res_vec, false);
  EXPECT_EQ(res_f16c, false);

  free(bracelet->faces);
  free(toyota->faces);
}

TEST(TestBVH, ToyotaTouchingBracelet) {
  std::cout << "ToyotaTouchingBracelet" << std::endl;
  std::vector<TriMesh>
      tm =
      tri_m_direct::load_meshes_to_TriMesh("../../Data/thicc/toyota-touching.obj");

  TriMesh* bracelet = &tm[0];
  TriMesh* toyota = &tm[1];

  // 5.4M triangles total
  std::cout << "Baseline" << std::endl;
  auto res_base = runBVH(bracelet, toyota);
  std::cout << "Vectorized" << std::endl;
  auto res_vec = runBVH_vec(bracelet, toyota);
  std::cout << "Vectorized + Quantized" << std::endl;
  auto res_f16c = runBVH_f16c(bracelet, toyota);
  EXPECT_EQ(res_base, true);
  EXPECT_EQ(res_vec, true);
  EXPECT_EQ(res_f16c, true);

  free(bracelet->faces);
  free(toyota->faces);
}

TEST(TestBVH, ToyotaHoveringAboveBraceletMultiRes) {
  std::cout << "ToyotaHoveringAboveBraceletMultiRes" << std::endl;
  std::vector<std::string> files{};
  files.emplace_back("../../Data/thicc/toyota.obj");
  files.emplace_back("../../Data/thicc/toyota-90.obj");
  files.emplace_back("../../Data/thicc/toyota-80.obj");
  files.emplace_back("../../Data/thicc/toyota-70.obj");
  files.emplace_back("../../Data/thicc/toyota-60.obj");
  files.emplace_back("../../Data/thicc/toyota-50.obj");
  files.emplace_back("../../Data/thicc/toyota-40.obj");
  files.emplace_back("../../Data/thicc/toyota-30.obj");
  files.emplace_back("../../Data/thicc/toyota-20.obj");
  files.emplace_back("../../Data/thicc/toyota-10.obj");

  for (size_t i = 0; i < files.size(); ++i) {
    std::cout << "Resolution: " << (100 - (10 * i)) << "%" << std::endl;
    std::vector<TriMesh>
        tm = tri_m_direct::load_meshes_to_TriMesh(files[i]);
    TriMesh* bracelet = &tm[0];
    TriMesh* toyota = &tm[1];

    std::cout << "Baseline:" << std::endl;
    auto res_base = runBVH(bracelet, toyota);
    std::cout << "Vectorized:" << std::endl;
    auto res_vec = runBVH_vec(bracelet, toyota);
    std::cout << "Vectorized + Quantized:" << std::endl;
    auto res_f16c = runBVH_f16c(bracelet, toyota);
    EXPECT_EQ(res_base, false);
    EXPECT_EQ(res_vec, false);
    EXPECT_EQ(res_f16c, false);

    free(bracelet->faces);
    free(toyota->faces);
  }
}

TEST(TestBVH, ToyotaTouchingBraceletMultiRes) {
  std::cout << "ToyotaTouchingBraceletMultiRes" << std::endl;
  std::vector<std::string> files{};
  files.emplace_back("../../Data/thicc/toyota-touching.obj");
  files.emplace_back("../../Data/thicc/toyota-touching-90.obj");
  files.emplace_back("../../Data/thicc/toyota-touching-80.obj");
  files.emplace_back("../../Data/thicc/toyota-touching-70.obj");
  files.emplace_back("../../Data/thicc/toyota-touching-60.obj");
  files.emplace_back("../../Data/thicc/toyota-touching-50.obj");
  files.emplace_back("../../Data/thicc/toyota-touching-40.obj");
  files.emplace_back("../../Data/thicc/toyota-touching-30.obj");
  files.emplace_back("../../Data/thicc/toyota-touching-20.obj");
  files.emplace_back("../../Data/thicc/toyota-touching-10.obj");

  for (size_t i = 0; i < files.size(); ++i) {
    std::cout << "Resolution: " << (100 - (10 * i)) << "%" << std::endl;
    std::vector<TriMesh>
        tm = tri_m_direct::load_meshes_to_TriMesh(files[i]);
    TriMesh* bracelet = &tm[0];
    TriMesh* toyota = &tm[1];

    std::cout << "Baseline:" << std::endl;
    auto res_base = runBVH(bracelet, toyota);
    std::cout << "Vectorized:" << std::endl;
    auto res_vec = runBVH_vec(bracelet, toyota);
    std::cout << "Vectorized + Quantized:" << std::endl;
    auto res_f16c = runBVH_f16c(bracelet, toyota);
    EXPECT_EQ(res_base, true);
    EXPECT_EQ(res_vec, true);
    EXPECT_EQ(res_f16c, true);

    free(bracelet->faces);
    free(toyota->faces);
  }
}

#endif

} // namespace test
