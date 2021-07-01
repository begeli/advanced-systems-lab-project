#include <string>
#include <cmath>

#include "test/utils/TestUtils.hpp"
#include "benchmarking/performance_tests.cpp"
#include "src/include/gjk.c"
#include "src/include/optimized_gjk.c"
#include "src/include/vectorized_gjk.c"
#include "src/include/gjk_soa.c"
#include "src/include/vectorized_gjk_soa.c"

/// Float accuracy: defining EPS < 1e-7 leads to test failures
#define EPS 1e-7
#define MEASURE

#define bench(s, o1, o2, n1, n2) \
std::cout << (s) << std::endl; \
auto cycles = rdtsc_gjk( &do_intersect3D, (o1), (o2), (n1), (n2)); \
std::cout << cycles << " cycles" << std::endl

/// The first 19 tests in this file are inspired by BourgondAries' tests for his Rust version of bGJK
/// The link to the file: https://github.com/BourgondAries/bgjk/blob/master/src/lib.rs
namespace test {

auto runGJK(const float (* o1)[3], const float (* o2)[3], int n1, int n2) -> int {  // NOLINT
  CHObject obj1{n1, o1};
  CHObject obj2{n2, o2};

  return do_intersect3D(&obj1, &obj2)
        & do_intersect3D_optimized(&obj1, &obj2)
        & do_intersect3D_optimized_inlined(&obj1, &obj2)
        & do_intersect3D_o_i_lu8(&obj1, &obj2)
        & do_intersect3D_o_i_branch_later(&obj1, &obj2)
        & do_intersect3D_vectorized(&obj1, &obj2)
        & do_intersect3D_vectorized_inlined(&obj1, &obj2);
}

auto runGJK_SoA(float** o1, float** o2, int n1, int n2) -> int {  // NOLINT
  CHObject_soa obj1{n1, o1};
  CHObject_soa obj2{n2, o2};

  return do_intersect3D_soa(&obj1, &obj2)
        & do_intersect3D_soa_v_slow_idx(&obj1, &obj2)
        & do_intersect3D_soa_v_slow_idx_inlined(&obj1, &obj2)
        & do_intersect3D_soa_vectorized(&obj1, &obj2)
        & do_intersect3D_soa_vectorized_lu2(&obj1, &obj2)
        & do_intersect3D_soa_vectorized_inlined(&obj1, &obj2);
}

TEST(TestBooleanGJK, Square1) {
  float o1[][3] = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {1.0, 1.0, 0.0}
  };
  float o2[][3] = {
      {-2.0, 0.0, 0.0},
      {-3.0, 0.0, 0.0},
      {-2.0, 1.0, 0.0},
      {-3.0, 1.0, 0.0}
  };

#ifdef MEASURE
  bench("Square1", o1, o2, 4, 4);
#endif
  auto res = runGJK(o1, o2, 4, 4);
  EXPECT_EQ(res, false);
}

TEST(TestBooleanGJK, Square1_SoA) {
  float** o1 = new float*[3];
  o1[0] = new float[4]{0.0, 1.0, 0.0, 1.0};
  o1[1] = new float[4]{0.0, 0.0, 1.0, 1.0};
  o1[2] = new float[4]{0.0, 0.0, 0.0, 0.0};

  float** o2 = new float*[3];
  o2[0] = new float[4]{-2.0, -3.0, -2.0, -3.0};
  o2[1] = new float[4]{0.0, 0.0, 1.0, 1.0};
  o2[2] = new float[4]{0.0, 0.0, 0.0, 0.0};

  auto res = runGJK_SoA(o1, o2, 4, 4);
  EXPECT_EQ(res, false);
}

TEST(TestBooleanGJK, ExactOverlap) {
  float o1[][3] = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {1.0, 1.0, 0.0}
  };
  float o2[][3] = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {1.0, 1.0, 0.0}
  };

#ifdef MEASURE
  bench("ExactOverlap", o1, o2, 4, 4);
#endif
  auto res = runGJK(o1, o2, 4, 4);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, ExactOverlap_SoA) {
  float** o1 = new float*[3];
  o1[0] = new float[4]{0.0, 1.0, 0.0, 1.0};
  o1[1] = new float[4]{0.0, 0.0, 1.0, 1.0};
  o1[2] = new float[4]{0.0, 0.0, 0.0, 0.0};

  float** o2 = new float*[3];
  o2[0] = new float[4]{0.0, 1.0, 0.0, 1.0};
  o2[1] = new float[4]{0.0, 0.0, 1.0, 1.0};
  o2[2] = new float[4]{0.0, 0.0, 0.0, 0.0};

  auto res = runGJK_SoA(o1, o2, 4, 4);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, LineOverlap) {
  float o1[][3] = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0}
  };
  float o2[][3] = {
      {0.5, 1.0, 0.0},
      {0.5, -1.0, 0.0}
  };

#ifdef MEASURE
  bench("LineOverlap", o1, o2, 2, 2);
#endif
  auto res = runGJK(o1, o2, 2, 2);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, LineOverlap_SoA) {
  float** o1 = new float*[3];
  o1[0] = new float[2]{0.0, 1.0};
  o1[1] = new float[2]{0.0, 0.0};
  o1[2] = new float[2]{0.0, 0.0};

  float** o2 = new float*[3];
  o2[0] = new float[2]{0.5, 0.5};
  o2[1] = new float[2]{1.0, -1.0};
  o2[2] = new float[2]{0.0, 0.0};

  auto res = runGJK_SoA(o1, o2, 2, 2);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, LineNonOverlap) {
  float o1[][3] = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0}
  };
  float o2[][3] = {
      {1.5, 1.0, 0.0},
      {1.5, -1.0, 0.0}
  };

#ifdef MEASURE
  bench("LineNonOverlap", o1, o2, 2, 2);
#endif
  auto res = runGJK(o1, o2, 2, 2);
  EXPECT_EQ(res, false);
}

TEST(TestBooleanGJK, LineNonOverlap_SoA) {
  float** o1 = new float*[3];
  o1[0] = new float[2]{0.0, 1.0};
  o1[1] = new float[2]{0.0, 0.0};
  o1[2] = new float[2]{0.0, 0.0};

  float** o2 = new float*[3];
  o2[0] = new float[2]{1.5, 1.5};
  o2[1] = new float[2]{1.0, -1.0};
  o2[2] = new float[2]{0.0, 0.0};

  auto res = runGJK_SoA(o1, o2, 2, 2);
  EXPECT_EQ(res, false);
}

TEST(TestBooleanGJK, SmolLinePointNonOverlap) {
  float o1[][3] = {
      {0.0, 0.0, 0.0},
      {0.01, 0.0, 0.0}
  };
  float o2[][3] = {
      {0.005, 0.0, 0.1}
  };

#ifdef MEASURE
  bench("SmolLinePointNonOverlap", o1, o2, 2, 1);
#endif
  auto res = runGJK(o1, o2, 2, 1);
  EXPECT_EQ(res, false);
}

TEST(TestBooleanGJK, SmolLinePointNonOverlap_SoA) {
  float** o1 = new float*[3];
  o1[0] = new float[2]{0.0, 0.01};
  o1[1] = new float[2]{0.0, 0.0};
  o1[2] = new float[2]{0.0, 0.0};

  float** o2 = new float*[3];
  o2[0] = new float[1]{0.005};
  o2[1] = new float[1]{0.0};
  o2[2] = new float[1]{0.1};

  auto res = runGJK_SoA(o1, o2, 2, 1);
  EXPECT_EQ(res, false);
}

TEST(TestBooleanGJK, SmolLinePointOverlap) {
  float o1[][3] = {
      {0.0, 0.0, 0.0},
      {0.01, 0.0, 0.0}
  };
  float o2[][3] = {
      {0.005, 0.0, 0.0}
  };

#ifdef MEASURE
  bench("SmolLinePointOverlap", o1, o2, 2, 1);
#endif
  auto res = runGJK(o1, o2, 2, 1);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, SmolLinePointOverlap_SoA) {
  float** o1 = new float*[3];
  o1[0] = new float[2]{0.0, 0.01};
  o1[1] = new float[2]{0.0, 0.0};
  o1[2] = new float[2]{0.0, 0.0};

  float** o2 = new float*[3];
  o2[0] = new float[1]{0.005};
  o2[1] = new float[1]{0.0};
  o2[2] = new float[1]{0.0};

  auto res = runGJK_SoA(o1, o2, 2, 1);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, PointOverlap) {
  float o1[][3] = {
      {0.5, 1.0, 0.0}
  };
  float o2[][3] = {
      {0.5, 1.0, 0.0}
  };

#ifdef MEASURE
  bench("PointOverlap", o1, o2, 1, 1);
#endif
  auto res = runGJK(o1, o2, 1, 1);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, PointOverlap_SoA) {
  float** o1 = new float*[3];
  o1[0] = new float[1]{0.5};
  o1[1] = new float[1]{1.0};
  o1[2] = new float[1]{0.0};

  float** o2 = new float*[3];
  o2[0] = new float[1]{0.5};
  o2[1] = new float[1]{1.0};
  o2[2] = new float[1]{0.0};

  auto res = runGJK_SoA(o1, o2, 1, 1);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, SquareSquareSideOverlap) {
  float o1[][3] = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {1.0, 1.0, 0.0}
  };
  float o2[][3] = {
      {1.0, 0.0, 0.0},
      {2.0, 0.0, 0.0},
      {1.0, 1.0, 0.0},
      {2.0, 1.0, 0.0}
  };

#ifdef MEASURE
  bench("SquareSquareSideOverlap", o1, o2, 4, 4);
#endif
  auto res = runGJK(o1, o2, 4, 4);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, SquareSquareSideOverlap_SoA) {
  float** o1 = new float*[3];
  o1[0] = new float[4]{0.0, 1.0, 0.0, 1.0};
  o1[1] = new float[4]{0.0, 0.0, 1.0, 1.0};
  o1[2] = new float[4]{0.0, 0.0, 0.0, 0.0};

  float** o2 = new float*[3];
  o2[0] = new float[4]{1.0, 2.0, 1.0, 2.0};
  o2[1] = new float[4]{0.0, 0.0, 1.0, 1.0};
  o2[2] = new float[4]{0.0, 0.0, 0.0, 0.0};

  auto res = runGJK_SoA(o1, o2, 4, 4);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, SquareSquarePointOverlap) {
  float o1[][3] = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {1.0, 1.0, 0.0}
  };
  float o2[][3] = {
      {1.0, 1.0, 0.0},
      {2.0, 1.0, 0.0},
      {1.0, 2.0, 0.0},
      {2.0, 2.0, 0.0}
  };

#ifdef MEASURE
  bench("SquareSquarePointOverlap", o1, o2, 4, 4);
#endif
  auto res = runGJK(o1, o2, 4, 4);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, SquareSquarePointOverlap_SoA) {
  float** o1 = new float*[3];
  o1[0] = new float[4]{0.0, 1.0, 0.0, 1.0};
  o1[1] = new float[4]{0.0, 0.0, 1.0, 1.0};
  o1[2] = new float[4]{0.0, 0.0, 0.0, 0.0};

  float** o2 = new float*[3];
  o2[0] = new float[4]{1.0, 2.0, 1.0, 2.0};
  o2[1] = new float[4]{1.0, 1.0, 2.0, 2.0};
  o2[2] = new float[4]{0.0, 0.0, 0.0, 0.0};

  auto res = runGJK_SoA(o1, o2, 4, 4);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, SquareSquareNonOverlap) {
  float o1[][3] = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {1.0, 1.0, 0.0}
  };
  float o2[][3] = {
      {1.0 + EPS, 0.0, 0.0},
      {2.0, 0.0, 0.0},
      {1.0 + EPS, 1.0, 0.0},
      {2.0, 1.0, 0.0}
  };

#ifdef MEASURE
  bench("SquareSquareNonOverlap", o1, o2, 4, 4);
#endif
  auto res = runGJK(o1, o2, 4, 4);
  EXPECT_EQ(res, false);
}

TEST(TestBooleanGJK, SquareSquareNonOverlap_SoA) {
  float** o1 = new float*[3];
  o1[0] = new float[4]{0.0, 1.0, 0.0, 1.0};
  o1[1] = new float[4]{0.0, 0.0, 1.0, 1.0};
  o1[2] = new float[4]{0.0, 0.0, 0.0, 0.0};

  float** o2 = new float*[3];
  o2[0] = new float[4]{1.0  + EPS, 2.0, 1.0 + EPS, 2.0};
  o2[1] = new float[4]{0.0, 0.0, 1.0, 1.0};
  o2[2] = new float[4]{0.0, 0.0, 0.0, 0.0};

  auto res = runGJK_SoA(o1, o2, 4, 4);
  EXPECT_EQ(res, false);
}

TEST(TestBooleanGJK, PolyPointOverlap) {
  float o1[][3] = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {1.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
      {1.0, 0.0, 1.0},
      {0.0, 1.0, 1.0},
      {1.0, 1.0, 1.0}
  };
  float o2[][3] = {
      {1.0, 1.0, 1.0},
      {2.0, 1.0, 1.0},
      {1.0, 2.0, 1.0},
      {2.0, 2.0, 1.0},
      {1.0, 1.0, 2.0},
      {2.0, 1.0, 2.0},
      {1.0, 2.0, 2.0},
      {2.0, 2.0, 2.0}
  };

#ifdef MEASURE
  bench("PolyPointOverlap", o1, o2, 8, 8);
#endif
  auto res = runGJK(o1, o2, 8, 8);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, PolyPointOverlap_SoA) {
  float** o1 = new float*[3];
  o1[0] = new float[8]{0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
  o1[1] = new float[8]{0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0};
  o1[2] = new float[8]{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};

  float** o2 = new float*[3];
  o2[0] = new float[8]{1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0};
  o2[1] = new float[8]{1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 2.0, 2.0};
  o2[2] = new float[8]{1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0};

  auto res = runGJK_SoA(o1, o2, 8, 8);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, PolyPointNonOverlap) {
  float o1[][3] = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {1.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
      {1.0, 0.0, 1.0},
      {0.0, 1.0, 1.0},
      {1.0, 1.0, 1.0}
  };
  float o2[][3] = {
      {1.0, 1.0, 1.0 + EPS},
      {2.0, 1.0, 1.0 + EPS},
      {1.0, 2.0, 1.0 + EPS},
      {2.0, 2.0, 1.0 + EPS},
      {1.0, 1.0, 2.0},
      {2.0, 1.0, 2.0},
      {1.0, 2.0, 2.0},
      {2.0, 2.0, 2.0}
  };

#ifdef MEASURE
  bench("PolyPointNonOverlap", o1, o2, 8, 8);
#endif
  auto res = runGJK(o1, o2, 8, 8);
  EXPECT_EQ(res, false);
}

TEST(TestBooleanGJK, PolyPointNonOverlap_SoA) {
  float** o1 = new float*[3];
  o1[0] = new float[8]{0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
  o1[1] = new float[8]{0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0};
  o1[2] = new float[8]{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};

  float** o2 = new float*[3];
  o2[0] = new float[8]{1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0};
  o2[1] = new float[8]{1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 2.0, 2.0};
  o2[2] = new float[8]{1.0 + EPS, 1.0 + EPS, 1.0 + EPS, 1.0 + EPS, 2.0, 2.0, 2.0, 2.0};

  auto res = runGJK_SoA(o1, o2, 8, 8);
  EXPECT_EQ(res, false);
}

TEST(TestBooleanGJK, PolyLineOverlap) {
  float o1[][3] = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {1.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
      {1.0, 0.0, 1.0},
      {0.0, 1.0, 1.0},
      {1.0, 1.0, 1.0}
  };
  float o2[][3] = {
      {1.0, 1.0, 0.0},
      {2.0, 1.0, 0.0},
      {1.0, 2.0, 0.0},
      {2.0, 2.0, 0.0},
      {1.0, 1.0, 1.0},
      {2.0, 1.0, 1.0},
      {1.0, 2.0, 1.0},
      {2.0, 2.0, 1.0}
  };

#ifdef MEASURE
  bench("PolyLineOverlap", o1, o2, 8, 8);
#endif
  auto res = runGJK(o1, o2, 8, 8);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, PolyLineOverlap_SoA) {
  float** o1 = new float*[3];
  o1[0] = new float[8]{0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
  o1[1] = new float[8]{0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0};
  o1[2] = new float[8]{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};

  float** o2 = new float*[3];
  o2[0] = new float[8]{1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0};
  o2[1] = new float[8]{1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 2.0, 2.0};
  o2[2] = new float[8]{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};

  auto res = runGJK_SoA(o1, o2, 8, 8);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, PolyProjectiveNonOverlap) {
  float o1[][3] = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {1.0, 1.0, 0.0},
      {1.0, 0.0, 1.0},
      {2.0, 0.0, 1.0},
      {1.0, 1.0, 1.0},
      {2.0, 1.0, 1.0}
  };
  float o2[][3] = {
      {1.1, 1.0, 0.0},
      {2.1, 1.0, 0.0},
      {1.1, 2.0, 0.0},
      {2.1, 2.0, 0.0},
      {2.1, 1.0, 1.0},
      {3.1, 1.0, 1.0},
      {2.1, 2.0, 1.0},
      {3.1, 2.0, 1.0}
  };

#ifdef MEASURE
  bench("PolyProjectiveNonOverlap", o1, o2, 8, 8);
#endif
  auto res = runGJK(o1, o2, 8, 8);
  EXPECT_EQ(res, false);
}

TEST(TestBooleanGJK, PolyProjectiveNonOverlap_SoA) { // Fails...
  float** o1 = new float*[3];
  o1[0] = new float[8]{0.0, 1.0, 0.0, 1.0, 1.0, 2.0, 1.0, 2.0};
  o1[1] = new float[8]{0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0};
  o1[2] = new float[8]{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};

  float** o2 = new float*[3];
  o2[0] = new float[8]{1.1, 2.1, 1.1, 2.1, 2.1, 3.1, 2.1, 3.1};
  o2[1] = new float[8]{1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 2.0, 2.0};
  o2[2] = new float[8]{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};

  auto res = runGJK_SoA(o1, o2, 8, 8);
  EXPECT_EQ(res, false);
}

TEST(TestBooleanGJK, PolyProjectiveOverlap) {
  float o1[][3] = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {1.0, 1.0, 0.0},
      {1.0, 0.0, 1.0},
      {2.0, 0.0, 1.0},
      {1.0, 1.0, 1.0},
      {2.0, 1.0, 1.0}
  };
  float o2[][3] = {
      {1.1, 1.0, 0.0},
      {2.1, 1.0, 0.0},
      {1.1, 2.0, 0.0},
      {2.1, 2.0, 0.0},
      {2.0, 1.0, 1.0},
      {3.1, 1.0, 1.0},
      {2.0, 2.0, 1.0},
      {3.1, 2.0, 1.0}
  };

#ifdef MEASURE
  bench("PolyProjectiveOverlap", o1, o2, 8, 8);
#endif
  auto res = runGJK(o1, o2, 8, 8);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, PolyProjectiveOverlap_SoA) {
  float** o1 = new float*[3];
  o1[0] = new float[8]{0.0, 1.0, 0.0, 1.0, 1.0, 2.0, 1.0, 2.0};
  o1[1] = new float[8]{0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0};
  o1[2] = new float[8]{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};

  float** o2 = new float*[3];
  o2[0] = new float[8]{1.1, 2.1, 1.1, 2.1, 2.0, 3.1, 2.0, 3.1};
  o2[1] = new float[8]{1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 2.0, 2.0};
  o2[2] = new float[8]{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};

  auto res = runGJK_SoA(o1, o2, 8, 8);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, ProceduralOverlap) {
  float o1[100][3];
  float o2[100][3];

  for (int i = 0; i < 100; ++i) {
    float radian = i * 2.0 * M_PI / 100.;
    o1[i][0] = cos(radian);
    o1[i][1] = sin(radian);
    o1[i][2] = 0.;

    o2[i][0] = cos(radian);
    o2[i][1] = sin(radian);
    o2[i][2] = 0.;
  }

#ifdef MEASURE
  bench("ProceduralOverlap", o1, o2, 100, 100);
#endif
  auto res = runGJK(o1, o2, 100, 100);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, ProceduralOverlap_SoA) {
  float** o1 = new float*[3];
  float** o2 = new float*[3];

  for (int i = 0; i < 3; i++) {
    o1[i] = new float[100];
    o2[i] = new float[100];
  }

  for (int i = 0; i < 100; ++i) {
    float radian = i * 2.0 * M_PI / 100.;
    o1[0][i] = cos(radian);
    o1[1][i] = sin(radian);
    o1[2][i] = 0.;

    o2[0][i] = cos(radian);
    o2[1][i] = sin(radian);
    o2[2][i] = 0.;
  }

  auto res = runGJK_SoA(o1, o2, 100, 100);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, ProceduralNonOverlap) {
  float o1[100][3];
  float o2[100][3];

  for (int i = 0; i < 100; ++i) {
    float radian = i * 2.0 * M_PI / 100.;
    o1[i][0] = cos(radian);
    o1[i][1] = sin(radian);
    o1[i][2] = 0.;

    o2[i][0] = cos(radian);
    o2[i][1] = sin(radian);
    o2[i][2] = EPS;
  }

#ifdef MEASURE
  bench("ProceduralNonOverlap", o1, o2, 100, 100);
#endif
  auto res = runGJK(o1, o2, 100, 100);
  EXPECT_EQ(res, false);
}

TEST(TestBooleanGJK, ProceduralNonOverlap_SoA) {
  float** o1 = new float*[3];
  float** o2 = new float*[3];

  for (int i = 0; i < 3; i++) {
    o1[i] = new float[100];
    o2[i] = new float[100];
  }

  for (int i = 0; i < 100; ++i) {
    float radian = i * 2.0 * M_PI / 100.;
    o1[0][i] = cos(radian);
    o1[1][i] = sin(radian);
    o1[2][i] = 0.;

    o2[0][i] = cos(radian);
    o2[1][i] = sin(radian);
    o2[2][i] = EPS;
  }

  auto res = runGJK_SoA(o1, o2, 100, 100);
  EXPECT_EQ(res, false);
}

TEST(TestBooleanGJK, ProceduralSectionOverlap) {
  float o1[100][3];
  float o2[100][3];

  for (int i = 0; i < 100; ++i) {
    float radian = i * 2.0 * M_PI / 100.;
    o1[i][0] = cos(radian);
    o1[i][1] = sin(radian);
    o1[i][2] = 0.;

    o2[i][0] = cos(radian) + 0.5;
    o2[i][1] = sin(radian);
    o2[i][2] = 0.;
  }

#ifdef MEASURE
  bench("ProceduralSectionOverlap", o1, o2, 100, 100);
#endif
  auto res = runGJK(o1, o2, 100, 100);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, ProceduralSectionOverlap_SoA) {
  float** o1 = new float*[3];
  float** o2 = new float*[3];

  for (int i = 0; i < 3; i++) {
    o1[i] = new float[100];
    o2[i] = new float[100];
  }

  for (int i = 0; i < 100; ++i) {
    float radian = i * 2.0 * M_PI / 100.;
    o1[0][i] = cos(radian);
    o1[1][i] = sin(radian);
    o1[2][i] = 0.;

    o2[0][i] = cos(radian) + 0.5;
    o2[1][i] = sin(radian);
    o2[2][i] = 0.;
  }


  auto res = runGJK_SoA(o1, o2, 100, 100);
  EXPECT_EQ(res, true);
}

TEST(TestBooleanGJK, ProceduralAwayNonOverlap) {
  float o1[100][3];
  float o2[100][3];

  for (int i = 0; i < 100; ++i) {
    float radian = i * 2.0 * M_PI / 100.;
    o1[i][0] = cos(radian);
    o1[i][1] = sin(radian);
    o1[i][2] = 0.;

    o2[i][0] = cos(radian) + 2. + 2. * EPS;
    o2[i][1] = sin(radian);
    o2[i][2] = 0.;
  }

#ifdef MEASURE
      bench("ProceduralAwayNonOverlap", o1, o2, 100, 100);
#endif
  auto res = runGJK(o1, o2, 100, 100);
  EXPECT_EQ(res, false);
}

TEST(TestBooleanGJK, ProceduralAwayNonOverlap_SoA) {
  float** o1 = new float*[3];
  float** o2 = new float*[3];

  for (int i = 0; i < 3; i++) {
    o1[i] = new float[100];
    o2[i] = new float[100];
  }

  for (int i = 0; i < 100; ++i) {
    float radian = i * 2.0 * M_PI / 100.;
    o1[0][i] = cos(radian);
    o1[1][i] = sin(radian);
    o1[2][i] = 0.;

    o2[0][i] = cos(radian) + 2. + 2. * EPS;
    o2[1][i] = sin(radian);
    o2[2][i] = 0.;
  }

  auto res = runGJK_SoA(o1, o2, 100, 100);
  EXPECT_EQ(res, false);
}

// TODO: test with more complex shapes
// TODO: test procedurally generated shapes
// TODO: make sure all cases in do_simplex3D_vectorized are covered by the tests.
// From here on out, the tests are not based off other work.

} // namespace test
