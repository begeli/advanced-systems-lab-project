#include "common.h"
#include "tsc_x86.h"

//#define CALIBRATE
#define CYCLES_REQUIRED 100000000

/**
 * This file holds various implementations (currently only straightforward gjk)
 * of the performance tests we might need.
 * */

int result = -1;  // global variable, so repeated runs don't get optimized.

/**
 * Measures the runtime of gjk function
 * */
myInt64 rdtsc_gjk(gjk_func func, const float (* o1)[3], const float (* o2)[3], int n1, int n2) {
    int i, num_runs;
    myInt64 cycles;
    myInt64 start;
    num_runs = 1;

    CHObject obj1{n1, o1};
    CHObject obj2{n2, o2};

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
      result = func(&obj1, &obj2);
    }
    cycles = stop_tsc(start);

    if (cycles >= CYCLES_REQUIRED) break;

    num_runs *= 2;
  }
#endif

    start = start_tsc();
    for (i = 0; i < num_runs; ++i) {
        result = func(&obj1, &obj2);
    }

    cycles = stop_tsc(start) / num_runs;
    return cycles;
}

myInt64 rdtsc_gjk_soa(gjk_func_soa func, float** o1, float** o2, int n1, int n2) {
  int i, num_runs;
  myInt64 cycles;
  myInt64 start;
  num_runs = 1;

  CHObject_soa obj1{n1, o1};
  CHObject_soa obj2{n2, o2};
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
      result = func(&obj1, &obj2);
    }
    cycles = stop_tsc(start);

    if (cycles >= CYCLES_REQUIRED) break;

    num_runs *= 2;
  }
#endif

  start = start_tsc();
  for (i = 0; i < num_runs; ++i) {
    result = func(&obj1, &obj2);
  }

  cycles = stop_tsc(start) / num_runs;
  return cycles;
}

