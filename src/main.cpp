#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>
#include "util/OBJ_Loader.h"
#include "benchmarking/add_functions.cpp"
#include "benchmarking/common.h"
#include "benchmarking/performance_tests.cpp"

#define SHIFT 2
#define MEASURE_SIZE_IMPACT

/**
 * Global vars, used to keep track of gjk functions
 * */
std::vector<gjk_func> userFuncs;
std::vector<std::string> funcNames;
int numFuncs = 0;

std::vector<gjk_func_soa> userFuncs_soa;
std::vector<std::string> funcNames_soa;
int numFuncs_soa = 0;

/**
 * Initialize the points
 * */
void kernel_base(float (* pts1)[3], float (* pts2)[3], int N, intersection_type type) {
  float delta = (type == PROCEDURAL_OVERLAP) ? 0.5 : SHIFT;
  unsigned seed = 4010595406; // FIX SEED, ELIMINATE RANDOMNESS BETWEEN DIFFERENT RUNS std::chrono::system_clock::now().time_since_epoch().count();

  std::mt19937 generator(seed);
  std::uniform_real_distribution<float> uniform01(0.0, 1.0);

  for (int i = 0; i < N; ++i) {
      float theta = 2 * M_PI * uniform01(generator);
      float phi = acos(1 - 2 * uniform01(generator));

      pts1[i][0] = sin(phi) * cos(theta);
      pts1[i][1] = sin(phi) * sin(theta);
      pts1[i][2] = cos(phi);

      pts2[i][0] = sin(phi) * cos(theta);
      pts2[i][1] = sin(phi) * sin(theta);
      pts2[i][2] = delta + cos(phi);
  }
}

/**
 * Initialize the points
 * */
void kernel_base_soa(float** pts1, float** pts2, int N, intersection_type type) {
  float delta = (type == PROCEDURAL_OVERLAP) ? 0.5 : SHIFT;
  unsigned seed = 4010595406; // FIX SEED, ELIMINATE RANDOMNESS BETWEEN DIFFERENT RUNS std::chrono::system_clock::now().time_since_epoch().count();

  std::mt19937 generator(seed);
  std::uniform_real_distribution<float> uniform01(0.0, 1.0);

  for (int i = 0; i < N; ++i) {
    float theta = 2 * M_PI * uniform01(generator);
    float phi = acos(1 - 2 * uniform01(generator));

    pts1[0][i] = sin(phi) * cos(theta);
    pts1[1][i] = sin(phi) * sin(theta);
    pts1[2][i] = cos(phi);

    pts2[0][i] = sin(phi) * cos(theta);
    pts2[1][i] = sin(phi) * sin(theta);
    pts2[2][i] = delta + cos(phi);
  }
}

/**
 * Registers a user function to be tested by the driver program. Registers a
 * string description of the function as well
 * */
void add_gjk_function(gjk_func f, const std::string& name, int type) {
  userFuncs.push_back(f);
  funcNames.emplace_back(name);

  numFuncs++;
}

void add_gjk_soa_function(gjk_func_soa f, const std::string& name, int type) {
  userFuncs_soa.push_back(f);
  funcNames_soa.emplace_back(name);

  numFuncs_soa++;
}

int main(int argc, char* argv[]) {
  std::cout << "Starting program. ";

  myInt64 perf;
  int i;

  register_functions();

  if (numFuncs == 0 && numFuncs_soa == 0){
      std::cout << std::endl;
      std::cout << "No functions registered - nothing for driver to do" << std::endl;
      std::cout << "Register functions by calling register_func(f, name)" << std::endl;
      std::cout << "in register_funcs()" << std::endl;

      return 0;
  }
  std::cout << numFuncs << " function(s) registered." << std::endl;

  // This section measures how the runtime is affected with increasing input size.
#ifdef MEASURE_SIZE_IMPACT
  int n_var;
  for (int j = 0; j < numFuncs; j++) {
    std::cout << "################ Performance Tests for " << funcNames[j] << " ###################" << std::endl;
    for (i = 3; i <= 26; i++) { // i = 19 reaches program stack's maximum size. (probably...)
      n_var = 1 << i;
      //float pts1_var[n_var][3];
      auto pts1_var = new float[n_var][3];
      auto pts2_var = new float[n_var][3];

      kernel_base(pts1_var, pts2_var, n_var, PROCEDURAL_OVERLAP);

      perf = rdtsc_gjk(userFuncs[j], pts1_var, pts2_var, n_var, n_var);
      std::cout << std::endl << "Running: " << funcNames[j] << " for input size n = 2^" << i << std::endl;
      std::cout << perf << " cycles" << std::endl;

      // Should I free arrays of size 3 separately?
      delete [] pts1_var;
      delete [] pts2_var;
    }
  }

  for (int j = 0; j < numFuncs_soa; j++) {
    for (int i = 3; i <= 26; i++) {
      int n_base = 1 << i;
      float** pts1_soa = new float*[3];
      float** pts2_soa = new float*[3];

      pts1_soa[0] = new float[n_base];
      pts1_soa[1] = new float[n_base];
      pts1_soa[2] = new float[n_base];

      pts2_soa[0] = new float[n_base];
      pts2_soa[1] = new float[n_base];
      pts2_soa[2] = new float[n_base];

      kernel_base_soa(pts1_soa, pts2_soa, n_base, PROCEDURAL_OVERLAP);

      perf = rdtsc_gjk_soa(userFuncs_soa[j], pts1_soa, pts2_soa, n_base, n_base);
      std::cout << std::endl << "Running: " << funcNames_soa[j] << " for input size n = 2^" << i << std::endl;
      std::cout << perf << " cycles" << std::endl;

      for (int k = 0; k < 3; k++) {
        delete [] pts1_soa[k];
        delete [] pts2_soa[k];
      }
      delete [] pts1_soa;
      delete [] pts2_soa;
    }
  }
#endif

  // Declare the the objects that will be used in the algorithm.
  const int n_base = (1 << 10);
  auto pts1 __attribute__ ((aligned (32))) = new float[n_base][3];
  auto pts2 __attribute__ ((aligned (32))) = new float[n_base][3];

  // Initialized the objects.
  kernel_base(pts1, pts2, n_base, PROCEDURAL_OVERLAP);

  // Declare object in SoA format.
  float** pts1_soa = new float*[3];
  float** pts2_soa = new float*[3];

  pts1_soa[0] = new float[n_base];
  pts1_soa[1] = new float[n_base];
  pts1_soa[2] = new float[n_base];

  pts2_soa[0] = new float[n_base];
  pts2_soa[1] = new float[n_base];
  pts2_soa[2] = new float[n_base];

  kernel_base_soa(pts1_soa, pts2_soa, n_base, PROCEDURAL_OVERLAP);

  for (i = 0; i < numFuncs; i++) {
    perf = rdtsc_gjk(userFuncs[i], pts1, pts2, n_base, n_base);
    std::cout << "Running: " << funcNames[i] << std::endl;
    std::cout << perf << " cycles" << std::endl;
  }

  for (i = 0; i < numFuncs_soa; i++) {
    perf = rdtsc_gjk_soa(userFuncs_soa[i], pts1_soa, pts2_soa, n_base, n_base);
    std::cout << "Running: " << funcNames_soa[i] << std::endl;
    std::cout << perf << " cycles" << std::endl;
  }

  // Delete the arrays
  delete [] pts1;
  delete [] pts2;
  for (int k = 0; k < 3; k++) {
    delete [] pts1_soa[k];
    delete [] pts2_soa[k];
  }
  delete [] pts1_soa;
  delete [] pts2_soa;

  return 0;
}