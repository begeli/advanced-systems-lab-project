#ifndef GJK_COMMON_H
#define GJK_COMMON_H

#include <string>
#include "../src/include/gjk.h"

/**
 * Might have to change the name of the file to something more appropriate.
 * The file is supposed to contain the type definitions and function signatures
 * that are used to facilitate benchmarking of various gjk versions.
 * /



/**
 * gjk_func is the general function pointer for all
 * of the optimized versions of our gjk algorithm.
 * */
typedef int(*gjk_func)(const CHObject*, const CHObject*);
typedef int(*gjk_func_soa)(const CHObject_soa*, const CHObject_soa*);

typedef enum {
    PROCEDURAL_OVERLAP,
    PROCEDURAL_NON_OVERLAP
} intersection_type;

/**
 * Adds gjk version to a list of function pointers.
 * */
void add_gjk_function(gjk_func f, const std::string& name, int);

/**
 * Adds gjk version to a list of function pointers.
 * */
void add_gjk_soa_function(gjk_func_soa f, const std::string& name, int);

/**
 * Initializes arguments to the methods that will be tested.
 * */
void kernel_base(float (*)[3], float (*)[3], int N, intersection_type);

/**
 * Initializes arguments to the methods that will be tested.
 * Object are initialized in SoA format.
 * */
void kernel_base_soa(float**, float**, int N, intersection_type);

/**
 * Allows multiple functions to be added. Might not be necessary
 * depending on how and where the optimized versions of gjk are implemented.
 * */
void register_functions();

#endif //GJK_COMMON_H
