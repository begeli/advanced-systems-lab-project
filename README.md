# Collision Detection: The Gilbert-Johnson-Keerthi Algorithm and Bounding Volume Hierarchies 

## Project Setup

Currently only importing into CLion is confirmed to work.
Running CMake directly (passing the desired C and C++ compiler flags explicitly) should probably also work.

## For CLion
Open the outer folder in CLion as a project. Setup any Profile under `Build, Execution, Deployment` > `CMake`, passing the desired Compiler Flags (such as target architecture) as CMake options.

For now the C compiler flags are: `-std=c11 -O3 -march=[MACHINE ARCHITECTURE] -mprefer-vector-width=512 -ffast-math`.

`[MACHINE ARCHITECTURE]` can be:
* `skylake-avx512` (for the rented `Intel Xeon D-2141I CPU` bare metal server)
* `icelake` (for Berke's `Intel Core i7-1068NG7` Laptop)

The C++ flags are analogous.

## Testing
### Using GoogleTest 
The `test` folder contains all tests.
Tests are compiled in C++ using the `GoogleTest` library. They test the GJK and BVH code, which was written and compiled in C.

The tests can also be run from CLion conveniently, like unit tests in Java with JUnit.
(Note: you will have to build the `gjk` target, before building and running the `tester` target).
