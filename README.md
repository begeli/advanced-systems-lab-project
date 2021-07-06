# Collision Detection: The Gilbert-Johnson-Keerthi Algorithm and Bounding Volume Hierarchies 

## Files and Folders
* `14_report.pdf`: This file is the final version of our report. It contains brief explanations of the algorithms of interest, detailed explanations of the optimization methods we applied and our experimental results.
* `src`: This folder contains full implementations of the algorithms and optimized versions.
* `test`: This folder contains our unit test cases for the algorithms. 
* `benchmarking`: This folder contains the code used to measure execution times of the various implementations.

## Project Setup

Currently only importing into CLion is confirmed to work.
Running CMake directly (passing the desired C and C++ compiler flags explicitly) should probably also work.

## For CLion
Open the outer folder in CLion as a project. Setup any Profile under `Build, Execution, Deployment` > `CMake`, passing the desired Compiler Flags (such as target architecture) as CMake options.

For now the C compiler flags are: `-std=c11 -O3 -march=[MACHINE ARCHITECTURE] -mprefer-vector-width=512 -ffast-math`.

`[MACHINE ARCHITECTURE]` can be:
* `skylake-avx512` (for the rented `Intel Xeon D-2141I CPU` bare metal server) for BVH
* `icelake` (`Intel Core i7-1068NG7`) for GJK

The C++ flags are analogous.

## Testing
### Using GoogleTest 
The `test` folder contains all tests.
Tests are compiled in C++ using the `GoogleTest` library. They test the GJK and BVH code, which was written and compiled in C.

The tests can also be run from CLion conveniently, like unit tests in Java with JUnit.
(Note: you will have to build the `gjk` target, before building and running the `tester` target).

## Contributors

* Mihai Zorca
* Berke Egeli
* Chris Müller
* Liam van der Poel
