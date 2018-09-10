# oscode: Solve oscillatory differential equations with RKWKB 

## Description

Restructured version of cRKWKB.

## Requirements

- Testing makes use of various routines of the NAG library, namely the
  self-contained library with shared linkage `libnagc_nag.so`. Therefore you
  will need to:
    - download and install the NAG C library from [here](https://www.nag.co.uk/content/downloads-nag-c-library-versions)
    - modify `INCDIR` and `SHARDIR` in the Makefile according to your
      installation directory.

## Files and dependencies

- `system.hpp`: defines an object that stores information about the system of
  equations being solved such as F, w, g, gradients, Hessians, etc.
- `integrators.hpp`: defines an object each for every method of numerical
  integration available, e.g. RKF, WKB. Each of these has a step method that
  takes one step with the given method. The results of the step are stored in a
  Step data structure, defined here.
- `solver.hpp`: defines a solution object.
- `typedefs.hpp`: for Eigen typedefs.
- All files live in the `namespace RKWKB` (apart from `test.cpp`) to avoid
  typedef clashes.

## Testing

- is done with [catch2](https://github.com/catchorg/Catch2)
- `test-airy.cpp` is used for testing an RKWKB solver on the Airy equation.
    - `make test-airy`, then `./test-airy`.

## TODO

- Write the dS, S, etc. functions for WKB orders other than 1
- Truncation error estimate in WKB steps
