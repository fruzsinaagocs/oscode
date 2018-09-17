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
    - extend/set the environment variable `LD_LIBRARY_PATH`
      according to [these](https://www.nag.co.uk/doc/inun/cl26/l6i1dl/un.html#example) instructions.


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
- `test/test-airy.cpp` is used for testing an RKWKB solver on the Airy equation.
    - `make test/test-airy`, then `./test/test-airy` with the optional argument `nag` or `rkwkb`. The latter will only run one of the solvers as specified.
    - This will append to `outputs/airy-<nag/rkwkb>.txt`, creating it if necessary.
    - `plot.py`'s functions can then be used to plot the results.
- `test/test-ms.cpp` test the solution of the Mukhanov--Sasaki equation.
    - `make test/test-ms`, then `./test/test-ms` with same optional argument as before.
    - This will append to `outputs/ms-<nag/rkwkb>.txt`.
    - `plot.py` contains useful functions to plot the output with.

## TODO

- Write the dS, S, etc. functions for WKB orders other than 1
- Warning/error message and termination of integration when troublesome point encountered
- Various error messages and not quiet failure
- Improve root finding to detect end of integration.
