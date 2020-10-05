.. title:: oscode (C++ interface)

============================
Using oscode's C++ interface
============================

.. sectnum:: 

about 
../api/library_root

Overview
--------

This documentation illustrates how one can use ``oscode`` via its C++ interface.
Usage of ``oscode`` involves

- defining an equation to solve,
- solving the equation,
- and extracting the solution and other statistics about the run.

The next sections will cover each of these. For a quick
reference of function arguments, see the `Quick Reference`_. 

Defining an equation
--------------------

The equations ``oscode`` can be used to solve are of the form 

.. math::

   \ddot{x}(t) + 2\gamma(t)\dot{x}(t) + \omega^2(t)x(t) = 0,

where :math:`x(t)`, :math:`\gamma(t)`, :math:`\omega(t)` can be complex. We will
call :math:`t` the independent variable, :math:`x` the dependent variable,
:math:`\omega(t)` the frequency term, and :math:`\gamma(t)` the friction or
first-derivative term. 

Defining an equation consists of the following:

- giving the frequency :math:`\omega(t)`,
- giving the first-derivative term :math:`\gamma(t)`,
- (implicitly) defining the range of integration, from :math:`t_i` to :math:`t_f`.

Defining the frequency and the first-derivative term can either be done by
giving them as **functions explicitly**, or by giving them as **time-series
values** on a grid of values of :math:`t`.

:math:`\omega` and :math:`\gamma` as explicit functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If :math:`\omega` and :math:`\gamma` are closed-form functions of time, then
define them as

.. code:: c
    
    #include "solver.hpp" // de_system, Solution defined in here 

    std::complex<double> g(double t){
        return 0.0;
    };
    
    std::complex<double> w(double t){
        return std::pow(9999,0.5)/(1.0 + t*t);
    };

Then feed them to the solver via the de_system class:

.. code:: c
    
    de_system sys(&w, &g);   
    Solution solution(sys, ...) // other arguments left out

:math:`\omega` and :math:`\gamma` as time series
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes :math:`\omega` and :math:`\gamma` will be results of numerical
integration, and they will have no closed-form functional form. In this case,
they can be specified on a grid, and ``oscode`` will perform linear
interpolation on the given grid to find their values at any timepoint. Because
of this, some important things to **note** are:

- For the sake of speed, ``oscode`` will assume the grid of timepoints :math:`\omega` and :math:`\gamma` are given on is **evenly spaced**. It is possible not to supply an evenly spaced grid, this will be discussed later.
- The timepoints grid needs to be **monotonically increasing/decreasing**.
- The timepoints grid needs to **include the range of integration** (:math:`t_i` and :math:`t_f`). 
- The grids for the timepoints, frequencies, and first-derivative terms have to be the **same size**.
- The speed/efficiency of the solver depends on how accurately it can carry out numerical integrals of the frequency and the first-derivative terms, therefore the **grid fineness** needs to be high enough. (Typically this means that linear interpolation gives a :math:`\omega(t)` value that is accurate to 1 part in :math:`10^{-9}` or so.)

To define the grids, use any array-like container (arrays, std::arrays,
std::vectors, Eigen::vectors are all accepted): 

.. code:: c
    
    #include "solver.hpp" // de_system, Solution defined in here 

    // Create a fine grid of timepoints and 
    // a grid of values for w, g
    N = 10000; 
    std::vector<double> ts(N);
    std::vector<std::complex<double>> ws(N), gs(N);
    
    // Fill up the grids
    for(int i=0; i<N; i++){
        ts[i] = i;
        ws[i] = std::sqrt(i);
        gs[i] = 0.0;
    }   

They can then be given to the solver again via the de_system class:

.. code:: c
    
    de_system sys(ts, ws, gs);   
    Solution solution(sys, ...) // other arguments left out


Often :math:`\omega` and :math:`\gamma` are much easier to perform linear
interpolation on once taken natural log of. This is what the optional ``islogw``
and ``islogg`` arguments of the overloaded ``de_system::de_system()``
constructor are for:

.. code:: c
    
    #include "solver.hpp" // de_system, Solution defined in here 

    // Create a fine grid of timepoints and 
    // a grid of values for w, g
    N = 10000; 
    std::vector<double> ts(N);
    std::vector<std::complex<double> logws(N), gs(N); // Note the log!
    
    // Fill up the grids
    for(int i=0; i<N; i++){
        ts[i] = i;
        logws[i] = 0.5*i;
        gs[i] = 0.0; // Will not be logged
    }   
    
    // We want to tell de_system that w has been taken natural log of, but g
    // hasn't. Therefore islogw=true, islogg=false:
    de_system sys(ts, logws, gs, true, false);
    Solution solution(sys, ... ) // other arguments left out


DIY interpolation
=================

For some problems, linear interpolation of :math:`\omega` and :math:`\gamma` (or
their natural logs) on an evenly spaced grid might simply not be enough, or the
user may want to carry out linear interpolation instead of letting ``oscode`` do
it for the sake of speed.

For example the user could carry out linear interpolation on an unevenly spaced
grid and feed :math:`\omega` and :math:`\gamma` as functions to ``de_system`` as
given below. Quadratic or other interpolation schemes can also be coded and used
like this.

.. code:: c

    std::complex<double> g(double t){
        int i;
        // Find index of element in ts closest to t from above
        i = std::distance(t.begin(), std::lower_bound(ts.begin(), ts.end(), t));
        std::complex<double> g0 = g[i-1];
        std::complex<double> g1 = g[i];
        return (g0+(g1-g0)*(t-ts[i-1])/(ts[i]-ts[i-1]));
    };


An example for wanting to do linear interpolation outside of ``oscode`` is
when ``Solution.solve()`` is ran in a loop, and for each iteration a large grid
of :math:`\omega` and :math:`\gamma` is required, depending on some parameter.
Instead of generating them over and over again, one could define them as
functions, making use of some underlying vectors that are independent of the
parameter we iterate over:

.. code:: c

    // A, B, and C are large std::vectors, same for each run
    // k is a parameter, different for each run
    // the grid of timepoints w, g are defined on starts at tstart, and is
    // evenly spaced with a spacing tinc.

    // tstart, tinc, A, B, C defined here

    std::complex<double> g(double t){
        int i;
        i=int((t-tstart)/tinc);
        std::complex<double> g0 = 0.5*(k*k*A[i] + 3.0 - B[i] + C[i]*k;
        std::complex<double> g1 = 0.5*(k*k*A[i+1] + 3.0 - B[i+1] + C[i+1]*k);
        return (g0+(g1-g0)*(t-tstart-tinc*i)/tinc);
    };



Solving an equation
-------------------

Once the equation to be solver has been defined as an instance of the
``de_system`` class, the following additional information is necessary to solve
it: 

- initial conditions, :math:`x(t_i)` and :math:`\dot{x}(t_f)`,
- the range of integration, from :math:`t_i` and :math:`t_f`,
- (optional) order of WKB approximation to use, ``order=3``,
- (optional) relative tolerance, ``rtol=1e-4``,
- (optional) absolute tolerance ``atol=0.0``,
- (optional) initial step ``h_0=1``,
- (optional) output file name ``full_output=""``,

**Note** the following about the optional arguments:

- ``rtol``, ``atol`` are tolerances on the local error. The global error in the solution is not guaranteed to stay below these values, but the error per step is. In the RK regime (not oscillatory solution), the global error will rise above the tolerance limits, but in the WKB regime, the global error usually stagnates.
- The initial step should be thought of as an initial estimate of what the first stepsize should be. The solver will determine the largest possible step within the given tolerance limit, and change ``h_0`` if necessary.
- The full output of ``solve()`` will be written to the filename contained in ``full_output``, if specified.  

Here's an example to illustrate usage of all of the above variables:

.. code:: c
    
    #include "solver.hpp" // de_system, Solution defined in here 

    // Define the system
    de_system sys(...) // For args see previous examples

    // Necessary parameters:
    // initial conditions
    std::complex<double> x0=std::complex<double>(1.0,1.0), dx0=0.0;
    // range of integration
    double ti=1.0, tf=100.0;
    
    // Optional parameters:
    // order of WKB approximation to use
    int order=2;
    // tolerances
    double rtol=2e-4, atol=0.0;
    // initial step
    double h0 = 0.5;
    // write the solution to a file
    std::string outfile="output.txt";

    Solution solution(sys, x0, dx0, ti, tf, order, rtol, atol, h0, outfile);
    // Solve the equation:
    solution.solve()

Here, we've also called the ``solve()`` method of the ``Solution`` class, to
carry out the integration. Now all information about the solution is in
``solution`` (and written to ``output.txt``).

Using the solution
------------------

Let's break down what ``solution`` contains (what ``Solution.solve()`` returns).
An instance of a ``Solution`` object is returned with the following attributes:

- ``times`` [std::list of double]: timepoints at which the solution was determined. These are **not** supplied by the user, rather they are internal steps that the solver has takes. The list starts with :math:`t_i` and ends with :math:`t_f`, these points are always guaranteed to be included.
- ``sol`` [std::list of std::complex<double>]: the solution at the timepoints specified in ``times``.
- ``dsol`` [std::list of std::complex<double>]: first derivative of the solution at timepoints specified in ``times``. 
- ``wkbs`` [std::list of int/bool]: types of steps takes at each timepoint in ``times``. **1** if the step was WKB, **0** if it was RK.  
- ``ssteps`` [int]: total number of accepted steps.  
- ``totsteps`` [int]: total number of attempted steps (accepted + rejected).  
- ``wkbsteps`` [int]: total number of successful WKB steps. 


Quick Reference
---------------

To construct a system, use the overloaded ``de_system`` constructor:

.. code:: c

    // For use with w, g as arrays
    template<typename X, typename Y, typename Z> de_system(const X &ts, const Y &ws, const Z &gs, bool isglogw=false, bool islogg=false);
    
    // For use with w, g as functions
    de_system(std::complex<double> (*w)(double), std::complex<double> (*g)(double));

To solve an equation, first build a ``Solution`` object with the constructor

.. code:: c

    Solution(de_system &de_sys, std::complex<double> x0, std::complex<double>
    dx0, double t_i, double t_f, int o=3, double r_tol=1e-4, double a_tol=0.0,
    double h_0=1, const char* full_output="");
   
And then to solve, simply call ``Solution``'s ``solve`` method
    
.. code:: c
    
    void solve();





