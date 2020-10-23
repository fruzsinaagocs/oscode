import sys
import os
import _pyoscode
import numpy

def solve(ts, ws, gs, ti, tf, x0, dx0, t_eval=[], logw=False, logg=False, order=3,
rtol=1e-4, atol=0.0, h=None, full_output="", even_grid=False, check_grid=False):
    """Solve a differential equation with the RKWKB method.
    
    Parameters
    ----------
    ts: numpy.ndarray [float] or list [float]
       An array of real numbers representing the values of the independent
       variable at which the frequency and friction term are evaluated. 

    ws: numpy.ndarray [complex] or list [complex]
       An array-like object of real or complex 
       numbers, representing the values of frequency w at the points given in ts.

    gs: numpy.ndarray [complex] or list [complex]
        An array-like object of real or complex numbers representing the values
        of the friction term g at the points given in ts.

    ti,tf: float
        Start and end of integration range.

    x0, dx0: complex
        Initial values of the dependent variable and its derivative.

    t_eval: numpy.ndarray [float] or list [float]
        An array of times where the solution is to be returned.

    logw, logg: boolean, optional
        If true, the array of frequencies and friction values, respectively, will be
        exponentiated (False, False by default).

    order: int, optional
        Order of WKB approximation to use, 3 (the highest value) by default.
    
    rtol, atol: float, optional
        Relative and absolute tolerance of the solver, 1e-4 and 0 by default.
        Note that atol at the moment is not implemented. 

    h: float, optional
        Size of the initial step, 1 by default.

    full_output: str , optional
        If given, the return dictionary will be written to a file with the supplied
        name.

    even_grid: boolean, optional
        False by default. Set this to True if the ts array is evenly spaced for
        faster interpolation.

    check_grid: boolean, optional
        False by default. If True, the fineness of the ws, gs grids will be
        checked based on how accurate linear interpolation would be on them, and
        a warning will be issued if this accuracy is deemed too low. It's a good
        idea to set this to True when solving an equation for the first time.
    
    Returns
    -------
    A dictionary with the following keywords and values:

        sol: list [complex]
            A list containing the solution evaluated at timepoints listed under
            the 't' keyword.

        dsol: list [complex]
            A list containint the first deriv ative of the solution evaluated at
            timepoints listed under the 't' keyword.

        t: list [float]
            Contains the values of the independent variable where the solver
            stepped, i.e. evaluated the solution at. This is not determined by
            the user, rather these are the internal steps the solver naturally
            takes when integrating.

        types: list [float]
            A list of True/False values corresponding to the step types the
            solver chose at the timepoints listed under the keyword 't'. If
            True, the step was WKB, and RK otherwise.

        x_eval: list [complex]
            Values of the solution at the points specified in t_eval.

        dx_eval: list [complex]
            Values of the derivative of the solution at the points specified in
            t_eval.
     
    """
    # Set direction of integration if initial stepsize, h, not given
    if h==None:
        h = numpy.sign(tf - ti)
        # Handle the case of ti = tf
        if h==0:
            h=1
    
    # Run oscode from module library
    resdict = _pyoscode.solve(ts, ws, gs, ti, tf, x0, dx0, t_eval=t_eval, logw=logw, logg=logg,
    order=order, rtol=rtol, atol=atol, h=h, full_output=full_output,
    even_grid=even_grid, check_grid=check_grid) 
    
    return resdict
