#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <string>
#include <functional>
#include <Eigen/Dense>
#include <iomanip>
#include <fstream>
#include <complex>
#include <iostream>
#include <numpy/arrayobject.h>
#include <cmath>
#include <array>
#include "system.hpp"
#include "solver.hpp"

// Docstrings for the module and its members
static char module_docstring[] = 
"This module provides an interface for solving oscillatory ordinary differential equations with the RKWKB method in C++.";

static char solve_docstring[] =
"Solve a differential equation of the form x(t)'' + 2*g(t)*x(t)' + w(t)^2*x(t) = 0 with the RKWKB method.\n\nParameters\n----------\nts : `numpy.ndarray` [`float`] or `list` [`float`]\n\tAn array of real numbers representing the values of the independent variable at which the frequency and friction term are evaluated. \nws : `numpy.ndarray` [`complex`] or `list` [`complex`]\n\tAn array-like object of real or complex numbers, representing the values of frequency w at the points given in ts.\ngs : `numpt.ndarray` [`complex`] or `list` [`complex`]\n\tAn array-like object of real or complex numbers representing the values of the friction term g at the points given in ts.\nti, tf : `float`\n\tStart and end of integration range.\nx0,dx0 : `complex`\n\tInitial values of the dependent variable and its derivative.\nlogw,logg : `boolean`, optional\n\tIf true, the array of frequencies and friction values, respectively, will be exponentiated (False, False by default).\norder : `int`, optional\n\tOrder of WKB approximation to use, 3 (the highest value) by default.\nrtol,atol : `float`, optional\n\tRelative and absolute tolerance of the solver, 1e-4 and 0 by default.\nh : `float, optional\n\tSize of the initial step, 1 by default.\nfull_output : `boolean`, optional\n\tIf true, the full solution will be written to a file.\n\nReturns\n----------\n";

// Declaration of the solver function
static PyObject *_solver_solve(PyObject *self, PyObject *args, PyObject *kwargs);

// Defining the members of this module
static PyMethodDef module_methods[] = {
    {"solve", (PyCFunction) _solver_solve, METH_VARARGS | METH_KEYWORDS, solve_docstring},
    {NULL, NULL, 0, NULL}
};

// Init function for the module
PyMODINIT_FUNC init_solver(void){
    
    PyObject *m = Py_InitModule3("_solver", module_methods, module_docstring);
    if(m==NULL)
        return;
    // Load numpy functionality
    import_array();
}

// Define constructor function
static PyObject *_solver_solve(PyObject *self, PyObject *args, PyObject *kwargs){

    int islogw=0,islogg=0,full_output=0,order=3,interp=1;
    double ti,tf,rtol=1e-4,atol=0.0,h0=1.0;
    std::complex<double> x0,dx0;
    PyObject *tsobj, *wsobj, *gsobj;
    // Define keywords
    static const char *kwlist[] =
    {"ts","ws","gs","ti","tf","x0","dx0","logw","logg","order","rtol","atol","h","full_output",NULL};

    // Interpret input arguments.
    if (!PyArg_ParseTupleAndKeywords(args,kwargs,"OOOddDD|iiidddi",const_cast<char**>(kwlist),&tsobj,&wsobj,&gsobj,&ti,&tf,&x0,&dx0,&islogw,&islogg,&order,&rtol,&atol,&h0,&full_output))
        return NULL;
    // Interpret input objects as numpy arrays
    PyObject *tsarray = PyArray_FROM_OTF(tsobj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *wsarray = PyArray_FROM_OTF(wsobj, NPY_CDOUBLE, NPY_IN_ARRAY);
    PyObject *gsarray = PyArray_FROM_OTF(gsobj, NPY_CDOUBLE, NPY_IN_ARRAY);
    // If that didn't work, throw an exception
    if(tsarray==NULL or wsarray==NULL or gsarray==NULL){
        Py_XDECREF(tsarray);    
        Py_XDECREF(wsarray);    
        Py_XDECREF(gsarray);    
        return NULL;
    }

    // Get pointers to the data as c++-types
    // TODO: find size of numpy arrays and check if consistent, then pass the
    // result of the check to system.
    // TODO: check overlap of integration range and supplied arrays
    double *ts = (double*)PyArray_DATA(tsarray);
    std::complex<double> *ws = (std::complex<double>*)PyArray_DATA(wsarray);
    std::complex<double> *gs = (std::complex<double>*)PyArray_DATA(gsarray);
   
    // Call the C++ functions to construct system and solve
    de_system sys = de_system(ts,ws,gs,islogw,islogg);
    Solution solution(sys,x0,dx0,ti,tf,order,rtol,atol,h0,full_output,interp);
    solution.solve();
    // Build output values
    std::list<std::complex<double>> sol,dsol;
    std::list<double> times;
    std::list<bool> wkbs;
    sol = solution.sol;
    dsol = solution.dsol;
    times = solution.times;
    wkbs = solution.wkbs;
    auto itd = dsol.begin();
    auto itt = times.begin();
    auto itwkb = wkbs.begin();
    int Nsol = sol.size();
    PyObject *pysol=PyList_New(Nsol),*pydsol=PyList_New(Nsol),*pytimes=PyList_New(Nsol),*pywkbs=PyList_New(Nsol),*retdict;
    Nsol = 0;
    for(auto it=sol.begin(); it!=sol.end(); ++it){
        PyList_SetItem(pysol,Nsol,Py_BuildValue("O",PyComplex_FromDoubles(std::real(*it), std::imag(*it)))); 
        PyList_SetItem(pydsol,Nsol,Py_BuildValue("O",PyComplex_FromDoubles(std::real(*itd), std::imag(*itd)))); 
        PyList_SetItem(pytimes,Nsol,Py_BuildValue("d",*itt)); 
        PyList_SetItem(pywkbs,Nsol,Py_BuildValue("i",*itwkb));
        ++itd; ++itt; ++itwkb, ++Nsol;
    };
    retdict = Py_BuildValue("{s:O,s:O,s:O,s:O}","sol",pysol,"dsol",pydsol,"t",pytimes,"types",pywkbs);

    // Clean up
    Py_DECREF(tsarray);
    Py_DECREF(wsarray);
    Py_DECREF(gsarray);
    Py_INCREF(Py_None);
    return retdict;
}
