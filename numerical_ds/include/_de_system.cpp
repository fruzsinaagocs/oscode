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

// Docstrings for the module and its members
static char module_docstring[] = 
"This module provides an interface for setting up a system of differential equations to be solved with the RKWKB method, using C++.";

static char system_docstring[] = 
"Set up a system of differential equations to be solved with the RKWKB method.";

// Declaration of the de_system constructor
static PyObject *_de_system_de_system(PyObject *self, PyObject *args);

// Defining the members of this module
static PyMethodDef module_methods[] = {
    {"de_system", _de_system_de_system, METH_VARARGS, system_docstring},
    {NULL, NULL, 0, NULL}
};

// Init function for the module
PyMODINIT_FUNC init_de_system(void){
    
    PyObject *m = Py_InitModule3("_de_system", module_methods, module_docstring);
    if(m==NULL)
        return;
    // Load numpy functionality
    import_array();
}

// Define constructor function
static PyObject *_de_system_de_system(PyObject *self, PyObject *args){

    int islogw, islogg;
    PyObject *tsobj, *wsobj, *gsobj;

    // Interpret input arguments.
    if (!PyArg_ParseTuple(args,"iiOOO",&islogw,&islogg,&tsobj,&wsobj,&gsobj))
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
    double *ts = (double*)PyArray_DATA(tsarray);
    std::complex<double> *ws = (std::complex<double>*)PyArray_DATA(wsarray);
    std::complex<double> *gs = (std::complex<double>*)PyArray_DATA(gsarray);
   
    // Call the C++ function
    de_system sys = de_system(ts,ws,gs,islogw,islogg);

    // Clean up
    Py_DECREF(tsarray);
    Py_DECREF(wsarray);
    Py_DECREF(gsarray);
    return Py_None;
}
