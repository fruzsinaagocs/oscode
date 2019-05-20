#pragma once
#include <array>
#include "_python.hpp"
#include <string>
#include <functional>
#include <Eigen/Dense>
#include <iomanip>
#include <fstream>
#include <complex>
#include <iostream>
#include <cmath>

/* Docstrings */
static char module_docstring[] = 
"pyoscode: this module provides an interface for oscode, for solving oscillatory ordinary differential equations with the RKWKB method.";

static char solve_docstring[] =
"Runs the solver";

/* Available functions */
static PyObject *_pyoscode_solve(PyObject *self, PyObject *args, PyObject *kwargs);

/* Module interface */
static PyMethodDef module_methods[] = {
    {"solve", (PyCFunction) _pyoscode_solve, METH_VARARGS | METH_KEYWORDS, solve_docstring},
    {NULL, NULL, 0, NULL}
};
 
#ifdef PYTHON3
static struct PyModuleDef _pyoscodemodule = {
    PyModuleDef_HEAD_INIT,
    "_pyoscode", 
    module_docstring,
    -1,
    module_methods
};
#endif
