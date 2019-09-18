#define PY_SSIZE_T_CLEAN
#include "_python.hpp"
#include "_pyoscode.hpp"
#include "system.hpp"
#include "solver.hpp"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>


/* Initialise the module */
#ifdef PYTHON3
PyMODINIT_FUNC PyInit__pyoscode(void){
    import_array();
    return PyModule_Create(&_pyoscodemodule);
}
#else 
PyMODINIT_FUNC init_pyoscode(void){
    PyObject *m = Py_InitModule3("_pyoscode", module_methods, module_docstring);
    if(m==NULL)
        return;
    // Load numpy functionality
    import_array();
}
#endif

/* Function to run the solver */ 
static PyObject *_pyoscode_solve(PyObject *self, PyObject *args, PyObject *kwargs){

    int islogw=0,islogg=0,order=3;
    const char* full_output="";
    double ti,tf,rtol=1e-4,atol=0.0,h0=1.0;
    std::complex<double> x0,dx0;
    PyObject *tsobj, *wsobj, *gsobj;
    // Define keywords
    static const char *kwlist[] =
    {"ts","ws","gs","ti","tf","x0","dx0","logw","logg","order","rtol","atol","h","full_output",NULL};

    // Interpret input arguments.
    if (!PyArg_ParseTupleAndKeywords(args,kwargs,"OOOddDD|iiiddds",const_cast<char**>(kwlist),&tsobj,&wsobj,&gsobj,&ti,&tf,&x0,&dx0,&islogw,&islogg,&order,&rtol,&atol,&h0,&full_output))
        return NULL;
    // Interpret input objects as numpy arrays
    PyObject *tsarray = PyArray_FROM_OTF(tsobj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *wsarray = PyArray_FROM_OTF(wsobj, NPY_CDOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *gsarray = PyArray_FROM_OTF(gsobj, NPY_CDOUBLE, NPY_ARRAY_IN_ARRAY);
    // If that didn't work, throw an exception
    if(tsarray==NULL or wsarray==NULL or gsarray==NULL){
        Py_XDECREF(tsarray);    
        Py_XDECREF(wsarray);    
        Py_XDECREF(gsarray);    
        return NULL;
    }

    // Get pointers to the data as c++-types
    PyArrayObject *tsarray_arr = reinterpret_cast<PyArrayObject*>(tsarray);
    PyArrayObject *wsarray_arr = reinterpret_cast<PyArrayObject*>(wsarray);
    PyArrayObject *gsarray_arr = reinterpret_cast<PyArrayObject*>(gsarray);
    double *ts = (double*)PyArray_DATA(tsarray_arr);
    std::complex<double> *ws = (std::complex<double>*)PyArray_DATA(wsarray_arr);
    std::complex<double> *gs = (std::complex<double>*)PyArray_DATA(gsarray_arr);
    
    // Get array sizes
    int tssize = (int)PyArray_SIZE(tsarray_arr);
    int wssize = (int)PyArray_SIZE(wsarray_arr);
    int gssize = (int)PyArray_SIZE(gsarray_arr);

    try{
        // Check 1: same size arrays, large enough arrays
        if(tssize != wssize or wssize != gssize or tssize != gssize)
            throw "The sizes of the input arrays (ts, ws, gs) do not match. Please supply arrays of the same size.";
        else if(tssize < 2)
            throw "The input arrays (ts, ws, gs) have to have at least size 2.";

        // Check 2: ts must be strictly monotonous
        double dir = ts[1] - ts[0];
        double diff;
        for(int i=1; i<tssize; i++){
            if(dir > 0)
                diff = ts[i] - ts[i-1]; 
            else
                diff = ts[i-1] - ts[i];
            if(diff < 0)
                throw "The time array (ts) must be strictly monotonous.";
        }

        // Check 3: integration limits ti, tf must lie inside ts
        if(ti < ts[0] or ti > ts[tssize-1])
            throw "The start of integration, ti, must lie inside the array ts.";
        else if(tf < ts[0] or tf > ts[tssize-1])
            throw "The end of integration, tf, must lie inside of the array ts.";
    }
    catch(const char* errormsg){
        PyErr_SetString(PyExc_TypeError, errormsg);
        // Clean up
        Py_DECREF(tsarray);
        Py_DECREF(wsarray);
        Py_DECREF(gsarray);
        return (PyObject *) NULL;
    }

    // Call the C++ functions to construct system and solve
    de_system sys = de_system(ts,ws,gs,islogw,islogg);
    Solution solution(sys,x0,dx0,ti,tf,order,rtol,atol,h0,full_output);
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
        Py_complex x_complex, dx_complex;
        x_complex.real = std::real(*it);
        x_complex.imag = std::imag(*it);
        dx_complex.real = std::real(*itd);
        dx_complex.imag = std::imag(*itd);
        PyList_SetItem(pysol,Nsol,Py_BuildValue("D", &x_complex)); 
        PyList_SetItem(pydsol,Nsol,Py_BuildValue("D",&dx_complex)); 
        PyList_SetItem(pytimes,Nsol,Py_BuildValue("d",*itt)); 
        PyList_SetItem(pywkbs,Nsol,Py_BuildValue("i",*itwkb));
        ++itd; ++itt; ++itwkb; ++Nsol;
    };
    retdict = Py_BuildValue("{s:O,s:O,s:O,s:O}","sol",pysol,"dsol",pydsol,"t",pytimes,"types",pywkbs);

    // Clean up
    Py_DECREF(tsarray);
    Py_DECREF(wsarray);
    Py_DECREF(gsarray);
    Py_DECREF(pysol);
    Py_DECREF(pydsol);
    Py_DECREF(pytimes);
    Py_DECREF(pywkbs);
    return retdict;
}
