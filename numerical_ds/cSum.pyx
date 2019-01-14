from numpy cimport ndarray
cimport numpy
cimport cython

@cython.boundscheck(False)
def cSum(ndarray[numpy.complex128_t, ndim=1] a not None):
    cdef Py_ssize_t i
    cdef Py_ssize_t n = a.shape[0]
    cdef double complex m = 0.0
    for i in range(n):
        m += a[i]
    return m 
 
def cSumd(ndarray[numpy.float64_t, ndim=1] a not None):
    cdef Py_ssize_t i
    cdef Py_ssize_t n = a.shape[0]
    cdef double m = 0.0
    for i in range(n):
        m += a[i]
    return m
