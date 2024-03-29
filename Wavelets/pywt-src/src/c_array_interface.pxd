# Copyright (c) 2006-2007 Filip Wasilewski <filip.wasilewski@gmail.com>
# See COPYING for license details.

# see http://numeric.scipy.org/array_interface.html
# see http://new.scipy.org/Wiki/Cookbook/ArrayStruct_and_Pyrex

# $Id: c_array_interface.pxd 81 2007-09-11 23:21:36Z filipw $

cimport c_python

cdef extern from "array_interface.h":
    
    ctypedef c_python.Py_intptr_t intp

    ctypedef struct PyGenericArrayInterface:
        int two                         # contains array interace version number (min. 2)
        int nd                          # number of dimensions
        char typekind                   # kind in array --- character code of typestr
        int itemsize                    # size of each element
        int flags                       # flags indicating how the data should be interpreted
        c_python.Py_intptr_t *shape     # A length-nd array of shape information
        c_python.Py_intptr_t *strides   # A length-nd array of stride information
        void *data                      # A pointer to the first element of the array
        c_python.PyObject* descr        # NULL or data-description -- must set ARR_HAS_DESCR flag

    ctypedef enum PyGenericArray_KINDS:
        PyArrayKind_BOOL
        PyArrayKind_INT
        PyArrayKind_UINT
        PyArrayKind_FLOAT
        PyArrayKind_COMPLEX
        PyArrayKind_STRING
        PyArrayKind_UNICODE
        PyArrayKind_OBJECT
        PyArrayKind_RECORD
        PyArrayKind_VOID
        PyArrayKind_BIT
        PyArrayKind_OTHER

    ctypedef enum PyGenericArray_FLAGS:
        GA_CONTIGUOUS
        GA_FORTRAN
        GA_ALIGNED
        GA_NOTSWAPPED
        GA_WRITEABLE
        GA_ARR_HAS_DESCR

    cdef double* PyArrayInterface_DATA_AS_FLOAT64_C_ARRAY(PyGenericArrayInterface* )
    cdef double* PyArrayInterface_DATA_AS_FLOAT64_C_ARRAY_RO(PyGenericArrayInterface* )

    cdef float* PyArrayInterface_DATA_AS_FLOAT32_C_ARRAY(PyGenericArrayInterface* )
    cdef float* PyArrayInterface_DATA_AS_FLOAT32_C_ARRAY_RO(PyGenericArrayInterface* )

    cdef void* PyArrayInterface_DATA_AS_FLOAT_C_ARRAY(PyGenericArrayInterface* )
    cdef void* PyArrayInterface_DATA_AS_FLOAT_C_ARRAY_RO(PyGenericArrayInterface* )

    cdef int PyArrayInterface_IS_C_ARRAY(PyGenericArrayInterface* )
    cdef int PyArrayInterface_IS_C_ARRAY_RO(PyGenericArrayInterface* )

    cdef int PyArrayInterface_CHECK(PyGenericArrayInterface* )
    cdef int PyArrayInterface_CHECK_1D(PyGenericArrayInterface* )
    cdef int PyArrayInterface_CHECK_2D(PyGenericArrayInterface* )

    int PyArrayInterface_TWO(PyGenericArrayInterface*)            
    int PyArrayInterface_ND(PyGenericArrayInterface*)             
    char PyArrayInterface_TYPEKIND(PyGenericArrayInterface*)     
    int PyArrayInterface_ITEMSIZE(PyGenericArrayInterface*)     
    int PyArrayInterface_FLAGS(PyGenericArrayInterface*)        
    c_python.Py_intptr_t * PyArrayInterface_SHAPES(PyGenericArrayInterface*)
    c_python.Py_intptr_t PyArrayInterface_SHAPE(PyGenericArrayInterface*, int)
    c_python.Py_intptr_t * PyArrayInterface_STRIDES(PyGenericArrayInterface*)
    c_python.Py_intptr_t PyArrayInterface_STRIDE(PyGenericArrayInterface*, int)
    void* PyArrayInterface_DATA(PyGenericArrayInterface*)
    c_python.PyObject* PyArrayInterface_DESCR(PyGenericArrayInterface*)
