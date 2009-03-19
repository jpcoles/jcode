# -*- coding: utf-8 -*-

# Copyright (c) 2006-2007 Filip Wasilewski <filip.wasilewski@gmail.com>
# See COPYING for license details.

# $Id: numerix.py 86 2007-10-29 22:36:37Z filipw $

"""
Thin wrapper for numeric modules. Modify this to use wavelets with libraries other than NumPy.

Provides efficient mathematical functions and array datatypes.
"""

from numpy import ndarray, array, asarray
from numpy import empty, zeros, linspace, arange
from numpy import intp, float64, float32
from numpy import transpose, concatenate
from numpy import cumsum, cos, diff, exp, sinc
from numpy import argmax, mean
from numpy import convolve
from numpy import where, less, greater
from numpy.fft import fft

default_dtype = float64

def as_float_array(source):
    if isinstance(source, ndarray) and (source.dtype == float64 or source.dtype == float32):
        return source
    return array(source, default_dtype)

def contiguous_array_from_any(source):
    raise DeprecationWarning()
    return contiguous_float64_array_from_any(source)

def contiguous_float64_array_from_any(source):
    return array(source, float64) # ensure contiguous

def contiguous_float32_array_from_any(source):
    return array(source, float32) # ensure contiguous

def astype(source, dtype):
    return asarray(source, dtype)

def float64_memory_buffer_object(size):
    return zeros((size,), float64)

def float32_memory_buffer_object(size):
    return zeros((size,), float32)

def is_array_type(arr, typ):
    return isinstance(arr, ndarray) and arr.dtype == typ

def keep(arr, keep_length):
    length = len(arr)
    if keep_length < length:
        left_bound = (length - keep_length) / 2
        return arr[left_bound:left_bound+keep_length]
    return arr

def integrate(arr, step):
    integral = cumsum(arr)
    integral *= step
    return integral
