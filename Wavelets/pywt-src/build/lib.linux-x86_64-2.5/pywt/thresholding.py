# -*- coding: utf-8 -*-

# Copyright (c) 2006-2007 Filip Wasilewski <filip.wasilewski@gmail.com>
# See COPYING for license details.

# $Id: thresholding.py 85 2007-09-22 14:49:38Z filipw $

"""Thresholding routines"""

__all__ = ['soft', 'hard', 'greater', 'less', 'zero', 'copy']

import numerix

def soft(data, value, substitute=0):
    mvalue = -value
    
    cond_less = numerix.less(data, value)
    cond_greater = numerix.greater(data, mvalue)

    data = numerix.where(cond_less & cond_greater, substitute, data)
    data = numerix.where(cond_less, data + value, data)
    data = numerix.where(cond_greater, data - value, data)
    
    return data

def hard(data, value, substitute=0):
    mvalue = -value
    
    cond = numerix.less(data, value)
    cond &= numerix.greater(data, mvalue)
    
    return numerix.where(cond, substitute, data)
    
def greater(data, value, substitute=0):
    return numerix.where(numerix.less(data, value), substitute, data)

def less(data, value, substitute=0):
    return numerix.where(numerix.greater(data, value), substitute, data)

def zero(data, *args):
    if isinstance(data, numerix.ndarray):
        return numerix.zeros(data.shape, data.dtype)
    return numerix.zeros(len(data))

def copy(data, *args):
    return numerix.array(data)
