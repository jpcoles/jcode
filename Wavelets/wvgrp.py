from copy import deepcopy
import numpy
from numpy import *
import pylab
from pywt import *
import random
from math import log
import pywt

def my_wavedec(data, wavelet, mode='sym', level=None):
    """
    Multilevel 1D Discrete Wavelet Transform of data.
    Returns coefficients list - [cAn, cDn, cDn-1, ..., cD2, cD1]

    data    - input data
    wavelet - wavelet to use (Wavelet object or name string)
    mode    - signal extension mode, see MODES
    level   - decomposition level. If level is None then it will be
              calculated using `dwt_max_level` function.
    """

    if not isinstance(wavelet, pywt.Wavelet):
        wavelet = pywt.Wavelet(wavelet)

    if level is None:
        level = pywt.dwt_max_level(len(data), wavelet.dec_len)
    elif level < 0:
        raise ValueError("Level value of %d is too low . Minimum level is 0." % level)

    coeffs_list = []

    a = data
    for i in xrange(level):
        a, d = pywt.dwt(a, wavelet, mode)
        d = list(d)
        d.reverse()
        coeffs_list.append(d)

    a = list(a)
    a.reverse()
    coeffs_list.append(a)
    #coeffs_list.reverse()

    return coeffs_list

def xfrm1d(L0):
    global N, level, w
    coeffs = my_wavedec(L0, w, mode=mode, level=level)
    ca = coeffs[0]
    #print coeffs[:0:-1]
    #for c in coeffs[:0:-1]:
    for c in coeffs[1:]:
        ca = numpy.concatenate((ca,c))
    return ca

def xfrm2d(L0):
    global N, level
    for i in xrange(N):
        L0[i] = xfrm1d(L0[i])
    return L0

w = Wavelet('db3')
mode = MODES.per
N = 2**8
P = zeros((N,N), 'double')

#pylab.gray()
for i in xrange(N):
    P[i][random.sample(xrange(N), N/8)] = 1

def load_image():
    global P,N
    I = pylab.imread("arnold.png")

    assert len(I) == N
    assert len(I[0]) == N

    for i in xrange(N):
        for j in xrange(N):
            P[i][j] = I[i][j][0]
    del I
    return P

def one():
    Pw = []
    for i,n in enumerate(xrange(0,log(N,2))):
        pylab.figure()
        Pw.append(deepcopy(P))
        print len(Pw)
        level = n
        print "n=",n
        Pw[i] = xfrm2d(Pw[i]).T
        Pw[i] = xfrm2d(Pw[i]).T
        #Pw[Pw==0.0] = 255
        pylab.matshow(Pw[i])

def two():
    global P
    eps = 1
    Pw = deepcopy(P)
    for i,n in enumerate(xrange(1,4)):
        X = pywt.wavedec2(Pw, 'sym5', level=n)
        Pw = X[1][0]
        #Pw[abs(Pw) < eps] = 0
        #Pw[abs(Pw) >= eps] = 1
        pylab.matshow(Pw)

print pywt.families()

pylab.gray()
P = load_image()
two()

pylab.show()
