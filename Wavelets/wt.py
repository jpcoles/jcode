import pywt, numpy, sys, random
from numpy import array, ones, zeros, linalg, eye, dot
from numpy.linalg import norm,inv
from scipy.linalg.iterative import cg, cgs
from scipy.fftpack import fft, ifft, fft2, ifft2, rfft, irfft, fftshift, ifftshift
from math import log, fabs
from copy import deepcopy

import pylab

N = 0
w = None
L0 = None
mode = None
level = None

def startup():
    global w, N, L0, mode, level

    print pywt.families()
    print pywt.wavelist('db')

    w = pywt.Wavelet('db6')
    mode = pywt.MODES.per

    print w
    print "vanishing_moments_psi:", w.vanishing_moments_psi
    print "vanishing_moments_phi:", w.vanishing_moments_phi

    N = 2**9
    print "max level = ", pywt.dwt_max_level(N, w.dec_len)
    L0 = numpy.zeros((N,N), 'double')
    if True:
        for i in xrange(0,N):
            L0[i][i-1], L0[i][i], L0[i-1][i] = (1., -2., 1.)
    else:
        for i in xrange(1,N):
            L0[i][i-1], L0[i][i], L0[i-1][i] = (1., -2., 1.)
        L0[0][0] = -2.

#L0[0][N-1], L0[0][0], L0[N-1][0] = (1, -2, 1)

#L0 = numpy.eye(N)

#for i in xrange(0,N):
        #L0[i] = [1,2,3,4,5,6,7,8]

    numpy.core.arrayprint.set_printoptions(threshold=N*N+1, linewidth=100000)

#print L0

#coeffs = pywt.wavedec2(L0, w) #, level=pywt.dwt_max_level(N, w.dec_len))
#print coeffs
#print pywt.waverec2(coeffs, w)
#
#coeffs = pywt.wavedec2(L0, w) #, level=pywt.dwt_max_level(N, w.dec_len))
#print coeffs

    print "max level = ", pywt.dwt_max_level(N, w.dec_len)
    #level = log(N, 2) /2 - 1 
#level = pywt.dwt_max_level(N, w.dec_len)


def cond(x,p=2):
    '''Compute the condition number of a matrix using
    the p-norm.'''

    return norm(x,p)*norm(inv(x),p)

def cond2(x):
    m = 1e10
    M = 1e-7
    evals = abs(numpy.linalg.eigvals(x))
    for v in evals:
        if 1e-7 < v < m: m = v
        if v > M: M = v
    c = M/m
    print m, M, c
    return c

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
        #d.reverse()
        coeffs_list.append(d)

    a = list(a)
    #a.reverse()
    coeffs_list.append(a)
    #coeffs_list.reverse()

    return coeffs_list


def my_waverec(data, wavelet, mode='sym', level=None):
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

    a = data[0:2]
    d = data[2:4]
    for i in xrange(int(level)):
        offs = pow(2,i+2)
        a = pywt.idwt(a, d, wavelet, mode)
        d = data[offs:2*offs]

    return a

def xfrm1d(L0):
    global N, level, w
    coeffs = my_wavedec(L0, w, mode=mode, level=level)
    ca = coeffs[0]
    #print coeffs[:0:-1]
    #for c in coeffs[:0:-1]:
    for c in coeffs[1:]:
        ca = numpy.concatenate((ca,c))
    return ca

def ixfrm1d(L0):
    global N, level, w
    return my_waverec(L0, w, mode=mode, level=level)

def xfrm2d(L0):
    global N, level
    for i in xrange(N):
        L0[i] = xfrm1d(L0[i])
    return L0

def ixfrm2d(L0):
    global N, level
    for i in xrange(N):
        L0[i] = ixfrm1d(L0[i])
    return L0

def test0():
    global L0, N

    L = deepcopy(L0)

    #xfrm(L)
    #ixfrm(L)
    #L[abs(L)<0.0001] = 0
    #print L
    #sys.exit(0)

    xfrm(L)
    L = L.T
    xfrm(L)
    L = L.T

    #L[L<0.01] = 0
    #print L[0]
    print L
    #print Li

    #Li= numpy.linalg.inv(L)
    #Li[Li<0.01] = 0
    #print  Li

    L[abs(L)<0.01] = 0
    print L
    Li= numpy.linalg.inv(L)
    Li[abs(Li)<0.01] = 0

    print "Li=", Li

    L = L.T
    ixfrm(L)
    #L[L<0.01] = 0
    #print L

    L = L.T
    ixfrm(L)

    L[abs(L)<0.01] = 0
    #print L[0]
    print L

def test1():
    global L0, N
    
    L = deepcopy(L0)

    xfrm(L)
    L = L.T
    xfrm(L)
    L = L.T

    evals = abs(numpy.linalg.eigvals(L))
    print max(evals), min(evals), max(evals)/min(evals)

def test2():
    global L0, N

    L = deepcopy(L0)
    rho = zeros(N, 'double')
    rho[0]  = 1.
    rho[N/2] = 1.

    print rho

    rho = xfrm1d(rho)

    xfrm2d(L)
    L = L.T
    xfrm2d(L)
    L = L.T

    x = linalg.solve(L, rho)
    #x[abs(x)<0.001] = 0

    x = ixfrm1d(x)
    print x
    F = []
    for i in xrange(len(x)-1):
        F.append(x[i+1] - x[i])

    print "F =", F
    print "*********************************"

def test3():
    global L0, N

    L = deepcopy(L0)

    rho = zeros(N, 'double')
    rho[0]  = 1.
    rho[N/2] = 1.

    print rho

    x = linalg.solve(L, rho)
    #x = linalg.solve(eye(N), rho)
    #x[abs(x)<0.001] = 0

    print x
    F = []
    for i in xrange(len(x)-1):
        F.append(x[i+1] - x[i])

    print "F =", F
    print "================================="
    
def test4():
    global L0, N

    L = deepcopy(L0)
    rho = zeros(N, 'double')
    rho[0]  = 1.
    rho[N/2] = 1.

    print rho

    rho = fft(rho)
    print "fft(rho) =", rho

    L = fft2(L)
    #print L

    x = linalg.solve(L, rho)
    print "x =", x
    #x[abs(x)<0.001] = 0

    x = ifft(x).real * N
    print "ifft(x) =", x
    F = []
    for i in xrange(len(x)-1):
        F.append(x[i+1] - x[i])

    print "F =", F
    print "--------------------------------"

def test5():
    global L0, N

    L = deepcopy(L0)
    rho = zeros(N, 'double')
    rho[0]  = 1.
    rho[N/2] = 1.

    print rho

    print fft(rho)
    rho = fftshift(fft(rho))
    print "fft(rho) =", rho

    L = fft(L).T
    L = fft(L).T

    L = fftshift(L)

    #print L

    x = linalg.solve(L, rho)
    print "x =", x
    #x[abs(x)<0.001] = 0

    x = ifftshift(ifft(x)).real * N
    print "ifft(x) =", x
    F = []
    for i in xrange(len(x)-1):
        F.append(x[i+1] - x[i])

    print "F =", F
    print "--------------------------------"

def test6():
    global L0, N

    L = deepcopy(L0)

    rho = zeros(N, 'double')
    rho[0]  = 1.
    rho[N/2] = 1.

    print rho

    x, info = cg(L, rho)
    #x = linalg.solve(eye(N), rho)
    #x[abs(x)<0.001] = 0

    print x
    F = []
    for i in xrange(len(x)-1):
        F.append(x[i+1] - x[i])

    print "F =", F
    print "================================="

def precondition2(L):
    m=0
    k = 1
    for i in range(log(N,2), 0, -1):
        k *= 2.0
        for j in range(0, 2**(i-1)):
            L[m] *= k*k
            m += 1
    L[N-1] *= k*k 

def precondition(L):
    m=0
    k = 1
    for i in range(log(N,2), 0, -1):
        k *= 2.0
        for j in range(0, 2**(i-1)):
            L[m] *= k
            m += 1
    L[N-1] *= k

def test7():
    global L0, N

    pylab.gray()
    pylab.matshow(L0)

    L = deepcopy(L0)
#   L = numpy.zeros((N,N), 'double')
#   for i in xrange(0,N):
#       for j in xrange(0,N):
#           if i == j: continue
#           L[i][j] = 1. / (i-j)

    P = zeros((N,N), 'double')
    m = 0
    k = 1.0
    for i in range(log(N,2), 0, -1):
        for j in range(0, 2**(i-1)):
            P[m][m] = k*k
            m += 1
        k *= 2.0
    P[N-1][N-1] = 1.0 * N*N

    rho = zeros(N, 'double')
    rho[random.sample(xrange(N), N/2)] = 1
    print rho
    rho = xfrm1d(rho)

    #print "cond(L):", cond2(L) 

    Lw = deepcopy(L)

    #Lw[numpy.abs(Lw) < 1e-5] = 0

    Lw = xfrm2d(Lw).T
    Lw = xfrm2d(Lw).T

    precondition2(Lw)


    #print "cond(Lw):", cond2(Lw) 

    PLw = Lw
    #PLw  = dot(P, Lw)
    #PLw  = dot(dot(P, Lw), P)
    #print numpy.max(PLw)
    #print numpy.min(PLw)
    #PLw[numpy.abs(PLw) <= 200.0] = 0.
    print "cond(PLwP):", cond2(PLw) 

    pylab.matshow(PLw)

    LI = linalg.inv(PLw)
    #print L
    #L[abs(L) < 0.000001] = 0
    L2 = zeros((N,N), 'double')
    #L2[numpy.abs(L) == 0] = 1.0
    L2[abs(LI) < 1e-7] = 0
    L2[abs(LI) >= 1e-7] = 1
    pylab.matshow(L2)
    pylab.show()

    #x, info = cg(L, rho, tol=1e-10)
    x, info = cgs(Lw, rho, tol=1e-7)
    x[abs(x)<0.001] = 0

    precondition(x)

    x = ixfrm1d(x)
    print x
    F = []
    for i in xrange(len(x)-1):
        F.append(x[i+1] - x[i])

    print "F =", F
    print "*********************************"

def test8():
    global L0, N

    L = deepcopy(L0)
    rho = zeros(N, 'double')

    rho[random.sample(xrange(N), N/2)] = 1

    print rho

    LI = linalg.inv(L)
    #print L
    #print LI
    #I = numpy.dot(L,LI)
    #I[abs(I)<0.001] = 0
    #print I

    t = numpy.greater(rho, 0)
    X = numpy.zeros((N,N))
    for i in xrange(N):
        X[0][i] = i
    print X
    LIC = numpy.compress(t, LI, 1)
    print LIC
    LIC = numpy.compress(t, LIC, 0)
    print LIC
    LICI = linalg.inv(LIC)
    print LICI
    
def test9():
    global L0, N

    #L = numpy.diag([-2]*N)
    L = deepcopy(L0)
    rho = zeros(N, 'double')

    rho[0] = 1
    rho[N/2] = 1
    #rho[random.sample(xrange(N), N/2)] = 1

    print rho
    print L

    LI = linalg.inv(L)
    #print L
    print LI
    I = numpy.dot(LI, rho)
    I[abs(I)<0.001] = 0
    print I

    print "===================================="

    L = numpy.diag([-2]*N)
    print L

    LI = linalg.inv(L)
    #print L
    print LI
    I = numpy.dot(LI, rho)
    I[abs(I)<0.001] = 0
    print I
    return

    t = numpy.greater(rho, 0)
    X = numpy.zeros((N,N))
    for i in xrange(N):
        X[0][i] = i
    print X
    LIC = numpy.compress(t, LI, 1)
    print LIC
    LIC = numpy.compress(t, LIC, 0)
    print LIC
    LICI = linalg.inv(LIC)
    print LICI

def test10():
    global L0, N

    #L = numpy.diag([-2]*N)
    L = deepcopy(L0)
    rho = zeros(N, 'double')

    rho[1] = 1
    rho[N/2] = 1
    #rho[random.sample(xrange(N), N/2)] = 1

    print rho
    print L

    LI = linalg.inv(L)
    #print L
    print LI
    I = numpy.dot(LI, rho)
    I[abs(I)<0.001] = 0
    print I

    print "===================================="

    t = numpy.greater(rho, 0)
    t = xrange(1,N-1)

    LC = numpy.compress(t, L, 1)
    print "LC", LC
    LC = numpy.compress(t, LC, 0)
    print "LC", LC
    LCI = linalg.inv(LC)
    print "LCI", LCI

    rho = numpy.compress(t, rho)
    print "rho", rho
    I = numpy.dot(LCI, rho)
    I[abs(I)<0.001] = 0
    print I

startup()
#test0()
#test1()
#test2()
#test3()
#test4()
#test5()
#test6()
test7()
#test8()
#test9()
#test10()

