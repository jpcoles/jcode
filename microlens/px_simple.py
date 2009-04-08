import sys
from numpy import empty, reshape, ones, ndindex, vectorize, apply_over_axes, apply_along_axis, array, repeat, fromfunction, zeros_like, clip
from numpy.random import standard_normal, poisson, uniform
from math import pow
import numpy
from numpy import hypot, sum, max, linspace, arange, arange, dot, matrix
from numpy.core.umath import exp, sqrt, log

from scipy.linalg import det
from scipy.special import hermite
from scipy.misc import factorial

import pylab
#from pymc import DiscreteUniform

# Constants
G           =       6.67300e-11             # m^3/kg/s^2
c2          =       2.99792458e8            # m/s
dLS         =       1.0
dL          =       1.0
dS          =       1.0
Msol        =       1.98892e30              # kg

num_trials  =       1000
gamma_tot   =       10000000                   # Total photon count
grid_width  =       2.0                     # arcsec
grid_rad    =       30                      # pixels
grid_size   =       2*grid_rad + 1          # pixels
grid_pos    =       (0.0, 0.0)              # arcsec
cell_size   =       grid_width / grid_size  # arcsec/pixel

# True values of the star
zT = (Msol, 0.025, complex(-0.75, 0.25))                 # Mass and position

##############################################################################

def deVaucouleurs(R, Ie=.40, m=4, Re=.50):
    bm = 2*m - 0.324
    d = hypot(R.real, R.imag)
    return Ie * exp(-bm * (d/Re)**.25)

def profile_exponential(R, Ie=.40, Rd=.50):
    d = hypot(R.real, R.imag)
    return Ie * exp(-d/Rd)

def profile_exponential2(R):
    d = hypot(R.real, R.imag)
    return (1+d)**-3

def raytrace(theta_E, z):
    """theta_E - Einstein radius of star
       z - position of star (complex)"""
    global theta
    return theta - (theta_E**2) * (theta-z) / hypot((theta-z).real, (theta-z).imag)

def mapSB(sb, mp):
    res = zeros_like(sb)
    for i in xrange(mp.shape[0]):
        for j in xrange(mp.shape[1]):
            r = mp[i,j] / cell_size
            if 0 <= r.real <= grid_size and 0 <= r.imag <= grid_size:
                res[i,j] = sb[int(r.real), int(r.imag)]

def normalize(sb):
    global gamma_tot
    sb *= gamma_tot / sum(sb)
    return sb

def random_sample(n):
    global zT
    for i in xrange(n):
        #x = uniform(low=zT[1].real-.2, high=zT[1].real+.2)
        #y = uniform(low=zT[1].imag-.2, high=zT[1].imag+.2)
        x,y = uniform(low=-1, high=1, size=2)
        yield 0.025, complex(x,y)

def rE_sample(min, max, steps):
    global zT

    step = (max-min) / steps
    for i in xrange(steps):
        yield i, (min + i*step), zT[2]

def B_n(n1,n2, beta):
    beta2 = beta**2
    K1    = 1.0/sqrt(2**n1 * sqrt(pi) * factorial(n1,1) * beta))
    K2    = 1.0/sqrt(2**n2 * sqrt(pi) * factorial(n2,1) * beta))
    H1    = hermite(n1)
    H2    = hermite(n2)
    return lambda x: (K1 * H1(x.real/beta) * exp(-x.real**2/(2*beta2))) * 
                     (K2 * H2(x.imag/beta) * exp(-x.imag**2/(2*beta2)))

def Bfuncs(n, beta):
    l = []
    for i in xrange(n):
        for j in xrange(n):
            l.append(B_n(i,j,beta))
    return l

##############################################################################

L = 1

sample = random_sample
surface_brightness = profile_exponential2

# Create the position grid
theta = empty((grid_size, grid_size), 'complex')
for i in xrange(grid_size):
    for j in xrange(grid_size):
        theta[j,i] = complex(i-grid_rad, -(j-grid_rad)) * cell_size
    

# Create the target galaxy and data (with noise)
galaxySB  = normalize(surface_brightness(raytrace(zT[1], zT[2])))
vaddnoise = vectorize(lambda x: poisson(lam=x))
data      = vaddnoise(galaxySB)
data[data < 1] = 1

sigma = sqrt(data)
data2 = grid_size**2

#re = DiscreteUniform('re', lower=0.015, upper=0.035)
#zs = 

# Now sample the space of rE,z and compute the liklihood of each point
probs = empty((num_trials, 4), 'float')
for i,rE,z in rE_sample(0.01, .22, num_trials):

    f = normalize(surface_brightness(raytrace(rE, z)))
    assert sum(f) - gamma_tot < 1e-4

    P     = empty((1,L), 'float')
    for l in xrange(L):
        P[l] = sum(
    P     = matrix(sum(f))
    C_inv = matrix(sum(f**2 / data))
    C     = C_inv.I
    prob  = log(det(C) ** (L/2.)) + P.T * C * P

    probs[i] = rE, z.real, z.imag, prob

#data = mapSB(galaxySB, beta)

# Normalize so that the maximum liklihood is 1
probs[:,3] -= max(probs[:,3])
probs[:,3]  = exp(probs[:,3])

if True:
    pylab.figure()
    pylab.plot(probs[:,0], probs[:,3])
    pylab.axvline(x=zT[1])
    pylab.xlabel('rE')
    pylab.ylabel('Liklihood')

if False:
    pylab.figure()
    pylab.imshow(noise)
    pylab.title('Noise')

if False:
    pylab.figure()
    pylab.imshow(galaxySB)
    pylab.title('Galaxy Surface Brightness')

if True:
    pylab.figure()
    pylab.imshow(data)
    pylab.title('data')

if False:
    pylab.figure()
    pylab.plot(probs)

if False:
    xi = linspace(-grid_width/2,grid_width/2, grid_size*10)
    yi = linspace(-grid_width/2,grid_width/2, grid_size*10)
    zi = pylab.griddata(xs,ys,probs, xi,yi)

    pylab.figure()
    pylab.contour(xi,yi,zi)
    pylab.colorbar()
    pylab.scatter(xs, ys, marker='o',c='b',s=5)
    pylab.xlim(-grid_width/2, grid_width/2)
    pylab.ylim(-grid_width/2, grid_width/2)


pylab.show()

