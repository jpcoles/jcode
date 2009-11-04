#
# ml.py - Weak microlensing simulator
#
# Copyright 2009 Jonathan Coles
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division
import sys, time
from math import pow, pi

import numpy
from numpy import empty, zeros, ndindex, vectorize, abs, logical_not
from numpy import sum, linspace, mat, abs, inner, asarray_chkfinite, ravel
from numpy.random import poisson, seed
from numpy.core.umath import exp, sqrt, log

from scipy.linalg import eig
from scipy.special import hermite, jn
from scipy.misc import factorial

# Use LaTeX to generate text for the plots
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('text', usetex=True)

import pylab

numpy.set_printoptions(threshold=numpy.nan)


#===============================================================================
# A variety of different surface brightness distribution functions.
#===============================================================================

def deVaucouleurs(R, Ie=1, m=4, Re=1.0):
    bm = 2*m - 0.324
    return Ie * exp(-bm * sqrt((abs(R)/Re)**2 + .056**2)**(1/m))

def profile_exponential(R, Ie=.40, Rd=.50):
    return Ie * exp(-abs(R)/Rd)

def profile_exponential2(R):
    return (1+abs(R))**-2

def profile_exponential3(R): # AfterEaster run1
    return (1 + sqrt(R.real**2 + 0.8*R.imag**2))**-2

def profile_exponential4(R): # AfterEaster run2
    return (1 + sqrt(0.8*R.real**2 + R.imag**2))**-2

def profile_exponential5(R):
    return (.2+abs(R))**-2

def profile_exponential5_2(R):
    dr = complex(0.03,0.03)
    return (.2+abs(R-dr))**-2 + (.2+abs(R+dr))**-2

#===============================================================================

def star_track_symmetric(n):
    yield 0, complex(0)

    if n > 2: 
        pos = linspace(.5, -.5, n-1)
    elif n == 2:
        pos = [0.0]
    else:
        return

    for i, p in enumerate(pos):
        yield i+1, complex(p, -closest_star_approach)

def star_track_nonsymmetric(n):
    yield 0, None #complex(0)

    for i in range(1, n):
        yield i, complex(i*cell_size*5, -closest_star_approach)

    #pos = linspace(2*cell_size, 10*cell_size, n-1)
    #for i, p in enumerate(pos):
        #yield i+1, complex(p, -closest_star_approach)

#===============================================================================

def model(m):
    """ Generate model data for all epochs.  Accepts a function that returns
    the surface brightness distribution.  This uses star_track() to get the
    position of the star at each epoch.  This will be normalized to the global
    value gamma_tot by normalize().
    """
    data = empty((nepochs, grid_size, grid_size), 'float')
    for t,z in star_track(nepochs):
        if t == 0:
            xx = raytrace()
            mask = 1
        else:
            xx = raytrace(rE_true, z)
            mask = star_mask(z)

        data[t] = normalize(m(xx) * mask)

    return data

def model_two_basis_functions():
    """ A test function that returns a model similar to model(), except
    that it uses the shapelet basis functions as the surface brightness
    and does not normalize.
    """
    data = empty((nepochs, grid_size, grid_size))
    beta2 = beta**2
    for t,z in star_track(nepochs):

        if t == 0:
            x = raytrace()
        else:
            x = raytrace(rE_true, z)

        n  = 0
        K1 = 1.0/sqrt(2**n * sqrt(pi) * factorial(n,1) * beta)
        H1 = hermite(n)

        data[t]  = (K1 * H1(x.real/beta) * exp(-x.real**2/(2*beta2))) * \
                   (K1 * H1(x.imag/beta) * exp(-x.imag**2/(2*beta2)))
#        data[t] *= 100

#        n  = 1
#        K1 = 1.0/sqrt(2**n * sqrt(pi) * factorial(n,1) * beta)
#        H1 = hermite(n)
#
#        data[t] += (K1 * H1(x.real/beta) * exp(-x.real**2/(2*beta2))) * \
#                   (K1 * H1(x.imag/beta) * exp(-x.imag**2/(2*beta2)))

    return data

#===============================================================================

def raytrace(theta_E=None, z=None):
    """ Raytrace an image plane theta (global) to a source plane assuming the
    presence of a star at z with Einstein radias theta_E. If theta_E and z are
    both none, just return the unlensed source plane (i.e., theta). z is
    complex number with real and imaginary parts corresponding to the x and y
    position, respectively.
    """
    if theta_E is None and z is None: return theta
    f = theta - (theta_E**2) * (theta-z) / (abs(theta-z)**2)
    return asarray_chkfinite(f)

def normalize(sb):
    """ Normalize a surface brightness distribution by gamma_tot (global). sb
    may be an n-dimensional array.
    """
    sb *= gamma_tot / sum(sb)
    return sb

#===============================================================================

def star_mask(z):
    return logical_not(abs(theta-z) < star_mask_size)

#===============================================================================

def demo_lensing(N):
    """Generate some plots of raytraced surface brightnesses."""
    surface_brightness = profile_exponential5
    pylab.figure()
    for t,z in star_track(nepochs):
        xx = normalize(surface_brightness(raytrace(rE_true, z)))
        pylab.subplot(N/2, N/2, t+1)
        pylab.imshow(xx.T)

    pylab.suptitle("Lensing demo  rE=%.3f" % rE_true)
    pylab.savefig("LensingDemo.png")

def save_basis_functions(N):
    """Generate N 2D shapelets and plot."""

    L = N**2

    beta2 = beta**2
    B     = empty((grid_size, grid_size)) # Don't want matrix behaviour here

    #---------------------------------------------------------------------------
    # Basis function constants, and hermite polynomials
    #---------------------------------------------------------------------------

    vals = [[n, 1.0/sqrt((2**n) * sqrt(pi) * factorial(n,1) * beta), hermite(n), 0,0] 
            for n in xrange(N)]
    expreal = exp(-theta.real**2/(2*beta2))
    expimag = exp(-theta.imag**2/(2*beta2))
    for n,K,H,v0,v1 in vals:
        vals[n][3] = K * H(theta.real/beta) * expreal
        vals[n][4] = K * H(theta.imag/beta) * expimag

    pylab.figure()
    l=0
    for v1 in vals:
        for v2 in vals:
            B = v1[3] * v2[4]
            pylab.subplot(N,N, l+1)
            pylab.axis('off')
            pylab.imshow(B.T)
            l += 1
    pylab.suptitle("Shapelets N=%i Beta=%.4f" % (N, beta))
    #pylab.savefig("B%i.png" % N)
    pylab.show()

def save_bessel_functions(N):
    """Generate N 2D shapelets and plot."""

    beta2 = beta**2
    B     = empty((grid_size, grid_size)) # Don't want matrix behaviour here

    #---------------------------------------------------------------------------
    # Basis function constants, and hermite polynomials
    #---------------------------------------------------------------------------

    vals = [[n, 1.0/sqrt((2**n) * sqrt(pi) * factorial(n,1) * beta), 0, 0,0] 
            for n in xrange(N)]
    expreal = exp(-theta.real**2/(2*beta2))
    expimag = exp(-theta.imag**2/(2*beta2))
    for n,K,H,_,_ in vals:
        vals[n][3] = K*jn(n,theta.real) * expreal
        vals[n][4] = K*jn(n,theta.imag) * expimag

    pylab.figure()
    l=0
    for v1 in vals:
        for v2 in vals:
            B = v1[3] * v2[4]
            pylab.subplot(N,N, l+1)
            pylab.axis('off')
            pylab.imshow(B.T)
            l += 1
    pylab.suptitle("Shapelets N=%i Beta=%.4f" % (N, beta))
    #pylab.savefig("B%i.png" % N)
    pylab.show()

# Note: This has been heavily tailored to generate specially labeled plots for
# the paper.
def plot_reconstructions(i, data, rE, L, C, P, N):
    """ Reconstruct the surface brightness distribution from the model and plot."""
    nepochs   = data.shape[0]
    grid_size = data.shape[1]

    fn = C*P
    #for n in xrange(N): print fn[n,0]

    for_paper = True

    for t,z in star_track(nepochs):
        f = zeros((grid_size, grid_size), 'float')

        if t == 0:
            mask = 1
        else:
            mask = star_mask(z)

        # Could probably use tensordot() here.
        for n in xrange(N): f += fn[n,0] * L[n,t]

        f       *= mask
        data[t] *= mask

        if t == 0: f0 = f

        n = 3
        m = nepochs
        m = 1
        cmap = pylab.cm.jet

        diff  = 100 * abs(data[t] - f) / data[t]
        diff1 = 100 * abs(data[t] - f0) / data[t]

        print 'diff shape is', diff.shape

#       pylab.figure()
#       pylab.hist(diff.flatten())

        if for_paper:

            print "PAPER INFO"
            print "rE = %g\n" % rE

            fig=pylab.figure(figsize=(8,3))
            #figax=pylab.gca()


            nr, nc = data[t].shape
            extent = [-0.5, nc-0.5, nr-0.5, -0.5]
            kw = {'extent': extent,
                  'origin': 'upper',
                  'interpolation': 'nearest',
                  'aspect': 'equal',
                  'cmap': cmap}

            ax=pylab.subplot(m, n, 1) # Original 
            im=pylab.imshow(data[t], **kw)
            ax.xaxis.set_major_locator(pylab.NullLocator())
            ax.yaxis.set_major_locator(pylab.NullLocator())
            if [t,i] in [ [0,14], [1,9] ]:
                ax.set_title('Original')

            if t == 0:
                pylab.ylabel('Unlensed observation')
            else:
                pylab.ylabel(r'Lensed, $\theta_{E,\mathrm{test}}=%.4f$' % rE)

            ax=pylab.subplot(m, n, 2) # Reconstructed 
            pylab.imshow(f, **kw)
            ax.xaxis.set_major_locator(pylab.NullLocator())
            ax.yaxis.set_major_locator(pylab.NullLocator())
            if [t,i] in [ [0,14], [1,9] ]:
                ax.set_title('Reconstructed')

            ax=pylab.subplot(m, n, 3) # % Error in difference 
            im=pylab.imshow(abs(data[t] - f), **kw)
            ax.xaxis.set_major_locator(pylab.NullLocator())
            ax.yaxis.set_major_locator(pylab.NullLocator())
            bnds = ax.get_position().bounds
            #cax = pylab.axes([0.95, 0.1, 0.04, 0.2])
            #print cax
            print [bnds[0]+bnds[2], bnds[1], .2, bnds[3]]
            #pylab.colorbar(mappable=im, cax=pylab.axes([bnds[0]+bnds[2]+0.01, .20, .02, .6]))
            pylab.colorbar(cax=pylab.axes([bnds[0]+bnds[2]+0.01, .20, .02, .6]))
            if [t,i] in [ [0,14], [1,9] ]:
                ax.set_title('Residual')

            pylab.savefig("recon.%i.%02i.%i.color.eps" % (run_id, i, t), 
                          bbox_inches='tight',
                          pad_inches=0);

            continue

        pylab.figure()
        pylab.matshow(diff)
        pylab.figure()
        pylab.matshow(data[t])
        pylab.colorbar()
        pylab.gray()

def run_sim(data, N0):
    """ Run the lensing simulation with the parameters provided in params.py.
    This consists of four steps:
    (1)  Choose different values for rE.
    (2)  Move the star across the sky and generate "observations". The first is unlensed.
    (3)  Compute the projection and covariance matrices.
    (4)  Computer the effective chi^2.

    data is a 3 dimensional array consisting of n 2D normalized surface brightness 
    distributions. Using shaplets as the basis function, the data is reconstructed
    using test masses of a lensing star. The marginalized likelihood for each
    test mass is returned.
    """

    N = N0**2       # Total number of basis functions

    nepochs   = data.shape[0]
    grid_size = data.shape[1]

    print "nepochs = %i grid_size=%i" % (nepochs, grid_size)
    assert nepochs > 0
    assert grid_size > 0

    #---------------------------------------------------------------------------
    #
    #---------------------------------------------------------------------------

    # The basis functions evaluated at the lensed positions
    L     = empty((N, nepochs, grid_size, grid_size), numpy.float64)
    # Copied, flattened version of the above divided by sigma^2
    Lt    = empty((N, nepochs*grid_size*grid_size), numpy.float64)

    # Projection of data on the model
    P     = mat(empty((N,1), numpy.float64))
    # Covariance matrix
    C_inv = mat(empty((N,N), numpy.float64))
    
    # Output array
    probs = empty((num_samples, 2), numpy.float64)

    #---------------------------------------------------------------------------
    # Basis function constants, and hermite polynomials.
    # Precompute the coefficients and setup the actual form of each hermite
    # function now, since they are constant over the run.
    #---------------------------------------------------------------------------

    #vals = [[n, 1.0/sqrt((2**n) * sqrt(pi) * factorial(n,1) * beta), hermite(n), 0,0] 
    vals = [[n, 1.0/sqrt( (2**n) * sqrt(pi) * factorial(n,1) ), hermite(n), 0,0] 
            for n in xrange(N0)]

    sqrt_data = sqrt(data)
    beta2     = beta**2


    #---------------------------------------------------------------------------
    # Now we start the simulation.
    #
    # (1) Choose different values for rE.
    #---------------------------------------------------------------------------


    for i,rE in enumerate(linspace(rE_sample[0], rE_sample[1], num_samples)):

        # XXX: Just for some specific plots
        #if i not in [14]: continue
        # XXX: Just for some specific plots


        #-----------------------------------------------------------------------
        # (2) Move the star across the sky 
        #-----------------------------------------------------------------------


        for t,z in star_track(nepochs):


            #-------------------------------------------------------------------
            # (2a) Generate an "observation". The first is unlensed.
            #-------------------------------------------------------------------
            if t == 0:
                print "%4i] rE=%f Epoch %i NO LENS" % (i, rE, t)
                xx = raytrace()
                mask = 1
            else:
                print "%4i] rE=%f Epoch %i @ %f,%f" % (i, rE, t, z.real,z.imag)
                xx = raytrace(rE, z)
                mask = star_mask(z)

            #-------------------------------------------------------------------
            # Basis function approximation
            #-------------------------------------------------------------------
            expreal = exp(-xx.real**2/(2*beta2))
            expimag = exp(-xx.imag**2/(2*beta2))
            for n,K,H,_,_ in vals:
                vals[n][3] = K * H(xx.real/beta) * expreal
                vals[n][4] = K * H(xx.imag/beta) * expimag

            n=0
            for _,_,_,b1,_ in vals:
                for _,_,_,_,b2 in vals:
                    L[n,t] = b1 * b2 / beta * mask
                    n += 1


        #-----------------------------------------------------------------------
        # (3) Compute the projection and covariance matrices.
        #-----------------------------------------------------------------------


        for n in xrange(N):
            sum(L[n], out=P[n])
            Lt[n] = (L[n] / sqrt_data).flatten()

        print "Building C_inv"
        C_inv = mat(inner(Lt,Lt))
        print "Done"
        C = C_inv.I

        #-----------------------------------------------------------------------
        # (3a) Optionally calculate the basis function coefficients and plot
        # the reconstruction of the image.
        #-----------------------------------------------------------------------
        if 0: plot_reconstructions(i,data,rE,L,C,P,N)


        #-----------------------------------------------------------------------
        # (4) Computer the effective chi^2. Note that we do not subtract the
        # third term (effectively the gamma_tot) here; we leave that for the 
        # plotting program to do. chi2 here will then be about equal to
        # gamma_tot.
        #-----------------------------------------------------------------------


        log_det = sum(log(eig(C, right=False).real))
        PCP     = P.T * C * P

        chi2 = log_det + PCP

        probs[i] = rE, chi2
        print "rE,chi2 =", rE, chi2

    return probs

################################################################################
################################################################################
################################################################################

if __name__ == "__main__":

    # The identify specifying which parameter set to use is given as the 
    # first argument.
    id = int(sys.argv[1])

    # Do this so that we can have params.py in the current directory, but have
    # ml.py in a different one.
    sys.path.insert(0,'.')

    #---------------------------------------------------------------------------
    # params.py is generated by run_ml.py. It contains all the parameter
    # combinations that will be explored. Each invocation of ml.py receives a
    # different id, which is passed to the get_params() function defined in
    # params.py to get the correct parameters.
    #---------------------------------------------------------------------------
    from params import get_params
    get_params = __import__('params', globals(), locals(), ['get_params'], 0).get_params

    run_id, nepochs, gamma_tot, \
    rE_true, closest_star_approach, rE_sample, num_samples = get_params(id)

    #---------------------------------------------------------------------------
    # Global parameters.
    #---------------------------------------------------------------------------

    beta           = .20                     # arcsec - Basis function normalization
    Nbases         = 20                      # sqrt(Number of basis functions)
    grid_phys      = 2.0                     # arcsec - Physical size across grid
                   
    grid_radius    = 35                      # pixels
    grid_size      = 2*grid_radius + 1       # pixels
                   
    cell_size      = grid_phys / grid_size   # arcsec/pixel
    star_mask_size = 1 * cell_size           # arcsec

    # A fixed seed is used so that all the invocations generate the same test data.
    seed(12)


    if 0:
        r = grid_phys*10
        X = linspace(-r, r, 2*(grid_radius+1)*10+1)
        Y = (deVaucouleurs(X))

        pylab.plot(X, normalize(deVaucouleurs(X)), label="deV")
        pylab.plot(X, normalize(profile_exponential5(X)), label="exp5")
        #pylab.plot(X, normalize(profile_exponential4(X)), label="exp4")
        pylab.plot(X, normalize(profile_exponential5_2(X)), label="exp52")
        pylab.legend()
        pylab.show()
        sys.exit(0)


    if rE_sample[0] == rE_sample[1]:
        print "WARNING: Sample range is one value"
        num_samples = 1

    print "cell_size = %.4f arcsec/pixel" % cell_size


    #---------------------------------------------------------------------------
    # Create the position grid
    # 
    # The grid has a center on (0,0) and the upper left corner is 
    # (-grid_radius, grid_radius).
    #---------------------------------------------------------------------------

    theta = empty((grid_size, grid_size), 'complex')
    for i,j in ndindex(grid_size,grid_size):
        theta[i,j] = complex(j-grid_radius, -(i-grid_radius)) * cell_size

    #save_basis_functions(Nbases)
    #save_bessel_functions(Nbases)

    #---------------------------------------------------------------------------
    # Select data set
    #---------------------------------------------------------------------------


    star_track = star_track_nonsymmetric

    data = model(deVaucouleurs)
    #data = model(profile_exponential5_2)
    #data = model(profile_exponential5)
    #data = model(profile_exponential2)
    #data = model(profile_exponential4)
    #data = model(profile_exponential3)
    #data = model_two_basis_functions()

    #---------------------------------------------------------------------------
    # Add some poissonian noise
    #---------------------------------------------------------------------------

    vaddnoise = vectorize(lambda x: poisson(lam=x))
    data = vaddnoise(normalize(data))
    data[data < 1] = 1

    if 0:
        for d in data:
            for t,z in star_track(nepochs):
                if t > 0: d *= star_mask(z)
            pylab.matshow(d)
            pylab.colorbar()
        pylab.show()
        sys.exit(0)


    if 0:
        d = data[0]
        for t,z in star_track(nepochs):
            if t > 0: d *= star_mask(z)
        pylab.matshow(d)
        #pylab.colorbar()
        pylab.show()
        sys.exit(0)
    #---------------------------------------------------------------------------
    #save_basis_functions(8)
    #demo_lensing(6)
    #sys.exit(0)
    #---------------------------------------------------------------------------


    #---------------------------------------------------------------------------
    # Run the simulation
    #---------------------------------------------------------------------------

    probs = run_sim(data, Nbases)

    sys.exit(0)


    #---------------------------------------------------------------------------
    # Write some output
    #---------------------------------------------------------------------------

    f = open('L%i.txt' % run_id, "w")
    print >>f, "# Generated on: %s" % time.asctime()
    print >>f, "# %19s %21s" % ("rE", "Effective \chi^2 (+ ~\gamma_tot)")
    for x,y in probs:
        print >>f, "%21.15e %21.15e" % (x,y)
    print >>f
    f.close()

    print "%i FINISHED" % run_id

    try:
        pylab.show()
    except:
        pass

