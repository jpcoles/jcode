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
from numpy import empty, zeros, ndindex, vectorize, abs
from numpy import sum, linspace, mat, abs, inner, asarray_chkfinite, ravel
from numpy.random import poisson, seed
from numpy.core.umath import exp, sqrt, log

from scipy.linalg import eig
from scipy.special import hermite
from scipy.misc import factorial

# Use LaTeX to generate text for the plots
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('text', usetex=True)

import pylab

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
    yield 0, complex(0)

    pos = linspace(0, 0.5, n-1)

    for i, p in enumerate(pos):
        yield i+1, complex(p, -closest_star_approach)

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
        else:
            xx = raytrace(rE_true, z)

        data[t] = normalize(m(xx))

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
# A variety of different surface brightness distribution functions.
#===============================================================================

def deVaucouleurs(R, Ie=1, m=4, Re=1.0):
    bm = 2*m - 0.324
    return Ie * exp(-bm * (abs(R)/Re)**(1/m))

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
    pylab.savefig("B%i.png" % N)

# Note: This has been heavily tailored to generate specially labeled plots for
# the paper.
def plot_reconstructions(i, data, rE, L, C, P, N):
    """ Reconstruct the surface brightness distribution from the model and plot."""
    nepochs   = data.shape[0]
    grid_size = data.shape[1]

    fn = C*P
    #for n in xrange(N): print fn[n,0]

    for t,z in star_track(nepochs):
        f = zeros((grid_size, grid_size), 'float')


        # Could probably use tensordot() here.
        for n in xrange(N): f += fn[n,0] * L[n,t]

        if t == 0: f0 = f

        n = 3
        m = nepochs
        m = 1
        cmap = pylab.cm.jet
        #cmap = pylab.cm.gray_r

        diff  = 100 * abs(data[t] - f) / data[t]
        diff1 = 100 * abs(data[t] - f0) / data[t]

        fig=pylab.figure(figsize=(8,3))
        figax=pylab.gca()

        ax=pylab.subplot(m, n, 1) # Original 
        pylab.imshow(data[t], cmap=cmap)
        ax.xaxis.set_major_locator(pylab.NullLocator())
        ax.yaxis.set_major_locator(pylab.NullLocator())
        if (t == 0 and i==14) or (t==1 and i == 9):
            ax.set_title('Original')

        if t == 0:
            pylab.ylabel('Unlensed observation')
        else:
            pylab.ylabel(r'Lensed, $\theta_{E,\mathrm{test}}=%.4f$' % rE)

        ax=pylab.subplot(m, n, 2) # Reconstructed 
        pylab.imshow(f, cmap=cmap)
        ax.xaxis.set_major_locator(pylab.NullLocator())
        ax.yaxis.set_major_locator(pylab.NullLocator())
        if (t == 0 and i==14) or (t==1 and i == 9):
            ax.set_title('Reconstructed')

        ax=pylab.subplot(m, n, 3) # % Error in difference 
        im=pylab.imshow(abs(data[t] - f), cmap=cmap)
        ax.xaxis.set_major_locator(pylab.NullLocator())
        ax.yaxis.set_major_locator(pylab.NullLocator())
        bnds = ax.get_position().bounds
        #cax = pylab.axes([0.95, 0.1, 0.04, 0.2])
        #print cax
        print [bnds[0]+bnds[2], bnds[1], .2, bnds[3]]
        pylab.colorbar(cax=pylab.axes([bnds[0]+bnds[2]+0.01, .20, .02, .6]))
        if (t == 0 and i==14) or (t==1 and i == 9):
            ax.set_title('Residual')

        #pylab.savefig("recon.%i.%02i.%i.png" % (run_id, i, t), bbox_inches='tight');
        pylab.savefig("recon.%i.%02i.%i.color.eps" % (run_id, i, t), 
                      bbox_inches='tight',
                      pad_inches=0);

        #pylab.show()

#       ax=pylab.subplot(nepochs, n, n*t + 4) # % Error in difference 
#       pylab.imshow(diff)
#       ax.xaxis.set_major_locator(pylab.NullLocator())
#       ax.yaxis.set_major_locator(pylab.NullLocator())
#       pylab.colorbar()
#       pylab.gray()

#       ax=pylab.subplot(nepochs, n, n*t + 5) # % Error in difference 
#       pylab.hist(ravel(diff), bins=40, normed=True, histtype='step')
#       ax.yaxis.set_major_locator(pylab.NullLocator())
#       ax.set_xlim(0,100)

#       ax=pylab.subplot(nepochs, n, n*t + 6) # % Error in difference 
#       pylab.imshow(diff1)
#       ax.xaxis.set_major_locator(pylab.NullLocator())
#       ax.yaxis.set_major_locator(pylab.NullLocator())
#       pylab.colorbar()
#       pylab.gray()

#       ax=pylab.subplot(nepochs, n, n*t + 7) # % Error in difference 
#       pylab.hist(ravel(diff1), bins=40, normed=True, histtype='step')
#       ax.yaxis.set_major_locator(pylab.NullLocator())
#       ax.set_xlim(0,100)

    title = r"""\
rE=%.4f (%.4f) sample %i/%i $\gamma_\mathrm{tot}$=%i
Original, Reconstructed, %% Err""" % \
(rE, rE_true, i+1, num_samples, gamma_tot)

    #pylab.suptitle(title)
    #pylab.savefig("recon.%i.%02i.png" % (run_id, i));
    #pylab.savefig("recon.%i.%02i.eps" % (run_id, i));



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
    L     = empty((N, nepochs, grid_size, grid_size))
    # Copied, flattened version of the above divided by sigma^2
    Lt    = empty((N, nepochs*grid_size*grid_size))

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

    vals = [[n, 1.0/sqrt((2**n) * sqrt(pi) * factorial(n,1) * beta), hermite(n), 0,0] 
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
            else:
                print "%4i] rE=%f Epoch %i @ %f,%f" % (i, rE, t, z.real,z.imag)
                xx = raytrace(rE, z)

            #-------------------------------------------------------------------
            # Basis function approximation
            #-------------------------------------------------------------------
            expreal = exp(-xx.real**2/(2*beta2))
            expimag = exp(-xx.imag**2/(2*beta2))
            for n,K,H,v0,v1 in vals:
                vals[n][3] = K * H(xx.real/beta) * expreal
                vals[n][4] = K * H(xx.imag/beta) * expimag

            n=0
            for v1 in vals:
                for v2 in vals:
                    L[n,t] = v1[3] * v2[4]
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
        # Optionally calculate the basis function coefficients and plot the
        # reconstruction of the image.
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
        print rE, chi2

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

    beta        =       .20                     # arcsec - Basis function normalization
    Nbases      =       25                      # sqrt(Number of basis functions)
    grid_phys   =       2.0                     # arcsec - Physical size across grid
    grid_radius =       60                      # pixels
    grid_size   =       2*grid_radius + 1       # pixels
    cell_size   =       grid_phys / grid_size   # arcsec/pixel

    # A fixed seed is used so that all the invocations generate the same test data.
    seed(0)


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

    #---------------------------------------------------------------------------
    # Select data set
    #---------------------------------------------------------------------------

    star_track = star_track_nonsymmetric

    data = model(deVaucouleurs)
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

    #---------------------------------------------------------------------------
    #save_basis_functions(8)
    #demo_lensing(6)
    #sys.exit(0)
    #---------------------------------------------------------------------------


    #---------------------------------------------------------------------------
    # Run the simulation
    #---------------------------------------------------------------------------

    probs = run_sim(data, Nbases)

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

    if False and len(probs) != 0:
        pylab.figure()
        pylab.plot(probs[:,0], gamma_tot - probs[:,1])
        pylab.axvline(x=rE_true)
        pylab.xlabel('rE')
        pylab.ylabel('Likelihood')
        title = "rE(%.4f)/closest approach(%.4f)=%.4f\nnepochs=%i gamma_tot=%ld" % \
            (rE_true, closest_star_approach, rE_true/closest_star_approach, nepochs,gamma_tot)

        pylab.title(title)
        pylab.savefig("L%i.png" % run_id);

