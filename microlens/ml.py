import sys
from math import pow, pi
from numpy import array, empty, zeros, ones, ndindex, vectorize, repeat, zeros_like, abs, argmax, argmin
from numpy import sum, max, linspace, arange, mat, abs, min, inner, ndenumerate, tensordot
from numpy import savetxt, asarray_chkfinite, copy, ravel, any, sort, cumsum
from numpy.random import poisson, seed
from numpy.core.umath import exp, sqrt, log

from scipy.linalg import det, eig
from scipy.special import hermite
from scipy.misc import factorial

import pylab

iter_params = ''

if True:

    #-------------------------------------------------------------------------

    from params import *

    #-------------------------------------------------------------------------

    beta        =       .20                     # arcsec - Basis function normalization
    Nbases      =       28                      # sqrt(Number of basis functions)
    grid_phys   =       2.0                     # arcsec - Physical size across grid
    grid_radius =       60                      # pixels
    grid_size   =       2*grid_radius + 1       # pixels
    cell_size   =       grid_phys / grid_size   # arcsec/pixel
    #grid_pos    =       (0.0, 0.0)              # arcsec

##############################################################################
##############################################################################
##############################################################################

if False:
    nepochs     =       6
    num_samples  =       20
    gamma_tot   =       100000000L              # Total photon count over all observations
    #gamma_tot   =       20000000L              # Total photon count over all observations
    rE_true     =       0.025                   # arcsec
    rE_sample   =       (0.02, 0.03)
    beta        =       .2                      # arcsec - Basis function normalization
    Nbases      =       20

    grid_phys   =       2.0                     # arcsec
    grid_radius =       60                      # pixels
    grid_size   =       2*grid_radius + 1       # pixels
    grid_pos    =       (0.0, 0.0)              # arcsec
    cell_size   =       grid_phys / grid_size   # arcsec/pixel

elif False:

    #-------------------------------------------------------------------------

    run_id      =       0
    nepochs     =       2
    gamma_tot   =       10000000L               # Total photon count over all observations
    num_samples  =       5
    rE_true     =       0.025                   # arcsec
    closest_star_approach = 0.2

    #-------------------------------------------------------------------------

    rE_sample   =       (0.010, 0.04)
    beta        =       .2                      # arcsec - Basis function normalization
    Nbases      =       20
    grid_phys   =       2.0                     # arcsec
    grid_radius =       60                      # pixels
    grid_size   =       2*grid_radius + 1       # pixels
    grid_pos    =       (0.0, 0.0)              # arcsec
    cell_size   =       grid_phys / grid_size   # arcsec/pixel

##############################################################################
##############################################################################
##############################################################################

#seed(0)
seed()

print "cell_size = %.4f arcsec/pixel" % cell_size

##############################################################################

def star_track(n):
#    pos = [complex(-100, -100),
#           complex(-.75, .3)]
#
#    #pos = [complex(-.25, .3)]
#
#    for i,p in enumerate(pos):
#        if i == n: break
#        yield i, p

    #yield 0, complex(-.75, -.5)
    yield 0, complex(0)

    if n > 2: 
        pos = linspace(.5, -.5, n-1)
    elif n == 2:
        pos = [0.0]
    else:
        return

    for i, p in enumerate(pos):
        yield i+1, complex(p, -closest_star_approach)

def model(m):
    #profsb = m(theta)
    #surface_brightness = lambda R: mapSB(profsb, R)
    data = empty((nepochs, grid_size, grid_size), 'float')
    for t,z in star_track(nepochs):
        if t == 0:
            xx = raytrace()
        else:
            xx = raytrace(rE_true, z)

        data[t] = normalize(m(xx))

        #data[t] = normalize(surface_brightness(xx))

    return data

def model_two_basis_functions():

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

##############################################################################

def deVaucouleurs(R, Ie=.40, m=4, Re=.50):
    bm = 2*m - 0.324
    return Ie * exp(-bm * (abs(R)/Re)**.25)

def profile_exponential(R, Ie=.40, Rd=.50):
    return Ie * exp(-abs(R)/Rd)

def profile_exponential2(R):
    return (1+abs(R))**-2

def profile_exponential3(R):
    return (1 + sqrt(R.real**2 + 0.8*R.imag**2))**-2

def raytrace(theta_E=None, z=None):
    """theta_E - Einstein radius of star
       z       - position of star (complex)"""
    if theta_E is None and z is None: return theta
    f = theta - (theta_E**2) * (theta-z) / abs(theta-z)
    return asarray_chkfinite(f)

def mapSB(sb, mp):
    res = zeros_like(sb)
    for (i,j),r in ndenumerate(mp):
        r /= cell_size
        x = int(r.real) + grid_radius
        y = int(r.imag) + grid_radius
        if 0 <= x <= grid_size and 0 <= y <= grid_size:
            res[i,j] = sb[x,y]
    return res

def normalize(sb):
    #global gamma_tot
    #print gamma_tot, sb
    sb *= float(gamma_tot) / sum(sb)
    return sb

##############################################################################

def demo_lensing(N):

    pylab.figure()
    for t,z in star_track(nepochs):
        xx = normalize(surface_brightness(raytrace(rE_true, z)))
        pylab.subplot(N/2, N/2, t+1)
        pylab.imshow(xx.T)

    pylab.suptitle("Lensing demo  rE=%.3f" % rE_true)
    pylab.savefig("LensingDemo.png")

def save_basis_functions(N=Nbases):
    #global theta, beta

    L = N**2

    beta2 = beta**2
    B     = empty((grid_size, grid_size)) # Don't want matrix behaviour here

    #-------------------------------------------------------------------------
    # Basis function constants, and hermite polynomials
    #-------------------------------------------------------------------------

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


def run_sim(data, N0=Nbases):

    N = N0**2

    nepochs   = data.shape[0]
    grid_size = data.shape[1]

    print "nepochs = %i grid_size=%i" % (nepochs, grid_size)
    assert nepochs > 0
    assert grid_size > 0

    #-------------------------------------------------------------------------
    #
    #-------------------------------------------------------------------------

    B     = empty((N, nepochs, grid_size, grid_size)) # Don't want matrix behaviour here
    Bt    = empty((N, nepochs*grid_size*grid_size)) # Don't want matrix behaviour here
    #Bt    = empty((N, nepochs, grid_size, grid_size)) # Don't want matrix behaviour here
    Br    = empty((N, nepochs * grid_size * grid_size)) # Don't want matrix behaviour here
    tmp   = empty((nepochs, grid_size, grid_size)) # Don't want matrix behaviour here
    P     = mat(empty((N,1), 'float'))
    C_inv = mat(empty((N,N), 'float'))
    probs = empty((num_samples, 2), 'float')
    beta2 = beta**2

    #-------------------------------------------------------------------------
    # Basis function constants, and hermite polynomials
    #-------------------------------------------------------------------------

    vals = [[n, 1.0/sqrt((2**n) * sqrt(pi) * factorial(n,1) * beta), hermite(n), 0,0] 
            for n in xrange(N0)]

    sqrt_data = sqrt(data)

    #-------------------------------------------------------------------------
    # Choose different values for rE
    #-------------------------------------------------------------------------

    for i,rE in enumerate(linspace(rE_sample[0], rE_sample[1], num_samples)):

        #---------------------------------------------------------------------
        # Move the star across the sky 
        #---------------------------------------------------------------------

        for t,z in star_track(nepochs):
            #-----------------------------------------------------------------
            # Generate an "observation". The first is unlensed.
            #-----------------------------------------------------------------
            if t == 0:
                print "%4i] rE=%f Epoch %i NO LENS" % (i, rE, t)
                xx = raytrace()
            else:
                print "%4i] rE=%f Epoch %i @ %f,%f" % (i, rE, t, z.real,z.imag)
                xx = raytrace(rE, z)

            #-----------------------------------------------------------------
            # Basis function approximation
            #-----------------------------------------------------------------
            expreal = exp(-xx.real**2/(2*beta2))
            expimag = exp(-xx.imag**2/(2*beta2))
            for n,K,H,v0,v1 in vals:
                vals[n][3] = K * H(xx.real/beta) * expreal
                vals[n][4] = K * H(xx.imag/beta) * expimag

            n=0
            for v1 in vals:
                for v2 in vals:
                    B[n,t] = v1[3] * v2[4]
                    n += 1

        #---------------------------------------------------------------------
        # 
        #---------------------------------------------------------------------

        for n in xrange(N):
            sum(B[n], out=P[n])
            Bt[n] = ravel(B[n] / sqrt_data)

        print "Building C_inv"

        C_inv = mat(tensordot(Bt,Bt,axes=(-1,-1)))

#       for k in xrange(N):
#           for n in xrange(k, N):
#               #print Bt[k]
#               #C_inv[n,k] = C_inv[k,n] = inner(ravel(Bt[k]), ravel(Bt[n]))
#               C_inv[n,k] = C_inv[k,n] = inner(Bt[k], Bt[n])

        print "Done"

        #---------------------------------------------------------------------
        # Invert the matrix
        #---------------------------------------------------------------------

        C = C_inv.I

        #---------------------------------------------------------------------
        # Optionally calculate the basis function coefficients and plot the
        # reconstruction of the image.
        #---------------------------------------------------------------------

        if True:
            fn = C*P
            pylab.figure()
            for t,z in star_track(nepochs):
                f = zeros((grid_size, grid_size), 'float')

                #for n in xrange(N): print fn[n,0]

                for n in xrange(N): f += fn[n,0] * B[n,t]

                diff = 100 * abs(data[t] - f) / data[t]

                ax=pylab.subplot(nepochs, 4, 4*t + 1) # Original 
                pylab.imshow(data[t])
                ax.xaxis.set_major_locator(pylab.NullLocator())
                ax.yaxis.set_major_locator(pylab.NullLocator())
                pylab.jet()

                ax=pylab.subplot(nepochs, 4, 4*t + 2) # Reconstructed 
                pylab.imshow(f)
                ax.xaxis.set_major_locator(pylab.NullLocator())
                ax.yaxis.set_major_locator(pylab.NullLocator())
                pylab.jet()

                ax=pylab.subplot(nepochs, 4, 4*t + 3) # % Error in difference 
                pylab.imshow(diff)
                ax.xaxis.set_major_locator(pylab.NullLocator())
                ax.yaxis.set_major_locator(pylab.NullLocator())
                pylab.colorbar()
                pylab.gray()

                ax=pylab.subplot(nepochs, 4, 4*t + 4) # % Error in difference 
                pylab.hist(ravel(diff), bins=40, normed=True, histtype='step')
                ax.yaxis.set_major_locator(pylab.NullLocator())
                ax.set_xlim(0,100)

            title = """\
rE=%.4f (%.4f) sample %i/%i gamma_tot=%i
Original, Reconstructed, %% Err""" % \
(rE, rE_true, i+1, num_samples, gamma_tot)

            pylab.suptitle(title)
            pylab.savefig("recon.%s.%i.png" % (run_id, i));

        #---------------------------------------------------------------------
        # 
        #---------------------------------------------------------------------


        #prob = (1/2.)*log(abs(det(C))) + P.T * C * P
#       print "1"
        LN   = (1/2.)*log(eig(C, right=False).real)
#       print L
#       print "2"
#       LF  = L.flatten()
#       print "3"
#       LFS = sort(LF)
#       print "4"
#       LCS = cumsum(LFS)
#       print "5"

        #print LCS
        #pylab.figure()
        #pylab.plot(LCS)
        #pylab.show()

        #print "6"
        LS  = sum(LN)
        #print "7"
        PCP = P.T * C * P

        prob = LS + PCP
        #prob = PCP
        #prob = (1/2.)*sum(log(eig(C, right=False).real)) + P.T * C * P

        probs[i] = rE, prob
        print rE, prob

    return probs

##############################################################################
##############################################################################
##############################################################################

if __name__ == "__main__":

    #-------------------------------------------------------------------------
    # Create the position grid
    # 
    # The grid has a center on (0,0) and the upper left corner is 
    # (-grid_radius, grid_radius).
    #-------------------------------------------------------------------------

    theta = empty((grid_size, grid_size), 'complex')
    for i,j in ndindex(grid_size,grid_size):
        theta[i,j] = complex(j-grid_radius, -(i-grid_radius)) * cell_size

    #-------------------------------------------------------------------------
    # Select data set
    #-------------------------------------------------------------------------

    data = model(profile_exponential3)
    #data = model_two_basis_functions()

    #-------------------------------------------------------------------------
    # Add some noise
    # This causes a lot of disagreement in the reconstructed images.
    #-------------------------------------------------------------------------

    vaddnoise = vectorize(lambda x: x)
    #vaddnoise = vectorize(lambda x: poisson(lam=x))
    data = vaddnoise(normalize(data))
    data[data < 1] = 1

    #-------------------------------------------------------------------------
    #save_basis_functions(8)
    #demo_lensing(6)
    #sys.exit(0)
    #-------------------------------------------------------------------------


    #-------------------------------------------------------------------------
    # Run the simulation
    #-------------------------------------------------------------------------

    probs = run_sim(data)

    f = open('L%s.txt' % run_id, "w")
    print >>f, '%s' % iter_params
    print >>f, "# %19s %21s" % ("rE", "Probability")
    for x,y in probs:
        print >>f, "%21.15e %21.15e" % (x,y)
    print >>f
    Mi = argmax(probs[:,1])
    mi = argmin(probs[:,1])
    print >>f, "%21.15e %21.15e" % (probs[0,0], probs[0,1])
    print >>f, "%21.15e %21.15e" % (probs[-1,0], probs[-1,1])
    print >>f, "%21.15e %21.15e" % (probs[mi,0], probs[mi,1])
    print >>f, "%21.15e %21.15e" % (probs[Mi,0], probs[Mi,1])
    f.close()

    #-------------------------------------------------------------------------
    # Write some output
    #-------------------------------------------------------------------------

    #data = mapSB(galaxySB, beta)

    # Normalize so that the maximum liklihood is 1
    #probs[:,1] -= max(probs[:,1])
    #probs[:,1]  = exp(probs[:,1])

    #print >>f, "nepochs=%i" % nepochs
    #print >>f, "gamma_tot=%ld" % gamma_tot
    #print >>f, "rE_true=%.4f" % rE_true
    #print >>f, "closest_star_approach=%.4f" % closest_star_approach
    #print >>f, "rE_sample=(%.4f, %.4f)" % (rE_sample[0], rE_sample[1])
    #print >>f, "num_samples=%i" % num_samples
    #print >>f, "\n"

    if True and len(probs) != 0:
        pylab.figure()
        pylab.plot(probs[:,0], gamma_tot - probs[:,1])
        pylab.axvline(x=rE_true)
        pylab.xlabel('rE')
        pylab.ylabel('Liklihood')
        title = "rE(%.4f)/closest approach(%.4f)=%.4f\nnepochs=%i gamma_tot=%ld" % \
            (rE_true, closest_star_approach, rE_true/closest_star_approach, nepochs,gamma_tot)

        pylab.title(title)
        pylab.savefig("L%s.png" % run_id);


    if False:
        pylab.figure()
        pylab.imshow(noise)
        pylab.title('Noise')

    if False:
        pylab.figure()
        pylab.imshow(galaxySB)
        pylab.title('Galaxy Surface Brightness')

    if False:
        pylab.figure()
        pylab.imshow(data)
        pylab.title('data')

    if False:
        pylab.figure()
        pylab.plot(probs)

    if False:
        xi = linspace(-grid_phys/2,grid_phys/2, grid_size*10)
        yi = linspace(-grid_phys/2,grid_phys/2, grid_size*10)
        zi = pylab.griddata(xs,ys,probs, xi,yi)

        pylab.figure()
        pylab.contour(xi,yi,zi)
        pylab.colorbar()
        pylab.scatter(xs, ys, marker='o',c='b',s=5)
        pylab.xlim(-grid_phys/2, grid_phys/2)
        pylab.ylim(-grid_phys/2, grid_phys/2)



##############################################################################
##############################################################################

    if False and len(probs) == 0:

        N = 20
        L = N**2

        img_beta = raytrace(rE_true, complex(.25, .25))
        #img_beta = theta

        P = mat(empty((L,1), 'float'))
        B = empty((L, grid_size, grid_size)) # Don't want matrix behaviour here
        beta = .2
        beta2 = beta**2

        vals = []
        for n in xrange(N):
            v = []
            v.append(1.0/sqrt(2**n * sqrt(pi) * factorial(n,1) * beta))
            v.append(hermite(n))
            vals.append(v)

        l=0
        for n1 in xrange(N):
            K1,H1 = vals[n1]
            for n2 in xrange(N):
                K2,H2 = vals[n2]
                x = img_beta
                B[l] = (K1 * H1(x.real/beta) * exp(-x.real**2/(2*beta2))) * \
                       (K2 * H2(x.imag/beta) * exp(-x.imag**2/(2*beta2)))
                P[l]  = sum(B[l])
                l += 1

        print "Here", data.shape, data[0].shape
        C_inv = mat(empty((L,L), 'float'))
        for k in xrange(L):
            for l in xrange(k, L):
                C_inv[l,k] = C_inv[k,l] = sum(B[k] * B[l] / data[0])

        print "There"

        C = C_inv.I
        print C.shape, P.shape
        fn = C*P
        print fn.shape, B.shape

        print C

        f = zeros((grid_size, grid_size), 'float')
        for l in xrange(L):
            f += fn[l,0] * B[l]


    #   for l in xrange(9):
    #       pylab.figure()
    #       pylab.imshow(B[l])
    #       pylab.title(str(l))


        print f.shape
        diff = galaxySB - f

        pylab.figure()
        pylab.imshow(galaxySB)
        pylab.title("Original")
        pylab.colorbar()
        pylab.figure()
        pylab.imshow(f)
        pylab.title("Reconstruction")
        pylab.colorbar()
        pylab.figure()
        pylab.imshow(diff)
        pylab.title("Difference")
        pylab.colorbar()

    #pylab.show()

    sys.exit(0)


##############################################################################
##############################################################################

