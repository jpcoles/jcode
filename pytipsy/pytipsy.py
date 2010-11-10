#
# pytipsy.py
#
# Read a TIPSY file into numpy arrays. Support record array views of the
# data and the ability to memmap the file for large data sets.
#
# Written by Jonathan Coles
# Institute for Theoretical Physics, University of Zurich
#
# This code is public domain.
#

import numpy as np

def load_tipsy(filename, fmt='auto', memmap=True, which='gds', merge=None):

    hdr_type  = np.dtype(
                [('time',    'f8'),      # Time
                 ('nBodies', 'u4'),      # Number of particles
                 ('nDim',    'u4'),      # Number of dimensions
                 ('nSph',    'u4'),      # Number of SPH particles
                 ('nDark',   'u4'),      # Number of dark matter particles
                 ('nStar',   'u4'),      # Number of star particles
                 ('',        'u4')])     # <padding>
    
    dark_type = np.dtype(
                [('m',       'f4'),      # Mass
                 ('r',       'f4', 3),   # Position
                 ('v',       'f4', 3),   # Velocity
                 ('eps',     'f4'),      # Softening
                 ('phi',     'f4')])     # Potential

    star_type = np.dtype(
                [('m',       'f4'),      # Mass
                 ('r',       'f4', 3),   # Position
                 ('v',       'f4', 3),   # Velocity
                 ('metals',  'f4'),      # Metals
                 ('tform',   'f4'),      # Formation time
                 ('eps',     'f4'),      # Softening
                 ('phi',     'f4')])     # Potential

    sph_type =  np.dtype(
                [('m',       'f4'),      # Mass
                 ('r',       'f4', 3),   # Position
                 ('v',       'f4', 3),   # Velocity
                 ('rho',     'f4'),      # Density
                 ('temp',    'f4'),      # Temperature
                 ('hsmooth', 'f4'),      # 
                 ('metals',  'f4'),      # Metals
                 ('phi',     'f4')])     # Potential

    def newarray_file(file, type, offs, count):
        file.seek(offs)
        return np.recarray(np.fromfile(file=file, dtype='uint8', count=count*type.itemsize), 
                           offset=offs, shape=count, dtype=type)

        #return np.recarray(buf=file, offset=offs, shape=count, dtype=type)
        #return np.rec.array(np.fromfile(file=file, dtype=type, count=count), copy=False)

    def newarray_mem(buf, type, offs, count):
        return np.recarray(buf=buf, offset=offs, shape=count, dtype=type)

    def newbyteorder(t):
        return [hdr_type.newbyteorder(t),
                sph_type.newbyteorder(t),
                dark_type.newbyteorder(t),
                star_type.newbyteorder(t)]

    if memmap:
        newarray = newarray_mem
        f = np.memmap(filename, mode='c')
    else:
        newarray = newarray_file
        f = open(filename, 'rb')

    if fmt == 'standard':      hdr_type,sph_type,dark_type,star_type = newbyteorder('B')
    if fmt == 'little-endian': hdr_type,sph_type,dark_type,star_type = newbyteorder('L')
    if fmt == 'native':        hdr_type,sph_type,dark_type,star_type = newbyteorder('N')

    hdr = newarray(f, hdr_type, 0, 1)[0]

    if fmt == 'auto':
        if hdr.ndim not in [2,3]:
            hdr_type,sph_type,dark_type,star_type = newbyteorder('S')
            hdr = newarray(f, hdr_type, 0, 1)[0]

    if hdr.nDim not in [2,3]:
        if memmap:
            del f
        else:
            f.close()
        raise IOError("Corrupt TIPSY file? Dimensions are not 2 or 3. Perhaps try fmt='standard' or 'native'")

    if hdr.nBodies != hdr.nSph + hdr.nDark + hdr.nStar:
        raise IOError("Corrupt TIPSY file? Particle counts are inconsistent with total")

    sph_bytes  = hdr.nSph  * sph_type.itemsize
    dark_bytes = hdr.nDark * dark_type.itemsize
    star_bytes = hdr.nStar * star_type.itemsize

    nSph  = hdr.nSph  if 'g' in which else 0
    nDark = hdr.nDark if 'd' in which else 0
    nStar = hdr.nStar if 's' in which else 0

    offs  = hdr_type.itemsize; sph   = newarray(f, sph_type,  offs, nSph)
    offs += sph_bytes;         dark  = newarray(f, dark_type, offs, nDark)
    offs += dark_bytes;        star  = newarray(f, star_type, offs, nStar)

    if not memmap: 
        f.close()
        f = None

    class T: pass
    t = T()
    t.hdr  = t.h = hdr
    t.sph  = t.g = sph
    t.dark = t.d = dark
    t.star = t.s = star
    t.file = t.f = f

    if merge is not None:
        lst = []
        if nSph:  lst.append(t.g[merge])
        if nDark: lst.append(t.d[merge])
        if nStar: lst.append(t.s[merge])
        t.p = np.hstack(lst).view(np.recarray)

    return t
    #hdr, sph, dark, star, f


def com(t):
    def findmin(p):
        mi = np.argmin(p.phi)
        return [p, mi, p[mi].phi]

    l = []
    if t.g.size: l.append(findmin(t.g))
    if t.d.size: l.append(findmin(t.d))
    if t.s.size: l.append(findmin(t.s))

    p, pot_min_i, phi = l[ np.argmin(map(lambda x: x[2], l)) ]

    rc = p[pot_min_i].r

    t.g.r -= rc
    t.d.r -= rc
    t.s.r -= rc

    print rc, pot_min_i, phi



if __name__ == '__main__':

    from pylab import figure, plot, show

    #h,g,d,s,f = load_tipsy('galaxy1.std', memmap=True, fmt='standard')
    #t = load_tipsy('galaxy1.std', memmap=False, fmt='auto')
    #h,g,d,s,f = load_tipsy('galaxy1.std', memmap=True, fmt='little-endian')
    #h,g,d,s = load_tipsy('galaxy1.std')

    #t = load_tipsy('galaxy1.std', memmap=False)
    t = load_tipsy('galaxy1.std', memmap=True, which='ds', merge=['m','r', 'tform'])
    #h,g,d,s,f = memmap_tipsy('galaxy1.std')
    #t = load_tipsy('/Users/jonathan/TIP2VTK/tip2vtk/dark.std', memmap=False, which='')

    #com(t)
    print t.p.size, t.hdr.nBodies, t.p.dtype

    #print 'Radial distances (DM)'
    #print np.sum(d.r**2, 1)
    #print d.r[0:10,0]**2
    #print np.sum(d.r**2, 1)
    #print d.m
    #print g.m
    #print s.m

    #print d.r[0,0], d.r[0,1]
    #print d.r[:,0]
    #print d.r[:,0]

    #print len(d[::].r[0])
    #plot(d.r[:,0], d.r[:,1], 'k,', alpha=0.3)
    #show()
        
