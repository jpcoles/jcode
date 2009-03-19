#!/usr/bin/python
#
# @file
# @brief Dump particles in a Tipsy file.
# @author Doug Potter
#

import sys, os
import re as Regex
from getopt import gnu_getopt
from pytipsy import ifTipsy, ofTipsy, TipsyHeader, TipsyDarkParticle, TipsyGasParticle, TipsyStarParticle, tipsypos

def Display0(o, split=False, precision=5):
    
    fmt = "%%+.%ie " % precision

    if isinstance(o, TipsyHeader):
        red = 1.0 / o.h_time - 1
        print "TIME: %f (redshift %f)" % (o.h_time, red)
        print "#: %i total (%i dimensions), %i dark, %i gas, %i star." % \
            (o.h_nBodies, o.h_nDims, o.h_nDark, o.h_nGas, o.h_nStar)
        print "     x            y             z     ",
        print "     Vx           Vy            Vz     ",
        if bSplit: print
        print "   Mass         phi        eps/Smooth ",
        print "   Metals        Temp           Rho"
    elif isinstance(o, TipsyDarkParticle) or isinstance(o, TipsyStarParticle):
        sys.stdout.write(fmt*6 % (o.pos[0], o.pos[1], o.pos[2],o.vel[0], o.vel[1], o.vel[2]))
        if split: print
        print fmt*3 % (o.mass, o.phi, o.eps)
    elif isinstance(o, TipsyGasParticle):
        sys.stdout.write(fmt*6 % (o.pos[0], o.pos[1], o.pos[2],o.vel[0], o.vel[1], o.vel[2]))
        if split: print
        print fmt*6 % (o.mass, o.phi, o.eps, o.metals, o.temp, o.rho)

def Usage(name):
    print "Usage:", name, "[-h1] <tipsy> [g|d|s]#[-[+]#] ..." 
    print "  -h,--help        Show help (this text)"
    print "  --standard,--std Tipsy standard file (default)"
    print "  --native,--nat   Tipsy native file"
    print "  -p,--precision=5 Set output precision"
    print "  -1,--single      Display each particle on a single line",
    print 
    print "  --markfile=file  Dump only particles from the mark file",
    print 
    print "  --tipsyout=file  Write a Tipsy file instead",
    print
    print
    print "  You must specify at least one ",
    print "particle to dump.  Particle numbers"
    print "  start at zero.  Specify 'g' for gas, 'd' for dark,",
    print " or 's' for star"
    print "  particles.  If you omit g, d or s, the ",
    print "absolute particle number is"
    print "  is used.  For example, "
    print
    print "    tdump in.std 0          ",
    print "Dump the first particle (particle 0)"
    print "    tdump in.std d0         ",
    print "Dump the first dark particle (particle 0)"
    print "    tdump in.std g0-10      ",
    print "Dump the first 11 gas particles"
    print "    tdump in.std d0 g0 s0   ",
    print "Dump the first dark, gas and star"
    print "    tdump in.std d100-101   ",
    print "Dump dark particles 100 and 101"
    print "    tdump in.std d100000-+10",
    print "Dump dark particles 100,000 to 100,010"

if __name__ == '__main__':

    bHelp  = False    #!< Help was requested on the command line
    bError = False    #!< An error occurred on the command line
    bSplit = True     #!< Split lines
    iPrecision=5
    ftype="standard"
    filename=None
    markname=None
    tipsyname=None

    inp = ifTipsy()
    out = ofTipsy()
    worklist = []


    args = []
    rest = []

    long_options = \
       ["help",
        "single",
        "standard",
        "std",
        "native",
        "nat",
        "markfile=",
        "tipsyout=",
        "precision="]

    args,rest = gnu_getopt( sys.argv[1:], "h1p:", long_options)
    for option, value in args:
        if option in ['-h', '--help']:
            bHelp = True
        elif option in ['-1', '--single']:
            bSplit = False
        elif option in ['--standard','--std']:
            ftype = "standard"
        elif option in ['--native','--nat']:
            ftype = "native"
        elif option == '--markfile':
            assert(len(value) != 0)
            markname = value
        elif option == '--tipsyout':
            assert len(value) != 0
            tipsyname = value
        elif option in ['-p', '--precision']:
            assert len(value) != 0
            iPrecision = int(value)
        else:
            bError = True

    if bHelp or bError:
	    Usage(sys.argv[0])
	    sys.exit(1)

    if len(rest) != 0:
        parts = rest[0].split(':', 1)
        if len(parts) == 1: 
            filename = parts[0]
        else:
            ftype, filename = parts
    else:
        print >>sys.stderr, "Missing tipsy input file"
        print >>sys.stderr, "Try", sys.argv[0], "--help"
        sys.exit(2)

    inp.open(filename, ftype)
    if not inp.is_open():
        print >>sys.stderr, "Unable to open Tipsy binary", filename
        sys.exit(2)

    # Read the header
    h = inp.read()

    if h.h_nBodies != h.h_nGas + h.h_nDark + h.h_nStar or not (1 <= h.h_nDims <= 3):
        print >>sys.stderr, filename, "has crazy dimensions"
        print >>sys.stderr, "(have you tried --std or --nat?)"
        sys.exit(2)


    re = Regex.compile("^([gsd])?(\d+)(-(\\+?)(\d+))?$")

    for arg in rest[1:]:
        m = re.search(arg)
        if m is not None:
            grps = m.groups()

            if grps[0] is None: 
                p = tipsypos.particle
            else:
                p = [tipsypos.gas, tipsypos.dark, tipsypos.star]['gds'.index(grps[0])]

            b = int(grps[1])
            if grps[4] is None: 
                e = b
            else:
	            e = int(grps[4])
	            if grps[3] is not None: e += b

            if p == tipsypos.gas:
                bError = b>=h.h_nGas or e>=h.h_nGas
            elif p == tipsypos.dark:
                bError = b>=h.h_nDark or e>=h.h_nDark
            elif p == tipsypos.star:
                bError = b>=h.h_nStar or e>=h.h_nStar
            elif p == tipsypos.particle:
                bError = b>=h.h_nBodies or e>=h.h_nBodies
            else:
                assert(False) # Invalid file section

            if bError:
                print >>sys.stderr, "Range", arg, "is outside the file"
                sys.exit(2)

            worklist.append([p,b,e])
        else:
            print >>sys.stderr, "Invalid particle selection:", arg
            print >>sys.stderr, "Try", sys.argv[0], "--help"
            sys.exit(2)

    if markname is not None:
        try:
            mf = open(markname, 'r')
            nd, ng, ns = mf.readline().split()
            if nd != h.h_nDark or ng != h.h_nGas or ns != h.h_nStar:
                print >>sys.stderr, markname, "does not match the tipsy file"
                sys.exit(2)

            for line in mf:
                i = int(line)
                worklist.append([tipsypos.particle,i-1,i-1])

        except IOError:
            print >>sys.stderr, "Unable to open", markname
            sys.exit(2)

    if len(worklist) == 0:
        worklist.append([tipsypos.particle,0,h.h_nBodies-1])

    if tipsyname is not None:
        try:
            statinfo = os.stat(tipsyname)
            print >>sys.stderr, "Output file exists, delete it first"
            sys.exit(3)
        except OSError:
            pass

        out.open(tipsyname)
        if not out.is_open():
            print >>sys.stderr, "Unable to create", tipsyname
            sys.exit(2)

    if out.is_open():
        Display = out.write
    else:
        Display = lambda x: Display0(x, bSplit, iPrecision)

    Display(h)

    h.h_nGas = h.h_nDark = h.h_nStar = h.h_nBodies = 0

    for section,begin,end in worklist:
        inp.seekg(tipsypos(section, begin))

        for o in xrange(begin, end+1):
            sect = inp.tellg().section()
            p = inp.read()
            if sect == tipsypos.gas:
                h.h_nGas += 1
            elif sect == tipsypos.dark:
                h.h_nDark += 1
            elif sect == tipsypos.star:
                h.h_nStar += 1

            Display(p)

    if out.is_open():
        # Rewrite the header
        h.h_nBodies = h.h_nGas + h.h_nDark + h.h_nStar
        out.seekp(tipsypos(tipsypos.header, 0))
        out.write(h)
        out.close()

    inp.close()

