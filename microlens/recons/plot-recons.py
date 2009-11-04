from __future__ import division, with_statement
from numpy import array, where, argwhere, sqrt, sum, logical_and
from matplotlib import rc
rc('font',**{'family':'serif',
             'serif':['Computer Modern Roman'],
             'size':28})
rc('text', usetex=True)
rc('xtick.major', **{'pad' : 10})
rc('ytick.major', **{'pad' : 10})
#rc('ytick', **{'major.pad' : 8})
              
from pylab import scatter, show, figure, title, xlabel, ylabel, subplot,\
                  subplots_adjust, ylim, xlim, savefig, gca


fmt = [['CNS NAME',           15,    5, None],
       ['# comp',              1,   21, None],
       ['LHS',                 4,   25, None],
       ['RA',                 10,   32, None],
       ['DEC',                 9,   43, None],
       ['position ref',        1,   53, None],
       ['proper motion "/yr',  6,   56, float],
       ['proper motion angle', 5,   63, float],
       ['pm ref',              1,   69, None],
       ['parallax',            7,   73, float],
       ['parallax err',        7,   81, float],
       ['Spectral',            4,   90, None],
       ['Spectral ref',        1,  104, None],
       ['V',                   6,  107, None],
       ['Mv',                  5,  116, None],
       ['mass',                4,  124, lambda x: float(x or 0)],
       ['notes',              17,  132, None],
       ['Common Name',        36,  152, None]]


if __name__ == "__main__":

    reconsSOA = []
    reconsAOS = {}
    for col,length,offs,func in fmt:
        reconsAOS[col] = []

    with open('RECONS-100-nearest-systems', 'r') as f:
    #with open('x', 'r') as f:
        for lineno,line in enumerate(f):
            line = line.rstrip()
            #print line
            if not 9 <= lineno <= 257: continue
            if len(line) == 0: continue
        
            l = {}
            for col,length,offs,func in fmt:
                val = line[offs:offs+length].strip()
                #print col, val
                #if val and func: val = func(val)
                if func: val = func(val)
                reconsAOS[col].append(val)
                l[col] = val
            reconsSOA.append(l)

    mass = array(reconsAOS['mass'])                 # Msun
    dist = 1/array(reconsAOS['parallax'])           # pc
    rE   = 90*sqrt(mass/dist)                       # mas
    ppm  = array(reconsAOS['proper motion "/yr'])   # arsec/yr

    w = mass > 0
    #w = logical_and(w, mass < 0.08)

    # PAPER PARAMS
    w0 = logical_and(w,
         logical_and(ppm > 0.5, 
         logical_and(rE > 10, mass < .5)))

    clr = [ [ (1.,0.,0.,1.), (0.,0.,0.,1.) ][x] for x in w0 ]

    print len(rE), len(rE[w])

    figure()
    scatter(dist[w], ppm[w], s=80,marker='+', edgecolors=clr)
    ylabel(r'Proper motion [arcsec/yr]', labelpad=20)
    xlabel(r'Distance [parsecs]', labelpad=20)
    ylim(ymin=0)
    xlim(xmin=0)
    subplots_adjust(bottom=0.20, top=0.95, left=0.20)
    savefig('recons-ppm.eps')

    figure()
    scatter(dist[w], mass[w], s=80,marker='+', edgecolors=clr)
    ylabel(r'Estimated mass [$M_\odot$]', labelpad=20)
    xlabel(r'Distance [parsecs]', labelpad=20)
    ylim(ymin=0)
    xlim(xmin=0)
    subplots_adjust(bottom=0.20, top=0.95, left=0.20)
    savefig('recons-mass.eps')

#   subplot(223)
#   scatter(dist[w], rE[w], s=80,marker='+')
#   ylabel(r'$\theta_E$ [milliarcsec]')
#   xlabel(r'Distance [parsecs]')
#   ylim(ymin=0)

    figure()
    scatter(dist[w], 2*rE[w]/1000*ppm[w], s=80,marker='+', edgecolors=clr)
    ylabel(r'$2\theta_E \times$ proper motion [arcsec$^2$/yr]', labelpad=20)
    xlabel(r'Distance [parsecs]', labelpad=20)
    ylim(ymin=0)
    xlim(xmin=0)
    subplots_adjust(bottom=0.20, top=0.95, left=0.20)
    savefig('recons-area.eps')

    figure()
    scatter(mass[w], 2*rE[w]/1000*ppm[w], s=80,marker='+', edgecolors=clr)
    ylabel(r'$2\theta_E \times$ proper motion [arcsec$^2$/yr]', labelpad=20)
    xlabel(r'Estimated mass [$M_\odot$]', labelpad=20)
    ylim(ymin=0)
    xlim(xmin=0)
    subplots_adjust(bottom=0.20, top=0.95, left=0.20)
    savefig('recons-area-v-mass.eps')


    print len(rE[w0])
    print sum(2*rE[w0]/1000*10*ppm[w0])

    show()
