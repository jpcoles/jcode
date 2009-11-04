#
# plot_Ls.py - Plot parameter likelihoods from ml.py
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
import sys
from collections import defaultdict
from numpy import array, amax, log2, sort, loadtxt, argmin, argmax, sqrt, pi, abs

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('text', usetex=True)
rc('font', size=12)
from matplotlib.pyplot import subplots_adjust
from matplotlib.cm import get_cmap

import pylab
from pylab import scatter, plot, hist, savefig, figure, subplot, colorbar, legend
from pylab import show, xlabel, ylabel, gray, xticks, yticks, xlim, ylim, axvline, axhline

sys.path.insert(0, '.')
from params import get_params, get_all_params

ps = get_all_params()

files = sys.argv[1:]

def mass(r):
    G = 6.67300e-11
    c = 299792458.0
    pc = 3.08568025e16
    return r**2 * pc * c**2 / (0.09**2 * 4*G)

def my_label_outer(ax):
    """
    set the visible property on ticklabels so xticklabels are
    visible only if the subplot is in the last row and yticklabels
    are visible only if the subplot is in the first column
    """
    lastrow = ax.is_last_row()
    lastcol = ax.is_last_col()
    for label in ax.get_xticklabels():
        label.set_visible(lastrow)

    for label in ax.get_yticklabels():
        label.set_visible(lastcol)

    if lastcol:
        ax.yaxis.set_ticks_position('right')
        ax.yaxis.set_label_position('right')

def plot7(data2):

    subn = len(ps['nepochs'])
    subm = len(ps['gamma_tot'])

    err = []

    gray()
    cm_gray = get_cmap()

    all_good  = [[], [], [], [], []]
    all_bad   = [[], [], [], [], []]
    all_ugly  = [[], [], [], [], []]

    rEs = defaultdict(list)

    cell_size = 2./71
    ys = ["%.3f"%abs(complex(2*cell_size,-p)) for p in ps['ps']]
    print p

    i = 0
    for epoch in ps['nepochs'][::-1]:
        for g in ps['gamma_tot']:
            i += 1
            #if epoch is None or g is None: continue
            data = data2[epoch][g]
            #if not data: continue

            p = subplot(subn, subm, i)
            my_label_outer(p)

#           if p.is_first_col() and p.is_last_row():
#               xlabel(r'$\theta_{E,\mathrm{true}}$ [arcsec]')
#               ylabel(r'$p$ [arcsec]')

            if p.is_first_row():
                p.xaxis.set_label_position('top')
                xlabel(r'%i' % g)

            if p.is_first_col():
                p.yaxis.set_label_position('left')
                ylabel(r'%i' % epoch)


            pylab.setp(p.get_xticklabels(), rotation=45)
            #pylab.setp(p.get_yticklabels(), rotation=45)
            good  = [[], [], [], [], []]
            bad   = [[], [], [], [], []]
            ugly  = [[], [], [], [], []]
            for rE_true,closest_star_approach,vals,_,_,_,M in data:

                #rE_err = 100 * abs(rE_true - vals[M][0]) / rE_true
                rE_err = 100 * abs(mass(rE_true) - mass(vals[M][0])) / mass(rE_true)

                if rE_err > 100: rE_err = 100
                x = rE_true
                y = closest_star_approach

                rEs[rE_true].append(rE_err)

                R = lambda x: 15*sqrt(x/(2*pi))
                maxR = R(100)

                if 0 < M < len(vals)-1:
                    print rE_err
                    if rE_err <= 20:
                        good[0].append(x)
                        good[1].append(y)
                        r = maxR - R(rE_err)
                        #good[2].append(r)
                        #good[2].append(int(rE_err*5))
                        if rE_err < 5:
                            good[3].append((0.,0.8,0.,1.))
                            good[2].append(90)
                        elif rE_err < 11:
                            good[3].append((1.,0.5,0.,1.))
                            good[2].append(50)
                        elif rE_err < 20:
                            good[3].append((1.,0.,0.,1.))
                            good[2].append(30)
                        good[4].append(rE_err)
                    elif rE_err:
                        bad[0].append(x)
                        bad[1].append(y)
                        r = maxR - R(rE_err)
                        bad[2].append(r)
                        #bad[2].append(int(rE_err*5))
                        bad[3].append((1.,1.,1.,1.))
                        bad[4].append(rE_err)
                    err.append(rE_err)
                    z = 1
                else:
                    ugly[0].append(x)
                    ugly[1].append(y)
                    ugly[2].append(int(360))
                    ugly[3].append((1.,1.,1.,1.))
                    ugly[4].append(100)
                    z = -1

            #if not good[0] and not bad[0]: continue

            if 0:
                if good[0]:
                    scatter(array(good[1])/array(good[0]), good[4], marker='o')
                if bad[0]:
                    scatter(array(bad[1])/array(bad[0]), bad[4], marker='x')
                ylim(0,50)
                xlim(0,60)
            else:

                if good[0]:
                    scatter(good[0], log2(good[1]), s=good[2], c=good[3], marker='o', zorder=1)
                if bad[0]:
                    scatter(bad[0], log2(bad[1]), s=bad[2], c=bad[3], marker='o', zorder=0, facecolors='none')

                #yticks(log2(ps['ps']), ps['ps'])
                yticks(log2(ps['ps']), ys)
                xticks(ps['rE_true'])

                xlim(0.00, 0.06)
                ylim(-5, 1)

            a,b,c,d,e = good
            all_good[0]+=a; all_good[1]+=b; all_good[2]+=c; all_good[3]+=d; all_good[4]+=e;
            a,b,c,d,e = bad
            all_bad[0]+=a; all_bad[1]+=b; all_bad[2]+=c; all_bad[3]+=d; all_bad[4]+=e;
            a,b,c,d,e = ugly
            all_ugly[0]+=a; all_ugly[1]+=b; all_ugly[2]+=c; all_ugly[3]+=d; all_ugly[4]+=e;

    #print 'min/max err:', min_err, max_err
    #subplots_adjust(wspace=0, hspace=0, bottom=0.12,left=0.16,right=0.94)
    subplots_adjust(wspace=0, hspace=0, bottom=0.13,left=0.09,right=0.85,top=0.89)
#   figure()
#   hist(err, bins=30)

    bins = [0,5,11,20,30,40,50,60,70,80,90,100]
    try:
        figure()
        hist(all_good[4]+all_bad[4]+all_ugly[4], 
             bins=bins, cumulative=True, histtype='step', normed=True)
        hist(all_good[4],
             bins=bins, cumulative=True, histtype='step', normed=True)
        ylim(0, 1)
        savefig('fullsim-t10-hist.eps')
    except:
        pass

#   figure()
#   hist(all_bad[4], bins=30)
#   figure()
#   hist(err,         bins=30, normed=True, cumulative=True, histtype='step', edgecolor='k')
#   hist(all_good[4], bins=30, normed=True, cumulative=True, histtype='step', edgecolor='g')
#   hist(all_bad[4],  bins=30, normed=True, cumulative=True, histtype='step', edgecolor='r')
#   ylim(0,1)
#   figure()
#   hist(all_good[4], bins=30, cumulative=True, histtype='step', edgecolor='g', hold=True)
#   hist(all_bad[4],  bins=30, cumulative=True, histtype='step', edgecolor='r')
#   #hist(all_ugly[4],  bins=30, cumulative=True, histtype='step', edgecolor='r')
#   axvline(15, ls=':', c='k')
#   xlabel(r'\% Relative error')
#   ylabel(r'Fraction of tests')

#   figure()
#   for i,k in enumerate(sort(rEs.keys())[::-1]):
#       v = rEs[k]
#       
#       p=subplot(1,5,i+1, title='%.4f'%k)
#       p.label_outer()
#       hist(v, bins=10, histtype='step', label='%.4f'%k, range=(0,50))
#       ylim(0,100)
#       xlim(0,100)
    #legend()

def plot8(data, label=None):
    rE_true,closest_star_approach,vals,b,e,minL,maxL = data
    xs = [ v[0] for v in vals ]
    ys = 3000000-array([ v[1] for v in vals ])
    plot(xs, ys, 'ko-', label=label)
    axvline(rE_true, ls='-', c='k')
    axvline(vals[maxL][0], ls=':', c='k')
    xlabel(r'$\theta_E$ [arcsec]', size=14)
    ylabel(r'Effective $\chi^2$', size=14)
    print vals[maxL][0], 3000000-vals[maxL][1], rE_true



#===============================================================================
# Read the input and generate the plots
#===============================================================================


data = []
data2 = defaultdict(lambda: defaultdict(list))

for j, file in enumerate(files):

    print file
    id = int(file[-8:-4])
    params = get_params(id)
    if not params: continue

    run_id, nepochs, gamma_tot, \
    rE_true, closest_star_approach, rE_sample, num_samples = params

    try:
    #if 1:
        vals = loadtxt(file)    
        minL = argmin(vals[:,1])
        maxL = argmax(vals[:,1])
        b,e = vals[0], vals[-1]

        d = [rE_true,closest_star_approach,vals,b,e,minL,maxL]
        data2[nepochs][gamma_tot].append(d)

        data.append([nepochs,gamma_tot,rE_true,closest_star_approach,vals,b,e,minL,maxL])
#        plot8(d, label=id)
        pass
    except:
        print "Skipping %s\n" % file


#print data2[2][3000000]
#print [ x[2] for x in data2[2][3000000] if x[0] == 0.125 and .02<x[1]< .04]
#f = figure(); plot8([ x for x in data2[2][3000000] if x[1] == 0.125 and .02<x[0]< .04][0])
#f = figure(); plot8([ x for x in data2[3][500000] if x[1] == 1.00 and .03<x[0]< .05][0])
#f = figure(); plot8([ x for x in data2[2][3000000] if x[1] == 0.25 and .04<x[0]< .06][0])

try:
    f = figure(); 
    plot8([ x for x in data2[2][3000000] if 0.110<x[1]<0.114 and .04<x[0]< .06][0])
    savefig('likelihood.eps')
    #f = figure(); plot8([ x for x in data2[2][3000000] if x[1] == 0.25 and .04<x[0]< .06][0])
    pass
except:
    pass

#savefig('likelihood.eps')
#show()
#sys.exit(0)

if 1:
    f = figure(); 
    plot7(data2)
    f=figure(f.number)
    f.text(.5, .95, 
        r'Number of photons $\gamma_\mathrm{tot}$',
        horizontalalignment='center', size=14)
    #f.text(.5, .92, 
    #    r'$\gamma_\mathrm{tot}=0.5,1,3,5,7,12\quad\times 10^6$ (left to right)', 
    #    horizontalalignment='center')
    f.text(.02, .5, 
        r'Number of observations $N_\mathrm{obs}$',
        verticalalignment='center', rotation='vertical', size=14)
    #f.text(.04, .5, 
    #    r'$N_\mathrm{obs}=5,4,3,2$ (bottom to top)', 
    #    verticalalignment='center', rotation='vertical')

    f.text(.5, .015, 
        r'Einstein radius $\theta_{E,\mathrm{true}}$ [arcsec]',
        horizontalalignment='center', size=14)

    f.text(.95, .5, 
        r'Impact parameter $p$ [arcsec]',
        verticalalignment='center', rotation=270, size=14)
        #verticalalignment='center', rotation='vertical')

    #               ylabel(r'$p$ [arcsec]')
    savefig('fullsim-t10.eps')

show()

