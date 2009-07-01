from __future__ import division
import sys
sys.path.append('.')
import pylab
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
#rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Sans serif']})
rc('text', usetex=True)
from pylab import scatter, show, xlabel, ylabel, figure, subplot, gca, hist, colorbar, gray, jet, hot, pink, xticks, yticks, xlim, ylim, axvline, legend, plot, savefig
from matplotlib.pyplot import subplots_adjust
from matplotlib.cm import get_cmap
from random import random
from math import log10, exp, log
from collections import defaultdict
from numpy import inf, array, amax, log2, sort

from params import get_params, get_all_params


ps = get_all_params()

files = sys.argv[1:]

colors = ['blue', 'green', 'red', 'magenta', 'cyan']

def gamma_color(v):
    return colors[ps['gamma_tot'].index(v)]
#   if v == 1200000L:  return 'red'
#   if v == 7000000L:  return 'green'
#   if v == 12000000L: return 'blue'
#   assert 0, 'Unrecognized photon count'

def gamma_offs(v):
    return [0,1,2][ps['gamma_tot'].index(v)]
#   if v == 1200000L:  return 0
#   if v == 7000000L:  return 1
#   if v == 12000000L: return 2
#   assert 0, 'Unrecognized photon count'

def gamma_size(v):
    return [40, 80, 160][ps['gamma_tot'].index(v)]
#   if v == 1200000L:  return 40
#   if v == 7000000L:  return 80
#   if v == 12000000L: return 160
#   assert 0, 'Unrecognized photon count'

def nepochs_color(v):
    return colors[ps['nepochs'].index(v)]
#   if v == 2: return 'red'
#   if v == 3: return 'green'
#   if v == 5: return 'blue'
#   if v == 8: return 'magenta'
#   assert 0, 'Unrecognized epoch count'

def nepochs_size(v):
    return 2**v

def mass(r):
    G = 6.67300e-11
    c = 299792458 
    pc = 3.08568025e16
    return r**2 * pc * c**2 / (0.09**2 * 4*G)

def plot1(data):

    xs,ys,cs,ss = [], [], [], []

    count = 0
    for nepochs,gamma_tot,rE_true,closest_star_approach,vals,b,e,m,M in data:

        rE_err = 100 * abs(rE_true - vals[M][0]) / rE_true
        xs.append(1000 * (rE_true / closest_star_approach) + nepochs)
        ys.append(rE_err)
        ss.append(2**nepochs)

        #if nepochs==2 and gamma_tot==1200000: print '*', M, rE_err

        if 0 < M < len(vals)-1:
            cs.append(gamma_color(gamma_tot))
        else:
            #ss.append(10)
            cs.append('black')
            print 'bad', closest_star_approach/rE_true
            count+=1

    print 'count=%i' % count

    figure()
    ylabel('% err')
    xlabel('rE_true / closest_star_approach')
    scatter(xs, ys, s=ss, c=cs, marker='o', alpha=0.5) #, facecolors='none')

def plot2(data):

    xs,ys,cs,ss = [], [], [], []

    for nepochs,gamma_tot,rE_true,closest_star_approach,vals,b,e,m,M in data:

        err = 100 * abs(rE_true - vals[M][0]) / rE_true
        #err = 100 * abs(gamma_tot - vals[M][1]) / gamma_tot
        #xs.append(1000 * closest_star_approach + nepochs)
        x = closest_star_approach
        x += random()*0.001-0.0005

        xs.append(x)
        ys.append(err)
        ss.append(1000*rE_true)
        #ss.append(80)
        #ss.append(2**nepochs)

        if 0 < M < len(vals)-1:
            cs.append(gamma_color(gamma_tot))
        else:
            cs.append('black')

    figure()
    ylabel('% err')
    xlabel('closest star approach')
    scatter(xs, ys, s=ss, c=cs, marker='o', alpha=0.5) #, facecolors='none')

def plot3a(data):

    xs,ys,cs,ss = [], [], [], []

    for nepochs,gamma_tot,rE_true,closest_star_approach,vals,b,e,m,M in data:

        rE_err = 100 * abs(rE_true - vals[M][0]) / rE_true
        x = rE_true
        x += random()*0.001-0.0005
        xs.append(x)
        ys.append(rE_err)
        ss.append(1000*closest_star_approach)
        #ss.append(2**nepochs)
        #ss.append(80)

        if 0 < M < len(vals)-1:
            #ss.append(nepochs_size(nepochs))
            cs.append(gamma_color(gamma_tot))
        else:
            cs.append('black')

    figure()
    ylabel('% err')
    xlabel('rE_true')
    scatter(xs, ys, s=ss, c=cs, marker='o', alpha=0.5) #, facecolors='none')

def plot3b(data):

    xs,ys,cs,ss = [], [], [], []

    for nepochs,gamma_tot,rE_true,closest_star_approach,vals,b,e,m,M in data:

        rE_err = 100 * abs(rE_true   - vals[M][0]) / rE_true
        lh_err = 100 * abs(gamma_tot - vals[M][1]) / gamma_tot # Likelihood error

        xs.append(1000 * rE_true + nepochs)
        ys.append(rE_err)
        #ss.append(gamma_size(gamma_tot))
        #ss.append(2**nepochs)
        ss.append(80)

        if 0 < M < len(vals)-1:
            #ss.append(nepochs_size(nepochs))
            cs.append(gamma_color(gamma_tot))
        else:
            cs.append('black')

    figure()
    ylabel('% err')
    xlabel('% rE_true err')
    scatter(xs, ys, s=ss, c=cs, marker='o', alpha=0.5) #, facecolors='none')

def plot4(data):

    xs,ys,cs,ss = [], [], [], []

    figure()
    for nepochs,gamma_tot,rE_true,closest_star_approach,vals,b,e,m,M in data:

        err = 100 * abs(rE_true - vals[M][0]) / rE_true
        #xs.append(10 * log(closest_star_approach / rE_true, 2) + nepochs)

        x = closest_star_approach / rE_true
        x += random()*0.5-0.25

        xs.append(x)
        ys.append(err)
        #ss.append(2**nepochs)
        ss.append(80)

        if 0 < M < len(vals)-1:
            cs.append(gamma_color(gamma_tot))
        else:
            cs.append('black')

    ylabel('% err')
    xlabel('closest_star_approach / rE_true')
    scatter(xs, ys, c=cs, marker='o', alpha=0.5, lw=1) #, facecolors='none')

def plot5(data):

    xs,ys,cs,ss = [], [], [], []

    for nepochs,gamma_tot,rE_true,closest_star_approach,vals,b,e,m,M in data:

        #err = 100 * abs(gamma_tot - vals[M][1]) / gamma_tot
        err = 100 * abs(rE_true - vals[M][0]) / rE_true
        xs.append(nepochs)
        ys.append(err)
        ss.append(100)

        if 0 < M < len(vals)-1:
            cs.append(gamma_color(gamma_tot))
        else:
            cs.append('black')

    figure()
    ylabel('% err')
    xlabel('nepochs')
    scatter(xs, ys, s=ss, c=cs, marker='o', alpha=0.5) #, facecolors='none')

def plot6(data):

    xs,ys,cs,ss = [], [], [], []

    for nepochs,gamma_tot,rE_true,closest_star_approach,vals,b,e,m,M in data:

        #err = 100 * abs(gamma_tot - vals[M][1]) / gamma_tot
        err = 100 * abs(rE_true - vals[M][0]) / rE_true
        xs.append(10*nepochs + gamma_offs(gamma_tot))
        ys.append(rE_true / closest_star_approach)
        ss.append(10*err + 10)

        if 0 < M < len(vals)-1:
            cs.append(gamma_color(gamma_tot))
        else:
            cs.append('black')

    figure()
    ylabel('rE_true / closest_star_approach')
    xlabel('closest_star_approach')
    scatter(xs, ys, s=ss, c=cs, marker='o', alpha=0.5) #, facecolors='none')

def plot7(data2):


    subn = len(ps['nepochs'])
    subm = len(ps['gamma_tot'])

    err = []

    gray()
    cm_gray = get_cmap()
    #pink()
    #cm_jet = get_cmap()

    all_good  = [[], [], [], [], []]
    all_bad   = [[], [], [], [], []]
    all_ugly  = [[], [], [], [], []]

    rEs = defaultdict(list)

    i = 0
    for epoch in ps['nepochs']:
        for g in ps['gamma_tot']:
            i += 1
            #if epoch is None or g is None: continue
            data = data2[epoch][g]
            if not data: continue

            p = subplot(subn, subm, i)
            if p.is_first_col() and p.is_last_row():
                xlabel(r'$r_{E,\mathrm{true}}$ [arcsec]')
                ylabel(r'$p$ [arcsec]')

            p.label_outer()
            pylab.setp(p.get_xticklabels(), rotation=45)
            #pylab.setp(p.get_yticklabels(), rotation=45)
            good  = [[], [], [], [], []]
            bad   = [[], [], [], [], []]
            ugly  = [[], [], [], [], []]
            for rE_true,closest_star_approach,vals,b,e,m,M in data:

                #rE_err = 100 * abs(rE_true - vals[M][0]) / rE_true
                rE_err = 100 * abs(mass(rE_true) - mass(vals[M][0])) / mass(rE_true)

                if rE_err > 100: rE_err = 100
                x = rE_true #+ random()*0.001-0.0005
                y = closest_star_approach #+ random()*0.001-0.0005

                rEs[rE_true].append(rE_err)

                if 0 < M < len(vals)-1:
                    #ss.append(nepochs_size(nepochs))
                    #cs.append(gamma_color(gamma_tot))
                    if rE_err <= 15:
                        good[0].append(x)
                        good[1].append(y)
                        good[2].append(int(rE_err*5))
                        if rE_err < 5:
                            good[3].append((0.,0.,0.,1.))
                        elif rE_err < 10:
                            good[3].append((0.3,0.3,0.3,1.))
                        else:
                            good[3].append((0.6,0.6,0.6,1.))
                        good[4].append(rE_err)
                    else:
                        bad[0].append(x)
                        bad[1].append(y)
                        bad[2].append(int(rE_err*5))
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
                    #cs.append(cm(1.))
                    #cs.append(cm_gray(1.))
                    #ss.append(90)
                    #cs.append(100)
                    z = -1
                    #cs.append('black')
                    #cs.append('black')

            #print cs
            #cm.set_over()
#           w = cs > 15
#           cs = array(cs)
#           cs[w] = cm_gray(cs[w] / 100.)
#           w = cs <= 15
#           cs[w] = cm_jet(cs[w] / 100.)

            #cs /= amax(cs)
            #if not good[0] and not bad[0]: continue

            if 0:
                if good[0]:
                    scatter(array(good[1])/array(good[0]), good[4], marker='o')
                if bad[0]:
                    scatter(array(bad[1])/array(bad[0]), bad[4], marker='x')
                ylim(0,50)
                xlim(0,60)
                #yticks(ps['ps'])
                #xticks(ps['rE_true'])
            else:

                if good[0]:
                    scatter(good[0], log2(good[1]), s=good[2], c=good[3], marker='o', zorder=1)
                if bad[0]:
                    scatter(bad[0], log2(bad[1]), s=bad[2], c=bad[3], marker='o', zorder=0, facecolors='none')
                #if ugly[0]:
                    #scatter(ugly[0], log2(ugly[1]), s=ugly[2], c=ugly[3], marker='o', zorder=-1, facecolors='none')
                yticks(log2(ps['ps']), ps['ps'])
                xticks(ps['rE_true'])

            a,b,c,d,e = good
            all_good[0]+=a; all_good[1]+=b; all_good[2]+=c; all_good[3]+=d; all_good[4]+=e;
            a,b,c,d,e = bad
            all_bad[0]+=a; all_bad[1]+=b; all_bad[2]+=c; all_bad[3]+=d; all_bad[4]+=e;
            a,b,c,d,e = ugly
            all_ugly[0]+=a; all_ugly[1]+=b; all_ugly[2]+=c; all_ugly[3]+=d; all_ugly[4]+=e;

            #scatter(xs, ys, s=ss, c=cs, marker='o', alpha=0.5) #, facecolors='none')
    #ylabel(r'$p$')
    #xlabel(r'$r_E$ true')

    #print 'min/max err:', min_err, max_err
    subplots_adjust(wspace=0, hspace=0, bottom=0.12,left=0.16,right=0.94)
#   figure()
#   hist(err, bins=30)
#   figure()
#   hist(all_good[4], bins=30)
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

def plot8(data):
    rE_true,closest_star_approach,vals,b,e,minL,maxL = data
    xs = [ v[0] for v in vals ]
    ys = 3000000-array([ v[1] for v in vals ])
    #ys = [ v[1] for v in vals ]
    plot(xs, ys, 'ko-')
    axvline(rE_true, ls='-', c='k')
    axvline(vals[maxL][0], ls=':', c='k')
    xlabel(r'$\theta_E$ [arcsec]')
    ylabel(r'Effective $\chi^2$')
    #ylabel(r'$\gamma$ $[+2.966\times10^6]$')
    print vals[maxL][0], rE_true

data = []
data2 = defaultdict(lambda: defaultdict(list))

if 0:
    for j, file in enumerate(files):
        print file
        f = open(file, 'r')
        found_header = False
        for l in f:
            l = l.strip()
            if len(l) == 0: 
                if found_header: 
                    break;
                else:
                    continue
            else:
                found_header = True
                    
            try:
                exec l 
                pass
            except:
                print "Reached end of code"
                break

for j, file in enumerate(files):
    print file
    id = int(file[-8:-4])

    f = open(file, 'r')
#   for l in f:
#       l = l.strip()
#       if len(l) == 0: break

    params = get_params(id)
    if not params: continue

    run_id, nepochs, gamma_tot, \
    rE_true, closest_star_approach, rE_sample, num_samples = params

    #if nepochs != 3: continue

    vals = []

    minL = [0,inf]
    maxL = [0,-inf]

    for i,l in enumerate(f):
        l = l.strip()
        if len(l) == 0: break;
        if l.startswith('#'): continue
        l = map(float, l.split())
        if l[1] < minL[1]: minL = [len(vals), l[1]]
        if l[1] > maxL[1]: maxL = [len(vals), l[1]]
        vals.append(l)

    minL = minL[0]
    maxL = maxL[0]
    b = vals[0]
    e = vals[-1]

    data2[nepochs][gamma_tot].append([rE_true,closest_star_approach,vals,b,e,minL,maxL])

    data.append([nepochs,gamma_tot,rE_true,closest_star_approach,vals,b,e,minL,maxL])

    f.close()

#plot1(data)
#plot2(data)
#plot3a(data)
#plot3b(data)
#plot4(data)
#plot5(data)
#plot6(data)
#print data2[2][3000000]
#print [ x[2] for x in data2[2][3000000] if x[0] == 0.125 and .02<x[1]< .04]
#f = figure(); plot8([ x for x in data2[2][3000000] if x[1] == 0.125 and .02<x[0]< .04][0])
f = figure(); plot8([ x for x in data2[2][3000000] if x[1] == 0.25 and .04<x[0]< .06][0])
savefig('likelihood.eps')
#show()
#sys.exit(0)

f = figure(); 
plot7(data2)
f=figure(f.number)
f.text(.5, .95, 
    r'Number of photons $\gamma_\mathrm{tot}$',
    horizontalalignment='center')
f.text(.5, .92, 
    r'$\gamma_\mathrm{tot}=0.5,1,3,5,7,12\quad\times 10^6$ (left to right)', 
    horizontalalignment='center')
f.text(.02, .5, 
    r'Number of epochs $N_\mathrm{obs}$',
    verticalalignment='center', rotation='vertical')
f.text(.04, .5, 
    r'$N_\mathrm{obs}=5,4,3,2$ (bottom to top)', 
    verticalalignment='center', rotation='vertical')
savefig('fullsim-t2.eps')
#colorbar()

show()

