import sys
from pylab import scatter, show, xlabel, ylabel, figure
from random import random
from math import log10, exp, log 

files = sys.argv[1:]

def gamma_color(v):
    if v == 1200000L:  return 'red'
    if v == 7000000L:  return 'green'
    if v == 12000000L: return 'blue'

def gamma_offs(v):
    if v == 1200000L:  return 0
    if v == 7000000L:  return 1
    if v == 12000000L: return 2

def gamma_size(v):
    if v == 1200000L:  return 40
    if v == 7000000L:  return 80
    if v == 12000000L: return 160

def nepochs_color(v):
    if v == 2: return 'red'
    if v == 3: return 'green'
    if v == 5: return 'blue'
    if v == 8: return 'magenta'

def nepochs_size(v):
    return 2**v

def plot1(data):

    xs,ys,cs,ss = [], [], [], []

    for nepochs,gamma_tot,rE_true,closest_star_approach,vals,b,e,m,M in data:

        if m[0] in [0, len(vals)-1]:
            ss.append(10)
            cs.append('black')
            assert False
        else:
            err = 100 * abs(gamma_tot - m[1]) / gamma_tot
            xs.append(100 * (rE_true / closest_star_approach) + nepochs)
            ys.append(err)
            ss.append(2**nepochs)
            cs.append(gamma_color(gamma_tot))

    figure()
    ylabel('% err')
    xlabel('rE_true / closest_star_approach')
    scatter(xs, ys, s=ss, c=cs, marker='o', alpha=0.5) #, facecolors='none')

def plot2(data):

    xs,ys,cs,ss = [], [], [], []

    for nepochs,gamma_tot,rE_true,closest_star_approach,vals,b,e,m,M in data:

        if m[0] in [0, len(vals)-1]:
            ss.append(10)
            cs.append('black')
            assert False
        else:
            err = 100 * abs(gamma_tot - m[1]) / gamma_tot
            xs.append(100 * closest_star_approach + nepochs)
            ys.append(err)
            ss.append(2**nepochs)
            cs.append(gamma_color(gamma_tot))

    figure()
    ylabel('% err')
    xlabel('closest star approach')
    scatter(xs, ys, s=ss, c=cs, marker='o', alpha=0.5) #, facecolors='none')

def plot3(data):

    xs,ys,cs,ss = [], [], [], []

    for nepochs,gamma_tot,rE_true,closest_star_approach,vals,b,e,m,M in data:

        rE_err = 100 * abs(rE_true   - vals[M][0]) / rE_true
        lh_err = 100 * abs(gamma_tot - vals[M][1]) / gamma_tot # Likelihood error

        xs.append(100 * rE_true + nepochs)
        #xs.append(100 * (rE_true / closest_star_approach) + nepochs)
        ys.append(rE_err)

        if 0 < M < len(vals)-1:
            ss.append(gamma_size(gamma_tot))
            #ss.append(nepochs_size(nepochs))
            cs.append(gamma_color(gamma_tot))
        else:
            ss.append(10)
            cs.append('black')

    figure()
    ylabel('% err')
    xlabel('rE_true')
    scatter(xs, ys, s=ss, c=cs, marker='o', alpha=0.5) #, facecolors='none')

def plot4(data):

    xs,ys,cs,ss = [], [], [], []

    for nepochs,gamma_tot,rE_true,closest_star_approach,vals,b,e,m,M in data:

        if m[0] in [0, len(vals)-1]:
            ss.append(10)
            cs.append('black')
            assert False
        else:
            err = 100 * abs(gamma_tot - m[1]) / gamma_tot
            xs.append(10 * log(closest_star_approach / rE_true, 2) + nepochs)
            ys.append(err)
            ss.append(2**nepochs)
            cs.append(gamma_color(gamma_tot))

    figure()
    ylabel('% err')
    xlabel('log_2(closest_star_approach / rE_true)')
    scatter(xs, ys, s=ss, c=cs, marker='o', alpha=0.5) #, facecolors='none')

def plot5(data):

    xs,ys,cs,ss = [], [], [], []

    for nepochs,gamma_tot,rE_true,closest_star_approach,vals,b,e,m,M in data:

        if m[0] in [0, len(vals)-1]:
            ss.append(10)
            cs.append('black')
            assert False
        else:
            err = 100 * abs(gamma_tot - m[1]) / gamma_tot
            xs.append(nepochs)
            ys.append(err)
            ss.append(100)
            cs.append(gamma_color(gamma_tot))

    figure()
    ylabel('% err')
    xlabel('nepochs')
    scatter(xs, ys, s=ss, c=cs, marker='o', alpha=0.5) #, facecolors='none')

def plot6(data):

    xs,ys,cs,ss = [], [], [], []

    for nepochs,gamma_tot,rE_true,closest_star_approach,vals,b,e,m,M in data:

        if m[0] in [0, len(vals)-1]:
            ss.append(10)
            cs.append('black')
            assert False
        else:
            err = 100 * abs(gamma_tot - m[1]) / gamma_tot
            xs.append(10*nepochs + gamma_offs(gamma_tot))
            ys.append(rE_true / closest_star_approach)
            ss.append(10*err + 10)
            cs.append(gamma_color(gamma_tot))

    figure()
    ylabel('rE_true / closest_star_approach')
    xlabel('closest_star_approach')
    scatter(xs, ys, s=ss, c=cs, marker='o', alpha=0.5) #, facecolors='none')

data = []

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

    vals = []

    m = [0,+1e30]
    M = [0,-1e30]

    for i,l in enumerate(f):
        l = l.strip()
        if len(l) == 0: break;
        if l.startswith('#'): continue
        l = map(float, l.split())
        if l[1] < m[1]: m = [len(vals), l[1]]
        if l[1] > M[1]: M = [len(vals), l[1]]
        vals.append(l)

    m = m[0]
    M = M[0]
    b = vals[0]
    e = vals[-1]

    data.append([nepochs,gamma_tot,rE_true,closest_star_approach,vals,b,e,m,M])

    f.close()

#plot1(data)
#plot2(data)
plot3(data)
#plot4(data)
#plot5(data)
#plot6(data)
show()

