import sys
from pylab import scatter, show
from math import fabs, log10
from random import random

Xs = [ [] for i in range(2) ]
Ys = [ [] for i in range(2) ]

csas = [0.2, 0.05, 0.01]

for fname in sys.argv[1:]:
    f = open(fname, "r")
    lines = f.readlines()
    nepochs   = int(lines[0].split("=")[1])
    gamma_tot = float(lines[1].split("=")[1])
    rE_true   = float(lines[2].split("=")[1])
    csa       = float(lines[3].split("=")[1]) # closest star approach

    i = [2,5].index(nepochs)
    #i = csas.index(csa) 
    xs = Xs[i]
    ys = Ys[i]

    if gamma_tot == 1e7:
        for l in lines[7:]:
            x,y = map(float, l.split())
            if y > 0.9:
                xs.append(csa / rE_true + (random() - 0.5) * 0.2)
                #xs.append(log10(gamma_tot) + (random() - 0.5) * 0.2)
                #xs.append(gamma_tot/nepochs)
                #xs.append(nepochs)
                #xs.append(nepochs)
                #xs.append(csa)
                ys.append(fabs(rE_true - x) / rE_true + (random() - 0.5) * 0.05)

    f.close()

print Xs, Ys
colors = ['r', 'g', 'b', 'c', 'm', 'y' ]
for i in range(len(Xs)):
    scatter(Xs[i], Ys[i], color=colors[i])
print len(xs)
show()

