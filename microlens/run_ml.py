#
# run_ml.py - Parameter generating program for ml.py
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

import sys
from subprocess import Popen
from random import shuffle

params = 'params.py'

if len(sys.argv) == 2:
    mldir = sys.argv[1]+'/'

if False:  # For paper (3)
    gamma_tot  = [500000L, 1000000L, 3000000L, 7000000L, 12000000L]
    nepochs    = [2, 3, 4, 5]
    rE_true    = [0.01, 0.02, 0.03, 0.04, 0.05]
    ps         = [0.0625, 0.125, 0.2500, 0.5, 1]
    num_samples = 30

if False:  # For paper (3) First referee revision
    gamma_tot  = [500000L, 1000000L, 3000000L, 7000000L, 12000000L]
    nepochs    = [2, 3, 4, 5]
    rE_true    = [0.01, 0.02, 0.03, 0.04, 0.05]
    ps         = [0.0625, 0.125, 0.2500, 0.5, 1]
    num_samples = 30

if False:  # For paper (3) 
    x = 2./71 * 10
    gamma_tot  = [500000L, 1000000L, 3000000L, 7000000L, 12000000L]
    nepochs    = [2, 3, 4, 5]
    rE_true    = [0.01, 0.02, 0.03, 0.04, 0.05]
    ps         = [x/4, x/2, x, 2*x, 4*x]
    num_samples = 30

if True:  # For paper. FINAL VALUES.
    x = 2./71
    gamma_tot  = [500000L, 1000000L, 3000000L, 7000000L, 12000000L]
    nepochs    = [2, 3, 4, 5]
    rE_true    = [0.01, 0.02, 0.03, 0.04, 0.05]
    ps         = [x * 2**i for i in range(1,6)]
    #ps         = [0.0625, 0.125, 0.2500, 0.5, 1]
    num_samples = 30


f = open(params, 'w')

print >>f, 'def get_all_params():'
print >>f, '''\
    return {'gamma_tot':%s,
            'nepochs':%s,
            'rE_true':%s,
            'ps':%s,
            'num_samples':%s}''' % (gamma_tot, nepochs, rE_true, ps, num_samples)
print >>f
print >>f, 'def get_params(id):'

xargs_output = []

iter    = 0
maxiter = len(gamma_tot) * len(nepochs) * len(rE_true) * len(ps)
for i,g in enumerate(gamma_tot):
    if not g: continue
    for j,e in enumerate(nepochs):
        if not e: continue
        #if g == 12000000L and e == 3: print 'qwerqwerqwer'; continue
        for k,rE in enumerate(rE_true):
            if not rE: continue
            #for l,p in enumerate(p_over_rE):
            for l,p in enumerate(ps):
                if not p: continue

                iter += 1
                id = (i+1) * 1000 + (j+1)*100 + (k+1)*10 + (l+1)*1

                iter_params = '''\
    if id==%i: return [%i,%i,%ld,%.4f,%.4f,(%.4f,%.4f),%i]''' \
                    % (id, id,e, g, rE, p, rE*.5, rE*1.5, num_samples)
                print >>f, iter_params

                xargs_output.append([mldir, id])
                print '%sml.py %i' % (mldir, id)
f.close()

#shuffle(xargs_output)

xargs_fp = open('_xargs', 'w')
for mldir, id in xargs_output:
    print >>xargs_fp, '%sml.py %i' % (mldir, id)
xargs_fp.close()

print 'Usage: xargs --verbose -P 22 -n 2 nice -19 python < _xargs'
