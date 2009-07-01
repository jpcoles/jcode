import sys
from subprocess import Popen

xargs = len(sys.argv) >= 2 and sys.argv[1] == 'xargs'
if xargs:
    if len(sys.argv) == 3:
        mldir = sys.argv[2]+'/'


if False: 
    gamma_tot  = [1000000L, 10000000L, 100000000L]
    nepochs    = [2,5]
    rE_true    = [0.035, 0.015, 0.025, 0.10]
    #closest_star_approach = [0.2, 0.05, 0.01]
    p_over_rE  = [2,4,8,16]
    num_samples = 10

if False: # SkiTrip
    gamma_tot  = [1200000L, 7000000L, 12000000L]
    nepochs    = [2, 3, 5, 8] # 1 unlensed observation + (n-1) lensed observations
    rE_true    = [0.125, 0.25, 0.375, 0.5]
    p_over_rE  = [2,4,8]
    num_samples = 30

if False:  # After Easter
    gamma_tot  = [1200000L, 7000000L, 12000000L]
    nepochs    = [2, 3, 5, 8] # 1 unlensed observation + (n-1) lensed observations
    rE_true    = [0.0125, 0.025, 0.0375, 0.05]
    p_over_rE  = [2,4,8]
    num_samples = 30

if False:  # After Easter2
    gamma_tot  = [3000000L, 7000000L, 12000000L]
    nepochs    = [2, 3, 5, 8] # 1 unlensed observation + (n-1) lensed observations
    rE_true    = [0.0125, 0.025, 0.0375, 0.05]
    p_over_rE  = [2,4,8]
    num_samples = 30

if False:  # After Easter3 - with raytrace fix
    gamma_tot  = [3000000L, 7000000L, 12000000L]
    nepochs    = [2, 3, 5, 8] # 1 unlensed observation + (n-1) lensed observations
    rE_true    = [0.0125, 0.025, 0.0375, 0.05]
    p_over_rE  = [2,4,8]
    num_samples = 30

if False:  # After Easter2 - Large rE
    gamma_tot  = [3000000L, 7000000L, 12000000L]
    nepochs    = [2, 3, 5, 8] # 1 unlensed observation + (n-1) lensed observations
    rE_true    = [0.125, 0.25, 0.375, 0.5]
    p_over_rE  = [2,4,8]
    num_samples = 30

if False:
    gamma_tot  = [3000000L]
    #gamma_tot  = [12000000L]
    nepochs    = [5] # 1 unlensed observation + (n-1) lensed observations
    rE_true    = [0.25]
    #rE_true    = [0.375]
    p_over_rE  = [2]
    num_samples = 1

if False:  # For paper (1)
    gamma_tot  = [500000L, 1000000L, 3000000L, 7000000L, 12000000L]
    nepochs    = [2, 3, 4, 5]
    rE_true    = [0.0125, 0.025, 0.0375, 0.05]
    ps         = [0.0250, 0.0375, 0.0500, 0.0750, 0.1000, 0.1125, 0.1500, 0.2000]
    #p_over_rE  = [2,3,4]
    num_samples = 30

if False:  # For paper (2) - testing
    gamma_tot  = [500000L, None, None, None, None]
    nepochs    = [None, None, None, 5]
    rE_true    = [None, None, 0.0375, None]
    #ps         = [0.0250, 0.0375, 0.0500, 0.0750, 0.1000, 0.1125, 0.1500, 1.0000]
    ps         = [None,None,None,None,None,None,None, .5000]
    #p_over_rE  = [2,3,4]
    num_samples = 30

if True:  # For paper (3)
    gamma_tot  = [500000L, 1000000L, 3000000L, 7000000L, 12000000L]
    nepochs    = [2, 3, 4, 5]
    rE_true    = [0.01, 0.02, 0.03, 0.04, 0.05]
    ps         = [0.0625, 0.125, 0.2500, 0.5, 1]
    num_samples = 30

if xargs: xargs_fp = open('_xargs', 'w')

params = 'params.py'
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
                if not xargs:
                    print '*******************************************************'
                    print 'Iteration %i of %i' % (iter, maxiter)
                    print '*******************************************************'

                #id = '%i.%i.%i.%i' % (i,j,k,l)

                id = (i+1) * 1000 + (j+1)*100 + (k+1)*10 + (l+1)*1

#               iter_params = '''\
#   if id==%i: return [%i,%i,%ld,%.4f,%.4f,(%.4f,%.4f),%i]''' \
#       % (id, id,e, g, rE, (p*rE), rE*.5, rE*1.5, num_samples)

                iter_params = '''\
    if id==%i: return [%i,%i,%ld,%.4f,%.4f,(%.4f,%.4f),%i]''' \
        % (id, id,e, g, rE, p, rE*.5, rE*1.5, num_samples)
                print >>f, iter_params

                if xargs:
                    print >>xargs_fp, '%sml.py %i' % (mldir, id)
                    print '%sml.py %i' % (mldir, id)
                else:
                    proc = Popen([sys.executable,'ml.py', params])
                    #assert(proc.wait() == 0)
                    proc.wait()

f.close()

if xargs: 
    xargs_fp.close()
    print 'Usage: xargs -P 32 -n 2 python25 < _xargs'
