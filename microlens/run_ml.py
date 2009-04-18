import sys
from subprocess import Popen

xargs = len(sys.argv) == 2 and sys.argv[1] == 'xargs'

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

if True:  # After Easter3 - with raytrace fix
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

if xargs: xargs_fp = open('_xargs', 'w')

iter    = 0
maxiter = len(gamma_tot) * len(nepochs) * len(rE_true) * len(p_over_rE)
for i,g in enumerate(gamma_tot):
    for j,e in enumerate(nepochs):
        for k,rE in enumerate(rE_true):
            for l,p in enumerate(p_over_rE):

                iter += 1
                if not xargs:
                    print '*******************************************************'
                    print 'Iteration %i of %i' % (iter, maxiter)
                    print '*******************************************************'

                id = '%i.%i.%i.%i' % (i,j,k,l)
                iter_params = '''
run_id="%s"
nepochs=%i
gamma_tot=%ld
rE_true=%.4f
closest_star_approach=%.4f
rE_sample=(%.4f, %.4f)
num_samples=%i
''' % (id, e, g, rE, (p*rE), rE*.5, rE*1.5, num_samples)
#''' % (id, e, g, rE, (p*rE), rE*1, rE*1, num_samples)
                    
                params = 'params.%s.py' % id
                f = open(params, 'w')
                print >>f, iter_params
                print >>f, "iter_params='''%s'''" % iter_params
                f.close()

                if xargs:
                    print >>xargs_fp, 'ml.py %s' % params
                    print 'ml.py %s' % params
                else:
                    proc = Popen([sys.executable,'ml.py', params])
                    #assert(proc.wait() == 0)
                    proc.wait()

if xargs: xargs_fp.close()
