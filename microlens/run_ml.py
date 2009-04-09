import sys
from subprocess import Popen

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

if False:
    gamma_tot  = [1200000L, 7000000L, 12000000L]
    nepochs    = [2, 3, 5, 8] # 1 unlensed observation + (n-1) lensed observations
    rE_true    = [0.005, 0.0125, 0.025, 0.0375, 0.05]
    p_over_rE  = [2,4,8]
    num_samples = 30

if True:
    gamma_tot  = [3000000L]
    #gamma_tot  = [12000000L]
    nepochs    = [5] # 1 unlensed observation + (n-1) lensed observations
    rE_true    = [0.25]
    #rE_true    = [0.375]
    p_over_rE  = [2]
    num_samples = 1

times   = 5
iter    = 0
maxiter = len(gamma_tot) * len(nepochs) * len(rE_true) * len(p_over_rE)
for i,g in enumerate(gamma_tot):
    for j,e in enumerate(nepochs):
        for k,rE in enumerate(rE_true):
            for l,p in enumerate(p_over_rE):

                iter += 1
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
''' % (id, e, g, rE, (p*rE), rE*1, rE*1, num_samples)
#''' % (id, e, g, rE, (p*rE), rE*.7, rE*1.6, num_samples)
                    
                f = open('params.py', 'w')
                print >>f, iter_params
                print >>f, "iter_params='''%s'''" % iter_params
                f.close()

                proc = Popen(['/usr/local/bin/python','ml.py'])
                #assert(proc.wait() == 0)
                proc.wait()

