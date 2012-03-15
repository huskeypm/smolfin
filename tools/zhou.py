
import gating

class empty:pass

problem  = empty()
## inputs
a = 1.
problem.sigma = a 
problem.D = 1.

# from Fig3 
problem.k_E_ss = 4.* problem.D * a
pa_o_pc = 0.01
problem.L = 5.* a 
problem.expnBVa = 10.**3
problem.expnBVr = 1.0     

problem.p_a = pa_o_pc / (1+pa_o_pc) # see notes, since pr = (1-pa)
problem.p_r = 1 - problem.p_a


k_cs_ss = gating.slow_ligandinduced(problem)
k_ss = gating.slow_ligandindep(problem)
k_if_ss =  gating.fast_ligandindep(problem)
k_if_ss2 = gating.fast_ligandinduced(problem)

# results 
print "slow/induced %f " % k_cs_ss  
print "fig:  %e " % (k_cs_ss / problem.k_E_ss)
print "slow/indep %f " % k_ss  
print "fig:  %e " % (k_ss / problem.k_E_ss)
print "fast/induced: %f " % k_if_ss  
print "fig:  %e " % (k_if_ss / problem.k_E_ss)
print "fast/indep: %f (same as induced)" % k_if_ss  
print "fig:  %e " % (k_if_ss / problem.k_E_ss)






