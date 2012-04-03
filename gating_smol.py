#
# Coded according to Barreda 2011 JCP ref 
#

import smol
import numpy as np
from params import * # must do this for class

class empty:pass

parms = smol.parms


# to make this general, we replace problem.L term with the 
# integral appopriate for the PMF we are considering

# TODO - need to generalize for arbitray PMFs
#        Mostly need to replace computeChanellTerm with call to InteriorProblem gu#        but with a modified PMF (Ve vs Va, for inst)
# returns kPMF [1/Ms]
def computeChannelTerm(problem,expnBV):
  #print problem.D
  #print problem.L
  #print problem.sigma
  #print expnBV
  invkPMF = problem.L/(problem.D*problem.sigma*expnBV)
  kPMF = 1/invkPMF
  kPMF = kPMF * parms.um3_to_M
  return kPMF

# estimate Veff using 3.17
def computeVeff(problem):
  problem.expnBVeff = problem.p_a * problem.expnBVa + problem.p_r * problem.expnBVr

# expnBV = exp(-BV) term 
def slow(problem, expnBV):
  
  kPMF = computeChannelTerm(problem,expnBV)
   
  invk_ss = 1/problem.p_a * (1/problem.k_E_ss + 1/kPMF)
  # print "p_a %f " % problem.p_a

  k_ss = float(1/invk_ss)
  return k_ss


# slow, lig induced 
# eqn 3.24a assuming kappa0=inf
def slow_ligandinduced(problem):
  k_cs_ss = slow(problem,problem.expnBVa)
  return(k_cs_ss)
  # print "slow/induced %f " % k_cs_ss  
  # print "fig:  %e " % (k_cs_ss / problem.k_E_ss)

# slow, lig induced 
# eqn 3.24a assuming kappa0=inf, but Va is replaced by Veff
def slow_ligandindep(problem):
  computeVeff(problem)
  k_ss = slow(problem,problem.expnBVeff)
  return(k_ss)
  # print "slow/indep %f " % k_ss  
  # print "fig:  %e " % (k_ss / problem.k_E_ss)


# fast, lig induced 
# eqn 3.24a assuming kappa0=inf
def fast(problem):
  computeVeff(problem)
  kPMF = computeChannelTerm(problem,problem.expnBVeff)
  #print "kPMF %e" % kPMF 
  invk_if_ss = (1/problem.k_E_ss + 1/kPMF)

  k_if_ss = float(1/invk_if_ss)
  return k_if_ss

def fast_ligandinduced(problem):
  k_if_ss = fast(problem)
  return(k_if_ss)
 
  # # print "fast/induced: %f " % k_if_ss  
  # print "fig:  %e " % (k_if_ss / problem.k_E_ss)

def fast_ligandindep(problem):
  k_if_ss = fast(problem)
  return(k_if_ss)
 
  # print "fast/indep: %f (same as induced)" % k_if_ss  
  # print "fig:  %e " % (k_if_ss / problem.k_E_ss)
 



# inputs
problem = empty()
# Values -far- from protein (effectively ligand independent)
problem.wa=1e9 # rate of forming absorbing state [1/s] 
problem.wr=1e9 # rate of forming reflective state [1/s]
problem.pa=0.5 # probability of being in absorbing state  
problem.Vr = 0 # POtential associated with reflective state [kcal/mol]
# Near protein (can compute from p_a(near) = p_a(far) exp(-(Va-Vr)/kT)
problem.Va = -1 # potential associated when binding site

def Run(problem,result=0):
# compute outer problem 
  if(result != 0 and hasattr(result,"kss")==0):   # TODO replace kss with more informative name 
    result = channel.Run(problem)
    # SHOULD BE DEFINED FROM CHANNELproblem.Va =-1 # DEFINE [kcal/mol] 
    #parms.L = 1 # TODO - only relevant for gated, binding channel. Probably need to adjust problem to reflect this. 

  # dependent var
  problem.L   = result.xL - result.x0
  problem.D = result.D_x[0] # diff constant at booundary [TODO] 
  problem.sigma   = result.Sigma_x[0] # diff constant at booundary [TODO] 
  problem.p_r = 1 - problem.p_a
  problem.expnBVa  = np.exp(-parms.beta * problem.Va)
  problem.expnBVr  = np.exp(-parms.beta * problem.Vr)
  problem.k_E_ss = result.kon # TODO keep consistent
  
  print "TODO: Need to est. characteristic diff time"
  parms.td = 10e-8 # [s] 
  problem.w = problem.wa + problem.wr
  invw = 1/problem.w
  
  # Fast 
  if(invw < parms.td):
    print "1/w %e << tf %e --> Fast gating" % (invw,parms.td)
    k_g  = fast_ligandindep(problem)

  # Slow 
  elif (invw > parms.td):
    print "1/w %e >> tf %e --> Slow gating" % (invw,parms.td)

    if(problem.Va < problem.Veff):
      print "In slow, induced gating regime" 
      k_g  = slow_ligandinduced(problem)

    else:
      print "In slow, indifferent gating regime" 
      k_g  = slow_ligandindep(problem)

  else :
    print "Intermediate regime not supported yet"
    quit()


  result.k_g = k_g  

  return result 


if __name__ == "__main__":
  msg="stochastic_smol.py -root <file>"

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  elif(sys.argv[1]=="-root"):
    root = sys.argv[2]
    problem.fileMesh = root+"_mesh.xml.gz"
    problem.fileSubdomains= root+"_subdomains.xml.gz"
    problem.filePotential= root+"_values.xml.gz"

    RunStochasticSmol(problem)


