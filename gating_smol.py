"""
----------------------------------
Smolfin - a numerical solver of the Smoluchowski equation for interesting geometries
Copyright (C) 2012 Peter Kekenes-Huskey, Ph.D., huskeypm@gmail.com

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA


----------------------------------------------------------------------------
"""
#
# Coded according to Barreda 2011 JCP ref 
#

# Using log10 of gating/characteristic diffusion time numbers to 
# determine fast versus slow gating 

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
  kPMF = kPMF * parms.um3_to_invM
  return kPMF

# estimate Veff using 3.17
def computeVeff(problem):
  problem.expnBVeff = problem.p_a * problem.expnBVa + problem.p_r * problem.expnBVr

# expnBV = exp(-BV) term 
def slow(problem, expnBV,results=0):
  #quit()
  
  if(results==0):
    computeVeff(problem)
    kPMF = computeChannelTerm(problem,expnBV)
  else:
    print "WARNING: Need to compute Vr and Va separately pass to function. Not supported right now"
    kPMF = results.kPMF
   
  invk_ss = 1/problem.p_a * (1/problem.k_E_ss + 1/kPMF)
  k_ss = float(1/invk_ss)
  print "k_E_ss %e kPMF (Veff or Va) %e p_a %f --> k_ss %e" % (problem.k_E_ss,kPMF,problem.p_a,k_ss)
  return k_ss


# slow, lig induced 
# eqn 3.24a assuming kappa0=inf
def slow_ligandinduced(problem,results=0):
  k_cs_ss = slow(problem,problem.expnBVa,results=results)
  print "slow/induced %e " % k_cs_ss  
  print "fig:  %e " % (k_cs_ss / problem.k_E_ss)
  return(k_cs_ss)

# slow, lig induced 
# eqn 3.24a assuming kappa0=inf, but Va is replaced by Veff
def slow_ligandindep(problem,results=0):
  k_ss = slow(problem,problem.expnBVeff,results=results)
  print "slow/indep %e " % k_ss  
  print "fig:  %e " % (k_ss / problem.k_E_ss)
  return(k_ss)


# fast, lig induced 
# eqn 3.24a assuming kappa0=inf
def fast(problem,results=0):
  if(results==0):
    computeVeff(problem)
    kPMF = computeChannelTerm(problem,problem.expnBVeff)
  else:
    kPMF = results.kPMF

  invk_if_ss = (1/problem.k_E_ss + 1/kPMF)
  k_if_ss = float(1/invk_if_ss)
  print "Fast: k_E_ss %e kPMF %e (Veff) --> k_if_ss %e" % (problem.k_E_ss,kPMF,k_if_ss)

  return k_if_ss

def fast_ligandinduced(problem,results=0):
  k_if_ss = fast(problem,results=results)
  return(k_if_ss)
 
  #print "fast/induced: %e " % k_if_ss  
  #print "fig:  %e " % (k_if_ss / problem.k_E_ss)

def fast_ligandindep(problem,results=0):
  k_if_ss = fast(problem,results=results)
  return(k_if_ss)
 
  #print "fast/indep: %e (same as induced)" % k_if_ss  
  #print "fig:  %e " % (k_if_ss / problem.k_E_ss)
 



# inputs
problem = empty()
# Values -far- from protein (effectively ligand independent)
problem.wa=1e9 # rate of forming absorbing state [1/s] 
problem.wr=1e9 # rate of forming reflective state [1/s]
problem.pa=0.5 # probability of being in absorbing state  
problem.Vr = 0 # POtential associated with reflective state [kcal/mol]
# Near protein (can compute from p_a(near) = p_a(far) exp(-(Va-Vr)/kT)
problem.Va = -1 # potential associated when binding site

def Run(problem,results=0):
  # compute outer problem 
  if(results!= 0):
    # WARNING: there is a subtlety in the distinction between Va and Veff, which was lost in 
    # PMF calculations, since the fluctuations between open/closed were fast and thus both contributed
    # to the potential of mean force. To do this right, I need to more carefullly separate the 
    # statistics between the two states and treat them as separate PMFs 
    #if(hasattr(result,"kss")==0):   # TODO replace kss with more informative name 
    if(hasattr(results,"kPMF")==0):   # 
      results = channel.Run(problem)
    else:
     # do nothing, use pmf values inside result 
      1


    problem.expnBVa  = "NOT USED"
    problem.expnBVr  = "NOT USED"
    problem.expnBVeff= "NOT USED"
    

  # assume Va,Vr are constant along channel (channel will be called later)  
  else:
    print "Not sure why I should have ended up here. Should always pass results"
    # here we assume that Va,Vr are constant over the 
    problem.expnBVa  = np.exp(-parms.beta * problem.Va)
    problem.expnBVr  = np.exp(-parms.beta * problem.Vr)
    problem.L   = result.xL - result.x0
    problem.D = result.D_x[0] # diff constant at booundary [TODO, as in reality this should vary] 
    problem.sigma   = result.Sigma_x[0] # diff constant at booundary [TODO, as in reality this should vary] 
    computeVeff(problem)
    quit()

  # dependent var
  problem.p_r = 1 - problem.p_a
  problem.k_E_ss = results.kon # TODO keep consistent
  
  parms.td = 2e-10 #  for 10 A , 780 D [um^2/s] using t = r^2/6D 
  problem.w = problem.wa + problem.wr
  invw = 1/problem.w
  
  # Fast 
  logdiff = np.log10(invw) - np.log10(parms.td)
  if(logdiff < -1):       
    print "1/w %e << tf %e --> Fast gating" % (invw,parms.td)
    k_g  = fast_ligandindep(problem,results=results)
    result.k_g = "%e" % k_g

  # Slow 
  elif (logdiff > 1):       
    print "1/w %e >> tf %e --> Slow gating" % (invw,parms.td)

    if(problem.Va < problem.Veff):
      print "In slow, induced gating regime" 
      k_g  = slow_ligandinduced(problem,results=results)
      result.k_g = "%e" % k_g

    else:
      print "In slow, indifferent gating regime" 
      k_g  = slow_ligandindep(problem,results=results)
      result.k_g = "%e" % k_g

  else :
    print "Intermediate regime not supported yet - Giving bounds"
    k_gf = fast_ligandindep(problem,results=results)
    k_gsi= slow_ligandinduced(problem,results=results)
    k_gs = slow_ligandindep(problem,results=results)   
    min = np.min([k_gf,k_gsi,k_gs])
    max = np.max([k_gf,k_gsi,k_gs])
    results.k_g = "%e-%e" % (min,max)


  return results 


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


