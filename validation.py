#
# Simple example that imports mesh from APBS
#
import smol 
import InteriorProblemMaster
from params import * 
import numpy as np
import gating_smol as gating 
from LinearInteriorProblem import *
from view import *


import serca 
import troponin

class empty:pass

problem = empty()
parms = params()

V0 = 5





# validationo of CONSTANT potential with gating 
def ValidateGatedChannel(kon_ext): 
  ## inputs
  a = 1.
  problem.sigma = a
  problem.D = 1.
  Vr = 0
  problem.expnBVr  = np.exp(-parms.beta * Vr)
  
  # from Fig3 
  # Zhou: kon_ext = 4 * problem.D * a 
  problem.k_E_ss = kon_ext
  pa_o_pc = 0.01
  problem.L = 5.* a
  # Zhou: Va =-4.766
  
  
  problem.p_a = pa_o_pc / (1+pa_o_pc) # see notes, since pr = (1-pa)
  problem.p_r = 1 - problem.p_a

  # define range of Vas to test - at some point should run into ligand-indep?
  num=11
  Vas = np.linspace(-V0,0,num)
  k_cs = np.zeros(num)
  k_ss = np.zeros(num)
  k_if = np.zeros(num)
  
  for i in range(num):
    problem.expnBVa  = np.exp(-parms.beta * Vas[i])
    k_cs_ss = gating.slow_ligandinduced(problem)
    k_ss_ = gating.slow_ligandindep(problem)
    k_if_ss = gating.fast_ligandindep(problem)
    # s.b. same as fast_ligandindep gating.fast_ligandinduced(problem)

    k_cs[i] = k_cs_ss
    k_ss[i] = k_ss_
    k_if[i] = k_if_ss

    # for debugging 
    #print "Va (exp) %f (%f) " % (Vas[i],problem.expnBVa)
    #print "slow/induced %f " % k_cs_ss
    #print "fig:  %e " % (k_cs_ss / problem.k_E_ss)
    #print "slow/indep %f " % k_ss_
    #print "fig:  %e " % (k_ss_/problem.k_E_ss)
    #print "fast/induced: %f " % k_if_ss
    #print "fig:  %e " % (k_if_ss / problem.k_E_ss)
    #print "fast/indep: %f (same as induced)" % k_if_ss
    #print "fig:  %e " % (k_if_ss / problem.k_E_ss)



  # plot 
  journalfig = JournalFig()
  plt.title("Gating ",fontsize=journalfig.fontSize)
  plt.ylabel("$k_G/k_E$ $[\\frac{1}{Ms}]$",fontsize=journalfig.fontSize)
  plt.xlabel('$V_a$ $[kcal/mol]$',fontsize=journalfig.fontSize)
  p1, = plt.plot(Vas,k_cs/problem.k_E_ss, 'k-', color='black')
  p2, = plt.plot(Vas,k_ss/problem.k_E_ss, 'k--', color='black')
  p3, = plt.plot(Vas,k_if/problem.k_E_ss, 'k.-', color='black')
  ax = journalfig.ax
  ax.set_yscale('log')
  #handles, labels = ax.get_legend_handles_labels()
  ax.legend([p1,p2,p3],["Slow,induced","Slow,indifferent","Fast"])
  
  
  F = plt.gcf()
  F.savefig("fig1d_sphere.png")

  return k_ss




# validation of linear potential problem 
def  ValidateLinearPotential(kon_ext):
  sigma = 1 #
  num=11
  V0s = np.linspace(-V0,V0,num)
  invkPMFs = np.zeros(num)

  for i in range(num):
    linearIntProb = LinearInteriorProblem(V0=V0s[i])
    result=InteriorProblemMaster.Run(problem,InteriorObject=linearIntProb)
    invkPMFs[i] = result.invkPMF
  
  #invkss = 1/chargedresult.kon + invPMFs 
  invkss = 1/kon_ext + invkPMFs 
  k_ss = 1/invkss
  k_E_ss = kon_ext

  # plot 
  journalfig = JournalFig()
  plt.title("Linear Potential",fontsize=journalfig.fontSize)
  plt.ylabel("$k_{ss}/k_E$ [1/Ms]",fontsize=journalfig.fontSize)
  plt.xlabel('$V_0$ [kcal/mol]',fontsize=journalfig.fontSize)
  plt.plot(V0s,k_ss/k_E_ss, 'k-', color='black')
  journalfig.ax.set_yscale('log')
  F = plt.gcf()
  F.savefig("fig1c_sphere.png")

  return k_ss

 
# all validation cases (Fig 1s)
def ValidationSphere(useStored=0):
  # params 
  root = "/home/huskeypm/scratch/validation/sphere/sphere"
  problem.fileMesh = root+"_mesh.xml.gz"
  problem.fileSubdomains= root+"_subdomains.xml.gz"
  R = 12.5# [Ang]
  R = R / 10000  # [um]
  q = -1  # qlig * qprot [1/C]

  ## No electro 
  # problem 
  problem.filePotential= "none"
  if(useStored==0):
    unchargedresult = smol.Run(problem,pvdFileName="sphere_nocharge.pvd")
    plotslice(problem,unchargedresult,title="No charge",fileName="fig1a_sphere.png")

  else:
    print "Using stored values DONT TRUST"
    unchargedresult = empty()
    unchargedresult.kon = 9.274149e+08
  
  # results 
  kon_analy = 4 * np.pi * R * parms.D * parms.um3_to_M
  print "kon_anal %e pred %e " % (kon_analy, unchargedresult.kon)


  ## w electro 
  # problem 
  problem.filePotential= root+"_values.xml.gz"
  if(useStored==0):
    chargedresult = smol.Run(problem,pvdFileName="sphere_charge.pvd")
    plotslice(problem,chargedresult,title="Charge",fileName="fig1b_sphere.png")
    print chargedresult.kon
  else:
    print "Using stored values DONT TRUST";
    chargedresult = empty()
    chargedresult.kon = 1.603999e+09

  # results 
  # from (18) of Song 
  parms.D = 780 # 7.8 e4 A^2/us --> [um/s]
  R = 8e-4 # 8 A --> 8e-4 um
  kon_elec  = 4 * np.pi * parms.D * q * parms.bulk_conc * parms.um3_to_M
  # assuming r2-->inf
  kon_elec  = np.exp(q/R) / (np.exp(q/R) - np.exp(0)) 
  print "VERIFY kon_anal %e pred %e " % (kon_elec, chargedresult.kon)

  ## linear potential 
  k_lp = ValidateLinearPotential(chargedresult.kon)

  ## gating 
   # fig3 barrade
  print "WARNING: assuming constant potential here - needs to be generalized"
  k_gs = ValidateGatedChannel(chargedresult.kon)


  ## results
  i =0 # where V0 = -5 
  scale = 1.0e9 # normalization in figures 
  print "Sphere & %3.1e & %3.e &  & %3.e & %3.e & NA \\\\" % (unchargedresult.kon/scale,chargedresult.kon/scale,k_lp[i]/scale,k_gs[i]/scale) 
  print "Sphere & %3.1e & %3.e &  & %3.e & %3.e & NA \\\\" % (unchargedresult.kon/1,chargedresult.kon/1,k_lp[i]/1,k_gs[i]/1) 


if __name__ == "__main__":
  msg="smol.py run "



  import sys
  if len(sys.argv) < 2:  
      raise RuntimeError(msg)

  if(sys.argv[1]=="dbg"):
    ValidateGatedChannel(4.)
    ValidateLinearPotential(4.)
    #ValidationTnC(useStored=1)
  elif(sys.argv[1]=="run"):
    ValidationSphere(useStored=0)
    #serca.Validation(useStored=0)
    #troponin.Validation(useStored=0)

  else:
    raise RuntimeError(msg)

