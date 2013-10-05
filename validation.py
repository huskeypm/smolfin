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
# Simple example that imports mesh from APBS
#
import smol 
import InteriorProblemMaster
from params import * 
import numpy as np
import gating_smol as gating 
from LinearInteriorProblem import *
from view import *

import sys
sys.path.append("example/serca/") 
sys.path.append("example/troponin/") 



class empty:pass

problem = empty()
import smol
parms = smol.parms

#parms = params()

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
  num=101
  V0 = 10
  V0s = np.linspace(-V0,V0,num)
  invkPMFs = np.zeros(num)
  kPMFs = np.zeros(num)

  for i in range(num):
    linearIntProb = LinearInteriorProblem(
      V0=V0s[i], 
      L = 10,
      channelR = 5,
      diff_const=parms.Dchannel
    )
    result=InteriorProblemMaster.Run(problem,InteriorObject=linearIntProb)
    invkPMFs[i] = result.invkPMF
    kPMFs[i] = result.kPMF
 

  #invkss = 1/chargedresult.kon + invPMFs 
  
  #invkss = 1/kon_ext + invkPMFs 
  #k_ss = 1/invkss
  k_ss = kon_ext * kPMFs / (kon_ext + kPMFs)
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

def ValidationChargedSphere(problem,root,useStored=1):
  ## w electro 
  # problem (do once using defaults)  
  problem.filePotential= root+"_values.xml.gz"
  smol.parms = parms
  if(useStored==0):
    chargedresult = smol.Run(problem,pvdFileName="sphere_charge.pvd")
    #plotslice(problem,chargedresult,title="Charge",fileName="fig1b_sphere.png")
  else:
    print "Using stored values DONT TRUST";
    chargedresult = empty()
    chargedresult.kon = 1.818060e+10


  # doa couple times using different values 
  qs = np.array([1.0,0.6, 0.0, -0.6, -1.0])
  kons = np.zeros(np.size(qs))
  msgs=[]
  #for q in qs:
  #  chargedresult = smol.Run(problem,pvdFileName="temp.pvd",q=q)
  #  msgs.append("q=%f kon=%e [1/Ms] %e [1/Mmin]" % (q,chargedresult.kon,chargedresult.kon*60))
#
#  for m in msgs:
#    print m

  # results 
  R = parms.Rsphere
  q = parms.qsphere
  D = parms.D
  # from (18) of Song 
  #parms.D = 780 # 7.8 e4 A^2/us --> [um/s]
  #R = 8e-4 # 8 A --> 8e-4 um
   # I don't like this description, Check Schulten notes
  c= 4207955.8011049731 # value need to match Song paper for D=780, R= 8 A, z*z = 1 TODO
  Q_debye = 7; # [A] See notes from 120405 I think, where Debye shows Q = 7 A for monovalent ions 
  qEff = Q_debye * q
  r_Ang = parms.um_to_Ang * R
  kon_elec = 4 * np.pi * D * qEff * (np.exp(qEff/r_Ang)/(np.exp(qEff/r_Ang)-1)) * c
  kon_pred = chargedresult.kon*60. # 60 [s/min]
  print kon_pred-8.354e11
  assert(np.abs(kon_pred-8.354e11) < 0.01e11)
  print "kon_anal_elec %e pred %e [1/M min]" % (kon_elec, kon_pred)

  return chargedresult

 
# all validation cases (Fig 1s)
def ValidationSphere(useStored=0):
  # params 
  root = "/home/huskeypm/scratch/validation/sphere_120824a/sphere"
  root = "./example/sphere/sphere"
  #root = "/home/huskeypm/scratch/validation/sphere/sphere"
  #root = "./sphere"
  problem.fileMesh = root+"_mesh.xml.gz"
  problem.fileSubdomains= root+"_subdomains.xml.gz"
  smol.parms = parms

  ## No electro 
  # problem 
  problem.filePotential= "none"
  if(useStored==0):
    unchargedresult = smol.Run(problem,pvdFileName="sphere_nocharge.pvd")
    #plotslice(problem,unchargedresult,title="No charge",fileName="fig1a_sphere.png")

  else:
    print "Using stored values DONT TRUST"
    unchargedresult = empty()
    unchargedresult.kon = 3.689228e+09

  
  # results 
  # Made this to resemble Song paper. NOTE: they use an 8 A sphere and 2 A ion, thus
  # giving an effective radius of 10 A. To reproduce their expression, uncomment #SONG below
  # NOTE: they solve their problem on a finite domain, whereas I assume an infinite domain,
  # hence their solution is not exactly kon=4piRD 
  R = parms.Rsphere
  #R = 10 * parms.Ang_to_um  
  #parms.D = 780

  #kon_analy = 4 * np.pi * R * parms.D * parms.um3_to_invM
  # for conceptual reasons I redefine kon in terms of surface area
  kon_analy = (4 * np.pi * R**2)/R * parms.D * parms.um3_to_invM
  print parms.D   
  print "kon_anal %e [1/Ms] %e [1/Mmin] pred %e [1/Ms]" % (kon_analy, kon_analy*60,\
     unchargedresult.kon)
  assert(np.abs(unchargedresult.kon-6.587e9) < 0.01e9)

  ## Charged case
  chargedresult = ValidationChargedSphere(problem,root,useStored=useStored) 
  #print "%e %e\n" % (unchargedresult.kon,chargedresult.kon)
   

  ## linear potential 
  k_lp = ValidateLinearPotential(chargedresult.kon)

  ## gating 
   # fig3 barrade
  print "WARNING: assuming constant potential here - needs to be generalized"
  k_gs = ValidateGatedChannel(chargedresult.kon)

  ## results
  i =0 # where V0 = -5 
  scale = 1.0e9 # normalization in figures 
  msg = []
  m = "Sphere(%1e) & %7.5f & %7.5f &  %7.5f & %e & NA \\\\" % (
    scale,
    unchargedresult.kon/scale,
    chargedresult.kon/scale,
    k_lp[i]/scale,
    k_gs[i]) 

  msg.append(m) 
  m = "Sphere & %3.1e & %3.1e &  %3.1e & %e & NA \\\\" % (unchargedresult.kon/1,chargedresult.kon/1,k_lp[i]/1,k_gs[i]) 
  msg.append(m) 
 
  return msg

if __name__ == "__main__":
  msg="""
\nPurpose:
  Run validation

Usage:\n"""
  msg+= __file__+" dbg/sphere/run "
  msg+="""

Notes:

"""




  import sys
  if len(sys.argv) < 2:  
      raise RuntimeError(msg)

  if(sys.argv[1]=="dbg"):
    ValidateGatedChannel(4.)
    ValidateLinearPotential(4.)
    #ValidationTnC(useStored=1)
  elif(sys.argv[1]=="sphere"):
    m1 = ValidationSphere(useStored=0)
    
  elif(sys.argv[1]=="run"):
    m1 = ValidationSphere(useStored=0)
    import serca 
    m2 = serca.Validation(useStored=0)
    import troponin
    m3 = troponin.Validation(useStored=0)
    for i in m1:
      print i 
    for i in m2:
      print i 
    for i in m3:
      print i 
    
  else:
    raise RuntimeError(msg)

