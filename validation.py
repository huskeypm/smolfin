#
# Simple example that imports mesh from APBS
#
import smol 
import InteriorProblemMaster
from params import * 
import numpy as np
import gating_smol as gating 
from LinearInteriorProblem import *

class empty:pass

problem = empty()
parms = params()

V0 = 5

def plotslice(problem,result,title="no title",fileName="slice.png"):
    meshcoor = problem.mesh.coordinates()
    
    # assuming molecule is within 50 of middle of grid 
    # want 500 points in each dir (resolution)
    range = 50
    numpt = 500
    incr = numpt/ (2 * range) 
    #(grid_x,grid_y,grid_z) = np.mgrid[0:0:1j,-range:range:(incr*1j),-range:range:(incr*1j)]
    (grid_x,grid_y,grid_z) = np.mgrid[0:0:1j,-range:range:(numpt*1j),-range:range:(numpt*1j)]
    from scipy.interpolate import griddata
    slice = griddata(meshcoor, result.up.vector(), (grid_x, grid_y,grid_z),method="linear")
    slice[np.isnan(slice)]=0
    x1 = np.linspace(-range,range,numpt)
    x2 = np.linspace(-range,range,numpt)
    X1,X2 = np.meshgrid(x1,x2)

    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(14,10))
    gFontSize = 15
    plt.title(title)
    plt.ylabel("y [\AA]",fontsize=gFontSize)
    plt.xlabel('z [\AA]',fontsize=gFontSize)
    plt.pcolormesh( X1.T,X2.T,slice.T, shading='flat' )
    F = plt.gcf()
    F.savefig(fileName)
    print "Plotted %s" % fileName



# validationo of CONSTANT potential with gating 
def ValidateGatedChannel(kon_ext): 
  ## inputs
  a = 1.
  problem.sigma = a
  problem.D = 1.
  Vr = 0
  problem.expnBVr  = np.exp(-parms.beta * Vr)
  
  # from Fig3 
  problem.k_E_ss = kon_ext
  pa_o_pc = 0.01
  problem.L = 5.* a
  
  
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

  # plot 
  import matplotlib.pyplot as plt
  fig = plt.figure(figsize=(14,10))
  gFontSize = 15
  plt.title("Gating ")
  plt.ylabel("$k_G$ [1/Ms]",fontsize=gFontSize)
  plt.xlabel('$V_a$ [kcal/mol]',fontsize=gFontSize)
  plt.plot(Vas,k_cs, 'k-', color='black')
  plt.plot(Vas,k_ss, 'k--', color='black')
  plt.plot(Vas,k_if, 'k.-', color='black')
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

  # plot 
  import matplotlib.pyplot as plt
  fig = plt.figure(figsize=(14,10))
  gFontSize = 25
  plt.title("Linear Potential")
  plt.ylabel("$k_{ss}$ [1/Ms]",fontsize=gFontSize)
  plt.xlabel('$V_0$ [kcal/mol]',fontsize=gFontSize)
  plt.plot(V0s,k_ss, 'k-', color='black')
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
  else:
    print "Using stored values DONT TRUST";
    chargedresult = empty()
    chargedresult.kon = 1.603999e+09

  # results 
  # from (18) of Song 
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
  print "Sphere & %3.1e & %3.e &  & %3.e & %3.e & NA \\\\" % (unchargedresult.kon,chargedresult.kon,k_lp[i],k_gs[i]) 


if __name__ == "__main__":
  msg="smol.py run "



  import sys
  if len(sys.argv) < 2:  
      raise RuntimeError(msg)

  if(sys.argv[1]=="run"):
    ValidationSphere(useStored=0)
    #ValidationTnC(useStored=1)

  else:
    raise RuntimeError(msg)

