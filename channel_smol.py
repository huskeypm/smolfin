
import smol
import interiorProblem
import numpy as np
from params import * # must do this for class

parms = params()
problem = smol.problem

class empty:pass

 
# divided into interior (PMF-driven) and exterior (diffusion driven) problems
# This script primarily deals with interior problem and passes exterior problem to smol 
def RunChannelSmol(problem): 
  ## inputs
  problem.x0=0
  problem.xL=25
  invKappa0 = 0;    # (intrinsic reaction rate)^-1   # invKappa=0 implies IRR is infinitely gast
 
  ## PMF domain 
  print "Assuming values for range of x values over which pmf is defined"
  interiorResult = interiorProblem.Run(problem)
  print "kpmf %f " % interiorResult.invkPMF

  # need to form (2.10) Berez using quantities from interior problem and exterior problem
  # WARNING - need to make sure area of binding tunnel interface on the exterior problem 
  # matches that at the interior problem
  # if we asume steady state, I think we can use (2.8) at x=L for g_1(L,t) in (2.10)
  g1_L0 = interiorResult.sigma_xL * np.exp(-parms.beta*interiorResult.V_xL)   # g1(x=L,t=0)
  # binding Channel Term in (2.10) [external PMF due to electrostatics handled elsewhere) 
  bindChannelTerm = g1_L0 / (interiorResult.sigma_xL * np.exp(-parms.beta*interiorResult.V_xL) )

  problem.bindChannelTerm = bindChannelTerm 


  ## exterior domain 
  exteriorResult = smol.Run(problem,exteriorProblem=1)
  # defined in 4.1 Berez 
  invkE = 1/exteriorResult.kon
  print "kE  %f" % exteriorResult.kon


  ## kon 
  # Berez 4.1, assuming that Kappa0 >> s(0) e(-beta V(0))
  invkss = invkE + invKappa0 + interiorResult.invkPMF
  result = empty()
  result.kss = 1/invkss
  print "kss %f" % result.kss
   
  return result 


# pmf - refers to output from WHAM 
if __name__ == "__main__":
  msg="channel_smol.py -root <file> <pmf>"

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  elif(sys.argv[1]=="-root"):
    root = sys.argv[2]
    problem.fileMesh = root+"_mesh.xml.gz"
    problem.fileSubdomains= root+"_subdomains.xml.gz"
    problem.filePotential= root+"_values.xml.gz"

    # provble.FilePMF="/home/huskeypm/sources/dolfin_smol/example/out.pmf"
    problem.filePMF = "/home/huskeypm/sources//dolfin_smol/example/pmf/out.pmf";
    #problem.filePMF = sys.argv[3]

    RunChannelSmol(problem)

