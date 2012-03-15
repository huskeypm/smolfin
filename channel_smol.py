
#
# This script computes kon within binding channel 
#

import smol
import InteriorProblemMaster
import numpy as np
from params import * # must do this for class

parms = params()
problem = smol.problem

class empty:pass

 
# divided into interior (PMF-driven) and exterior (diffusion driven) problems
# This script primarily deals with interior problem and passes exterior problem to smol 
# useDefault - see 'InteriorProblemMaster.Run()'
def Run(problem,boundaries=0,pvdFileName="up.pvd",useDefault=1,results=0): 

  invKappa0 = 0;    # (intrinsic reaction rate)^-1   # invKappa=0 implies IRR is infinitely gast
  ## xterior domain 
  interiorResult = InteriorProblemMaster.Run(problem)

  ## exterior domain 
  if(results==0):
    exteriorResult = smol.Run(problem,boundaries=boundaries,pvdFileName=pvdFileName)
  else:
    print "Stored result used instead of copmpuyting smol"
    exteriorResult = results

  # defined in 4.1 Berez 
  invkE = 1/exteriorResult.kon
  print "kE  %f" % exteriorResult.kon


  ## kon 
  # Berez 4.1, assuming that Kappa0 >> s(0) e(-beta V(0))
  invkss = invkE + invKappa0 + interiorResult.invkPMF
  
  # technical should be copied to a new object (TODO)
  interiorResult.kon = exteriorResult.kon
  interiorResult.kss = 1/invkss
  print "kss %f" % interiorResult.kss
  print "xL %f " % interiorResult.xL
   
  return interiorResult 


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

    Run(problem)

