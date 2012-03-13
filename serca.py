
import smol
import interiorProblem
import numpy as np
from params import * # must do this for class
import channel_smol as channel


## boundaries 
#from TroponinBoundary import *
import SercaBoundary as sercaboundaries

#bound = TroponinBoundary()
parms = params()
problem = smol.problem

class empty:pass
boundaries = empty()


# pmf - refers to output from WHAM 
if __name__ == "__main__":
  msg="channel_smol.py run "

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  elif(sys.argv[1]=="run"):
    root = "/home/huskeypm/scratch/serca_mole/test2"
    problem.fileMesh = root+"_mesh.xml.gz"
    problem.fileSubdomains= root+"_subdomains.xml.gz"
    problem.filePotential= root+"_values.xml.gz"
    problem.filePotential= "none"

    # provble.FilePMF="/home/huskeypm/sources/dolfin_smol/example/out.pmf"
    problem.filePMF = "/home/huskeypm/sources//dolfin_smol/example/pmf/out.pmf";
    #problem.filePMF = sys.argv[3]

    # borrow values from testboundaries
    sercaboundaries.activeSiteLoc = np.array([-1.4566563606e+01,    3.3411479950e+01, -150])
    sercaboundaries.activeSiteR   = 10.0
    sercaboundaries.topZ = 0.0


     
    boundaries.activeSite = sercaboundaries.ActiveSite()
    boundaries.bulkBoundary = sercaboundaries.BulkBoundary()
    #smol.Run(problem, boundaries=boundaries)
    channel.RunChannelSmol(problem)

