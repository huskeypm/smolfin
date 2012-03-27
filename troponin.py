
import smol
import InteriorProblemMaster
import numpy as np
from params import * # must do this for class
import channel_smol as channel


## boundaries 
#from TroponinBoundary import *
import TroponinBoundary as troponinboundaries

#bound = TroponinBoundary()
parms = params()
problem = smol.problem

class empty:pass
boundaries = empty()

def TnCIsolated(problem,useStored=0):
    if (useStored==1):
      print "In debugging mode so skipping computation of smol. WARNING: results are wrong anyway!"
      results = empty()
      results.kon = 1 # for debugging 
    else:
      results = 0

    # exterior 
    root = "/home/huskeypm/scratch/validation/tnc_isolated/tnc_isolated"
    problem.fileMesh = root+"_mesh.xml.gz"
    problem.fileSubdomains= root+"_subdomains.xml.gz"

    # interior 
    problem.filePMF = root+".pmf"
    problem.sigma = 4  * np.pi * 10**2  # 10 Ang^2
    problem.D = parms.D
    problem.x0 = 9;
    problem.xL = 25;

    # wo electro
    problem.filePotential="none"
    uncharged = channel.Run(problem,pvdFileName="tnc_isolated_uncharged.pvd",results=results)

    # w electro 
    problem.filePotential= root+"_values.xml.gz"
    #print "Skipping electro for now"
    #problem.filePotential="none"

    charged = channel.Run(problem,pvdFileName="tnc_isolated_charged.pvd",results=results)
 
    return(uncharged,charged)


def TnCTroponin(problem,boundaries=0,useStored=0):
    if (useStored==1):
      print "In debugging mode so skipping computation of smol. WARNING: results are wrong anyway!"
      results = empty()
      results.kon = 1 # for debugging 
      return results
    else:
      results = 0

    root = "/home/huskeypm/scratch/validation/troponin/troponin"
    problem.fileMesh = root+"_mesh.xml.gz"
    problem.fileSubdomains= root+"_subdomains.xml.gz"
    problem.filePotential= root+"_values.xml.gz"
    print "Skipping electro for now"
    problem.filePotential="none"
    results = smol.Run(problem,boundaries=boundaries,pvdFileName="troponin.pvd")

    return (results)


def Validation(useStored=0):

    ## isolated
    (unchargedresultsisolated,chargedresultsisolated) = TnCIsolated(problem,useStored=useStored)

    ## troponin
    # see 120210_troponin.tex
    troponinboundaries.activeSiteLoc = np.array([-30.345, 39.371,216.75])       
    troponinboundaries.activeSiteR   = 10.0
    troponinboundaries.outerR = 690 # based on paraview

    boundaries.activeSite = troponinboundaries.ActiveSite()
    boundaries.bulkBoundary = troponinboundaries.BulkBoundary()

    resultstroponin = TnCTroponin(problem,boundaries=boundaries,useStored=useStored)
    #resultstroponin = TnCTroponin(problem,useStored=useStored)

    print "TnC   & %3.1e & %3.1e  & %3.1e & NA   & %3.1e \\\\" % (
      unchargedresultsisolated.kon,
      chargedresultsisolated.kon,
      chargedresultsisolated.kss,
      resultstroponin.kon)



# pmf - refers to output from WHAM 
if __name__ == "__main__":
  msg="troponin.py run "

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  elif(sys.argv[1]=="run"):
      Validation()

