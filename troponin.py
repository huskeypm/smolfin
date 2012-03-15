
import smol
import InteriorProblemMaster
import numpy as np
from params import * # must do this for class
import channel_smol as channel


## boundaries 
#from TroponinBoundary import *
import TroponinBoundary as tncboundaries

#bound = TroponinBoundary()
parms = params()
problem = smol.problem

class empty:pass
boundaries = empty()

def TnCIsolated(problem):
    root = "/home/huskeypm/scratch/validation/tnc_isolated/tnc_isolated"
    problem.fileMesh = root+"_mesh.xml.gz"
    problem.fileSubdomains= root+"_subdomains.xml.gz"

    # wo electro
    problem.filePotential="none"
    uncharged = channel.Run(problem,pvdFileName="tnc_isolated_uncharged.pvd")

    # w electro 
    problem.filePotential= root+"_values.xml.gz"
    print "Skipping electro for now"
    problem.filePotential="none"

    charged = channel.Run(problem,pvdFileName="tnc_isolated_charged.pvd")
 
    return(uncharged,charged)


def TnCTroponin(problem,boundaries=0):
    root = "/home/huskeypm/scratch/validation/troponin/troponin"
    problem.fileMesh = root+"_mesh.xml.gz"
    problem.fileSubdomains= root+"_subdomains.xml.gz"
    problem.filePotential= root+"_values.xml.gz"
    print "Skipping electro for now"
    problem.filePotential="none"
 

    results = smol.Run(problem,boundaries=boundaries,pvdFileName="troponin.pvd")

    return (results)


# pmf - refers to output from WHAM 
if __name__ == "__main__":
  msg="troponin.py run "

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  elif(sys.argv[1]=="run"):

    # provble.FilePMF="/home/huskeypm/sources/dolfin_smol/example/out.pmf"
    print "WARNING: Need to use correct PMF file"
    problem.filePMF = "/home/huskeypm/sources//dolfin_smol/example/pmf/out.pmf";
    #problem.filePMF = sys.argv[3]

     
    boundaries.activeSite = tncboundaries.ActiveSite()
    boundaries.bulkBoundary = tncboundaries.BulkBoundary()

    (unchargedresultsisolated,chargedresultsisolated) = TnCIsolated(problem)
    #resultstroponin = TnCTroponin(problem,boundaries=boundaries)
    resultstroponin = TnCTroponin(problem)

    print "TnC   & %3.1e & %3.e  & %3.e & NA   & %3.e \\\\" % (
      unchargedresultsisolated.kon,
      chargedresultsisolated.kon,
      chargedresultsisolated.kss,
      resultstroponin.kon)

