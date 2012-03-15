
import smol
import InteriorProblemMaster
import numpy as np
from params import * # must do this for class
import channel_smol as channel
import gating_smol as gating


## boundaries 
#from TroponinBoundary import *
import SercaBoundary as sercaboundaries

#bound = TroponinBoundary()
parms = params()
problem = smol.problem

class empty:pass
boundaries = empty()

def Setup():
    # smol 
    root = "/home/huskeypm/scratch/validation/serca/serca"
    print "WARNING: SERCA mesh is empyt?????"
    problem.fileMesh = root+"_mesh.xml.gz"
    problem.fileSubdomains= root+"_subdomains.xml.gz"
    problem.filePotential= root+"_values.xml.gz"
  
    # pmf 
    problem.filePMF = "/home/huskeypm/sources//dolfin_smol/example/pmf/out.pmf";

    # gating 
    problem.wa=1e9 # rate of forming absorbing state [1/s] 
    problem.wr=1e9 # rate of forming reflective state [1/s]
    problem.p_a=0.5 # probability of being in absorbing state 
    problem.Va = -5 # DEFINE [kcal/mol]
    print "WARNING: for now assuming Va = %f" % problem.Va
    problem.Vr = 0 # DEFINE [kcal/mol]
    print "WARNING: for now assuming Vr = %f" % problem.Vr


    # borrow values from testboundaries
    sercaboundaries.activeSiteLoc = np.array([-1.4566563606e+01,    3.3411479950e+01, -150])
    sercaboundaries.activeSiteR   = 10.0
    sercaboundaries.topZ = 0.0

     
    boundaries.activeSite = sercaboundaries.ActiveSite()
    boundaries.bulkBoundary = sercaboundaries.BulkBoundary()

    return root


def Run(problem,boundaries=0,pvdFileName="up.pvd"):
    print "In debugging mode so skipping computation of smol. WARNING: results are wrong anyway!"
    results = empty()
    results.kon = 1 # for debugging 

    #channelresults = channel.Run(problem,boundaries=boundaries,pvdFileName=pvdFileName) #,results=results)
    channelresults = channel.Run(problem,pvdFileName=pvdFileName,results=results)
    print "WARNING: need to pass channel info to gating (right now using computeChannelTerm"
    gatingresults = gating.Run(problem,result=channelresults)
    
    return gatingresults 


# pmf - refers to output from WHAM 
if __name__ == "__main__":
  msg="serca.py run "

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  elif(sys.argv[1]=="run"):
    root =Setup()


    # wo electro 
    #smol.Run(problem, boundaries=boundaries)
    problem.filePotential= "none"
    #unchargedresults = Run(problem,boundaries=boundaries,pvdFileName="serca_uncharged.pvd")
    unchargedresult = Run(problem,pvdFileName="serca_uncharged.pvd")

    # w electro 
    problem.filePotential= root+"_values.xml.gz"
    print "Skipping electro for now"
    problem.filePotential="none"
    #chargedresults = Run(problem,boundaries=boundaries,pvdFileName="serca_charged.pvd")
    chargedresult = Run(problem,pvdFileName="serca_charged.pvd")

    print "SERCA & %3.1e & %3.e &  & %3.e & %3.e & NA \\\\" % (
      unchargedresult.kon,
      chargedresult.kon,
      chargedresult.kss,
      chargedresult.k_g)

