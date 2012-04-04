
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
#$parms = params()
import smol
parms = smol.parms





problem = smol.problem

class empty:pass
boundaries = empty()

def Setup():
    # smol 
    root = "/home/huskeypm/scratch/validation/serca/serca"
    problem.fileMesh = root+"_mesh.xml.gz"
    problem.fileSubdomains= root+"_subdomains.xml.gz"
    problem.filePotential= root+"_values.xml.gz"
  
    # pmf 
    problem.filePMF = root+".pmf"
    problem.diameter = 10 # AA
    problem.channelR = problem.diameter
    problem.D = parms.D
    problem.x0 = 2;
    problem.xL = 25;

    # gating 
    problem.wa=0.67e9 # rate of forming absorbing state [1/s] 
    problem.wr=3.23e9 # rate of forming reflective state [1/s]
    problem.p_a=0.67e9 / 3.23e9; # probabilityy of being in absorbing state (?far from protein) 
    problem.Va = -1.74 # Potential felt near open state ([kcal/mol]
    problem.Vr = 0 # DEFINE [kcal/mol]
    print "Assuming Vr = %f" % problem.Vr


    # borrow values from testboundaries
    # compare this with value used in APBS calcs
    #sercaboundaries.activeSiteLoc = np.array([-1.4566563606e+01,    3.3411479950e+01, -150])
    # active site in pqr file - gridshift (see 120306_mesing.tx)
    # IN SERCABOUNDARY sercaboundaries.activeSiteLoc = np.array([20.730,-27.038 , 8.204 - 113])
    # IN SERCABOUNDARY sercaboundaries.activeSiteR   = 10.0
    # IN SERCABOUNDARY sercaboundaries.topZ = 0.0

     
    boundaries.activeSite = sercaboundaries.ActiveSite()
    boundaries.bulkBoundary = sercaboundaries.BulkBoundary()
    #molecularBoundary = bound.MolecularBoundary() # i think we can ignore, since zero anyway

    return root


def Run(problem,boundaries=0,pvdFileName="up.pvd",useStored=0):
    if (useStored==1):
      print "In debugging mode so skipping computation of smol. WARNING: results are wrong anyway!"
      results = empty()
      results.kon = 1e9 # for debugging 
    else:
      results = 0

    channelresults = channel.Run(problem,boundaries=boundaries,pvdFileName=pvdFileName,results=results)
    #channelresults = channel.Run(problem,pvdFileName=pvdFileName ,results=results)
    print "WARNING: need to pass channel info to gating (right now using computeChannelTerm"
    gatingresults = gating.Run(problem,result=channelresults)
    
    return gatingresults 

# validation example for SERCA 
def Validation(useStored=0):
    root = Setup()

    ## wo electro 
    #smol.Run(problem, boundaries=boundaries)
    problem.filePotential= "none"
    #unchargedresults = Run(problem,boundaries=boundaries,pvdFileName="serca_uncharged.pvd")
    unchargedresult = Run(problem,pvdFileName=root+"_uncharged.pvd",useStored=useStored)

    ## w electro 
    problem.filePotential= root+"_values.xml.gz"
    #print "Skipping electro for now"
    #problem.filePotential="none"
    chargedresult = Run(problem,boundaries=boundaries,pvdFileName=root+"_charged.pvd",useStored=useStored)
    #chargedresult = Run(problem,pvdFileName="serca_charged.pvd",useStored=useStored)

    ## with uncharged lipids
    unchgRoot  = "/home/huskeypm/scratch/validation/serca_no_lipid_chg/serca_no_lipid_chg"
    problem.filePotential= unchgRoot+"_values.xml.gz"
    nolipidchargedresult = Run(problem,boundaries=boundaries,pvdFileName=unchgRoot+".pvd",useStored=useStored)

    msg=[]
    m = "SERCA & %3.1e & %3.1e & %3.1e & %3.1e & %3.1e \\\\" % (
      unchargedresult.kon,
      chargedresult.kon,
      chargedresult.kss,
      chargedresult.k_g,
      nolipidchargedresult.kon)
    msg.append(m)

    scale = 1.e9
    m = "SERCA(%1e) & %3.1e & %3.1e &  & %3.1e & %3.1e & %3.1e \\\\" % (
      scale,
      unchargedresult.kon/scale,
      chargedresult.kon/scale,
      chargedresult.kss/scale,
      chargedresult.k_g/scale,
      nolipidchargedresult.kon/scale)
    msg.append(m)

    return msg


# pmf - refers to output from WHAM 
if __name__ == "__main__":
  msg="serca.py run "

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  elif(sys.argv[1]=="run"):
    Validation(useStored=0)


