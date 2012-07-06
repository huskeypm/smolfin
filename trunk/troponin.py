
import smol
import InteriorProblemMaster
import numpy as np
from params import * # must do this for class
import channel_smol as channel


## boundaries 
#from TroponinBoundary import *
import TroponinBoundary as troponinboundaries

#bound = TroponinBoundary()
#parms = params()
import smol
parms = smol.parms

problem = smol.problem

class empty:pass
boundaries = empty()

def TnCIsolated(problem,useStored=0):
    if (useStored==1):
      print "In debugging mode so skipping computation of smol. WARNING: results are wrong anyway!"
      results = empty()
      results.kon = 2.2e10 # for debugging 
    else:
      results = 0

    # exterior 
    root = "/home/huskeypm/scratch/validation/tnc_isolated/tnc_isolated"
    problem.fileMesh = root+"_mesh.xml.gz"
    problem.fileSubdomains= root+"_subdomains.xml.gz"

    # interior 
    problem.filePMF = root+".pmf"
    print "%s" %("NEED TO USE DIFF. SIGMA DEFINITION FOR TNC - DONT TRUST RESULTS UNTIL THEN") 
    problem.channelR = 5  # 10 Ang^2
    problem.D = parms.D
    problem.Dchannel = parms.Dchannel
    problem.x0 = 2;   # In binding site (Dist [AA] between E76 and Ca)
    problem.x0 = 3.6;   # In binding site (Dist [AA] between E76 and Ca)
    problem.xL = 4.3; # Just before rapid increase in PMF (5.5) 

    # wo electro
    problem.filePotential="none"
    uncharged = channel.Run(problem,pvdFileName=root+"_uncharged.pvd",results=results)

    # w electro 
    problem.filePotential= root+"_values.xml.gz"
    #print "Skipping electro for now"
    #problem.filePotential="none"

    charged = channel.Run(problem,pvdFileName=root+"_charged.pvd",results=results)
 
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
    noElectro=0


    #chg 
    results= smol.Run(problem,boundaries=boundaries,pvdFileName=root+".pvd")
    # unchg
    problem.filePotential="none"
    results_unchg = smol.Run(problem,boundaries=boundaries,pvdFileName=root+"_unchg.pvd")

    return (results)


def Validation(useStored=0):

    ## isolated
    (unchargedresultsisolated,chargedresultsisolated) = TnCIsolated(problem,useStored=useStored)

    ## troponin
    # see 120210_troponin.tex
    # in TroponinBoundary troponinboundaries.activeSiteLoc = np.array([-30.345, 39.371,216.75])       
    # in TroponinBoundary troponinboundaries.activeSiteR   = 10.0
    # in TroponinBoundary troponinboundaries.outerR = 450 # based on paraview, plus a little less 

    boundaries.activeSite = troponinboundaries.ActiveSite()
    boundaries.bulkBoundary = troponinboundaries.BulkBoundary()
    #molecularBoundary = bound.MolecularBoundary() # i think we can ignore, since zero anyway

    resultstroponin = TnCTroponin(problem,boundaries=boundaries,useStored=useStored)
    #resultstroponin = TnCTroponin(problem,useStored=useStored)

    scale = 1.e9
    msg = []

    m = "TnC(%1e)   & %7.5f & %7.5f  & %7.5f & NA   & %7.5f \\\\" % (
      scale,
      unchargedresultsisolated.kon/scale,
      chargedresultsisolated.kon/scale,
      chargedresultsisolated.kss/scale,
      resultstroponin.kon/scale)
    msg.append(m) 
    m = "TnC   & %3.1e & %3.1e  & %3.1e & NA   & %3.1e \\\\" % (
      unchargedresultsisolated.kon,
      chargedresultsisolated.kon,
      chargedresultsisolated.kss,
      resultstroponin.kon)
    msg.append(m) 
    return msg
 
    



# pmf - refers to output from WHAM 
if __name__ == "__main__":
  msg="troponin.py run "

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  elif(sys.argv[1]=="run"):
      Validation()

