import numpy as np
import smoltodolfin as std
import mcsftodolfin as mtd
import smol
from sys import *
from view import plotslicegeneral



def read_param(paramFilename,apbsFilenames):

  lines    =open(paramFilename).readlines()

  for line in lines:
    spl = line.split()
    if(spl[0] == "D"):
      smol.parms.D = np.float(spl[1])

    elif (spl[0] == "q"):
      smol.parms.q = np.float(spl[1])

    elif (spl[0] == "dx"):
      apbsFilenames.append(spl[1])
      wElectro=1

    else:
     msg = spl[0],"Not understood"
     raise RuntimeError(msg)


  #print dir(smol.parms)
  #print apbsFilenames
  #print wElectro 
  #quit()
   

    


if __name__ == "__main__":

    # ~/localTemp/NBCR/smol/apbs_fe/potential-0.dx
    apbsFilenames=[]
    # ~/localTemp/NBCR/smol/gamer/p.pqr.output.out.m
    mcsfFilename = "none"
    paramFilename = "none"
    wElectro=0; # do not include electrostatic conrtibution  




    import sys
    if len(sys.argv) <  2:
        msg ="""
\n
Purpose:
  Runs smol based on .m molecular mesh

Usage: 
  smol_wrapper.py -mcsf file.m <-param param.txt> <-tarball dxfiles.tgz>

  -mcsf  -  molecular mesh with outer, molecular and active sites marked (5,4, and 1)
  -param - optional parameter file listing .dx files from apbs 
  -tarball - zipped tarball of .dx files 

"""
        raise RuntimeError(msg)

    # parse args
    for i in range(len(sys.argv)):
      if(sys.argv[i] == '-mcsf'):
        mcsfFilename = sys.argv[i+1]

      if(sys.argv[i] == '-param'):
        paramFilename = sys.argv[i+1]

      if(sys.argv[i] == '-tarball'):
        tarball = sys.argv[i+1]
        # tarball is handled by opal, so don't care here 

    # read in param file 
    read_param(paramFilename,apbsFilenames)


    # get problem obhect 
    problem = smol.problem
    problem.root = mcsfFilename.replace(".m","")
    if(len(apbsFilenames)>0):
      wElectro=1 


    # convert mcsf mesh to dolfin format
    debug=0
    if(debug==0):

      # mcsf conversion 
      problem.mesh= mtd.read_and_mark(mcsfFilename)

      # electrostatic potential conversions
      if(wElectro==1):
      #meshFilename = problem.root+"_mesh.xml.gz"
	#      from dolfin import Mesh
	##      problem.mesh = Mesh(meshFilename)
        std.convertAllAPBS(problem, apbsFilenames,writePotentialOnly=1) 

      # filenames  
      root = problem.root 
      problem.fileMesh = root+"_mesh.xml.gz"
      problem.fileSubdomains= root+"_subdomains.xml.gz"

      if(wElectro==1):
        print "Running with electrostatics"
        problem.filePotential= root+"_values.xml.gz"
   
      else:
        print "Running without electrostatics"

      # run 
      smol.Run(problem)

    else:
      1











    


