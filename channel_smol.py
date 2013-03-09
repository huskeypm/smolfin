"
----------------------------------
Smolfin - a numerical solver of the Smoluchowski equation for interesting geometries
Copyright (C) 2012 Peter Kekenes-Huskey, Ph.D., huskeypm@gmail.com

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA


----------------------------------------------------------------------------
"

#
# This script computes kon within binding channel 
#

import smol
import InteriorProblemMaster
import numpy as np
from params import * # must do this for class

#parms = params()
parms = smol.parms
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

    if(exteriorResult.kon < 1):
      print "WARNING: Kon caclc failed"
      quit()
      exteriorResult.kon =1 

  else:
    print "Stored result used instead of copmpuyting smol"
    exteriorResult = results

  print "kPMF %e" % interiorResult.kPMF

  # defined in 4.1 Berez 
  invkE = 1/exteriorResult.kon
  #print "kE  %e" % exteriorResult.kon

  ## kon 
  # Berez 4.1, assuming that Kappa0 >> s(0) e(-beta V(0))
  invkss = invkE + invKappa0 + interiorResult.invkPMF
  
  # technical should be copied to a new object (TODO)
  interiorResult.kon = exteriorResult.kon
  interiorResult.kss = 1/invkss
  print "kss %e" % interiorResult.kss
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

