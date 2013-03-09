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
import smol
import numpy as np
from params import * # must do this for class

class empty:pass


#parms = params()
import smol
parms = smol.parms

problem = smol.problem


# I think I should be using 
# - 3.11d of Barreda for lig. indep gating 
# - 3.24a/b for lig dep gating 
# for in between these regimes SOmehow need to compute kE generally, but should be able to evoke steady state somewhere and assume fast/slow gating  limits 

# inputs
parms.Va = 1 # DEFINE [kcal/mol] 
parms.Veff=1 # DEFINE [kcal/mol]
parms.wa=1 # rate of forming absorbing state [1/s] 
parms.wr=1 # rate of forming reflective state [1/s]
parms.area=1 # i think this should come fmor our model 
#parms.L = 1 # TODO - only relevant for gated, binding channel. Probably need to adjust problem to reflect this. 

def RunStochasticSmol(problem):
# compute outer problem 
  exteriorResult = smol.Run(problem)
# give lig induced, fast/slow gatining 

  # using equations 3.24a & b 
  # assuming kappa_0 --> inf, so 2nd term is 0
  # assuming L-> 0 (e.g. gate is at mouth)
  kE = exteriorResult.kon
  pa = wa/(wa+wr)
  # 3.24a reduces to  1/kss = pa*1/kE
  kssslow = kE/pa;
  print "Slow gating %f " % kssslow
 
  # 3.24b reduces to 1/kss = 1/kE  
  kssfast = kE
  print "Fast gating %f " % kssfast
  



if __name__ == "__main__":
  msg="stochastic_smol.py -root <file>"

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  elif(sys.argv[1]=="-root"):
    root = sys.argv[2]
    problem.fileMesh = root+"_mesh.xml.gz"
    problem.fileSubdomains= root+"_subdomains.xml.gz"
    problem.filePotential= root+"_values.xml.gz"

    RunStochasticSmol(problem)


