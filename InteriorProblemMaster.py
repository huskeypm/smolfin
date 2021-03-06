"""
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
"""

#
# This script sets up relevant terms for PMF term describining binding channel 
#

import numpy as np
from params import *
from DefaultInteriorProblem import *

#parms = params()
import smol
parms = smol.parms




class empty:pass


# Based on integral term in (4.1) Berez ref 
# If InteriorObject!=0, be sure object has defined
#   problem.D_x      - diffusion constant (x) 
#   problem.Sigma_x  - area of pore perpendicular to x (x) 
#   problem.x = x    - x values (x) 
#   problem.V_x      - PMF values 
#
def Run(problem,InteriorObject=0):
  if(InteriorObject==0):
      interProb = DefaultInteriorProblem(problem.filePMF,channelR=problem.channelR,diff_const=problem.Dchannel,x0=problem.x0,xL=problem.xL)
  else:
      interProb = InteriorObject
  
  # init 
  # TODO can we do this from inside class? 
  interProb.D_x = interProb.DiffusionConst()
  interProb.Sigma_x = interProb.Sigma()


  # prepare integrand 
  #print "V_x"
  #print interProb.V_x
  boltz = np.exp(parms.beta * interProb.V_x) 
  interProb.boltz = boltz
  #print "Boltz"
  #print interProb.boltz
  #print interProb.D_x   
  #print interProb.Sigma_x
  integ = interProb.boltz / (interProb.D_x * interProb.Sigma_x)
  #print "integrand"
  #print integ
  
  # integrate 
  from numpy import trapz 
  # trapz(y,x)
  results = empty()
  # for consistency, defining kpmf as 1/IntegralTerm in 4.1
  integral = trapz(integ,interProb.x)
  # 120401 - verified w maxima
  #print "%f boltzz only " % trapz(np.exp(parms.beta * interProb.V_x),interProb.x)
  #integral = trapz(np.ones(np.size(interProb.x)),interProb.x)
  #print "x " 
  #print interProb.x
  #print "Interior integral %e " % integral

  kpmf = 1/integral * parms.um3_to_invM
  #print "rescaling into units of Ms %e" % kpmf 
  results.invkPMF= 1/kpmf     


  
  # provide Boltzman distribution at boundary x=L (channel mouth location)
  # note: be sure that this is a boltzman fact, as is the quantity in LoadPMF
  results.prob_x = np.exp(-parms.beta * interProb.V_x)
  xL = len(interProb.x)-1
  results.temp = integral
  results.kPMF = kpmf
  results.xL = interProb.x[ xL ]
  results.x0 = interProb.x[ 0 ]
  results.prob_xL= results.prob_x[ xL ]
  results.sigma_xL= interProb.Sigma_x[ xL ]
  results.D_xL= interProb.D_x[ xL ]
  results.V_xL= interProb.V_x[ xL ]
  results.D_x = interProb.D_x
  results.Sigma_x = interProb.Sigma_x


  #print "quit for debug"
  #quit()
 
  return (results)
  
# Interior domain problem 
def RunWRONG(problem,InteriorObject=0):
  ## inputs
 
  ## PMF domain 
  interiorResult = Run(problem,InteriorObject=InteriorObject)
  #print "kpmf %f " % interiorResult.invkPMF

  # need to form (2.10) Berez using quantities from interior problem and exterior problem
  # WARNING - need to make sure area of binding tunnel interface on the exterior problem 
  # matches that at the interior problem
  # if we asume steady state, I think we can use (2.8) at x=L for g_1(L,t) in (2.10)
  g1_L0 = interiorResult.sigma_xL * np.exp(-parms.beta*interiorResult.V_xL)   # g1(x=L,t=0)
  # binding Channel Term in (2.10) [external PMF due to electrostatics handled elsewhere) 
  bindChannelTerm = g1_L0 / (interiorResult.sigma_xL * np.exp(-parms.beta*interiorResult.V_xL) )

  problem.bindChannelTerm = bindChannelTerm 



if __name__ == "__main__":
    problem = empty() 

    # inputs 
    import sys
    if len(sys.argv) != 2:
        raise SetuptimeError("expected blahlah file as second argument")
    # pmf values
    #filePMF ="/net/home/huskeypm/localTemp/serca/120222/out.pmf"
    #filePMF ="/home/huskeypm/sources/dolfin_smol/example/out.pmf"
    problem.filePMF = sys.argv[1]

    Run(problem)  

