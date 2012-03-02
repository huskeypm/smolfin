
import numpy as np
from params import *

parms = params()


class empty:pass

# for now assuming the diffusion constant is constant
def DiffusionConst(x):
  diff_const = 1.0
  diff_x = np.ones( np.size(x)) * diff_const
  return(diff_x)

# area along reaction coordinate 
def Sigma(x):
  # R = 1 [nm]
  R = 1
  area = 4 * np.pi * R*R
  sigma_x = np.ones( np.size(x)) * area 
  return(sigma_x)

## PMF 
def PMFTerm(xrange,problem):
  # get data
  dat = np.loadtxt(problem.filePMF);
  x = dat[:,0]
  y = dat[:,1]

  # relevant ranges 
  inds = (x > xrange[0]) &  (x < xrange[1])
  inds = inds[:] 
  V_x = y[inds]
  x = x[inds]
  
  problem.x = x
  problem.V_x = V_x

  boltz = np.exp(parms.beta * V_x) 
  problem.boltz = boltz

  
# Based on integral term in (4.1) Berez ref 
def Run(problem):
  # prepare integrand 
  xrange = np.array([problem.x0,problem.xL])
  PMFTerm(xrange,problem)
  
  problem.D_x = DiffusionConst(problem.x)
  problem.Sigma_x = Sigma(problem.x)
  
  integ = problem.boltz / (problem.D_x * problem.Sigma_x)
  
  # integrate 
  from numpy import trapz 
  # trapz(y,x)
  results = empty()
  # for consistency, defining kpmf as 1/IntegralTerm in 4.1
  integral = trapz(integ,problem.x)
  results.invkPMF= integral
  print "Interior integral %f " % integral
  
  # provide Boltzman distribution at boundary x=L (channel mouth location)
  # note: be sure that this is a boltzman fact, as is the quantity in PMFTerm
  results.prob_x = np.exp(-parms.beta * problem.V_x)
  xL = len(problem.x)-1
  results.prob_xL= results.prob_x[ xL ]
  results.sigma_xL= problem.Sigma_x[ xL ]
  results.D_xL= problem.D_x[ xL ]
  results.V_xL= problem.V_x[ xL ]
 
  return (results)
  

if __name__ == "__main__":
    problem = empty() 

    # inputs 
    import sys
    if len(sys.argv) != 2:
        raise RuntimeError("expected blahlah file as second argument")
    # pmf values
    #filePMF ="/net/home/huskeypm/localTemp/serca/120222/out.pmf"
    #filePMF ="/home/huskeypm/sources/dolfin_smol/example/out.pmf"
    # range over which pmf is defined (allows discarding of points)
    # where active site is located
    problem.x0 = 0
    # where interfac is located
    problem.xL = 25        
    problem.filePMF = sys.argv[1]

    Run(problem)  

