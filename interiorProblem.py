
import numpy as np

beta = 0.6 # 1/kT

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
  area = 4 * pi * R^R
  sigma_x = np.ones( np.size(x)) * area 
  return(sigma_x)

## PMF 
def PMFTerm(xrange,problem):
  # get data
  dat = numpy.loadtxt(problem.filePMF);
  x = dat[:,0]
  y = dat[:,1]

  # relevant ranges 
  inds = (x > xrange[0]) &  (x < xrange[1])
  inds = inds[:] 
  V_x = y[inds]
  x = x[inds]
  
  problem.x = x
  problem.V_x = V_x

  boltz = exp(beta * V_x) 
  problem.boltz = boltz

  
# Based on integral term in (4.1) Berez ref 
def Run(problem):
  # prepare integrand 
  xrange = np.array([problem.x0,problem.xL])
  PMFTerm(xrange,problem)
  
  problem.D_x = DiffusionConst(x)
  problem.Sigma_x = Sigma(x)
  
  integ = problem.boltz / (problem.D_x * problem.Sigma_x)
  
  # integrate 
  from numpy import trapz 
  # trapz(y,x)
  results = empty()
  # for consistency, defining kpmf as 1/IntegralTerm in 4.1
  integral = trapz(integ,x)
  results.invkPMF= integral
  
  # provide Boltzman distribution at boundary x=L (channel mouth location)
  # note: be sure that this is a boltzman fact, as is the quantity in PMFTerm
  results.prob_x = exp(-beta * problem.V_x)
  results.prob_xL= results.prob_x[ xL ]
  results.sigma_xL= problem.Sigma_x[ xL ]
  results.diff_xL= problem.Diff_x[ xL ]
 
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

