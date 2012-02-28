
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
def PMFTerm(fileIn,xrange,problem):
  # get data
  dat = numpy.loadtxt(fileIn);
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

  
def InteriorProblem(fileIn,x0,xL):
  # prepare integrand 
  problem = empty() 
  xrange = np.array([x0,xL])
  PMFTerm(fileIn,xrange,problem)
  
  problem.D_x = DiffusionConst(x)
  problem.Sigma_x = Sigma(x)
  
  integ = problem.boltz / (problem.D_x * problem.Sigma_x)
  
  # integrate 
  from numpy import trapz 
  # trapz(y,x)
  results = empty()
  results.area = trapz(integ,x)
  
  # provide Boltzman distribution at boundary x=L (channel mouth location)
  # note: be sure that this is a boltzman fact, as is the quantity in PMFTerm
  results.prob_x = exp(-beta * problem.V_x)
  results.prob_xL= results.prob_x[ xL ]
  results.sigma_xL= problem.Sigma_x[ xL ]
  results.diff_xL= problem.Diff_x[ xL ]
 
  return (problem,results)
  

if __name__ == "__main__":
    # inputs 
    # pmf values
    #fileIn="/net/home/huskeypm/localTemp/serca/120222/out.pmf"
    fileIn="/home/huskeypm/sources/dolfin_smol/example/out.pmf"
    # range over which pmf is defined (allows discarding of points)
    # where active site is located
    x0 = 0
    # where interfac is located
    xL = 25        
    import sys
    if len(sys.argv) != 2:
        raise RuntimeError("expected blahlah file as second argument")
    filename = sys.argv[1]
    InteriorProblem(fileIn,x0,xL)

