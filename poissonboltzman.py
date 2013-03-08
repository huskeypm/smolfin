from dolfin import *
import numpy as np
#import modelparameters

#
# Revisions
# Todo
# -validate all aspects of code 
# -implement log function 


msg="""
Purpose:
  Tempate for solving (non)-linear poisson boltzmann equation 

Usage:
  poissonboltzmann.py run 
  
Notes: 
  Guaranteed to be wrong for right now!!

Author:
  the computational scientist formally known as pete 

"""

class params:
  dim = 2  # 2d mesh 
  molRad = 0.5 # radius of molecule [A]
  molMarker = 2 # marker for molecular boundary 
  eps = 0.001  # epsilson for 'error'
  domRad = 5.  # radius of domain [A] (kind of, since square)
  domMarker = 3  # marker for domain boundary 
  z = 1.       # unit charge 
  ec = 1. 	# electric constant [e]
  epsilonExterior = 80. # dielectric constant in exterior []
  kappa = 1/10.  # Inverse dybye length [1/A]
  center = np.zeros(dim)      
  kT = 0.69	# energy [kcal/mol]
  beta = 1/kT

  # ion
  ionC = 1. # ion conc [M]
  ionRad=1. # ion radius [A]

  # modes
  mode = "linear"# linear, nonlinear, finitesize 


class molecularBoundary(SubDomain):
  def inside(self,x,on_boundary):
    result = np.linalg.norm(x-params.center) < (params.eps+params.molRad)
    result = result and on_boundary
    return result      

class domainBoundary(SubDomain):
  def inside(self,x,on_boundary):
    result = np.linalg.norm(x-params.center) > (-params.eps+params.domRad)
    result = result and on_boundary
    return result      


# Create mesh and define function space
def doit(filename):
 
  # Create mesh and define function space
  debug=0
  if(debug):   
    mesh = UnitCube(8,8,8)
  else:
    mesh = Mesh(filename)

  doPB(mesh)

def doPB(mesh):
  V = FunctionSpace(mesh, "Lagrange", 1)

  ## Define boundary conditions
  subdomains = MeshFunction("uint",mesh,params.dim-1)
  bcs = []

  # define BC on molecular boundary 
  boundary = molecularBoundary()
  boundary.mark(subdomains,params.molMarker)
  # PKH: should maybe use sympy later 
  # Analytical solution for linearized PBE for atom of radius R and charge q
  # see Eq 5.1 [1]	
  #M. Holst, N. Baker, and F. Wang, JCC, v21, n15,pp1319-1342, 2000
  q = params.ec * params.z
  k = params.kappa/np.sqrt(params.epsilonExterior)
  u0 = Expression("q/epsilon*R*(1-k*R/(1+k*a))",\
    q=q,k=k,epsilon=params.epsilonExterior,\
    R=params.molRad,a=params.molRad)    
  bcs.append(DirichletBC(V, u0, subdomains,params.molMarker))

  # define BC on domainboundary (potential should be zero at boundary) 
  boundary = domainBoundary()
  boundary.mark(subdomains,params.domMarker)
  u0 = Constant(0.0)
  bcs.append(DirichletBC(V, u0, subdomains,params.domMarker))
  
  
  ## Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)
  
  # exterior
  # LHS  
  form = params.epsilonExterior*inner(grad(u), grad(v))*dx

  # eps*grad(u)grad(v) = kappa^2 uv
  if(params.mode=="linear"):
    form += -1*params.kappa*params.kappa *u*v*dx
  # eps*grad(u)grad(v) = kappa^2 sinh(u)*v
  elif(params.mode=="nonlinear"):
    print "Using nonlinear code"
    def sinh(x=0):
      print "Craptastic implementation. Need to use new version of fenics w sinh" 
      print "Craptastic implementation. Need to Johan's sympy model module"           
      return (exp(x)-exp(-x))/2.

    import ufl
    namespace = ufl.__dict__.copy()
    namespace["sinh"] = sinh
    ufl_exp = eval(str(sinh()), namespace, {"x":u})

    # Might need to cmopile this as UFL code, see pg 196 of dolfin manual
    # https://answers.launchpad.net/dolfin/+question/199782
    # WRONG form += -1*params.kappa*params.kappa *Expression("sinh(u)",u=u)*v*dx
    # WRONG form += -1*params.kappa*params.kappa *f(u)*v*dx

  # eq (7) notes.pdf
  elif(params.mode=="finite"):
    print "WARNING: hyperbolic sine not recognized"
    arg = params.z*params.beta*params.ec*u
    a= params.ionRad
    phi0 = 2*a*a*a*params.ionC     
    rhs1 = 8 *np.pi * params.z*params.ec*params.ionC
    rhs1 *= 1/params.epsilonExterior
    num = np.sinh(arg) 
    denom = (1-phi0+phi0*np.cosh(arg))
    rhs1 *= num/denom
    form += -1*rhs1*v*dx
  else:
    raise RuntimeError("What did you just tell me to do???")
  
  # Compute solution
  print "Solving %s form of PBE" % params.mode
  x = Function(V)
  solve(lhs(form)==rhs(form), x, bcs)
  # F = a
  # problem = NonlinearVariationalProblem(F, x, bcs=bc, J=J)
  # solver = NonlinearVariationalSolver(problem)
  #solver.parameters["linear_solver"] = "gmres"
  #solver.parameters["preconditioner"] = "ilu"
  #solver.solve()


  File("out.pvd") << x
  
  
  ### EXTRAS
  # use for interior problem later
  if(0): 
    # extract boundary values from previous solution and use in a new solution
    bc2 = DirichletBC(V,x,boundary)
    x2 = Function(V)
    solve(a==L,x2,bc2)

  return x

#sphere
if __name__ == "__main__":

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= "sphere2d.xml"
  if(len(sys.argv)==3):
    print "arg"

  for i,arg in enumerate(sys.argv):
    if(arg=="-nonlinear"):
      params.mode="nonlinear"




  doit(fileIn)


