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
from dolfin import *
import numpy as np
#import modelparameters

#
# Revisions
# Todo
# -validate all aspects of code 
#  The sinh expression initially comes out as zero, so it never looks like it is used in the form expression 
# -implement log function 
# -wholeDom prob. I think is behaving reasonably (need to visualize the potential in the exterior domain at a different dynamic range [0-10] than the interior domain)


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


print "WARNING: This code is no-where close to being validated (use at your own risk)"

class params:
  dim = 2  # 2d mesh 
  molRad = 1.5 # radius of molecule [A]
  molMarker = 2 # marker for molecular boundary 
  epsError = 0.001  # epsilson for 'error'
  domRad = 10.  # radius of domain [A] (kind of, since square)
  domMarker = 3  # marker for domain boundary 
  
  z = 1.       # unit charge 
  ec = 8.854187817e-12 # electric constant [C^2/(Jm)]
  M_TO_ANG = 1e-10
  J_TO_KCAL = 0.000239005736
  ec = ec / (M_TO_ANG * J_TO_KCAL) # ec [C^2/kcal A]
  ec = 1.0
  epsilonExterior = 80. # dielectric constant in exterior []
  center = np.zeros(dim)      
  kT = 0.59	# energy [kcal/mol]

  # ion
  ionC = .150 # ion conc [M]
  ionRad=2. # ion radius [A]

  # modes
  mode = "linear"# linear, nonlinear, finitesize 
  mode = "nonlinear"
  mode = "finite"    


  #kappa = 1/10.  # Inverse dybye length [1/A] (should be determined by ionic strengths directly) 
  beta = 1/kT
  ikappa  = 0.304 / np.sqrt(ionC) # Debye length [A], "Intermoplecular and surface forces, Israelachvili" 
  kappa = 1/ikappa


class molecularBoundary(SubDomain):
  def inside(self,x,on_boundary):
    result = np.linalg.norm(x-params.center) < (params.epsError+params.molRad)
    result = result and on_boundary
    return result      

class domainBoundary(SubDomain):
  def inside(self,x,on_boundary):
    result = np.linalg.norm(x-params.center) > (-params.epsError+params.domRad)
    result = result and on_boundary
    return result      

# exterior domain 
class OmegaOutside(SubDomain):
    def inside(self, x, on_boundary):
      result = np.linalg.norm(x-params.center) > (-params.epsError+params.molRad)
      #print "outside"
      #print result 
      return result 

# interior domain 
class OmegaInside(SubDomain):
    def inside(self, x, on_boundary):
      result = np.linalg.norm(x-params.center) <=(-params.epsError+params.molRad)
      #print x
      #print params.center
      #print "inside"
      #print result 
      return result 

## special functions 
def sinh(x=0):
 print "Craptastic implementation. Need to use new version of fenics w sinh/Johans sympy model module"  
 return (exp(x)-exp(-x))/2.

def cosh(x=0):
 print "Craptastic implementation. Need to use new version of fenics w sinh/Johans sympy model module"  
 return (exp(x)+exp(-x))/2.


# Create mesh and define function space
def doit(filename):
 
  # Create mesh and define function space
  debug=0
  if(debug):   
    mesh = UnitCube(8,8,8)
  else:
    mesh = Mesh(filename)

  doPB(mesh)

def doWholeDomainPB():
  scale = 2*params.domRad 
  ## mesh/functions 
  square=0 
  if square:
    mesh = UnitSquare(50,50)
    mesh.coordinates()[:] = scale * mesh.coordinates()
  else:
    # domain assigment doesn't seem to work correctly (fails on 'choose' funtion 
    mesh = Mesh("sphere_2d_entire.xml") 

  dim = 2
  V = FunctionSpace(mesh, "DG", 0)
  eps = Function(V)
  kappa= Function(V)


  ## point source 
  print "Point source is not functional"
  params.center = scale*np.array([0.5,0.5])
  p = Point(params.center[0],params.center[1])  
  #c = PointSource(V,p,magnitude=1) 
  q1 = 1. # [e]
  c = PointSource(V,p,q1)                 

  ## boundary/domain conditions 
  subdomains = MeshFunction('uint', mesh, dim)
  subdomainBoundary = MeshFunction('uint', mesh, dim-1)
  bcs = []
  boundary = domainBoundary()
  boundary.mark(subdomainBoundary,params.domMarker)
  u0 = Constant(0.0)
  bcs.append(DirichletBC(V, u0, subdomainBoundary,params.domMarker))

  ## discontinuous epsilon 
  # http://fenicsproject.org/documentation/tutorial/materials.html#working-with-two-subdomains
  subdomain0 = OmegaOutside()
  subdomain0.mark(subdomains, 0)
  subdomain1 = OmegaInside()
  subdomain1.mark(subdomains, 1)
 
  eps_values = [params.epsilonExterior, 1.]  # values of k in the two subdomains
  kappa_values = [params.epsilonExterior*params.kappa**2, params.epsError]  # values of k in the two subdomains
  #for cell_no in range(len(subdomains.array())):
  #  subdomain_no = subdomains.array()[cell_no]
  # eps.vector()[cell_no] = eps_values[subdomain_no]
  help = np.asarray(subdomains.array(), dtype=np.int32)
  eps.vector()[:] = np.choose(help, eps_values)

  kappa.vector()[:] = np.choose(help, kappa_values)
  #xkappa = eps


  ## Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)
  form = -1*eps*inner(grad(u), grad(v))*dx
#  form += eps *v*dx

  # eps*grad(u)grad(v) = kappa^2 uv
  form += -1*kappa *u*v*dx

  d= 1.0
  factor = 4 * np.pi * params.ec*params.ec / (params.kT)
  norm = 1/(d * np.sqrt(np.pi*2))
  form += Expression("factor*norm *exp(-( (x[0]-xC)*(x[0]-xC) + (x[1]-yC)*(x[1]-yC))/(2*d*d))",\
                      xC=params.center[0],yC=params.center[1],d=d,norm=norm,factor=factor)*v*dx
  A = lhs(form)
  b = rhs(form)

  # apply point charge 
   # not correct 
  #c.apply(b)
  

  x = Function(V)
  solve(A==b, x, bcs)
  # F = a
  # problem = NonlinearVariationalProblem(F, x, bcs=bc, J=J)
  # solver = NonlinearVariationalSolver(problem)
  #solver.parameters["linear_solver"] = "gmres"
  #solver.parameters["preconditioner"] = "ilu"
  #solver.solve()


  File("out.pvd") << x
  

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

  ## define form functions 
  import ufl
  namespace = ufl.__dict__.copy()
  
  # exterior
  # LHS  
  form = -1*params.epsilonExterior*inner(grad(u), grad(v))*dx

  # eps*grad(u)grad(v) = kappa^2 uv
  if(params.mode=="linear"):
    form += -1*params.epsilonExterior*params.kappa*params.kappa *u*v*dx
  # eps*grad(u)grad(v) = kappa^2 sinh(u)*v
  elif(params.mode=="nonlinear"):
    print "Idiot, your weak form is non-linear in u!"
    quit()
    print "Using nonlinear code"
    namespace["sinh"] = sinh
    sinh_expr = eval(str(-1*params.epsilonExterior*params.kappa*params.kappa * sinh()), namespace, {"x":u})
    print "sinh_expr %s" % sinh_expr

    print "WARNING: not sure if ufl exp is reflected in the form (need to compare linear w non-linear"
    form += sinh_expr*v*dx
    print form 

    # Might need to cmopile this as UFL code, see pg 196 of dolfin manual
    # https://answers.launchpad.net/dolfin/+question/199782
    # WRONG form += -1*params.kappa*params.kappa *Expression("sinh(u)",u=u)*v*dx
    # WRONG form += -1*params.kappa*params.kappa *f(u)*v*dx

  # eq (7) notes.pdf
  elif(params.mode=="finite"):
    print "Idiot, your weak form is non-linear in u!"
    quit()
    print "WARNING: hyperbolic sine not recognized"
    arg = params.z*params.beta*params.ec*u

    sinh_expr = eval(str(sinh()), namespace, {"x":arg})
    print "sinh_expr %s" % sinh_expr
    cosh_expr = eval(str(cosh()), namespace, {"x":arg})
    print "cosh_expr %s" % cosh_expr

    a= params.ionRad
    phi0 = 2*a*a*a*params.ionC     
    rhs1 = 8 *np.pi * params.z*params.ec*params.ionC
    rhs1 *= 1/params.epsilonExterior
    num = sinh_expr      
    denom = (1-phi0+phi0*cosh_expr)      
    rhs1 *= num/denom
    form += -1*rhs1*v*dx
    print form 
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
    if(arg=="-wholedom"):
      doWholeDomainPB()
      quit()


   



  doit(fileIn)


