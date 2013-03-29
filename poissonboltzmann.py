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
import sys
sys.path.append("/home/huskeypm/sources/jay/")
import poissonboltzmannSpecial as pbs
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
  poissonboltzmann.py -linear/-nonlinear/-finite/-wholedomain
  
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
  domRad = 5.  # radius of domain [A] (kind of, since square)
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
  #mode = "nonlinear"
  #mode = "finite"    


  #kappa = 1/10.  # Inverse dybye length [1/A] (should be determined by ionic strengths directly) 
  beta = 1/kT
  ikappa  = 0.304 / np.sqrt(ionC) # Debye length [A], "Intermoplecular and surface forces, Israelachvili" 
  kappa = 1/ikappa

# share params 
pbs.params = params

class molecularBoundary(SubDomain):
  def inside(self,x,on_boundary):
    result = np.linalg.norm(x-params.center) < (params.epsError+params.molRad)
    result = result and on_boundary
    return result      

class domainBoundary(SubDomain):
  def inside(self,x,on_boundary):
    result = np.linalg.norm(x-params.center) > (-params.epsError+params.domRad)
    result = result and on_boundary
    #print x
    #print result
    return result      

# Create mesh and define function space
def doOuterDomainPB(filename):
 
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
    form = pbs.Nonlinear()

  # eq (7) notes.pdf
  elif(params.mode=="finite"):
    form = pbs.Finite()

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

  fileIn= "./example/sphere/sphere2d.xml"
  mode = "outer"
  if(len(sys.argv)==3):
    print "arg"

  for i,arg in enumerate(sys.argv):
    if(arg=="-nonlinear"):
      params.mode="nonlinear"
    if(arg=="-finite"):
      params.mode="finite"
    if(arg=="-wholedomain"):
      mode = "wholedomain"
      fileIn= "./example/sphere/sphere_2d_entire.xml"
    if(arg=="-linear"):
      params.mode = "linear"
      fileIn= "./example/sphere/sphere2d.xml"


  if(mode=="wholedomain"):
    pbs.domainBoundary = domainBoundary
    pbs.molecularBoundary = molecularBoundary
    pbs.doWholeDomainPB(fileIn)
  else: 
    doOuterDomainPB(fileIn)


