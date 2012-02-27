#
# Simple example that imports mesh from APBS
# NOTE:
#	1. mesh does not yet correspond to molecule
#	2. have not yet coded Smol equation correctly
#
# ToDo:
#	1. Complete Born ion example
#	2. Compute flux around boundary to get kon
#
from dolfin import *
import numpy as np
import Sphere
#import Molecule
from view import *
from params import * # must do this for class

class empty:pass

parms  = params()
problem = empty()
results = empty()
 
## PDE terms

# PMF term in Smoluchowski equation
# F = del W(r)
# [1]	Y. Song, Y. Zhang, T. Shen, C. L. Bajaj, J. A. McCammon, and N. A. Baker, Finite element solution of the steady-state Smoluchowski equation for rate constant calculations.
def SmolPMF(problem,V,psi):

    pmf = parms.valence * psi
   
    # Not needed, given how I formulated the PDE 
    #dpmf = grad(pmf) # presumably this is differentiating wrt xyz (or see pg 309 in Dolfin manual)

    problem.pmf = pmf


# This works
def Test():

    mesh = UnitCube(32,32,32)

    # Function space
    V = FunctionSpace(mesh, "CG", 1)

    # The potential
    psi = Function(V)
    psi.vector()[:] = 1.0 # WRONG WRONG WRONG 

    Vv = VectorFunctionSpace(mesh,"CG",1) # need Vector, not scalar function space 
    # Need to know the spatial dimension to compute the shape of derivatives.
    Jp = project(D*grad(psi),Vv)
 



## get kon
# See eqn (4) of Notes
# following example on pg 619
# if sphere, validate against analytical result (see 111128_todo)
def ComputeKon(problem,results):
 
    # SOLVED: # ERROR: ufl.log.UFLException: Shape mismatch.
    Vv = VectorFunctionSpace(problem.mesh,"CG",1) # need Vector, not scalar function space 

    # Need to know the spatial dimension to compute the shape of derivatives.
    print "Check on negative sighns in expr"
    # don't need to reproject Jp = project(D * intfact * grad(invintfact * up),Vv)
    Jp = parms.D * grad(results.up) + parms.beta * results.up * grad(problem.pmf)

    boundary_flux_terms = assemble(dot(Jp, tetrahedron.n)*ds(parms.active_site_marker),
                                exterior_facet_domains = problem.subdomains)
    print "boundary fl %f" % boundary_flux_terms
    kon = boundary_flux_terms / float(parms.bulk_conc);

    print "My kon is %f" % (kon)

    results.kon = kon
    results.Jp  = Jp   

# load in APBS example, which has geometry and potential
# but no boundary
def MeshWPotential(fileMesh,fileSubdomains,filePotential):
  ## load data
  # coordinates
  mesh = Mesh(fileMesh);


  # Function space
  V = FunctionSpace(mesh, "CG", 1)

  # load subdomains 
  subdomains = MeshFunction("uint", mesh, fileSubdomains) 
  bc0 = DirichletBC(V,Constant(parms.active_site_absorb),subdomains,parms.active_site_marker)
  bc1 = DirichletBC(V,Constant(parms.bulk_conc),subdomains,parms.outer_boundary_marker)
  bcs = [bc0,bc1]
  #PrintBoundary(mesh,bc0,file="file1.pvd") 
  #PrintBoundary(mesh,bc1,file="file2.pvd") 
  #quit()

  # apply file values to fuction space
  if(filePotential!="none"):
    psi = Function(V,filePotential);
    psi.vector()[:]=0;
  else:
    psi = Function(V)
    psi.vector()[:] = 0.0
 


  # alternatively I need to interpolate the grid from APBS 
  # see email from Johand around 1128

  #problem = empty()
  problem.subdomains = subdomains
  problem.mesh=mesh
  problem.psi = psi
  problem.bcs = bcs
  problem.V   = V
  return problem


# boundary is given as input, here we compute potential (for sphere)
# given input boundary. Note: normally youd want to use APBS to get
# the electrostatic potential
def MeshNoPotential(debug=0):

  ## load data
  # data from ~/localTemp/NBCR/smol/born_ex
  fileMesh = "example/sphere/p.pqr.output.all_mesh.xml.gz"
  mesh = Mesh(fileMesh)
  

  # Function space
  V = FunctionSpace(mesh, "CG", 1)


  ## load markers
  useMarkers=1
  # Need to pull these from Gamer output at some point 

  ## define Dirichlet boundary
  bcs=1
  if(useMarkers==1):
    bcs=1
  else:
    # PKH: (Found no facets matching domain for boundary condition.) <-- probably from Dirichlet, since Neumann BC expressed in weak form of PDE
    bc_active = DirichletBC(V, Constant(parms.active_site_absorb), Sphere.DirichletActiveSite())
    bc_bulk = DirichletBC(V, Constant(parms.bulk_conc), Sphere.DirichletBulkBoundary())
    
    bcs = [bc_active, bc_bulk]

  if(debug==0):
      return mesh

  ## Define Neumann boundary
  if(useMarkers==1):
    1
  else:      
    # PKH: Verify that Neumann BC is implemented correctly
    subdomains = MeshFunction("uint",mesh,2)
    subdomain = Sphere.NeumannMolecularBoundary()
    subdomain.mark(subdomains,molecular_boundary_marker)

    # PKH: do I need to mark BulkBoundary as well in order to use assemble call later?
    # NOTE: is this for sure where I should be computing kon? Seems like the flux along this boundary s.b. zero!
    subdomainOuter = Sphere.DirichletBulkBoundary()
    subdomainOuter.mark(subdomains,outer_boundary_marker)

  haveAPBS=0
  if(haveAPBS):
    #psi = GetPotential()
    1
  ## compute potential
  else:
    # PKH Is this the best way of computing Sphere potential over mesh coordinates? 
    x = mesh.coordinates() 
    psi = interpolate(Sphere.Potential(x),V)
    # PKH: I need to check potential, since code x-plodes 
    #psi.vector()[:]=0;

  problem.subdomains = subdomains
  problem.mesh=mesh
  problem.psi = psi
  problem.bcs = bcs
  problem.V   = V
  return problem

def PDEPart(problem): # h,psi,bcs,V):

    # Compute W, dW from psi
    SmolPMF(problem,problem.V,problem.psi)
    V = problem.V

    # The solution function
    u = Function(V)

    # Test function
    v = TestFunction(V)


    # The diffusion part of PDE
    # Recasting as integration factor (see Eqn (5) in Notes)
    intfact = exp(- parms.beta * problem.pmf)
    invintfact = 1/intfact;
    # Create weak-form integrand (see eqn (6) in NOtes)
    # also refer ti Zhou eqn 2,3 uin 2011
    # NOTE: this is the u that satisfies Eqn (6), not the traditional Smol eqn
    # form of the PDE
    #  (no time dependence,so only consider del u del v term)
    F = parms.D * intfact*inner(grad(u), grad(v))*dx

    # apply Neumann cond 
    # PKH - check w Fenics example, since it seems like ds() no longer needed?
    # (because I already marked the subdomain)
    # PKH where is 'subdomain' used that I defined earlier? 
    # since noflux --> boundary==0, don't need to include this 
    #f_molecular_boundary = noflux_molecular_boundary
    #F += f_molecular_boundary*v*ds(molecular_boundary_marker)


    # Solve the problem
    # prob. is linear in u, so technically don't need to use a non-linear solver....
    solve(F==0, u, problem.bcs)

    # Project the solution
    # Return projection of given expression *v* onto the finite element space *V*
    # Solved for u that satisfies Eqn 6, so obtain u we want by transformation (See Eqn (7))
    up = project(intfact*u)

    File("solution.pvd") << up
    #plot(up, interactive=True)
 
    results.up = up
    #problem.pmf= pmf # shouldn't go here 
    return results

## Domain

def Run(fileMesh,fileSubdomains,filePotential="NONE"):


  # get stuff 
  problem = MeshWPotential(fileMesh,fileSubdomains,filePotential)

  # solve PDE
  results= PDEPart(problem)
  
  # print solution
  File("up.pvd") << results.up

  # compute something
  ComputeKon(problem, results)    



def Debug():
  # was fileMesh  ="example/molecule/potential-0_mesh.xml.gz"
  #fileMesh  ="example/molecule/p.pqr.output.out.mesh.xml.gz"
  fileMesh  ="example/molecule/p.pqr.output.out_mesh.xml.gz"
  # was potential="example/molecule/potential-0_values.xml.gz"
  filePotential="example/molecule/p.pqr.output.out_values.xml.gz"
  fileSubdomains="example/molecule/p.pqr.output.out_subdomains.xml.gz"

  # if loading from APBS
  apbs =1
  gamer= 0
  test = 0
  if 0:
    1
  ## NOT SUPPORTED YET 
  if(apbs==1):
    problem = MeshWPotential(fileMesh,fileSubdomains,filePotential)

  # sphere example
  elif(gamer==1):
    problem = MeshNoPotential()
 
  elif(test==1):
    # ignore me for testing PKH
    Test()
    quit()

  else:
    raise RuntimeError("Bug Pete to write usable code. You request dumbfounded me")


  problem.mesh=mesh
  problem.psi = psi
  problem.bcs = bcs
  problem.V   = V

  # solve PDE
  results = PDEPart(problem) # h,psi,bcs,V)
  
  # print solution
  File("up.pvd") << up

  # compute somethingnn
  ComputeKon(problem,results) 




def SimpleTest():
  fileMesh = "example/sphere/p.pqr.output.all_mesh.xml.gz"
  mesh = Mesh(fileMesh)
  V = FunctionSpace(mesh, "CG", 1)

  #bc_test = DirichletBC(V, Constant(2), DirichletActiveSite())
  #bc_test = DirichletBC(V, Constant(2), DirichletBulkBoundary())
  bc_test = DirichletBC(V, Constant(2), Sphere.NeumannMolecularBoundary())

  
  from view import PrintBoundary   
  PrintBoundary(mesh,bc_test)
  
  return mesh

#
# PDE
# 0 = del D [del - 
# F(r) = - del U(r)
# Eqn 3 of pg 4 in  notetaker doc
# but, we use integration factor instead to get Eqn (4)
# intfact = exp(-beta*pmf)
# 0 = del . (D intfact del (1/intfact p)
#


if __name__ == "__main__":
  msg="smol.py <test> or smol.py <mesh.gz> <subdomains.gz> <values.gz> or smol.py -root <file>"

  import sys
  if len(sys.argv) < 2:  
      raise RuntimeError(msg)

  if(sys.argv[1]=="test"):
    print "In testing mode"
    Debug()

  elif(sys.argv[1]=="-root"):
    root = sys.argv[2]
    fileMesh = root+"_mesh.xml.gz"
    fileSubdomains= root+"_subdomains.xml.gz"
    filePotential= root+"_values.xml.gz"
    import os.path
    if(os.path.isfile(filePotential)==0):
       filePotential= "none";
       print "Didn't see electrostatic potential..."

    Run(fileMesh,fileSubdomains,filePotential="none")

  elif(len(sys.argv)==3):
    print "In run mode"
    fileMesh = sys.argv[1]
    fileSubdomains= sys.argv[2]
    Run(fileMesh,fileSubdomains,filePotential="none")

  elif(len(sys.argv)==4):
    print "In run mode"
    fileMesh = sys.argv[1]
    fileSubdomains= sys.argv[2]
    filePotential= sys.argv[3]
    Run(fileMesh,fileSubdomains,filePotential)

  else:
    raise RuntimeError(msg)

