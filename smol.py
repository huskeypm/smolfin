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
from params import *
from view import *

class empty:pass
 
beta = 1/0.693 # 1/kT, kcal/mol
D = 1 # Diffusion const.
## PDE terms

# PMF term in Smoluchowski equation
# F = del W(r)
# [1]	Y. Song, Y. Zhang, T. Shen, C. L. Bajaj, J. A. McCammon, and N. A. Baker, Finite element solution of the steady-state Smoluchowski equation for rate constant calculations.
def SmolPMF(V,psi):
    valence = Constant(2)

    pmf = valence * psi
   
    # Not needed, given how I formulated the PDE 
    #dpmf = grad(pmf) # presumably this is differentiating wrt xyz (or see pg 309 in Dolfin manual)

    return pmf #,dpmf


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
def ComputeKon(mesh,intfact,invintfact,results,V,problem):
    c0 = bulk_conc;
 
    # SOLVED: # ERROR: ufl.log.UFLException: Shape mismatch.
    Vv = VectorFunctionSpace(mesh,"CG",1) # need Vector, not scalar function space 
# PKH  - why 1? why define for entire mesh? 
    # Need to know the spatial dimension to compute the shape of derivatives.
    print "Check on negative sighns in expr"
    # don't need to reproject Jp = project(D * intfact * grad(invintfact * up),Vv)
    Jp = D * grad(results.up) + beta * results.up * grad(results.pmf)
    #print "removed part, replace!!"
    #Jp = project(grad(up),Vv)

    # WHY?? subdomains = MeshFunction("uint", mesh, mesh.topology().dim() - 1)


    # SOLVED: ERROR: ufl.log.UFLException: Dot product requires non-scalar arguments, got arguments with ranks 0 and 1.
    # ???
    #boundary_flux_terms = assemble(dot(Jp, tetrahedron.n)*ds(outer_boundary_marker),
    #n = FacetNormal(mesh)
    #boundary_flux_terms = assemble(dot(Jp, tetrahedron.n)*ds(outer_boundary_marker),
    #                            exterior_facet_domains = subdomains)
    print "NEED TO FITURE OUT WHY I CANT INTEGRATE OVER ACTIVE SITE!!"
    boundary_flux_terms = assemble(dot(Jp, tetrahedron.n)*ds(active_site_marker),
                                exterior_facet_domains = problem.subdomains)
    print "boundary fl %f" % boundary_flux_terms
    kon = boundary_flux_terms / float(c0);

    print "My kon is %f" % (kon)

    return (kon,Jp) 


# TODO:
def GetPotential():
  # read in apbs mesh/values (johan did this)
  meshAPBS;
  valuesAPBS;

  # read in/pass in mesh from Gamer (different node locations)
  meshGamer; 
  V ; # meshGamer's functionspace  

  #PKH  APBS mesh differs from Gamer mesh, therefore need to interpolate APBS values onto Gamer mesh.
  # maybe do this: http://stackoverflow.com/questions/7701669/interpolate-large-irregular-grid-onto-another-irregular-grid-in-python
  v = scipy.interpolate(meshAPBS,valuesAPBS,meshGamer,V)
  

  1

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
  bc0 = DirichletBC(V,Constant(active_site_absorb),subdomains,active_site_marker)
  bc1 = DirichletBC(V,Constant(bulk_conc),subdomains,outer_boundary_marker)
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

  problem = empty()
  problem.subdomains = subdomains
  return mesh,psi,bcs,V,problem


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
    bc_active = DirichletBC(V, Constant(active_site_absorb), Sphere.DirichletActiveSite())
    bc_bulk = DirichletBC(V, Constant(bulk_conc), Sphere.DirichletBulkBoundary())
    
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
    psi = GetPotential()

  ## compute potential
  else:
    # PKH Is this the best way of computing Sphere potential over mesh coordinates? 
    x = mesh.coordinates() 
    psi = interpolate(Sphere.Potential(x),V)
    # PKH: I need to check potential, since code x-plodes 
    #psi.vector()[:]=0;

  return mesh, psi,bcs,V

def PDEPart(mesh,psi,bcs,V):

    # Compute W, dW from psi
    pmf = SmolPMF(V,psi)

    # The solution function
    u = Function(V)

    # Test function
    v = TestFunction(V)


    # The diffusion part of PDE
    # Recasting as integration factor (see Eqn (5) in Notes)
    intfact = exp(- beta * pmf)
    invintfact = 1/intfact;
    # Create weak-form integrand (see eqn (6) in NOtes)
    # also refer ti Zhou eqn 2,3 uin 2011
    # NOTE: this is the u that satisfies Eqn (6), not the traditional Smol eqn
    # form of the PDE
    #  (no time dependence,so only consider del u del v term)
    F = D * intfact*inner(grad(u), grad(v))*dx

    # apply Neumann cond 
    # PKH - check w Fenics example, since it seems like ds() no longer needed?
    # (because I already marked the subdomain)
    # PKH where is 'subdomain' used that I defined earlier? 
    # since noflux --> boundary==0, don't need to include this 
    #f_molecular_boundary = noflux_molecular_boundary
    #F += f_molecular_boundary*v*ds(molecular_boundary_marker)


    # Solve the problem
    solve(F==0, u, bcs)

    # Project the solution
    # Return projection of given expression *v* onto the finite element space *V*
    # Solved for u that satisfies Eqn 6, so obtain u we want by transformation (See Eqn (7))
    up = project(intfact*u)

    File("solution.pvd") << up
    #plot(up, interactive=True)
 
    results = empty()
    results.up = up
    results.pmf= pmf
    return intfact,invintfact,results

## Domain

def Run(fileMesh,fileSubdomains,filePotential="NONE"):


  # get stuff 
  mesh, psi,bcs,V,problem = MeshWPotential(fileMesh,fileSubdomains,filePotential)

  # solve PDE
  intfact,invintfact,results= PDEPart(mesh,psi,bcs,V)
  
  # print solution
  File("up.pvd") << results.up

  # compute something
  kon = ComputeKon(mesh,intfact,invintfact,results,V,problem) 



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
    mesh, psi,bcs,V,problem = MeshWPotential(fileMesh,fileSubdomains,filePotential)

  # sphere example
  elif(gamer==1):
    mesh, psi,bcs, V = MeshNoPotential()
 
  elif(test==1):
    # ignore me for testing PKH
    Test()
    quit()

  else:
    raise RuntimeError("Bug Pete to write usable code. You request dumbfounded me")



  # solve PDE
  print "FAIL!!!"
  quit()
  intfact,invintfact,up= PDEPart(mesh,psi,bcs,V)
  
  # print solution
  File("up.pvd") << up

  # compute something
  (kon,Jp) = ComputeKon(mesh,intfact,invintfact,up,V,problem) 
  return Jp




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
# 0 = del D [del - beta F(r)] p(r)
# F(r) = - del U(r)
# Eqn 3 of pg 4 in  notetaker doc
# but, we use integration factor instead to get Eqn (4)
# intfact = exp(-beta*pmf)
# 0 = del . (D intfact del (1/intfact p)
#


if __name__ == "__main__":
  msg="smol.py <test> or smol.py <mesh.gz> <subdomains.gz> <values.gz>"

  import sys
  if len(sys.argv) < 2:  
      raise RuntimeError(msg)

  if(sys.argv[1]=="test"):
    print "In testing mode"
    Jp=Debug()

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

