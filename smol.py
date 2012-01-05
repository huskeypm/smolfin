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
import MoleculeBoundaries
 
## VARIABLES
active_site_absorb = 0
bulk_conc = 1.0
D = 1.0 # diffusion constant

# temporary
temp_outerR = 5.0
temp_innerR = 1.0
temp_siteZ  = 0.0


# markers
outer_boundary_marker = 4 # verify 
molecular_boundary_marker = 5 # verify 


## Dirichlet Boundaries (loaded separately now)




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
    D = 1 
    Jp = project(D*grad(psi),Vv)
 



## get kon
# See eqn (4) of Notes
# following example on pg 619
# if sphere, validate against analytical result (see 111128_todo)
def ComputeKon(mesh,intfact,invintfact,up,V):
    c0 = bulk_conc;
 
    # SOLVED: # ERROR: ufl.log.UFLException: Shape mismatch.
    Vv = VectorFunctionSpace(mesh,"CG",1) # need Vector, not scalar function space 
    # Need to know the spatial dimension to compute the shape of derivatives.
    # PKH - something about invintfact is incompataible with 'up'
    Jp = project(D * intfact * grad(invintfact * up),Vv)

    subdomains = MeshFunction("uint", mesh, 2) # I called this earlier, ok to do again?


    # SOLVED: ERROR: ufl.log.UFLException: Dot product requires non-scalar arguments, got arguments with ranks 0 and 1.
    # ???
    #boundary_flux_terms = assemble(dot(Jp, tetrahedron.n)*ds(outer_boundary_marker),
    #n = FacetNormal(mesh)
    boundary_flux_terms = assemble(dot(Jp, tetrahedron.n)*ds(outer_boundary_marker),
                                exterior_facet_domains = subdomains)
    kon = boundary_flux_terms / c0;

    print "My kon is %f" % (kon)

    return kon 


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
def MeshWPotential():
  ## load data
  # coordinates
  # was fileMesh  ="example/molecule/potential-0_mesh.xml.gz"
  fileMesh  ="example/molecule/p.pqr.output.out.mesh.xml.gz"
  mesh = Mesh(fileMesh);

  # electrostatic potential
  # was potential="example/molecule/potential-0_values.xml.gz"
  filePotential="example/molecule/p.pqr.output.out.values.xml.gz"


  # Function space
  V = FunctionSpace(mesh, "CG", 1)

  ## define boundary
  define_activesite(mesh)
  define_bulkboundary(mesh)

  bc_active = DirichletBC(V, Constant(active_site_absorb), MoleculeBoundaries.DirichletActiveSite())
  bc_bulk = DirichletBC(V, Constant(bulk_conc), MoleculeBoundaries.DirichletBulkBoundary())
  # Neumann done later 
  bcs = [bc_active, bc_bulk]
  
  print "need to add subdomain for neuman marker"
  quit()

  # apply file values to fuction space
  psi = Function(V,potential);

  # alternatively I need to interpolate the grid from APBS 
  # see email from Johand around 1128

  return mesh,psi,bcs,V


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
  if(useMarkers==1):
    1
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
    psi.vector()[:]=0;

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
    beta = 1/0.693 # 1/kT, kcal/mol
    intfact = exp(- beta * pmf)
    invintfact = 1/intfact;
    # Create weak-form integrand (see eqn (6) in NOtes)
    # also refer ti Zhou eqn 2,3 uin 2011
    # NOTE: this is the u that satisfies Eqn (6), not the traditional Smol eqn
    # form of the PDE
    #  (no time dependence,so only consider del u del v term)
    F = D * intfact*inner(grad(u), grad(v))*dx

    # apply Neumann cond 
    f_molecular_boundary = Constant(0) # or some function...

    # PKH - check w Fenics example, since it seems like ds() no longer needed?
    # (because I already marked the subdomain)
    # PKH where is 'subdomain' used that I defined earlier? 
    F += f_molecular_boundary*v*ds(molecular_boundary_marker)


    # Solve the problem
    # ERROR:   Process 0: *** Warning: Found no facets matching domain for boundary condition.
    solve(F==0, u, bcs)

    # Project the solution
    # Return projection of given expression *v* onto the finite element space *V*
    # Solved for u that satisfies Eqn 6, so obtain u we want by transformation (See Eqn (7))
    up = project(intfact*u)

    File("solution.pvd") << up
    #plot(up, interactive=True)
 
    return intfact,invintfact,up

## Domain


def SimpleTest():
  fileMesh = "example/sphere/p.pqr.output.all_mesh.xml.gz"
  mesh = Mesh(fileMesh)
  V = FunctionSpace(mesh, "CG", 1)

  #bc_test = DirichletBC(V, Constant(2), DirichletActiveSite())
  #bc_test = DirichletBC(V, Constant(2), DirichletBulkBoundary())
  bc_test = DirichletBC(V, Constant(2), Sphere.NeumannMolecularBoundary())
  
  marked = Function(V)
  bc_test.apply(marked.vector())

  #plot(u, interactive=TRUE)
  plot(marked, interactive=1)
  interactive()

  File("marked.pvd") << marked
  
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

  # if loading from APBS
  apbs = 0
  gamer= 1
  test = 0
  if 0:
    1
  ## NOT SUPPORTED YET 
  if(apbs==1):
    mesh, psi,bcs,V = MeshWPotential()

  # sphere example
  elif(gamer==1):
    mesh, psi,bcs, V = MeshNoPotential()
 
  elif(test==1):
    # ignore me for testing PKH
    Test()
    quit()

  else:
    print "Bug Pete to write usable code. You request dumbfounded me"
    quit()



  # solve PDE
  intfact,invintfact,up= PDEPart(mesh,psi,bcs,V)

  # compute something
  kon = ComputeKon(mesh,intfact,invintfact,up,V) 


