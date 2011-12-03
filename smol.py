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


## Dirichlet Boundaries
def define_activesite(mesh):
    coor = mesh.coordinates()

    # find maximum of box and trace out small region
    ztop = np.max(coor[:,2])

    midpoint = np.mean(coor,axis=0)
    xymidpoint = midpoint[0:2]
    range=np.max(coor,axis=0) - np.min(coor,axis=0)
    xywidth = 0.5 * np.min(range) * 1/5.0;

def define_bulkboundary(mesh):
    coor = mesh.coordinates()

    # find maximum of box and trace out small region
    zbottom = np.min(coor[:,2])


# assuming active site is a small radius at top of boundary
def active_site(x):

    xydist = np.linalg.norm(xymidpoint-x[0:2])
    isInXY = xydist < xywidth
    isInZ  = x[2] > (ztop - DOLFIN_EPS)

    return isInXY * isInZ

def bulk_boundary(x):

    return x[2] < (zbottom + DOLFIN_EPS)

#def active_site_marked(x):
class DirichletActiveSite(SubDomain):
  def inside(self,x,on_boundary):
    # load in boundaries

    # assign '4' to Dirichlet outer boundary (i think)
    # assign '5' to reflective BC (verify)
    # assign '1' to active site (verify)

    # temporary
    r = np.linalg.norm(x)
    isOnR = r < temp_innerR + DOLFIN_EPS
    z = x[2]
    isOnSite = abs(z-temp_siteZ) < (2+ DOLFIN_EPS)
    return isOnR*isOnSite
  
    # check for marker
    #marker = 0
    #return marker

class DirichletBulkBoundary(SubDomain):
  def inside(self,x,on_boundary):
    # load in boundaries

    # assign '4' to Dirichlet outer boundary (i think)
    # assign '5' to reflective BC (verify)
    # assign '1' to active site (verify)

    r = np.linalg.norm(x)
    isOnOuterR = r > temp_outerR - DOLFIN_EPS
    return isOnOuterR

  
    # check for marker
    #marker = 0
    #return marker


class NeumannMolecularBoundary(SubDomain):
  def inside(self,x,on_boundary):

    # load in boundaries
  
    # assign '4' to Dirichlet outer boundary (i think)
    # assign '5' to reflective BC (verify)
    # assign '1' to active site (verify)
  
    # temporary
    r = np.linalg.norm(x)
    isOnR = r < temp_innerR + DOLFIN_EPS
    z = x[2]
    isNotOnSite = abs(z-temp_siteZ) > (2+ DOLFIN_EPS)
    return isOnR*isNotOnSite
  
  
  
    # check for marker
    #marker = 0
    #return marker


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

## get kon
# See eqn (4) of Notes
# following example on pg 619
# if sphere, validate against analytical result (see 111128_todo)
def ComputeKon(intfact,invintfact,up,V):
    c0 = bulk_conc;
 
    # ERROR: ufl.log.UFLException: Shape mismatch.
    # grad term is not  working, but it seems like its output should be the same size as input arguments?
    #Jp = project(D * intfact * grad(invintfact * up),V)
    Jp = project(D * intfact * (invintfact * up),V)

    subdomains = MeshFunction("uint", mesh, 2) # I called this earlier, ok to do again?


    # ERROR: ufl.log.UFLException: Dot product requires non-scalar arguments, got arguments with ranks 0 and 1.
    # ???
    #boundary_flux_terms = assemble(dot(Jp, tetrahedron.n)*ds(outer_boundary_marker),
    n = FacetNormal(mesh)
    boundary_flux_terms = assemble(dot(Jp, n)*ds(outer_boundary_marker),
                                exterior_facet_domains = subdomains)
    kon = boundary_flux_term / c0;

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

  #PKH  interpolate values onto my mesh?
  # maybe do this: http://stackoverflow.com/questions/7701669/interpolate-large-irregular-grid-onto-another-irregular-grid-in-python
  v = scipy.interpolate(meshAPBS,valuesAPBS,meshGamer,V)
  

  1

# given a point in xyz, computes electrostatic potential due to a sphere
# psi(r) = Q/(4 * pi e e0 (1+kappa R)) 1/r exp(-kappa(r-R))
# See equations (1) and (2) on pg 3 in Notetaker 2011-11-10
def SpherePotential(x):
  z = 2;
  e = 1; # replace w real value
  Q = z * e;
  ee = 1; # ee0, replace with real permitivities
  kappa = 1; # replace w Debye length
  bigR = 1; # replace w sphere size
  pi = 3.14;

  # assuming sphere is centered at 0,0,0
  r0 = np.linalg.norm(x)

  h = Expression("Q / (4 * pi * ee * (1 + kappa * bigR)) * 1/r * exp(-kappa * (r-bigR))",
                 Q = Q,   
                 pi = pi,
                 ee = ee, 
                 kappa = kappa,
                 bigR = bigR,  
                 r=r0       
                 );

  return h

# load in APBS example, which has geometry and potential
# but no boundary
def MeshNoBoundaryWPotential():
  ## load data
  # coordinates
  fileMesh  ="example/potential-0_mesh.xml.gz"
  # electrostatic potential
  potential="example/potential-0_values.xml.gz"
  mesh = Mesh(fileMesh);

  # Function space
  V = FunctionSpace(mesh, "CG", 1)

  ## define boundary
  define_activesite(mesh)
  define_bulkboundary(mesh)

  bc_active = DirichletBC(V, Constant(active_site_absorb), DirichletActiveSite())
  bc_bulk = DirichletBC(V, Constant(bulk_conc), DirichletBulkBoundary())
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
def MeshWBoundaryNoPotential():

  ## load data
  fileMesh = "example/p.pqr.output.all_mesh.xml.gz"
  mesh = Mesh(fileMesh)

  # Function space
  V = FunctionSpace(mesh, "CG", 1)


  ## load markers
  # Need to pull these from Gamer output at some point 

  ## define Dirichlet boundary
  bc_active = DirichletBC(V, Constant(active_site_absorb), DirichletActiveSite())
  bc_bulk = DirichletBC(V, Constant(bulk_conc), DirichletBulkBoundary())
  bcs = [bc_active, bc_bulk]

  ## Define Neumann boundary
  # PKH: not sure if this is right (added to PDE later) 
  subdomains = MeshFunction("uint",mesh,2)
  subdomain = NeumannMolecularBoundary()
  subdomain.mark(subdomains,molecular_boundary_marker)

  # PKH: do I need to mark BulkBoundary as well in order to use assemble call later?
  subdomainOuter = DirichletBulkBoundary()
  subdomainOuter.mark(subdomains,outer_boundary_marker)

  haveAPBS=0
  if(haveAPBS):
    psi = GetPotential()

  ## compute potential
  else:
    # PKH best way of doing this? 
    x = mesh.coordinates() 
    psi = interpolate(SpherePotential(x),V)
    # PKH: I need to check potential, since code x-plodes 
    psi.vector()[:]=0;

  return mesh, psi,bcs,V

def PDEPart():

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
    # it seems like I already marked the subdomain, so this shouldnt be
    # needed. That said, I need to see how to apply this f_mol. function
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
 
    return intfact,invintfact,up,V

## Domain

#
# PDE
# 0 = del D [del - beta F(r)] p(r)
# F(r) = - del U(r)
# Eqn 3 of pg 4 in  notetaker doc
# but, we use integration factor instead to get Eqn (4)
# intfact = exp(-beta*pmf)
# 0 = del . (D intfact del (1/intfact p)
#


# if loading from APBS
apbs = 0
gamer= 1
if(0):
  1
## NOT SUPPORTED YET 
##if(apbs==1):
##  mesh, psi,bcs,V = MeshNoBoundaryWPotential()

# sphere example
elif(gamer==1):
  mesh, psi,bcs, V = MeshWBoundaryNoPotential()

else:
  print "Dont understand"
  quit()




intfact,invintfact,up,V = PDEPart()
ComputeKon(intfact,invintfact,up,V) 
