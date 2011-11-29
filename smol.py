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

## VARIABLES 
active_site_absorb = 0
bulk_conc = 1.0

# temporary
temp_outerR = 5.0
temp_innerR = 1.0
temp_siteZ  = 0.0

## Dirichlet Boundaries
def define_activesite(mesh):
    coor = mesh.coordinates()

    # find maximum of box and trace out small region 
    ztop = np.max(coor[:,2])

    midpoint = np.mean(coor,axis=0)
    xymidpoint = midpoint[0:2]
    range=np.max(coor,axis=0) - np.min(coor,axis=0)
    xywidth = 0.5 * np.min(range) * 1/5.0;

def define_bulkboundary(mesh)
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

def active_site_marked(x):
  # load in boundaries 

  # assign '4' to Dirichlet outer boundary (i think) 
  # assign '5' to reflective BC (verify)
  # assign '1' to active site (verify)

  # temporary
  r = np.linalg.norm(x)
  isOnR = r < temp_innerR + DOLFIN_EPS 
  isOnSite = abs(z-temp_siteZ) < (2+ DOLFIN_EPS)
  return isOnR*isOnSite

  # check for marker
  marker = 0
  return marker

def molecular_boundary_marked(x):
  # load in boundaries 

  # assign '4' to Dirichlet outer boundary (i think) 
  # assign '5' to reflective BC (verify)
  # assign '1' to active site (verify)

  # temporary
  r = np.linalg.norm(x)
  isOnR = r < temp_innerR + DOLFIN_EPS 
  isNotOnSite = abs(z-temp_siteZ) > (2+ DOLFIN_EPS)
  return isOnR*isNotOneSite
   


  # check for marker
  marker = 0
  return marker

def bulk_boundary_marked(x):
  # load in boundaries 

  # assign '4' to Dirichlet outer boundary (i think) 
  # assign '5' to reflective BC (verify)
  # assign '1' to active site (verify)

  r = np.linalg.norm(x)
  isOnOuterR = r > temp_outerR - DOLFIN_EPS 
  return isOnOuterR


  # check for marker
  marker = 0
  return marker


## PDE terms

# W = del_x Psi * q
# rewrite, since needs to be exponentiated (might avoid by using 'beta' term below)
# [1]	Y. Song, Y. Zhang, T. Shen, C. L. Bajaj, J. A. McCammon, and N. A. Baker, “Finite element solution of the steady-state Smoluchowski equation for rate constant calculations.,” Biophysical Journal, vol. 86, no. 4, pp. 2017–2029, Apr. 2004.
def SmolPMF(V,psi):
    valence = Constant(2)

    pmf = valence * psi
    dpmf = grad(pmf) # presumably this is differentiating wrt xyz (or see pg 309 in Dolfin manual) 


# given a point in xyz, computes electrostatic potential due to a sphere
# psi(r) = Q/(4 * pi e e0 (1+kappa R)) 1/r exp(-kappa(r-R))
def SpherePotential(x):
  z = 2;
  e = 1; # replace w real value 
  Q = z * e;
  ee0 = 1; # replace with real permitivities
  kappa = 1; # replace w Debye length 
  bigR = 1; # replace w sphere size 
  pi = 3.14;

  # assuming sphere is centered at 0,0,0
  r = np.linalg.norm(x)

  h = Expression(4 * pi*ee0*(1+kappa*bigR))*1/r*exp(-kappa(r-R),r=0.0);
  h.r = r;
  # PKH - how do i return value 

  return h

# load in APBS example, which has geometry and potential
# but no boundary 
def MeshNoBoundaryWPotential()
  ## load data 
  # coordinates
  fileMesh  ="example/potential-0_mesh.xml.gz"
  # electrostatic potential 
  potential="example/potential-0_values.xml.gz"
  mesh = Mesh(fileMesh);

  ## define boundary 
  define_activesite(mesh)
  define_bulkboundary(mesh)

  bc_active = DirichletBC(V, Constant(active_site_absorb), active_site)
  bc_molecule=NeumannBC(V,Constant(0),molecular_boundary)
  bc_bulk = DirichletBC(V, Constant(bulk_conc), bulk_boundary)

  # apply file values to fuction space 
  psi = Function(V,potential);

  return mesh,psi


# boundary is given as input, here we compute potential (for sphere) 
# given input boundary. Note: normally youd want to use APBS to get
# the electrostatic potential 
def MeshWBoundaryNoPotential()
  ## load data 
  fileMesh = "example/p.pqr.output.all_mesh.xml.gz" 
  mesh = Mesh(fileMesh)
  # load markers 

  ## define boundary 
  bc_active = DirichletBC(V, Constant(active_site_absorb), active_site_marked)
  bc_molecule=NeumannBC(V,Constant(0),molecular_boundary_marked)
  bc_bulk = DirichletBC(V, Constant(bulk_conc), bulk_boundary_marked)

  ## compute potential 
  v= Function(V)
  psi = SpherePotential(v) 
 
  return mesh, psi

## Domain
# if loading from APBS
apbs = 0
gamer= 1
if(apbs==1):
  mesh, psi = MeshNoBoundaryWPotential()

# sphere example 
elif(gamer==1):
  mesh, psi = MeshWBoundaryNoPotential()

else:
  print "Dont understand"
  quit()

# Function space
V = FunctionSpace(mesh, "CG", 1)


# Boundaries
bcs = [bc_active, bc_molecule,bc_bulk]

# Compute dW 
SmolPMF(V,psi)

# The solution function
u = Function(V)

# Test function
v = TestFunction(V)

# The diffusion part of PDE
beta = 1 # Johan defines this as beta = exp(-valence*psi) to arrive at 
         # exponentiated form of PDE 
# refer ti Zhou eqn 2,3 uin 2011
F = beta*inner(grad(u), grad(v))*dx

# The PMF part of PDE (verify) 
beta_nrg = 0.693 # 1/kT, kcal/mol
pmfterm = beta_nrg * u * dpmf 
F = inner(pmfterm(u),v)*dx   # i think this is wrong 

# Solve the problem
solve(F==0, u, bcs)

# Project the solution
# Return projection of given expression *v* onto the finite element space *V*
# e.g. 'un-exponentiate'
# what is this (see pg 21 of Fenics manual)?
up = project(beta*u)

File("solution.pvd") << up
plot(up, interactive=True)


## get kon 
# following example on pg 619 
#boundary_flux = assemble( inner(surf_normal,Jp) ,boundary_bulk)  # verify, might need interior bondary 
Jp = Expression(D * del u .....) ???????
domain = MeshFunction("uint", mesh, 2)
boundary_flux_term  = assemble( Jp*ds, exterior_facet_domains = domain)
kon = boundary_flux_term / c0;

# if sphere, validate against analytical result (see 111128_todo)



