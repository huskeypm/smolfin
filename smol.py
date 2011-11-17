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

## PDE terms

# W = del_x Psi * q
# rewrite, since needs to be exponentiated
def smolforceterm(V,fileValues):
    psi = Function(V,fileValues);
    valence = Constant(2)
    beta = exp(-valence*psi)


# Domain
#mesh = SomeMesh()
fileMesh  ="example/potential-0_mesh.xml.gz"
fileValues="example/potential-0_values.xml.gz"
mesh = Mesh(fileMesh);
# Function space
V = FunctionSpace(mesh, "CG", 1)


# Boundaries
define_activesite(mesh)
define_bulkboundary(mesh)
bulk_conc = 1.0
bc_active = DirichletBC(V, Constant(0), active_site)
bc_bulk = DirichletBC(V, Constant(bulk_conc), bulk_boundary)

bcs = [bc_active, bc_bulk]

# The potential
#psi = Function(V)
#psi.vector()[:] = read_potential_values_from_somewhere()
#psi.vector()[:] = 1.0 # WRONG WRONG WRONG 


smolforceterm(V,fileValues)

# The solution function
u = Function(V)

# Test function
v = TestFunction(V)

# The form
F = beta*inner(grad(u), grad(v))*dx

# Solve the problem
solve(F==0, u, bcs)

# Project the solution
# Return projection of given expression *v* onto the finite element space *V*
up = project(beta*u)

File("solution.pvd") << up

plot(up, interactive=True)


