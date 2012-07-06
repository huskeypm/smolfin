
from dolfin import *

# Domain
#mesh = SomeMesh()
mesh = UnitCube(32,32,32)

# Function space
V = FunctionSpace(mesh, "CG", 1)

# Dirichlet Boundaries
def active_site(x):
    return x[0] < DOLFIN_EPS 

def bulk_boundary(x):
    return x[0] > 1.0 - DOLFIN_EPS

#active_site = SomeSubdomain()
#bulk_boundary = AnotherSubdomain()
bulk_conc = 1.0
bc_active = DirichletBC(V, Constant(0), active_site)
bc_bulk = DirichletBC(V, Constant(bulk_conc), bulk_boundary)

bcs = [bc_active, bc_bulk]

# The potential
psi = Function(V)
#psi.vector()[:] = read_potential_values_from_somewhere()
psi.vector()[:] = 1.0 # WRONG WRONG WRONG 
valence = Constant(2)
beta = exp(-valence*psi)

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


