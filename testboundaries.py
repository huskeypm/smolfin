# this script gives an easy test of boundary conditions defined in a
# boundary object

class empty:pass



from dolfin import *
def SercaTest():
  # WORKS def u0_boundary(x, on_boundary): 
  # WORKS   return on_boundary
  problem = empty()
  
  # Create mesh and define function space
  mesh = UnitCube(10,10,10)
  
  ## boundaries 
  import SercaBoundary as bound
  problem.mesh = mesh
  problem.bound=bound

  return (problem)
  

def TroponinTest():
  problem = empty()
  
  # Create mesh and define function space
  mesh = UnitCube(10,10,10)
  
  ## boundaries 
  import TroponinBoundary as bound

  problem.mesh = mesh
  problem.bound=bound

  return (problem)

def TroponinActual():
  root = "/home/huskeypm/scratch/troponin/marked"
  problem.fileMesh = root+"_mesh.xml.gz"

  mesh = Mesh(problem.fileMesh);

  ## boundaries 
  import TroponinBoundary as bound

  # modify active site 
  # from marked.m 13893    1     -5.1399998665e+00    -6.1920005798e+01     1.3070001221e+02
 
  # mody outer
  #    3252    5     -4.0252001953e+02    -2.6274002075e+02    -1.5098001099e+02


  problem.mesh = mesh
  problem.bound=bound

  return (problem)


def TestBoundaries(mode):

  if (mode="sercatest"):
    problem = SercaTest()
  elif (mode="tnctest"):
    problem = TroponinTest()
  elif (mode="tnc"):
    problem = TroponinActual()

  # apply 
  mesh = problem.mesh
  bound=problem.bound

  V = FunctionSpace(mesh, "Lagrange", 1)
  problem.V    = V
  
  # active site 
  activeSite = bound.ActiveSite()
  bc0 = DirichletBC(V,Constant(1),activeSite)
  marked1 = Function(V)
  bc0.apply(marked1.vector())
  
  
  # active site 
  bulkBoundary = bound.BulkBoundary()
  bc0 = DirichletBC(V,Constant(2),bulkBoundary)
  marked2 = Function(V)
  bc0.apply(marked2.vector())
  
  # active site 
  molecularBoundary = bound.MolecularBoundary()
  bc0 = DirichletBC(V,Constant(4),molecularBoundary)
  marked4 = Function(V)
  bc0.apply(marked4.vector())
  
  marked = marked1
  marked.vector()[:]= marked1.vector()[:] + marked2.vector()[:] + marked4.vector()[:]
  #plot(marked,interactive=0)
  File("test.pvd") << marked



########
TestBoundaries()




