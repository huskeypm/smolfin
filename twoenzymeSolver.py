import sys
from dolfin import *
from dolfin import nabla_grad as grad 
import smol
import view
import numpy as np


markerUnmarked=8
markerSubstrateEnzyme=1
markerProductEnzyme=4
markerOuter=5
concAbsorb= Constant(0.0)
concBulk = Constant(1.0) 

# use average area on boundary for dirichlet condition 
useAverageArea=0

class empty:pass

def DirichletSubstrate(problem):
  V = problem.V
  subdomains = problem.subdomains

  # assign BC
  bc1 = DirichletBC(V,concAbsorb,subdomains,markerSubstrateEnzyme)
  bc2 = DirichletBC(V,concBulk,subdomains,markerOuter)
  #Neumann: dcdn = 0 for markerProductEnzyme
  bc = [bc1,bc2]
  return bc

def DirichletProduct(problem):
  V = problem.V
  mesh = problem.mesh
  subdomains = problem.subdomains

  if(useAverageArea==1):
    cS_ProductEnzymeSurface = assemble(
      problem.cS*ds(markerProductEnzyme),exterior_facet_domains=subdomains)
    area= assemble(
      Constant(1)*ds(markerProductEnzyme),exterior_facet_domains=subdomains,mesh=mesh)
    cS_ProductEnzymeSurface /= area
#    print "cS_ProductEnzymeSurface %f " % cS_ProductEnzymeSurface
#    print "area %f" % area
    bc1 = DirichletBC(V,1-cS_ProductEnzymeSurface,subdomains,markerProductEnzyme)

    # testing 

  # Instead of a BC that has a single value along the entire 
  # boundary, I need to apply the BC based on the local concentration at 
  # each point along the boundary, e.g.
  else:
    vcS = Function(problem.V)
    len = (np.shape(problem.cS.vector()[:]))[0]
    vcS.vector()[:] = np.ones(len) -  problem.cS.vector()[:]
    bc1 = DirichletBC(V,vcS,subdomains,markerProductEnzyme)


    # Test dirichlet
    view.PrintBoundary(problem.mesh,bc1,file="testdir")


  bc2 = DirichletBC(V,Constant(0.),subdomains,markerOuter)
  #Neumann: dcdn = 0 for markerSubstrateEnzyme
  bc = [bc1,bc2]
  return bc

def DirichletIntermediate(problem):
  V = problem.V
  subdomains = problem.subdomains

  # use surface averaged area
  useAverageArea=0
  if(useAverageArea==1):
    cP_SubstrateEnzymeSurface = assemble(
      problem.cP*ds(markerSubstrateEnzyme),exterior_facet_domains=subdomains)
    area= assemble(
      Constant(1.)*ds(markerSubstrateEnzyme),exterior_facet_domains=subdomains)
    cP_SubstrateEnzymeSurface /= area
    #print "cP_SubstrateEnzymeSurface %f " % cP_SubstrateEnzymeSurface
    bc1 = DirichletBC(V,1-cP_SubstrateEnzymeSurface,subdomains,markerSubstrateEnzyme)

  # use exact values along surface as dirichlet
  else:
    vcP = Function(problem.V)
    len = (np.shape(problem.cP.vector()[:]))[0]
    vcP.vector()[:] = np.ones(len) -  problem.cP.vector()[:]
    bc1 = DirichletBC(V,vcP,subdomains,markerSubstrateEnzyme)

  # assign BC
  bc2 = DirichletBC(V,Constant(0.),subdomains,markerOuter)
  bc3 = DirichletBC(V,concAbsorb,subdomains,markerProductEnzyme)
  bc = [bc1,bc2,bc3]
  return bc
  


# solve PDE subject to differing boundary conditions
# dirichlet:
#   substrate
#   product
#   intermediate
def solvecase(problem, case="unk"):

  if(case=='substrate'):
    bcs = DirichletSubstrate(problem)

  elif(case=='product'):
    bcs = DirichletProduct(problem)

  elif(case=='intermediate'):
    bcs = DirichletIntermediate(problem)

  else:
    print "not understood"
    quit()

  #print "Commented out what we should be running - need to debug"
  problem.bcs = bcs
  results = smol.SolveSteadyState(problem,twoEnzymeVer=1)


  if(case=='substrate'):
    problem.cS = Function(problem.V) 
    problem.cS.vector()[:] = results.up.vector()[:]

  elif(case=='product'):
    problem.cP = Function(problem.V) 
    problem.cP.vector()[:] = results.up.vector()[:]

  elif(case=='intermediate'):
    problem.cI = Function(problem.V) 
    problem.cI.vector()[:] = results.up.vector()[:]

  # write 
  File(case+".pvd") << results.up


def solveProb(problem):


  # function space 
  problem.V = FunctionSpace(problem.mesh,"CG",1)

  # 1. solve for substrate distribution at steady state
  # 2. solve for product distribution at steady state
  # 3. solve for intermediate distribution at steady state
  solvecase(problem,case="substrate")
  solvecase(problem,case="product")
  solvecase(problem,case="intermediate")


  # TODO I think this is the quantity we want, e.g. how fast the product leaves the 
  # product enzyme 
  import smol
  results = empty()
  kout = smol.ComputeKon(problem,results,
    useSolutionVector=1,solutionVector = problem.cS,subdomainMarker=markerSubstrateEnzyme)
  print "kout S : %e " % kout

  kout = smol.ComputeKon(problem,results,
    useSolutionVector=1,solutionVector=problem.cI,subdomainMarker=markerProductEnzyme)
  print "kout I : %e " % kout

  kout = -1 * smol.ComputeKon(problem,results,
    useSolutionVector=1,solutionVector=problem.cP,subdomainMarker=markerProductEnzyme)
  print "kout P : %e " % kout

  problem.kout = kout
  
  return kout

def doit(problem):
  smol.LoadFiles(problem)

  #V = FunctionSpace(problem.mesh,"CG",1)
  #bc1 = DirichletBC(V, Constant(1), problem.subdomains,markerSubstrateEnzyme)
  #view.PrintBoundary(problem.mesh,bc1,file="domain1")
  #bc2 = DirichletBC(V, Constant(1), problem.subdomains,markerProductEnzyme)
  #view.PrintBoundary(problem.mesh,bc2,file="domain2")
  #bcO = DirichletBC(V, Constant(1), problem.subdomains,markerOuter)
  #view.PrintBoundary(problem.mesh,bcO,file="sphereO")
  #quit()

  solveProb(problem)



if __name__ == "__main__":
  msg="""
\nPurpose: 
  To simulate time-dep, reaction diff equation for the S->I->P reaction

Usage:
  .py _mesh.xml _subdomains.xml

"""

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  problem = empty()
  problem.fileMesh = sys.argv[1]
  problem.fileSubdomains = sys.argv[2]
  if(len(sys.argv)==3):
    1



  doit(problem)


