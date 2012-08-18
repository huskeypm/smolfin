import sys
from dolfin import *
from dolfin import nabla_grad as grad 
import smol


markerSubstrateEnzyme=1
markerProductEnzyme=4
markerOuter=5
concAbsorb= Constant(0.0)
concBulk = Constant(1.0) 

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
  cS_ProductEnzymeSurface = assemble(
    problem.cS*ds(markerProductEnzyme),exterior_facet_domains=subdomains)

  print "ERROR: Area calls are not working. Using fake value"
  area= assemble(
    Constant(1)*ds,exterior_facet_domains=subdomains,mesh=mesh)
    #Constant(1)*ds(markerProductEnzyme),exterior_facet_domains=subdomains,mesh=mesh)
  area = 300
  cS_ProductEnzymeSurface /= area
  print "cS_ProductEnzymeSurface %f " % cS_ProductEnzymeSurface
  print "area %f" % area
  
  print "Implementation NOT correct. Placeholder for now"
  # TODO: Instead of a BC that has a single value along the entire 
  # boundary, I need to apply the BC based on the local concentration at 
  # each point along the boundary, e.g.
  # uSProdEnz = uS( at product enzyme ) 
  # uPProdEnz =  1 - uSProdEnz
  # A = assemble('Function',mesh)
  # A[ where(x == product enzyme ) ] = uPProdEnz
  # assign BC

  
  print "WARNING: Not supposed to be using an average here"
  bc1 = DirichletBC(V,1-cS_ProductEnzymeSurface,subdomains,markerProductEnzyme)
  #PKH120816bc1 = DirichletBC(V,1-problem.cS,subdomains,markerProductEnzyme)
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
    #print "Using fake value for area"
    #area = 300
    cP_SubstrateEnzymeSurface /= area
    print "cP_SubstrateEnzymeSurface %f " % cP_SubstrateEnzymeSurface
    bc1 = DirichletBC(V,1-cP_SubstrateEnzymeSurface,subdomains,markerSubstrateEnzyme)

  # use exact values along surface as dirichlet
  else:
    print "Using exact values"
    #bc1 = DirichletBC(V,1.-problem.cP.vector()[:],subdomains,markerSubstrateEnzyme)
    bc1 = DirichletBC(V,problem.cP,subdomains,markerSubstrateEnzyme)


  # assign BC
  # PKH120816 
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
  print "kout P : %e " % kout

  kout = smol.ComputeKon(problem,results,
    useSolutionVector=1,solutionVector=problem.cI,subdomainMarker=markerProductEnzyme)
  print "kout I : %e " % kout

  kout = -1 * smol.ComputeKon(problem,results,
    useSolutionVector=1,solutionVector=problem.cP,subdomainMarker=markerProductEnzyme)
  print "kout S : %e " % kout

  problem.kout = kout
  
  return kout

def doit(problem):
  smol.LoadFiles(problem)
  solveProb(problem)



if __name__ == "__main__":
  msg="Purpose: Two simulate time-dep, reaction diff equation for the S->I->P reaction"
  msg=msg+"Usage: "
  msg=msg+"twoenzyme.py _mesh.xml _subdomains.xml"
  msg=msg+"Notes:"
  remap = "none"



  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  problem = empty()
  problem.fileMesh = sys.argv[1]
  problem.fileSubdomains = sys.argv[2]
  if(len(sys.argv)==3):
    print "arg"



  doit(problem)


