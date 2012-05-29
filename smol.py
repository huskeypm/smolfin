#
# Simple example that imports mesh from APBS
#
from dolfin import *
import numpy as np
#import Sphere
#import Molecule
from view import *
from params import * # must do this for class
# from smol import * 

class empty:pass

class problem:
  filePotential="none"
  fileMesh="none"
  fileSubdomains="none"

parms  = params()
problem = problem()
results = empty()
 
## PDE terms

# PMF term in Smoluchowski equation
# F = del W(r)
# [1]	Y. Song, Y. Zhang, T. Shen, C. L. Bajaj, J. A. McCammon, and N. A. Baker, Finite element solution of the steady-state Smoluchowski equation for rate constant calculations.
def ElectrostaticPMF(problem,psi,q="useparams"):

    problem.pmf = Function(problem.V)
    pmf = np.zeros( psi.vector().size() ) 

    if(q=="useparams"):
      pmf[:] = parms.valence * psi.vector()[:]
    else:
      print "Using q=%f for the ligand" % q
      pmf[:] = q * psi.vector()[:]

    # Sanity check, otherwise need to make grid bigger to get 0 at boundary 
    amin = np.min(pmf)
    amax = np.max(pmf)
    print "Potential range: %f - %f " %  (amin,amax)
    if(amin * amax > 0):
      print "WARNING: Your potential does not seem to cross through 0 (reguirement at outer boundary"

    
    # Not needed, givn how I formulated the PDE 
    #dpmf = grad(pmf) # presumably this is differentiating wrt xyz (or see pg 309 in Dolfin manual)

    # get too small
    THRESH=20
    idx = (np.where(pmf < -THRESH))[0]
    if(np.size(idx)>0): 
      print "Encountered %d values with very small PMF values/consider thresholding" % np.size(idx)

    # to large 
    idx = (np.where(pmf > THRESH))[0]
    if(np.size(idx)>0): 
      print "Encountered %d values with very high PMF values/consider thresholding" % np.size(idx)
    
    problem.pmf.vector()[:] = pmf

   

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
def ComputeKon(problem,results,subdomainMarker=-1,useSolutionVector=0,solutionVector=-1):
    ## smol defaults
    if(subdomainMarker==-1):
      subdomainMarker = parms.active_site_marker

    if(useSolutionVector==0):
      solutionVector=results.up
   
    # SOLVED: # ERROR: ufl.log.UFLException: Shape mismatch.
    #Vv = VectorFunctionSpace(problem.mesh,"CG",1) # need Vector, not scalar function space 

    # Need to know the spatial dimension to compute the shape of derivatives.
    # don't need to reproject Jp = project(D * intfact * grad(invintfact * up),Vv)
    Jp = parms.D * grad(solutionVector) + parms.beta * solutionVector * grad(problem.pmf)

    boundary_flux_terms = assemble(dot(Jp, tetrahedron.n)*ds(subdomainMarker),
                                exterior_facet_domains = problem.subdomains)
    kon = boundary_flux_terms / float(parms.bulk_conc);

    #print "Warning: using correction until i get units right. "
    #print "Based on smol for sphere w/o charge. "
    #print "Is something wrong w scaling by bulk conc?"
    kon = kon * parms.correction
    print "kon: %e [1/Ms] %e [1/M min]" % (kon, kon*60)

    results.kon = kon
    results.Jp  = Jp   
 
    return kon

# from (2.10) of Berez 
# problem s.b. globally defined, so we don't have to pass it in separately
def InterfaceFunction(x,on_boundary):
  # NOTE: need to figure out how to grab appropriate pmf value from mesh, since interfaceFunction passes in coordinate value 
  # should probably ask Dolfin community instead of pestering Johan!!
  g_xt = exp(-problem.pmf(x) * beta) * problem.bindingChannelTerm
  return g_xt
  

# load in APBS example, which has geometry and potential
# but no boundary
# ExteriorProblem means that we have two PDEs separated by an interface at the 'binding channel MOUTH ' in the 3D mesh. 
# ALthough we label this as an active site in the mesh, the active site is technically at the 'end' of the binding channel 
# and is handled in the interior problem. If exteriorProblem is defined, then we pass in the boundary condition at the 
# interface via interiorProblem  
import sys
import os.path
def ProblemDefinition(problem,boundaries=0):
  ## load data
  # coordinates
  if(os.path.exists(problem.fileMesh)==0):
    msg = "fileMesh %s does not exist" % problem.fileMesh
    raise RuntimeError(msg)
  else:
    mesh = Mesh(problem.fileMesh);

  # load subdomains 
  if(os.path.exists(problem.fileSubdomains)==0):
    msg = "fileSubdomains %s does not exist" % problem.fileSubdomains
    raise RuntimeError(msg)
  else:
    subdomains = MeshFunction("uint", mesh, problem.fileSubdomains) 

  # Load and apply electrostatic potential values to mesh 
  V = FunctionSpace(mesh, "CG", 1)
  if(problem.filePotential!="none"):
    if(os.path.exists(problem.filePotential) == 0):
      msg = "filePotential %s does not exist" % problem.filePotential 
      raise RuntimeError(msg)

    else:
      psi = Function(V,problem.filePotential);
    #psi.vector()[:]=0;
  else:
    psi = Function(V)
    psi.vector()[:] = 0.0


  ## assign init cond??????

  ## assign BC 
  bc0 = DirichletBC(V,Constant(parms.active_site_absorb),subdomains,parms.active_site_marker)
  bc1 = DirichletBC(V,Constant(parms.bulk_conc),subdomains,parms.outer_boundary_marker)

  # override with definition 
  if(boundaries!=0):
    print "Overriding marked boundaries with user definition (test with testboundaries.py)"
    bc0 = DirichletBC(V,Constant(parms.active_site_absorb),boundaries.activeSite)
    bc1 = DirichletBC(V,Constant(parms.bulk_conc),boundaries.bulkBoundary)
 
  bcs = [bc0,bc1]
 
  ## Package problem 
  problem.subdomains = subdomains
  problem.mesh=mesh
  problem.psi = psi
  problem.bcs = bcs
  problem.V   = V

  return problem



def SolveSteadyState(problem,pvdFileName="up.pvd",\
      q="useparams",twoEnzymeVer=0): 
  

    # Compute W, dW from psi
    if(twoEnzymeVer==0):
      ElectrostaticPMF(problem,problem.psi,q=q)
    else:
      print "Overriding PMF"
      problem.pmf = Function(V) 
      problem.pmf.vector()[:] = 0.

    V = problem.V

    # The solution function
    if(twoEnzymeVer==0):
      u = Function(V)
    else:
      print "WARNING: havent the slightest why I can't use same call"
      u = TrialFunction(V)

    # Test function
    v = TestFunction(V)

    import numpy as np 
    # The diffusion part of PDE
    # Recasting as integration factor (see Eqn (5) in Notes)
    intfact    =    exp(-parms.beta * problem.pmf)
    intfact_np = np.exp(-parms.beta * problem.pmf.vector()[:])

    # Create weak-form integrand (see eqn (6) in NOtes)
    # NOTE: this is the u that satisfies Eqn (6), not the traditional Smol eqn
    # form of the PDE
    #  (no time dependence,so only consider del u del v term)

    if(twoEnzymeVer==0):
      F = inner(grad(u),grad(v))*dx
      L = 0
    else:
      F = inner(grad(u),grad(v))*dx(1)
      beta = Constant(0.)
      L = inner(beta,v)*ds  # (1)


    # apply Neumann cond 
    # since noflux --> boundary==0, don't need to include this?
    #f_molecular_boundary = noflux_molecular_boundary
    #F += f_molecular_boundary*v*ds(molecular_boundary_marker)

    # Solve the problem
    # prob. is linear in u, so technically don't need to use a non-linear solver....
    if(twoEnzymeVer==0):
      solve(F==L, u, problem.bcs)
      x = u
    else:
      x=Function(V)
      solve(F==L, x, problem.bcs)


    # PKH: do I need to apply the initial condition? 
    # u_1 = params.bulk_conc * exp(-beta * pmf)
    # u_1.assign(u) # following pg 50 of Logg 

    

    # Project the solution
    # Return projection of given expression *v* onto the finite element space *V*
    # Solved for u that satisfies Eqn 6, so obtain u we want by transformation (See Eqn (7))
    up = project(intfact*x)
    print "Unprojected Solution range: %f - %f " %  (min(x.vector()),max(x.vector()))
    # gives non-negative values 
    #print "Projecting on difference basis"
    #up = project(intfact*u,FunctionSpace(problem.mesh,"DG",0))
    # overriding project with simply numpy op (to avoid issues w basis)
    #up.vector()[:] *= np.exp(-potential/kT) 
    up.vector()[:] = x.vector()[:] * intfact_np[:]
    results.up = up

    # debug
    print "Projected solution range: %f - %f " %  (min(up.vector()),max(up.vector()))

    ## print solution
    # File("solution.pvd") << up
    #plot(up, interactive=True)
    print "Printing solution to %s" % pvdFileName
    File(pvdFileName) << results.up

    #from view import plotslice
    #plotslice(problem,results,title="no title",fileName="slice.png")


 
    #problem.pmf= pmf # shouldn't go here 
    return results

## Domain

# boundaries - override markers used in mesh 
def Run(problem,boundaries=0,pvdFileName="up.pvd",q="useparams"):

  # get stuff 
  problem = ProblemDefinition(problem,boundaries=boundaries)

  # solve PDE
  results= SolveSteadyState(problem,pvdFileName=pvdFileName,q=q)
  
  # compute something
  ComputeKon(problem, results)    

  return results 



def Debug():
  # was fileMesh  ="example/molecule/potential-0_mesh.xml.gz"
  #fileMesh  ="example/molecule/p.pqr.output.out.mesh.xml.gz"
  fileMesh  ="example/molecule/p.pqr.output.out_mesh.xml.gz"
  # was potential="example/molecule/potential-0_values.xml.gz"
  filePotential="example/molecule/p.pqr.output.out_values.xml.gz"
  fileSubdomains="example/molecule/p.pqr.output.out_subdomains.xml.gz"

  # if loading from APBS
  apbs =1
  gamer= 0 # should never use this 
  test = 0 # should never use this 
  if 0:
    1
  ## NOT SUPPORTED YET 
  if(apbs==1):
    ProblemDefinition(problem)

  # sphere example
  elif(gamer==1):
    problem = MeshNoPotential() # in old.py 
 
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
  results = SolveSteadyState(problem) # h,psi,bcs,V)
  
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
    problem.fileMesh = root+"_mesh.xml.gz"
    problem.fileSubdomains= root+"_subdomains.xml.gz"
    problem.filePotential= root+"_values.xml.gz"
    import os.path
    if(os.path.isfile(problem.filePotential)==0):
       problem.filePotential= "none";
       print "Didn't see electrostatic potential..."

    Run(problem)

  elif(len(sys.argv)==3):
    print "In run mode"
    problem.fileMesh = sys.argv[1]
    problem.fileSubdomains= sys.argv[2]
    Run(problem)

  elif(len(sys.argv)==4):
    print "In run mode"
    problem.fileMesh = sys.argv[1]
    problem.fileSubdomains= sys.argv[2]
    problem.filePotential= sys.argv[3]
    Run(problem)
  
  else:
    raise RuntimeError(msg)

