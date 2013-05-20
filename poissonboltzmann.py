"""
----------------------------------
Smolfin - a numerical solver of the Smoluchowski equation for interesting geometries
Copyright (C) 2012 Peter Kekenes-Huskey, Ph.D., huskeypm@gmail.com

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA


----------------------------------------------------------------------------
"""
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from dolfin import *
import numpy as np
import sys
sys.path.append("/home/huskeypm/sources/jay/")
#import poissonboltzmannSpecial as pbs
class empty():pass
pbs = empty()
#import modelparameters

#
# Revisions
# Todo
# -validate all aspects of code 
#  The sinh expression initially comes out as zero, so it never looks like it is used in the form expression 
# -implement log function 
# -wholeDom prob. I think is behaving reasonably (need to visualize the potential in the exterior domain at a different dynamic range [0-10] than the interior domain)


msg="""
Purpose:
  Tempate for solving (non)-linear poisson boltzmann equation 

Usage:
  poissonboltzmann.py -linear/-nonlinear/-finite/-wholedomain/-validation
  
Notes: 
  Guaranteed to be wrong for right now!!

Author:
  the computational scientist formally known as pete 



"""
## interpolation has been validated 
# Verified that Expression is being appropriately interpolated to mesh (3D sphere, Gamer output, 1D line and analytical)


minr = 12.
maxr = 22.

class params:

  ## Parameters that are generally fixed
  # Standard units: A, mV, kcal, M

  # ec / 4 pi eps = 14.3840 [V Angstoms] 
  # --> eco4pieps = 14.384e3 [mV Angstroms]
  eco4pieps = 14.384e3 # conversion factor [mV Angstroms]
  #eco4pieps = 14.384 # [V Ang]
  kT = 0.59	# energy [kcal/mol]
  beta = 1/kT
  ec = 8.854187817e-12 # electric constant [C^2/(Jm)]
  M_TO_ANG = 1e-10
  J_TO_KCAL = 0.000239005736
  ec = ec / (M_TO_ANG * J_TO_KCAL) # ec [C^2/kcal A]
  epsilonExterior = 80. # dielectric constant in exterior []

  ## Domain-specific parameters 
  dim = 2 # grid dimensions 
  center = np.zeros(dim)      
  #dim = 2  # 2d mesh 
  molRad=12.5 # radius of internal boundary in sphere
  #a=10.0 # radius of internal boundary in sphere
  molMarker = 2 # marker for molecular boundary 
  domMarker = 3  # marker for domain boundary 
  epsError = 0.001  # epsilson for 'error'
  epsError = 3.000  # epsilson for 'error' in selecting boundary
  domRad = 5. * molRad # radius of domain [A] (kind of, since square)
  
  ## System-specific parameters
  z = -3.       # unit charge 

  # ion
  ionC = .150 # ion conc [M]
  ionRad=2. # ion radius [A]
  kappa = 0.328 * np.sqrt(ionC) # inverse debye length for monovalent specieds [1/A]
  ikappa  = 1/kappa # Debye length [A], "Intermoplecular and surface forces, Israelachvili" 

  ## solver modes
  mode = "linear"# linear, nonlinear, finitesize 
  #mode = "nonlinear"
  #mode = "finite"    



## expression
# Assumes charge is centered at origin
def DebyeHuckelExpr():
  # e^{ka}/(1+ka) Z ec / (4 pi eps eps0) * e^{-kr}/r
  prefac=params.z * np.exp(params.kappa*params.molRad) * params.eco4pieps 
  prefac/=params.epsilonExterior*(1+params.kappa*params.molRad)
  # f is in units [mV]
  f=Expression("prefac*exp(-k*sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]))/sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])",prefac=prefac,k=params.kappa)
  # validated using a=10.0, ionC=0.1, Z=5 
  return f

# share params 
pbs.params = params

class molecularBoundary(SubDomain):
  def inside(self,x,on_boundary):
    result = np.linalg.norm(x-params.center) < (params.epsError+params.molRad)
    result = result and on_boundary
#    if result:
#      print np.linalg.norm(x)
#      print x
    return result      

class domainBoundary(SubDomain):
  def inside(self,x,on_boundary):
    result = np.linalg.norm(x-params.center) > (-params.epsError+params.domRad)
    result = result and on_boundary
    #print x
    #print result
    return result      


def SolvePoissonBoltzmann(mesh,meshType="dolfin"):
  #
  print "Assuming your sphere is centered at 000"
  params.dim = np.shape(mesh.coordinates())[1]
  params.center = np.zeros(params.dim)      


  V = FunctionSpace(mesh, "Lagrange", 1)

  ## Define boundary conditions
  subdomains = MeshFunction("uint",mesh,params.dim-1)
  bcs = []

  # define BC on molecular boundary 
  boundary = molecularBoundary()
  boundary.mark(subdomains,params.molMarker)
  # PKH: should maybe use sympy later 
  # Analytical solution for linearized PBE for atom of radius R and charge q
  # see Eq 5.1 [1]	
  #M. Holst, N. Baker, and F. Wang, JCC, v21, n15,pp1319-1342, 2000
  #q = params.ec * params.z
  #k = params.kappa/np.sqrt(params.epsilonExterior)
  #u0 = Expression("q/epsilon*R*(1-k*R/(1+k*a))",\
  #  q=q,k=k,epsilon=params.epsilonExterior,\
  #  R=params.molRad,a=params.molRad)    
  f = DebyeHuckelExpr()
  bcs.append(DirichletBC(V, f, subdomains,params.molMarker))

  # test 
  m = Function(V)
  #bc = DirichletBC(V,Constant(1.),boundary)
  #bc.apply(m.vector())
  bcs[0].apply(m.vector())
  File("test.pvd") << m


  # define BC on domainboundary (potential should be zero at boundary) 
  boundary = domainBoundary()
  boundary.mark(subdomains,params.domMarker)
  bcs.append(DirichletBC(V, f, subdomains,params.domMarker))

  
  
  ## Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)

  ## define form functions 
  import ufl
  namespace = ufl.__dict__.copy()
  
  # exterior
  # LHS  
  if(meshType=="dolfin"): 
    form = -1*params.epsilonExterior*inner(grad(u), grad(v))*dx
  else:
    form = -1*params.epsilonExterior*inner(grad(u), grad(v))*dx(1)

  # eps*grad(u)grad(v) = kappa^2 uv
  if(params.mode=="linear"):
    if(meshType=="dolfin"): 
      form += -1*params.epsilonExterior*params.kappa*params.kappa *u*v*dx
    else:
      form += -1*params.epsilonExterior*params.kappa*params.kappa *u*v*dx(1)
  # eps*grad(u)grad(v) = kappa^2 sinh(u)*v
  elif(params.mode=="nonlinear"):
    form = pbs.Nonlinear()

  # eq (7) notes.pdf
  elif(params.mode=="finite"):
    form = pbs.Finite()

  else:
    raise RuntimeError("What did you just tell me to do???")
  
  # Compute solution
  print "Solving %s form of PBE" % params.mode
  x = Function(V)
  solve(lhs(form)==rhs(form), x, bcs)


  File("out.pvd") << x
  
  return (V,x)




# Mostly to check that interpolations are correct
def ValidateDebyeHuckel(meshFile):
  f = DebyeHuckelExpr()

  ## 3D mesh interpolation
  #mesh = Mesh("sphere_mesh.xml.gz")
  mesh = UnitSphere(50)
  mesh.coordinates()[:] = 50*mesh.coordinates()[:]
  V = FunctionSpace(mesh,"CG",1)
  phi = interpolate(f,V)
  #(gx,gy,gz) = np.mgrid[0:0:1j,0:0:1j,-255:255:500j]
  (gx,gy,gz) = np.mgrid[0:0:1j,0:0:1j,minr:maxr:100j]
  interp = griddata(mesh.coordinates(),phi.vector(),(gx,gy,gz))
  interp[np.isnan(interp)]=0
  interp = np.reshape(interp,100)
  gz = np.reshape(gz,100)
  
  ## GAMER 
  meshg = Mesh(meshFile)                
  Vg = FunctionSpace(meshg,"CG",1)
  phig = interpolate(f,Vg)
  #(gx,gy,gz) = np.mgrid[0:0:1j,0:0:1j,-255:255:500j]
  (gxg,gyg,gzg) = np.mgrid[0:0:1j,0:0:1j,minr:maxr:100j]
  interpg = griddata(meshg.coordinates(),phig.vector(),(gxg,gyg,gzg))
  interpg[np.isnan(interpg)]=0
  interpg = np.reshape(interpg,100)
  gzg = np.reshape(gzg,100)
  
  
  
  ## 1D mesh
  mesh1 = UnitInterval(150) # 0..1
  m=mesh1.coordinates()
  mesh1.coordinates()[:]=m*(maxr-minr) + minr
  V1 = FunctionSpace(mesh1,"CG",1)
  phi1 = interpolate(f,V1)
  gz1= np.mgrid[minr:maxr:100j]
  interp1 = griddata(mesh1.coordinates(),phi1.vector(),gz1)
  
  ## analytical
  print "WARNING: not sure what Expression assumes if just one argument is passed in"
  def vec(r):
    return f(r)
  vecf = np.vectorize(vec)
  
  # double check
  # agrees with Expression 
  def func(r,prefac,k):
    x = np.array([r,0,0])
    p =prefac*np.exp(-k*np.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]))/sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
    return p
  
  
  ## compare 
  plt.figure()
  plt.plot(gz,interp,"r", label="3d UnitSphere")
  plt.plot(gzg,interpg,"r.", label="3d gamer")
  plt.plot(gz1,interp1,"g",label="1d UnitInterval")
  plt.plot(gz1,vecf(gz1),"b",label="analytical")
  plt.xlim([minr,maxr])
  #plt.legend(['a','b','c','d'])
  plt.legend(loc=2)                         
  title = "Debye-Huckel potential for %4.1f [A] sphere, 1/k= %4.1f [A], z = %d" %\
    (params.molRad, params.ikappa, params.z)
  plt.title(title)
  plt.xlabel("r [A]")
  plt.ylabel("potential [mV]")
  plt.gcf().savefig("debyehuckel.png")

  
def ValidatePoissonBoltzmann(meshFile):

  ## GAMER 
  mesh = Mesh(meshFile)               
  (V,x) = SolvePoissonBoltzmann(mesh)

  ## validation of solution 
  # interpolate solution to line 
  (gx,gy,gz) = np.mgrid[0:0:1j,0:0:1j,minr:maxr:100j]
  interp = griddata(mesh.coordinates(),x.vector(),(gx,gy,gz))
  interp[np.isnan(interp)]=0
  interp = np.reshape(interp,100)
  gz = np.reshape(gz,100)
  
  # analytical 
  f = DebyeHuckelExpr()
  phi = interpolate(f,V)
  interpv = griddata(mesh.coordinates(),phi.vector(),(gx,gy,gz))
  interpv[np.isnan(interpv)]=0
  interpv = np.reshape(interpv,100)

  ## compare 
  plt.figure()
  plt.plot(gz,interp,"r", label="Numerical")
  plt.plot(gz,interpv,"b",label="analytical")
  #plt.plot(gz1,vecf(gz1),"b",label="analytical")
  plt.xlim([minr,maxr])
  #plt.legend(['a','b','c','d'])
  plt.legend(loc=2)                         
  title = "Debye-Huckel potential for %4.1f [A] sphere, 1/k= %4.1f [A], z = %d" %\
    (params.molRad, params.ikappa, params.z)
  plt.title(title)
  plt.xlabel("r [A]")
  plt.ylabel("potential [mV]")
  plt.gcf().savefig("debyehuckel_solve.png")

# Performs validation against DebyeHuckel equation 
def Validations():
  params.dim = 3
  params.center = np.zeros(params.dim)
  # 15 ensures all points at boundary r=12.5 are found   
  # 200 ensures all points at boundary r=225 are found   
#  params.innerBound = 15
#  params.outerBound = 220
  meshFile = "./example/sphere/sphere_mesh.xml.gz"
  ValidateDebyeHuckel(meshFile)
  ValidatePoissonBoltzmann(meshFile)


  

#sphere
if __name__ == "__main__":

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= "./example/sphere/sphere2d.xml"
  mode = "outer"
  if(len(sys.argv)==3):
    print "arg"

  for i,arg in enumerate(sys.argv):
    if(arg=="-nonlinear"):
      params.mode="nonlinear"
    if(arg=="-finite"):
      params.mode="finite"
    if(arg=="-wholedomain"):
      mode = "wholedomain"
      fileIn= "./example/sphere/sphere_2d_entire.xml"
    if(arg=="-linear"):
      params.mode = "linear"
      params.molRad = 1.5
      params.domRad = 5. 
      #fileIn= "./example/sphere/sphere2d.xml"
      fileIn= "./example/2d/volFrac_0.27.xml"
    if(arg=="-validation"):
      Validations()
      quit()


  if(mode=="wholedomain"):
    pbs.domainBoundary = domainBoundary
    pbs.molecularBoundary = molecularBoundary
    pbs.doWholeDomainPB(fileIn)
  else: 
    params.epsError = 1.000  # epsilson for 'error' in selecting boundary
    params.molRad=12.5 # actual radius is 
    params.domRad=15.
    mesh = Mesh(fileIn)
    SolvePoissonBoltzmann(mesh)


