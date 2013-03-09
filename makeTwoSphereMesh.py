"
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
"
#
# Revisions
#       10.08.10 inception
#

from dolfin import *
from meshmanipulate import *
import numpy as np
import twoenzymeSolver as tes
import view

res = 3 # sphere resolution 
sphereRadius = 5
dist = 20

# used to define boundary for spheres
margPct = 0.15
minPct = 1-margPct
maxPct = 1+margPct


# per Johan's suggeston for PETSC
#parameters.linear_algebra_backend="Epetra"
#parameters.linear_algebra_backend="uBLAS"  

class empty: pass

class Domain():
  def __init__(self,radius, location=np.array([0,0,0])):
    self.radius = float(radius)
    self.location = location # np.array([[-10,0,0]])

class GenericBoundary(SubDomain):
  def inside(self, x, on_boundary):
    #print "GB: %f %f %f " % (x[0],x[1],x[2])
    return on_boundary

class OutsideBoundary(SubDomain):
  def inside(self, x, on_boundary):
    dist = np.linalg.norm(x) # assuming zero centered 
    state = (dist  > (minPct*self.radius) and on_boundary)
    ##if(on_boundary):
    #print "OB: %f %f %d " % (dist, self.radius, state) 
    #if(dist > 10 and dist < 12):
    #print "OBa: %f %f %f " % (x[0],x[1],x[2])
    return state

class InsideBoundary(SubDomain):
#  def __init__(self):
#    self.midPoint = [0,0,0]
#    self.radius = 0

  def inside(self, x, on_boundary):
    dist = np.linalg.norm(x-self.midPoint)
    state = (dist  < (maxPct*self.radius) and on_boundary)
    #if(on_boundary):
    #print "IB: %f %f %d " % (dist, self.radius, state) 
    return state


def SimpleTest(problem): 
  mesh = problem.mesh
  subdomains = problem.subdomains

  V = FunctionSpace(mesh,"Lagrange",1)
  u = TrialFunction(V)
  v = TestFunction(V)

  ## mark boundaries 
  bc1 = DirichletBC(V, Constant(0), subdomains,tes.markerSubstrateEnzyme)
  # neumann implicit 
  bcO = DirichletBC(V, Constant(1), subdomains,tes.markerOuter)
  bcs = [bc1,bcO]
  
  Dii  = Constant(1.)
  beta = Constant(1.)
  form = inner(Dii*grad(u), grad(v))*dx(1)
  # source term indicative of flux betwen compartments 
  form += inner(beta,v)*dx(1)
  A = lhs(form)
  L = rhs(form)
  
  x = Function(V)
  solve(A==L,x,bcs)
  File("soln.pvd") << x

def markBoundaries(problem):

  mesh = problem.mesh
  #mesh = UnitCube(5,5,5)

  ## define boundaries 
  from dolfin import MeshFunction
  #OLDsubdomains = MeshFunction("uint", mesh, 2)
  subdomains = FacetFunction("uint", mesh, 0)
  problem.subdomains = subdomains 

  # mark everything with default
  generic = GenericBoundary()
  generic.mark(subdomains,tes.markerUnmarked)
  
  # sphere 1
  insideBoundary1 = InsideBoundary()
  insideBoundary1.midPoint = problem.domain1.location #@[-5,0,0]
  insideBoundary1.radius = problem.domain1.radius #5.1
  #print insideBoundary1.midPoint 
  #print insideBoundary1.radius
  insideBoundary1.mark(subdomains,tes.markerSubstrateEnzyme)
  
  # sphere 2
  insideBoundary2 = InsideBoundary()
  insideBoundary2.midPoint = problem.domain2.location #@[-5,0,0]
  insideBoundary2.radius = problem.domain2.radius #5.1
  #print insideBoundary2.midPoint 
  #print insideBoundary2.radius
  insideBoundary2.mark(subdomains,tes.markerProductEnzyme)
  
  # outside boundary 
  outsideBoundary = OutsideBoundary()
  outsideBoundary.radius = problem.sphereOuter.radius
  #print outsideBoundary.radius
  outsideBoundary.mark(subdomains,tes.markerOuter)

  # test
  File("subdoms.pvd") << subdomains

  ## mark boundaries 
  test = 1 
  if(test==1):
    V = FunctionSpace(mesh,"CG",1)
    bc1 = DirichletBC(V, Constant(1), subdomains,tes.markerSubstrateEnzyme)
    n1=view.PrintBoundary(problem.mesh,bc1,file="domain1")#,dbg=1)
    bc2 = DirichletBC(V, Constant(2), subdomains,tes.markerProductEnzyme)
    n2=view.PrintBoundary(problem.mesh,bc2,file="domain2")#,dbg=1)
    bcO = DirichletBC(V, Constant(4), subdomains,tes.markerOuter)
    nO=view.PrintBoundary(problem.mesh,bcO,file="sphereO")#,dbg=1)
    bcu = DirichletBC(V, Constant(8), subdomains,tes.markerUnmarked)
    nUnmarked = view.PrintBoundary(problem.mesh,bcu,file="unmarked",dbg=1)
    if(nUnmarked > 0):
      raise RuntimeError("SOmethign fishy happened with meshes. DYING") 
    if(n1!=n2 or n1!=nO):
      print "WARNING: number of vertices on two spheres differ -dbl check"
  #bc1 = DirichletBC(V, Constant(0), subdomains,tes.markerSubstrateEnzyme)
  ## neumann implicit 
  #bcO = DirichletBC(V, Constant(1), subdomains,tes.markerOuter)
  #bcs = [bc1,bcO]

def  pdbMeshMaker(pdbFile,dbg=0):
  from gamer import read_PDB_molsurf
  print "WARNING: need to implement this in makeCubicMesh series as well"
  sm = read_PDB_molsurf(pdbFile)
  print("vertices:", sm.num_vertices, "simplices:", sm.num_faces)

  if(dbg=="testsmol"):
    return sm 

  
  # Improve molecular surface
  iterMax=40
  iter = 0
  while not sm.smooth(20, 140, 1, True):
      iter=iter+1
      if(iter>iterMax):
        break 
      pass

  sm.coarse_flat()
  sm.coarse_dense(rate=2.0, numiter=iterMax)
  sm.smooth(20, 140, iterMax, True)
  sm.normal_smooth()
  sm.remove_unconnected_patches(10)
  for i in range(sm.num_faces):
    sm.face(i).m = 1

  return sm 


def make_twopdb_mesh(problem,dbg=0,dist=20,pdb1=-1,pdb2=-1):
  from gamer import SurfaceMesh, GemMesh

  print "WARNING: orrow makePDBMesh from makeCubicMeshPDB"

  inner1 = pdbMeshMaker(pdb1,dbg=dbg)
  if(dbg=="testmeshgen"):
    print "ERROR: REMOVING SERCOND MESH FOR NOW FIX"
    inner2=inner1
  else:
    inner2 = pdbMeshMaker(pdb2,dbg=dbg)


  domain1 = Domain(15,[-dist/2.,0,0])
  domain2 = Domain(15,[dist/2.,0,0])
  domain1.mesh = inner1
  domain2.mesh = inner2

  problem.domain1 = domain1
  problem.domain2 = domain2

  domain1.mesh = inner1
  domain2.mesh = inner2

  meshName = make_twodomain_mesh(problem,dbg=dbg,dist=dist)
  return meshName 


def make_sphere_mesh(problem,dbg=0,dist=20):
  from gamer import SurfaceMesh, GemMesh    

  print "Using sphereRadius %f and distance %f" % (sphereRadius,dist)

  # '5' just refers to the quality of the sphere
  inner1 = SurfaceMesh(res)
  inner2 = SurfaceMesh(res)

  domain1 = Domain(sphereRadius,[-dist/2.,0,0])
  scaleGeom(inner1,domain1.radius)

  domain2 = Domain(sphereRadius,[dist/2.,0,0])
  scaleGeom(inner2,domain2.radius)

  problem.domain1 = domain1
  problem.domain2 = domain2

  domain1.mesh = inner1
  domain2.mesh = inner2

  meshName = make_twodomain_mesh(problem,dbg=dbg,dist=dist)
  return meshName 

def make_twodomain_mesh(problem,dbg=0,dist=20):
  from gamer import SurfaceMesh, GemMesh    
  print "Using sphere separation of %f" % dist 
  base = "temp"
  meshName = base+"_mesh.xml"

  ## outer sphere
  # should base this number on distance/radii of spheres
  sphereOuter = Domain(500)
  outer  = SurfaceMesh(res)
  scaleGeom(outer,sphereOuter.radius)
  problem.sphereOuter = sphereOuter

  ## inner meshes 
  domain1 = problem.domain1
  domain2 = problem.domain2
  inner1 = domain1.mesh
  inner2 = domain2.mesh
  
  
  transvec = empty()
  transvec.x=domain1.location[0]
  transvec.y=domain1.location[1]
  transvec.z=domain1.location[2]
  translate(inner1,transvec)
  inner1.as_hole = True
  
  if(dbg=="testmeshgen"):
    print "ERROR: ignoreing '2'"
    gem_mesh = GemMesh([outer, inner1])
    gem_mesh.write_dolfin(meshName)
    quit()
  else:
    transvec.x=domain2.location[0]
    transvec.y=domain2.location[1]
    transvec.z=domain2.location[2]
    translate(inner2,transvec)
    inner2.as_hole = True
  

  # mesh stuff 
  if(dbg!="testsmol"):
    gem_mesh = GemMesh([outer, inner1,inner2])
    gem_mesh.write_off(base+".off")
    gem_mesh.write_mcsf(base+".m")
    gem_mesh.write_dolfin(meshName)
  #import gamerutil
  #gamerutil.write_dolfin(gem_mesh,base)
  
  print "Note: Printed mesh uses marker '1' for volume, so integrate over dx(1)"
  
  return meshName

# taken from mcsftodolfin
#def mark_neumann_facets(mesh, cellmarkers):
#    print "reintegrate into mcsftodolfin"
#    #from dolfin import cells, facets,mesh
#    from dolfin import cells, facets, MeshFunction, File
#    md = mesh.domains()
#    facet_markers = md.markers(2)
#
#    # 2 for top dim 2 (facets)
#    for j, cell in enumerate(cells(mesh)):
#       for i, face in enumerate(facets(cell)):
#           marker = int(cellmarkers[j,i])
#           facet_markers.set_value(cell.index(), i, marker)
#
#    subdomains = MeshFunction("uint", mesh, facet_markers)
#    print "Num faces at active_site_marker:", \
#          (subdomains.array()==parms.active_site_marker).sum()



def doit(dist=20,meshType="sphere",pdb1=-1,pdb2=-1,dbg=0,solve=1):
  problem=empty()

  # gamer 
  if(meshType=="sphere"):
    meshName = make_sphere_mesh(problem,dbg=dbg,dist=dist)
  elif(meshType=="pdb"):
    meshName = make_twopdb_mesh(problem,dbg=dbg,dist=dist,pdb1=pdb1,pdb2=pdb2)
  else:
    raise RuntimeError("Arg, don't understnad" )

  

  # load via dolfin and mark
  from dolfin import Mesh, FacetFunction
  problem.mesh = Mesh(meshName)
  markBoundaries(problem) 

  # write outputs 
  filename="twosphere"
  File(filename+"_mesh.xml.gz") << problem.mesh
  File(filename+"_subdomains.xml.gz") << problem.subdomains

  ## solve 
  if(solve==1):
    tes.solveProb(problem)
  elif(solve==-1):
    SimpleTest(problem)

  #if(test==1):
  #  from view import PrintBoundary
  #  # surface w no marker
  #  V = FunctionSpace(mesh,"CG",1)
  #  #bc0 = DirichletBC(V, Constant(1),subdomains,parms.unmarked_marker )
  #  #boundary = InsideBoundary()
  #  bc0 = DirichletBC(V, Constant(1),boundary)
  #  PrintBoundary(mesh, bc0,file="boundary")

  # write 
  #import mcsftodolfin as mtd
  #quit()
  #mtd.write_dolfin_files("test",mesh)

import sys

if __name__ == "__main__":
  msg=\
"""Purpose: Makes two-sphere system for dolfin

Usage: 
  script.py <-dist dist> [-sphereRadius val] [-pdb] [-noSolve/-simpleSolve] 

    where dist = distance in [A] 
    sphereRadius: define radius for spheres
    pdb: make mesh for pdbs ]hardcoded for now] 
    simpleSolve: simple diff pde test

Notes:
"""
  remap = "none"
  pdbMode = 0
  meshType = "sphere"
  solve = 1
  dbg=0


  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  #if(len(sys.argv)==2):
  #  dist = float(sys.argv[1])

  for i,arg in enumerate(sys.argv):
    if(arg=="-dist"):
      dist=float(sys.argv[i+1])
    if(arg=="-sphereRadius"):
      sphereRadius=float(sys.argv[i+1])
    if(arg=="-pdb"):
      pdbMode = 1
    if(arg=="-noSolve"):
      solve = 0
    if(arg=="-simpleSolve"):
      solve = -1
    if(arg=="-skip"):
      dbg = "testsmol"

  
  if(pdbMode==0):
    doit(dist=dist,meshType=meshType,dbg=dbg,solve=solve)        
  else:
    root = "/home/huskeypm/scratch/120806/"
    doit(dist=dist,meshType="pdb",\
      pdb1=root+"1O27_Ac.pdb",pdb2=root+"2C2T_Ac.pdb",dbg=dbg,solve=solve)        


