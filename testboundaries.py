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
# this script gives an easy test of boundary conditions defined in a
# boundary object

class empty:pass


problem = empty()

from dolfin import *
import numpy as np    
def SercaTest():
  # WORKS def u0_boundary(x, on_boundary): 
  # WORKS   return on_boundary
  
  # Create mesh and define function space
  mesh = UnitCube(10,10,10)
  
  ## boundaries 
  import SercaBoundary as bound
  problem.mesh = mesh
  problem.bound=bound

  return (problem)
  
def SercaActual():
  root = "/home/huskeypm/scratch/serca_mole/test2"
  problem.fileMesh = root+"_mesh.xml.gz"

  mesh = Mesh(problem.fileMesh);

  ## boundaries 
  import SercaBoundary as bound

  # modify active site 
  # from test.m       12458    2     -1.4566563606e+01     3.3411479950e+01     2.7281160355e+01
  #bound.activeSiteLoc = np.array([-1.4566563606e+01,    3.3411479950e+01,    2.7281160355e+01])
  # HACK 
  #bound.activeSiteLoc = np.array([-1.4566563606e+01,    3.3411479950e+01, -150])
  #bound.activeSiteR   = 10.0
  
 
  # mody outer
  #         12000    2      2.0618418884e+02    -1.4678872681e+02     1.4827406311e+02
  # NOT CORRECT? bound.topZ = 1.4827406311e+02
  #bound.topZ = 10.0 
  #HACKS 
  #bound.topZ = 0.0 


  problem.mesh = mesh
  problem.bound=bound

  return (problem)

def TroponinTest():
  
  # Create mesh and define function space
  mesh = UnitCube(10,10,10)
  
  ## boundaries 
  import TroponinBoundary as bound

  problem.mesh = mesh
  problem.bound=bound

  return (problem)

def TroponinActual():
  root = "/home/huskeypm/scratch/validation/troponin/troponin"
  problem.fileMesh = root+"_mesh.xml.gz"

  mesh = Mesh(problem.fileMesh);

  ## troponinboundariesaries 
  import TroponinBoundary as troponinboundaries

  # modify active site 
  # from marked.m 13893    1     -5.1399998665e+00    -6.1920005798e+01     1.3070001221e+02
  #troponinboundaries.activeSiteLoc = np.array([-5.1399998665e+00,   -6.1920005798e+01,    1.3070001221e+02])
  #troponinboundaries.activeSiteR   = 10.0
  
 
  # mody outer
  #    3252    5     -4.0252001953e+02    -2.6274002075e+02    -1.5098001099e+02
  #outerPoint = np.array([-4.0252001953e+02,   -2.6274002075e+02,   -1.5098001099e+02])
  #troponinboundaries.outerR = np.linalg.norm(troponinboundaries.dominantCentr[ [0,1] ])
  #print troponinboundaries.dominantCentr[ [0,1] ]
  #print "outerR: %f" % troponinboundaries.outerR
  #troponinboundaries.outerR = 690 # based on paraview

  # IN TROPONINBOUNDARY troponinboundaries.activeSiteLoc = np.array([-30.345, 39.371,216.75])
  # IN TROPONINBOUNDARY troponinboundaries.activeSiteR   = 10.0
  # IN TROPONINBOUNDARY troponinboundaries.outerR = 480 # based on paraview



  problem.mesh = mesh
  problem.bound=troponinboundaries

  return (problem)


def TestBoundaries(mode):

  if (mode=="sercatest"):
    problem = SercaTest()
  elif (mode=="tnctest"):
    problem = TroponinTest()
  elif (mode=="tnc"):
    problem = TroponinActual()
  elif (mode=="serca"):
    problem = SercaActual()

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
  nEle =  np.size(marked1.vector().array())
  print "Active site: (%d/%d)" % ( np.sum(marked1.vector().array()), nEle) 
  
  
  # bulk        
  bulkBoundary = bound.BulkBoundary()
  bc0 = DirichletBC(V,Constant(2),bulkBoundary)
  marked2 = Function(V)
  bc0.apply(marked2.vector())
  print "Bulk: (%d/%d)" % ( np.sum(marked2.vector().array())/2, nEle) 
  
  # molecular    
  molecularBoundary = bound.MolecularBoundary()
  bc0 = DirichletBC(V,Constant(4),molecularBoundary)
  marked4 = Function(V)
  bc0.apply(marked4.vector())
  print "Molecule: (%d/%d)" % ( np.sum(marked4.vector().array()/4), nEle) 
  
  marked = marked1
  marked.vector()[:]= marked1.vector()[:] + marked2.vector()[:] + marked4.vector()[:]
  #plot(marked,interactive=0)
  File(mode+".pvd") << marked



########
if __name__ == "__main__":
  msg="testboundaries.py <serca,sercatest,tnc,tnctest>"
  

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  TestBoundaries(mode=sys.argv[1])

#TestBoundaries(mode="tnc")
#TestBoundaries(mode="serca")




