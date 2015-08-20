
from dolfin import *
import numpy as np
class BoundaryDefinitions(object):
  def __init__(self, 
               # active site is defined as region within this 'sphere'
               activeSiteLoc = np.array([0,1,0]),
               activeSiteR = 20, # active site radius
               # bulk conc is assumed at outer region
               dominantAxis = np.array([0,0,1]),#?
               dominantCentr= np.array([0,0,0]),#?
               coordIdx = [0,1], # could determine these by grabbing 0 indices form domAxis
               outerR = 100 # based on paraview, plus a little less                
               ):

    self.activeSiteLoc =  activeSiteLoc
    self.activeSiteR  = activeSiteR 
    # bulk conc is assumed at outer region
    self.dominantAxis  = dominantAxis 
    self.dominantCentr = dominantCentr 
    self.coordIdx  = coordIdx 
    self.outerR  = outerR    

  def IsNearActiveSite(self,x):
        r = np.linalg.norm(self.activeSiteLoc-x)
        isNearActiveSite = r < (self.activeSiteR + DOLFIN_EPS)
        #if isNearActiveSite:
        #  print x
        #  print "r %f activeSiteR %f isNearActiveSite %d" % (r,self.activeSiteR,isNearActiveSite)

        return isNearActiveSite

  # assumes protein is centered at 0,0,0
  def IsNearOuter(self,x):
        #r = np.linalg.norm(self.dominantCentr[self.coordIdx] - x[self.coordIdx])
        r = np.linalg.norm(x) # assumes centered at 000
        isOnOuter= r  > (self.outerR - DOLFIN_EPS)
        #oif isOnOuter:
        #if 0:
        #  print x
        #  print self.dominantCentr[self.coordIdx] - x[self.coordIdx] 
        #  print "r %f R %f [%f,%f] is %d" % (r,self.outerR,x[0],x[1],isOnOuter)

        return isOnOuter

class GenericBoundary(SubDomain):
  def inside(self, x, on_boundary):
    #print "GB: %f %f %f " % (x[0],x[1],x[2])
    return on_boundary

# marks active site s region near active site location 
# marks active site s region near active site location 
class ActiveSite(SubDomain):
  def inside(self,x,on_boundary):
    isNearActiveSite = self.boundaryDefinitions.IsNearActiveSite(x)
    return isNearActiveSite and on_boundary

# marks outersphere 
class BulkBoundary(SubDomain):
  def inside(self,x,on_boundary):
    isOnOuter= self.boundaryDefinitions.IsNearOuter(x)                
    return isOnOuter and on_boundary

# opposite of DirichletActiveSite 
class MolecularBoundary(SubDomain):
  def inside(self,x,on_boundary):
    isNotNearActiveSite = (self.boundaryDefinitions.IsNearActiveSite(x)==False)
    isNotOnOuter= (self.boundaryDefinitions.IsNearOuter(x)==False)                
    isMol = isNotNearActiveSite and isNotOnOuter and on_boundary
#    if (on_boundary & isNotOnOuter):
#      print IsNearActiveSite(x)
#      print IsNearOuter(x)
#      print x
#      print "as %d t %d tot %d " % (isNotNearActiveSite,isNotOnOuter,isMol)




    return isMol

        
        
        

import smol
import numpy as np
from params import * # must do this for class???

import smol
parms = smol.parms
problem = smol.problem


# In[20]:

class empty:pass
    
def SetBoundaries(  
    problem,
    activeSiteLoc = np.array([24, 8,-5]),   
    activeSiteR = 5, 
    outerR=200
    ):
  # This is a bit convoluted, but can generalize later   
  boundaryDefinitions = BoundaryDefinitions(
    #activeSiteLoc = np.array([-30.345, 39.371,216.75]),
    activeSiteLoc =activeSiteLoc,
    activeSiteR   = activeSiteR   ,
    outerR=outerR
  )

  problem.mesh = Mesh(problem.fileMesh)
  mesh = problem.mesh
  subdomains = FacetFunction("uint", mesh, 0)
  problem.subdomains = subdomains
  # mark everything with default
  generic = GenericBoundary()
  generic.mark(subdomains,0) 

  activeSite = ActiveSite(); 
  activeSite.boundaryDefinitions = boundaryDefinitions
  activeSite.mark(subdomains,parms.active_site_marker) 

  bulkBoundary = BulkBoundary(); 
  bulkBoundary.boundaryDefinitions = boundaryDefinitions
  bulkBoundary.mark(subdomains,parms.outer_boundary_marker)

  molecularBoundary = MolecularBoundary(); # non-active region?
  molecularBoundary.boundaryDefinitions = boundaryDefinitions
  molecularBoundary.mark(subdomains,parms.molecular_boundary_marker)

  boundaries = empty()
  boundaries.activeSite = activeSite
  boundaries.molecularBoundary = molecularBoundary
  boundaries.bulkBoundary = bulkBoundary
  problem.boundary = boundaries

  # on PKH
  V= FunctionSpace(mesh,"Lagrange",1)

  # test
  File("subdoms.pvd") << subdomains
  File(problem.fileSubdomains) << problem.subdomains
  noElectro=0






