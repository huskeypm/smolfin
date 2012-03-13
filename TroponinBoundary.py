# For troponin C system

from dolfin import *
import numpy as np

## VARIABLES

# active site is defined as region within this 'sphere'
activeSiteLoc = np.array([0,1,0])
activeSiteR   = 0.5

# bulk conc is assumed at outer region
dominantAxis = np.array([0,0,1])
dominantCentr= np.array([0,0,0])
coordIdx = [0,1] # could determine these by grabbing 0 indices form domAxis
outerR = 0.5

def IsNearActiveSite(x):
    r = np.linalg.norm(activeSiteLoc-x)
    isNearActiveSite = r < (activeSiteR + DOLFIN_EPS)
    #print x
    #print "r %f activeSiteR %f isNearActiveSite %d" % (r,activeSiteR,isNearActiveSite)
    return isNearActiveSite

def IsNearOuter(x):
    r = np.linalg.norm(dominantCentr[coordIdx] - x[coordIdx])
    isOnOuter= r  > (outerR - DOLFIN_EPS)
    #print x
    #print dominantCentr[coordIdx] - x[coordIdx] 
    #print "r %f R %f [%f,%f] is %d" % (r,outerR,x[0],x[1],isOnOuter)
    return isOnOuter

# marks active site s region near active site location 
# marks active site s region near active site location 
class ActiveSite(SubDomain):
  def inside(self,x,on_boundary):
    isNearActiveSite = IsNearActiveSite(x)
    return isNearActiveSite * on_boundary
  
# marks outersphere 
class BulkBoundary(SubDomain):
  def inside(self,x,on_boundary):
    isOnOuter= IsNearOuter(x)                
    return isOnOuter * on_boundary

# opposite of DirichletActiveSite 
class MolecularBoundary(SubDomain):
  def inside(self,x,on_boundary):
    isNotNearActiveSite = (IsNearActiveSite(x)==0)
    isNotOnOuter= (IsNearOuter(x)==0)                
    isMol = isNotNearActiveSite*isNotOnOuter*on_boundary
    if (on_boundary & isNotOnOuter):
      print IsNearActiveSite(x)
      print IsNearOuter(x)
      print x
      print "as %d t %d tot %d " % (isNotNearActiveSite,isNotOnOuter,isMol)
      


    return isMol
  
