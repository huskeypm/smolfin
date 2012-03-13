# for SERCA

from dolfin import *
import numpy as np

## VARIABLES

# active site is defined as region within this 'sphere'
activeSiteLoc = np.array([0,1,0])
activeSiteR   = 0.5

# bulk conc is assumed at 'top' of box defined by z
topZ = 1.0

def IsNearActiveSite(x):
    r = np.linalg.norm(activeSiteLoc-x)
    isNearActiveSite = r < (activeSiteR + DOLFIN_EPS)
    #print x
    #print activeSiteLoc
    #print "r %f activeSiteR %f isNearActiveSite %d" % (r,activeSiteR,isNearActiveSite)

    return isNearActiveSite

def IsNearTop(x):
    isOnTop= x[2] > (topZ - DOLFIN_EPS)
    #print x
    #print "z %f topZ %f is %d" % (x[2],topZ,isOnTop)
    return isOnTop

# marks active site s region near active site location 
# marks active site s region near active site location 
class ActiveSite(SubDomain):
  def inside(self,x,on_boundary):
    isNearActiveSite = IsNearActiveSite(x)
    return isNearActiveSite * on_boundary
  
# marks outersphere 
class BulkBoundary(SubDomain):
  def inside(self,x,on_boundary):
    isOnTop= IsNearTop(x)                
    return isOnTop * on_boundary

# opposite of DirichletActiveSite 
class MolecularBoundary(SubDomain):
  def inside(self,x,on_boundary):
    isNotNearActiveSite = (IsNearActiveSite(x)==0)
    isNotOnTop= (IsNearTop(x)==0)                
    #print x
    #print "as %d t %d " % (isNotNearActiveSite,isNotOnTop)
    return isNotNearActiveSite*isNotOnTop*on_boundary
  
