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
# For troponin C system

from dolfin import *
import numpy as np

## VARIABLES

## DEFAULT 
# active site is defined as region within this 'sphere'
activeSiteLoc = np.array([0,1,0])
activeSiteR   = 0.5

# bulk conc is assumed at outer region
dominantAxis = np.array([0,0,1])
dominantCentr= np.array([0,0,0])
coordIdx = [0,1] # could determine these by grabbing 0 indices form domAxis
outerR = 0.5


## USER
activeSiteLoc = np.array([-30.345, 39.371,216.75])
activeSiteR   = 10.0
outerR = 450 # based on paraview, plus a little less 





def IsNearActiveSite(x):
    r = np.linalg.norm(activeSiteLoc-x)
    isNearActiveSite = r < (activeSiteR + DOLFIN_EPS)
#    print x
#    print "r %f activeSiteR %f isNearActiveSite %d" % (r,activeSiteR,isNearActiveSite)

    return isNearActiveSite

def IsNearOuter(x):
    r = np.linalg.norm(dominantCentr[coordIdx] - x[coordIdx])
    isOnOuter= r  > (outerR - DOLFIN_EPS)
#    print x
#    print dominantCentr[coordIdx] - x[coordIdx] 
#    print "r %f R %f [%f,%f] is %d" % (r,outerR,x[0],x[1],isOnOuter)

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
#    if (on_boundary & isNotOnOuter):
#      print IsNearActiveSite(x)
#      print IsNearOuter(x)
#      print x
#      print "as %d t %d tot %d " % (isNotNearActiveSite,isNotOnOuter,isMol)
      
      


    return isMol
  
