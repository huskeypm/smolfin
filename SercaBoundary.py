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
# for SERCA

from dolfin import *
import numpy as np

## VARIABLES

## DEFAULT 
# active site is defined as region within this 'sphere'
activeSiteLoc = np.array([0,1,0])
activeSiteR   = 0.5

# bulk conc is assumed at 'top' of box defined by z
topZ = 1.0

## USER 
activeSiteLoc = np.array([20.730,-27.038 , 8.204 - 113])
activeSiteR   = 10.0
topZ = 0.0


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
  
