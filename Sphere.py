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

from dolfin import *
import numpy as np

## VARIABLES

# temporary
temp_outerR = 5.0
temp_innerR = 10.0
temp_siteZ  = 0.0
temp_diam   = 10 

# given a point in xyz, computes electrostatic potential due to a sphere
# psi(r) = Q/(4 * pi e e0 (1+kappa R)) 1/r exp(-kappa(r-R))
# See equations (1) and (2) on pg 3 in Notetaker 2011-11-10
def SpherePotential(x):
  z = 2;
  e = 1; # replace w real value
  Q = z * e;
  ee = 1; # ee0, replace with real permitivities
  kappa = 1; # replace w Debye length
  bigR = 1; # replace w sphere size
  pi = 3.14;

  # assuming sphere is centered at 0,0,0
  r0 = np.linalg.norm(x)

  h = Expression("Q / (4 * pi * ee * (1 + kappa * bigR)) * 1/r * exp(-kappa * (r-bigR))",
                 Q = Q,   
                 pi = pi,
                 ee = ee, 
                 kappa = kappa,
                 bigR = bigR,  
                 r=r0       
                 );

  return h

# marks inner sphere of example (hopefully, not entirely sure)
class DirichletActiveSite(SubDomain):
  def inside(self,x,on_boundary):
    # load in boundaries

    # assign '4' to Dirichlet outer boundary (i think)
    # assign '5' to reflective BC (verify)
    # assign '1' to active site (verify)

    # temporary
    r = np.linalg.norm(x)
    #if(r <temp_innerR + DOLFIN_EPS):
    #    print "r %f" % r
      
    isOnR = r < temp_innerR + DOLFIN_EPS
    z = x[2]
    isOnSite = abs(z-temp_siteZ) < (temp_diam+ DOLFIN_EPS)
    #isOnSite=1
    return isOnR*isOnSite
  
    # check for marker
    #marker = 0
    #return marker

# marks outersphere 
class DirichletBulkBoundary(SubDomain):
  def inside(self,x,on_boundary):
    # load in boundaries

    # assign '4' to Dirichlet outer boundary (i think)
    # assign '5' to reflective BC (verify)
    # assign '1' to active site (verify)

    r = np.linalg.norm(x)
    isOnOuterR = r > temp_outerR - DOLFIN_EPS
    return isOnOuterR

  
    # check for marker
    #marker = 0
    #return marker

# opposite of DirichletActiveSite 
class NeumannMolecularBoundary(SubDomain):
  def inside(self,x,on_boundary):

    # load in boundaries
  
    # assign '4' to Dirichlet outer boundary (i think)
    # assign '5' to reflective BC (verify)
    # assign '1' to active site (verify)
  
    # temporary
    r = np.linalg.norm(x)
    isOnR = r < temp_innerR + DOLFIN_EPS
    z = x[2]
    isNotOnSite = abs(z-temp_siteZ) > (temp_diam + DOLFIN_EPS) 
    return isOnR*isNotOnSite
  
  
  
    # check for marker
    #marker = 0
    #return marker
