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

% inert plane for now 

% check to see if APBS can give electrostatic potential in half-plane
% otherwise might need to compute for subsection, then extrapolate the potential to an infinite plane


%gamer: create half-plane 
% workaround 
% - create globular molecule in Gamer
% - slice molecular in blender
% - take boolean of rectangle minus half-molecule 
% - mark: inert surface, molecule boundary, active site, 'far Dirichlet' and nearest walls (s.b. like myocyte problem)   	
% - tetrehderalize (and go to town....)

% unless I am mistaken, I believe I can use all of the same machinery as the spherical problem,
% since both find kon by integrating the stdy-state flux over the reactive site. 
% although accord to Berez 2011 eqn 2.11, I need a fully absorbing cond at the tunnel entrance 

% do steady state solution 
ke = smol()









