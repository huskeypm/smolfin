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









