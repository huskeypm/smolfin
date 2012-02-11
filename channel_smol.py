
from smol import *

active_site_radius = 1

# linear potential from 
def linear_potential(x):
  F = -1; # minimum energy 
  z = x[2];
  pmf = -(F/a) * z 

  if (z<-a):
    pmf = F 
  elsif (z>0):
    pmf = 0


  return pmf

# define 'area' (x,y) along channel as a function of z
# assume constant for now 
def sigma_channel(x)
  area = 4 * 3.14 * (active_site_radius^2)
  sigma = Constant(area)
  return sigma  


# input 1D (for now) PMF for channel 
def PMF_channel():
  # 
  linear_potential()
  
  
 


### problem definition 
# divided into interior (PMF-driven) and exterior (diffusion driven) problems
# This script primarily deals with interior problem and passes exterior problem to smol 

## PMF domain 
mesh_channel = UnitLine(10)
PMF_channel(mesh_channel)

## init cond
# g(r,0) = exp(-beta pmf(r))
# Since we assume the PMF is 0 outside of the channel,
# we let g(r>channel,0) = 1


# solve 1D PDE for g1 (eq 2.7 Berez)
# base this on simple 1D Fenics example 
#dg1/dt = d/dx D sigma(x) exp(-beta * pmf(x) ) d/dx g1 / (sigma(x) exp(-beta*pmf(x) 

# solve for g at t--> inf (stdy state)
kchannel = g(0)



## Diffusion domain (handled by Smol) 
# Cube of size a with reactive patch at x<b, z=0
# b = active_site_radius
mesh_domain = UnitCube(10,10,10)


## define boundaries
# Dirichlet 
# p(z=a) = 1
# Neumann (put into weak form) 
# dp(x>b,z<a) / dn = 0 
# dp(x<b,z=0) / dn = inf 


smolproblem <-- this stuff
ksurface = smol(smolproblem)


## k overall
# from eqn 2.9 of Berez
k = ksurface * kchannel






kon = smol()

k = Function(PMF_channel,kon)

