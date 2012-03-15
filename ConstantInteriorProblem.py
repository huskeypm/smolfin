
#
# This script sets up relevant terms for PMF term describining binding channel 
#

import numpy as np
from params import *

parms = params()

class empty:pass

def DefineConstantPMF(V0=0,L=1):
   x0 = 0
   xL = L
   incr = 10
   x = np.linspace(x0,xL,incr)
   V_x = np.ones(np.size(x)) * V0
   #print x
   #print V_x
    
   return (x,V_x)

#
# defining interior PMF problem with Constant potential 
#
class ConstantInteriorProblem:
   # __init__ is a class constructor
   # __****__ is usually a special class method.
   # L - length of pore 
   def __init__(self, sigmaR=1,diff_const=1,L=1,V0=0):
     # These values are created
     # when the class is instantiated.
     self.sigmaR = sigmaR
     self.diff_const = diff_const
     (self.x,self.V_x) = DefineConstantPMF(V0=V0,L=L)

   
   # area along reaction coordinate for cylinder
   def Sigma(self):
     # R = 1 [nm]
     area = 4 * np.pi * self.sigmaR*self.sigmaR
     sigma_x = np.ones( np.size(self.x)) * area 
     return(sigma_x)

   # for now assuming the diffusion constant is constant
   def DiffusionConst(self):
     diff_x = np.ones( np.size(self.x)) * self.diff_const
     return(diff_x)

