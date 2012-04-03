
#
# This script sets up relevant terms for PMF term describining binding channel 
#

import numpy as np
from params import *

#parms = params()
import smol
parms = smol.parms


class empty:pass

def DefineLinearPMF(V0=0,L=1):
   x0 = 0
   xL = L
   incr = 10
   incr = 1000
#   incr = 2
   x = np.linspace(x0,xL,incr)
   # Note - this is the opposite of the potential in Berez
   # Here:  V_x(x=0) = V0, V_x(x=L) = 0 
   V_x = V0 * (L -x ) / L
   #print x
   #print V_x
    
   return (x,V_x)

#
# defining interior PMF problem with Linear potential 
#
class LinearInteriorProblem:
   # __init__ is a class constructor
   # __****__ is usually a special class method.
   # L - length of pore 
   def __init__(self, channelR=1,diff_const=1,L=1,V0=0):
     # These values are created
     # when the class is instantiated.
     self.channelR = channelR
     self.diff_const = diff_const
     (self.x,self.V_x) = DefineLinearPMF(V0=V0,L=L)

   
   # area along reaction coordinate for cylinder
   def Sigma(self):
     # R = 1 [nm]
     area = 4 * np.pi * self.channelR*self.channelR
     sigma_x = np.ones( np.size(self.x)) * area 
     return(sigma_x)

   # for now assuming the diffusion constant is constant
   def DiffusionConst(self):
     diff_x = np.ones( np.size(self.x)) * self.diff_const
     return(diff_x)

