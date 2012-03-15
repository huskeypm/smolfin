
#
# This script sets up relevant terms for PMF term describining binding channel 
#

import numpy as np
from params import *

parms = params()

class empty:pass

def LoadPMF(filePMF,x0=0,xL=0):
   # get data
   dat = np.loadtxt(filePMF);
   x = dat[:,0]
   y = dat[:,1]
 
   # use relevant ranges 
   # xrange - optional - keep points falling in this range 
   if(x0 != xL): # check if they're defined 
     print "Using xrange to trim PMF (%f,%f)" % (x0,xL)
     xrange = np.array([x0,xL])
     inds = (x > xrange[0]) &  (x < xrange[1])
     inds = inds[:] 
     V_x = y[inds]
     x = x[inds]
   else:
     V_x = y
     x   = x
   
   return(x,V_x)


#
# General class for defining interior PMF problem 
#
class DefaultInteriorProblem:
   # __init__ is a class constructor
   # __****__ is usually a special class method.
   # pmfFileName name of pmfFile
   # x0 smallest x value to be read
   # xL largest x value to be read 
   def __init__(self, pmfFileName,sigmaR=1,diff_const=1,x0=0,xL=25):
     # These values are created
     # when the class is instantiated.
     self.pmfFileName = pmfFileName
     self.sigmaR = sigmaR
     self.diff_const = diff_const
     (self.x,self.V_x) = LoadPMF(pmfFileName,x0=0,xL=25)

   
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

