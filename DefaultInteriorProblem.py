
#
# This script sets up relevant terms for PMF term describining binding channel 
#

import numpy as np
from params import *

import smol
parms = smol.parms
#parms = params()

class empty:pass

def LoadPMF(filePMF,x0=0,xL=0):
   # get data
   print filePMF
   dat = np.loadtxt(filePMF);
   x = dat[:,0]
   y = dat[:,1]
   #print x
   #print y
 
   # use relevant ranges 
   # xrange - optional - keep points falling in this range 
   if(x0 != xL): # check if they're defined 
     print "Using xrange to trim PMF (%f,%f)" % (x0,xL)
     xrange = np.array([x0,xL])
     inds = (x > xrange[0]) &  (x < xrange[1])
     inds = inds[:] 
     V_x = y[inds]
     x = x[inds]
     #print V_x
     #print x
   else:
     V_x = y
     x   = x


   # Shift so that V(L) = 0
   idx_xL = np.size(x)-1
   V_x = V_x - V_x[idx_xL]
   #print V_x


   return(x,   V_x)


#
# General class for defining interior PMF problem 
#
class DefaultInteriorProblem:
   # __init__ is a class constructor
   # __****__ is usually a special class method.
   # pmfFileName name of pmfFile
   # radius [Ang]
   # diff const [um^2/s]
   # x0 smallest x value to be read
   # x0 [Ang]
   # xL largest x value to be read 
   def __init__(self, pmfFileName,channelR=1,diff_const=1,x0=0,xL=25):
     # These values are created
     # when the class is instantiated.
     self.pmfFileName = pmfFileName
     self.channelR = channelR * parms.Ang_to_um 
     print self.channelR
     self.diff_const = diff_const
     (x,self.V_x) = LoadPMF(pmfFileName,x0=x0,xL=xL)

     self.x = x * parms.Ang_to_um
   
  

     # from 2012-01-04 notes 
     # Validated 120329
     testMode = 0
     if(testMode==1):
       self.channelR = 5. * parms.Ang_to_um # [Ang]
       L = 10. * parms.Ang_to_um
       self.x = np.linspace(0,L,10)
       self.V_x = np.ones( np.size(self.x))
       # no pmf 
       self.V_x[:] = 0.0
       # higher pmf
       self.V_x[:] = 10.0
       self.V_x[:] = 3.0
       # lower pmf
       #self.V_x[:] = -10.0
       #self.V_x = self.x * -10.0
       
       

   
   # area along reaction coordinate for cylinder
   def Sigma(self):
     # R = 1 [nm]
     area = 4 * np.pi * self.channelR**2
     print "R %e A %e " % (self.channelR,area)
     sigma_x = np.ones( np.size(self.x)) * area 
     return(sigma_x)

   # for now assuming the diffusion constant is constant
   def DiffusionConst(self):
     diff_x = np.ones( np.size(self.x)) * self.diff_const
     return(diff_x)

