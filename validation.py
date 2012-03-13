#
# Simple example that imports mesh from APBS
#
import smol 

class empty:pass
 
def ValidationPaper():
def ValidationZhou():
  # params 
  root = "/home/huskeypm/scratch/validation/sphere/"
  problem.fileMesh = root+"_mesh.xml.gz"
  problem.fileSubdomains= root+"_subdomains.xml.gz"
  R = 12.5# [Ang]
  R = R / 10000  # [um]
  q = -1  # qlig * qprot [1/C]

  ## No electro 
  # problem 
  problem.filePotential= "none"
  result = smol.Run(problem)
  
  # results 
  kon_analy = 4 * np.pi * R * parms.D * parms.um3_to_M
  print "kon_anal %e pred %e " % (kon_analy, result.kon)


  ## w electro 
  # problem 
  problem.filePotential= "_values.xml.gz"
  result = smol.Run(problem)

  # results 
  # from (18) of Song 
  kon_elec  = 4 * np.pi * parms.D * q * parms.bulk_conc * parms.um3_to_M
  # assuming r2-->inf
  kon_elec  = exp(q/R) / (exp(q/R) - exp(0)) 
  print "VERIFY kon_anal %e pred %e " % (kon_elec, result.kon)

  ## linear potential 
  sigma = 1 #
  L = np.array([1,2,3,4])
  bVo = -3 
  print "Not correct vals. Need to recompute using Zhou's formula"
  kon_linear = np.array([0.2,0.22,0.245,0.28])  # from Fig 2 Barreda (not exact!)
  # pass in array of V values: beta V0 = -3
  x = np.linspace(0,L[0],10)
  V = ones(np.size(x)) * bV0


  ## gating 
   # fig3 barrade
  exps = np.linspace(-6,6,1)
  omegams = 10**exps

   

#
# PDE
# 0 = del D [del - 
# F(r) = - del U(r)
# Eqn 3 of pg 4 in  notetaker doc
# but, we use integration factor instead to get Eqn (4)
# intfact = exp(-beta*pmf)
# 0 = del . (D intfact del (1/intfact p)
#


if __name__ == "__main__":
  msg="smol.py run "



  import sys
  if len(sys.argv) < 2:  
      raise RuntimeError(msg)

  if(sys.argv[1]=="run"):
    Validation()

  else:
    raise RuntimeError(msg)

