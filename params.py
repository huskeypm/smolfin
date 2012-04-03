## VARIABLES
from dolfin import Constant


class params:
  def __init__(self):
    # CONVERTSION
    self.Ang_to_um  = 1.0e-4  # Ang into um
    self.Ang2_to_um2 = 1.0e-8  # Ang^2 into um^2
    self.um3_to_M = 602214150.0  # um^3 into M (1 cubic micron = 1.0 x 10-15 liter, times Avogadro's number)
  # markers
    self.unmarked_marker = 6; # Default marker - was '1' 
    self.active_site_marker = 1; # verified (red) 
    self.outer_boundary_marker = 5 # Verified on 120104_update.tex (green)
    self.molecular_boundary_marker = 4 # Verified (blue) 

    # Values 
    self.noflux_molecular_boundary = Constant(0) # generally won't change this 
    self.active_site_absorb = Constant(0)
    self.bulk_conc = Constant(1.0)
    self.D = 390.0 # diffusion constant, [um^2/s]
    #self.D = 100.0 # diffusion constant, [um^2/s]
    #print "WARNING: changed D temporarily"
    self.beta = 1/0.693 # 1/kT, [kcal/mol]
    self.valence = 2.0   # why do I need a constant? TODO  Constant(2)

    # temporary
    self.temp_outerR = 5.0
    self.temp_innerR = 1.0
    self.temp_siteZ  = 0.0

