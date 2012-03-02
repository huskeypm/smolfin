## VARIABLES
from dolfin import Constant


class params:
  def __init__(self):
    # CONVERTSION
    self.um3_to_M = 6*(10**8)  # um^3 into M # TODO Recompute exatcly 
  # markers
    self.unmarked_marker = 1; # Default marker 
    self.active_site_marker = 1; # verified (red) 
    self.outer_boundary_marker = 5 # Verified on 120104_update.tex (green)
    self.molecular_boundary_marker = 4 # Verified (blue) 

    # Values 
    self.noflux_molecular_boundary = Constant(0) # generally won't change this 
    self.active_site_absorb = Constant(0)
    self.bulk_conc = Constant(1.0)
    self.D = 100.0 # diffusion constant, [um^2/s]
    self.beta = 1/0.693 # 1/kT, [kcal/mol]
    self.valence = Constant(2)

    # temporary
    self.temp_outerR = 5.0
    self.temp_innerR = 1.0
    self.temp_siteZ  = 0.0

