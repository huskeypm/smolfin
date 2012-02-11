## VARIABLES
from dolfin import Constant

# markers
active_site_marker = 1; # verified 
outer_boundary_marker = 5 # Verified on 120104_update.tex
molecular_boundary_marker = 4 # Verified 

# Values 
active_site_absorb = Constant(0)
bulk_conc = Constant(1.0)
noflux_molecular_boundary = Constant(0)
D = 1.0 # diffusion constant

# temporary
temp_outerR = 5.0
temp_innerR = 1.0
temp_siteZ  = 0.0

