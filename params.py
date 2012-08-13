## VARIABLES
from dolfin import Constant



class params:
  def __init__(self):
    # CONVERTSION
    self.Ang_to_um  = 1.0e-4  # Ang into um
    self.um_to_Ang  = 1/self.Ang_to_um # Ang into um
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
    self.D = 780.0 # diffusion constant, [um^2/s]
    self.Dchannel = 0.6*self.D # diffusion constant, [um^2/s]
    #self.D = 100.0 # diffusion constant, [um^2/s]
    self.beta = 1/0.693 # 1/kT, [kcal/mol]
    self.valence = 2.0   # why do I need a constant? TODO  Constant(2)

    # validation 
    self.Rsphere = 12.5e-4 # 8 A --> 8e-4 um
    self.qsphere = 1  # 8 A --> 8e-4 um

    # temporary
    self.temp_outerR = 5.0
    self.temp_innerR = 1.0
    self.temp_siteZ  = 0.0

    # DETERMINE BY SETTING D=780 and valence=0 AND COMPARING AGAINST PREV. VERSION OF SMOL 
    # TO MATCH SMOL self.correction = -61568.52   
    #self.correction = -61568.52   
    # TO MATCH ANAL 
    self.correction =-67454.54
    #print "WARNING: CHANGING PARAMETERS FOR DEBUG"
    #self.D = 780;
    #self.valence = 1.0
    #self.valence = -1.0
    #self.valence = 0.0

    ## report 
    import socket
    import dolfin
    self.hostname=socket.gethostname()
    self.dolfinver = dolfin.__version__
    print "Running on %s using dolfin version %s" % (self.hostname,self.dolfinver)     

