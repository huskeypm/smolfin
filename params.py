"""
----------------------------------
Smolfin - a numerical solver of the Smoluchowski equation for interesting geometries
Copyright (C) 2012 Peter Kekenes-Huskey, Ph.D., huskeypm@gmail.com

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA


----------------------------------------------------------------------------
"""
## VARIABLES
from dolfin import Constant



class params:
  def __init__(self):
    # CONVERTSION
    self.Ang_to_um  = 1.0e-4  # Ang into um
    self.um_to_Ang  = 1/self.Ang_to_um # Ang into um
    self.Ang2_to_um2 = 1.0e-8  # Ang^2 into um^2
    self.um3_to_invM = 602214150.0  # um^3 into inverse M (1 cubic micron = 1.0 x 10-15 liter, times Avogadro's number)
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
    # Rsphere is based on the output from mcsftodolfin.py for sphere.m 
    self.Rsphere = 12.25 * self.Ang_to_um # based on 10.25 Ang sphere and 2 Ang ion 
    self.qsphere = 1  # 8 A --> 8e-4 um

    # temporary
    self.temp_outerR = 5.0
    self.temp_innerR = 1.0
    self.temp_siteZ  = 0.0

    ## report 
    import socket
    import dolfin
    self.hostname=socket.gethostname()
    self.dolfinver = dolfin.__version__
    print "Running on %s using dolfin version %s" % (self.hostname,self.dolfinver)     

