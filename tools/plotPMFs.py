"
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
"
from view import *
import numpy as np
files=[]

files.append("/home/huskeypm/scratch/validation/serca/serca.pmf")
dat0 = np.loadtxt(files[0]);
# for serca only to shift cood
newx =-0.9533 * dat0[:,0] + 25.789
dat0[:,0] = newx
# shifting so aligned at x=0
dat0[:,0] = dat0[:,0] - 3.184
# shifting so aligned at y=0
dat0[:,1] = dat0[:,1] -  2.0479


files.append("/home/huskeypm/scratch/validation/tnc_isolated/tnc_isolated.pmf")
dat1 = np.loadtxt(files[1]);
# shifting so aligned at x=0
dat1[:,0] = dat1[:,0] - 2.1
# shifting so aligned at y=0
dat1[:,1] = dat1[:,1] -  9.45  

journalfig = JournalFig()
plt.title("Potentials of mean force",fontsize=journalfig.fontSize)
plt.ylabel("PMF [kcal/mol]",fontsize=journalfig.fontSize)
plt.xlabel('Active Site Distance, $\\vec x$ [$\AA$]',fontsize=journalfig.fontSize)
p1, = plt.plot(dat0[:,0],dat0[:,1], 'k-', color='black')
p2, = plt.plot(dat1[:,0],dat1[:,1], 'k--', color='black')
ax = journalfig.ax
#handles, labels = ax.get_legend_handles_labels()
ax.legend([p1,p2],["SERCA","TnC"])

plt.xlim((0,10))

F = plt.gcf()
F.savefig("fig_pmfs.png")


