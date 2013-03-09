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
#> For test geom
def define_activesite(mesh):
    coor = mesh.coordinates()

    # find maximum of box and trace out small region
    ztop = np.max(coor[:,2])

    midpoint = np.mean(coor,axis=0)
    xymidpoint = midpoint[0:2]
    range=np.max(coor,axis=0) - np.min(coor,axis=0)
    xywidth = 0.5 * np.min(range) * 1/5.0;

def define_bulkboundary(mesh):
    coor = mesh.coordinates()

    # find maximum of box and trace out small region
    zbottom = np.min(coor[:,2])


# assuming active site is a small radius at top of boundary
def active_site(x):

    xydist = np.linalg.norm(xymidpoint-x[0:2])
    isInXY = xydist < xywidth
    isInZ  = x[2] > (ztop - DOLFIN_EPS)

    return isInXY * isInZ

def bulk_boundary(x):

    return x[2] < (zbottom + DOLFIN_EPS)
