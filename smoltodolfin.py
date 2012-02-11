import numpy
from apbstodolfin import *
from mcsftodolfin import *

def write_smol_files(filename, mesh,u): # , u,vertexmarkers,facetmarkers):
    from dolfin import File

    File(filename+"_mesh.xml.gz") << mesh
    File(filename+"_values.xml.gz") << u
    #File(filename+"_dirichmarkers.xml.gz") << vertexmarkers
    #File(filename+"_neumannmarkers.xml.gz") << facetmarkers

# mark vertices with cell markers
#def assign_vertex_markers(coordinates,cells,markers):
#    markedvertices = np.zeros( np.size(coordinates[:,0]))
#    # i believe the first three markers define the facet, while the
#    # fourth is the normal, therefore we extract the vertex ids for these ((NOT ALWAYS)!!)
#    for i in length(cells):
#        if(marker[i]!=0):
#            markedvertices[ cells[i,0:2] ] = marker[i]
#            
#            
#    return markedvertices
#

def read_fd_apbs_file(apbsfilename):
    # read in data 
    import sys, random
    from potentialextraction import  Vgrid
    file = open(apbsfilename, "r")
    dims = (None, None, None)
    spac = (None, None, None)
    origin = (None, None, None)
    n = None
    data = []
    vgrid = Vgrid(dims, spac, origin, data)
    vgrid.readOpenDX(file)
    file.close()

    # assign to numpy ar
    nx = vgrid.dims[0]
    ny = vgrid.dims[1]
    nz = vgrid.dims[2]
    hx = vgrid.spac[0]
    hy = vgrid.spac[1]
    hz = vgrid.spac[2]
    xmin = vgrid.origin[0]
    ymin = vgrid.origin[1]
    zmin = vgrid.origin[2]
    xlen = float(nx-1)*hx
    ylen = float(ny-1)*hy
    zlen = float(nz-1)*hz
    xmax = xmin + xlen
    ymax = ymin + ylen
    zmax = zmin + zlen
    xcent = xmin + 0.5*xlen
    ycent = ymin + 0.5*ylen
    zcent = zmin + 0.5*zlen


    coordinates = []
    values = []

    # read coords to arrays
    for i in range(1, (nx-1)):
        x = xmin + float(i)*hx

        for j in range(1, (ny-1)):
            y = ymin + float(j)*hy

            for k in range(1, (nz-1)):
                z = zmin + float(k)*hz
                coordinates.append([x,y,z])
                values.append(vgrid.value((x, y, z)))

    print "Read %d coords " % len(coordinates)

    return(coordinates,values)



def do_read_write(mcsffilename,apbsfilename,skipAPBS=0):

    # read apbs finite difference mesh 
    acoordinates,avalues = read_fd_apbs_file(apbsfilename);

    #read gamer
    #mcoordinates, mcells, mmarkers,mvertmarkers= read_mcsf_file(mcsffilename)
    mesh,mcoordinates = read_and_mark(mcsffilename)

    # have markers for each cell, now need to mark vertices accordingly
    #markedvertices = assign_vertex_markers(mcoordinates, mcells, mmarkers)

    # interpolate apbs values onto mcoordinates grid
    if(skipAPBS!=1):
        from scipy.interpolate import griddata
        mvalues = griddata(acoordinates, avalues,(mcoordinates),method='linear')
    else:
        mvalues = np.zeros( np.size(mcoordinates[:,1]))      

    # create mesh, etc
    #mesh = generate_dolfin_mesh(mcoordinates, mcells)
    values = generate_dolfin_function(mesh, mvalues)


    # save
    write_smol_files(mcsffilename.replace(".m", ""), mesh, values) #,markedvertices,mmarkers)


if __name__ == "__main__":

    # ~/localTemp/NBCR/smol/apbs_fe/potential-0.dx
    apbsfilename = "example/molecule/potential-0.dx"
    # ~/localTemp/NBCR/smol/gamer/p.pqr.output.out.m
    mcsffilename = "example/molecule/p.pqr.output.out.m"

    import sys
    if len(sys.argv) != 3:
        raise RuntimeError("expected an 1) mcsf and 2) apbs file")
    mcsffilename = sys.argv[1]
    apbsfilename = sys.argv[2]


    do_read_write(mcsffilename,apbsfilename)
    #do_read_write(mcsffilename,apbsfilename,skipAPBS=1)

