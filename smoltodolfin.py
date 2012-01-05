import numpy
from apbstodolfin import *
from mcsftodolfin import *

def write_smol_files(filename, mesh, u,markers):
    from dolfin import File

    File(filename+"_mesh.xml.gz") << mesh
    File(filename+"_values.xml.gz") << u
    File(filename+"_markers.xml.gz") << markers

# mark vertices with cell markers
def assign_vertex_markers(coordinates,cells,markers):
    markedvertices = np.zeros( np.size(coordinates(:,0)))
    # i believe the first three markers define the facet, while the
    # fourth is the normal, therefore we extract the vertex ids for these
    for i in length(cells):
        if(marker[i]!=0):
            markedvertices[ cells[i,0:2] ] = marker[i]
            
            
    return markedvertices

def do_read_write(mcsffilename,apbsfilename,skipAPBS=0):

    # read apbs
    acoordinates, acells, avalues = read_apbs_file(apbsfilename)
    #read gamer
    mcoordinates, mcells, mmarkers= read_mcsf_file(mcsffilename)

    # have markers for each cell, now need to mark vertices accordingly
    markedvertices = assign_vertex_markers(mcoor, mcells, mmarkers)

    # interpolate apbs values onto mcoordinates grid
    if(skipAPBS!=1):
        from scipy.interpolate import griddata
        mvalues = griddata(acoordinates, avalues,(mcoordinates),method='linear')
    else:
        mvalues = np.zeros( np.size(mcoordinates[:,1]))      

    # create mesh, etc
    mesh = generate_dolfin_mesh(mcoordinates, mcells)
    values = generate_dolfin_function(mesh, mvalues)


    # save
    write_smol_files(mcsffilename.replace(".m", ""), mesh, values,markedvertices)


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        raise RuntimeError("expected an 1) mcsf and 2) apbs file")
    mcsffilename = sys.argv[1]
    apbsfilename = sys.argv[2]

    # ~/localTemp/NBCR/smol/apbs_fe/potential-0.dx
    #apbsfilename = "example/molecule/potential-0.dx"
    # ~/localTemp/NBCR/smol/gamer/p.pqr.output.out.m
    #mcsffilename = "example/molecule/p.pqr.output.out.m"

    #do_read_write(mcsffilename,apbsfilename)
    do_read_write(mcsffilename,apbsfilename,skipAPBS=1)

