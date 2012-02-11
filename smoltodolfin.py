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


def do_read_write(mcsffilename,apbsfilename,skipAPBS=0):

    # read apbs finite difference mesh 
#    acoordinates,avalues = read_fd_apbs_file(apbsfilename);

    #read gamer
    #mcoordinates, mcells, mmarkers,mvertmarkers= read_mcsf_file(mcsffilename)
    mesh,mcoordinates = read_and_mark(mcsffilename,nomark=1)
    quit()

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

