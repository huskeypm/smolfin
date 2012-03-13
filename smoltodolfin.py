#
# Revisions
# 120312 Removing Nans in interpolated mesh 
#
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


# interpolate Finite difference APBS mesh onto finite element mesh
# replace nan with 0.0
def interpAPBS(mesh,apbs):
    from scipy.interpolate import griddata
    mcoordinates = mesh.coordinates()
    mvalues = griddata(apbs.coordinates, apbs.values,(mcoordinates),method='linear')
    print "Interpolated potential values [kT/e]: min (%e) max (%e) " % (min(mvalues),max(mvalues))
    #print mvalues
    bad = np.where(np.isnan(mvalues))        
    bad = bad[0]

    if(np.size(bad) > 0):
      print "Found %d bad entries, which suggest APBS grid does not overlap with FE grid." % np.size(bad)
      print "Replacing these regions with 0.0 for the time being"
      mvalues[bad] = 0.0 
      #print mvalues

    return mvalues

def do_read_write(mcsffilename,apbsfilename,skipAPBS=0):

    # read apbs finite difference mesh 
    acoordinates,avalues = read_fd_apbs_file(apbsfilename);
    # TODO goes in read_fd_apbs
    acoordinates = np.reshape( acoordinates,[np.size(avalues), 3]) 
    apbs = empty()
    apbs.coordinates = acoordinates
    apbs.values= avalues

    # debuging interpolation 
    # seems to be correct, 120312
    #minval = min(avalues)
    #print "min(%e) is at " % minval
    #ind = [i for i, v in enumerate(avalues) if v == minval]
    #loc = acoordinates[ ind[0], : ]
    #print "(%3.1f,%3.1f,%3.1f)" % (loc[0],loc[1],loc[2])
    
    


    #read gamer
    #mcoordinates, mcells, mmarkers,mvertmarkers= read_mcsf_file(mcsffilename)
    #mesh= read_and_mark(mcsffilename,nomark=1)
    mesh= read_and_mark(mcsffilename)
    mcoordinates = mesh.coordinates()

    # have markers for each cell, now need to mark vertices accordingly
    #markedvertices = assign_vertex_markers(mcoordinates, mcells, mmarkers)

    # interpolate apbs values onto mcoordinates grid
    if(skipAPBS!=1):
        mvalues = interpAPBS(mesh,apbs)
        print "Interpolated potential values [kT/e]: min (%e) max (%e) " % (min(mvalues),max(mvalues))
        #print np.isnan(mvalues)
        #print "WARNING: need to verify interpolation is correct (almost sure it is not correct)"
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

    print "WARNING: .M files are repositioned. not sure if APBS is in perfect alignment"

    do_read_write(mcsffilename,apbsfilename)
    #do_read_write(mcsffilename,apbsfilename,skipAPBS=1)

