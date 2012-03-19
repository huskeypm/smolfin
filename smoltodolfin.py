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
    File(filename+"_values.pvd") << u
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
# mvalues - add in potential to preexisting mesh. Allows multiple FD grids to be applied
# debug val - replaces potential with scalar to test interpolation 
def interpAPBS(mesh,apbs,usesubset=0,mvalues=0,debugValues=0):
    from scipy.interpolate import griddata

    # I suspect that very large/small data points are affecting the interpolation   
    # I am trying to rescale the potential here to clip any points beyond a certain range 
    pot = np.array(apbs.values)
    print "Min/Max of potential are (%e,%e)" % (
        np.min(pot),
        np.max(pot))

    # get too small
    THRESH=2 
    idx = (np.where(pot < -THRESH))[0]
    if(np.size(idx)>0):
      print "Encountered %d values with very small potential valus. Clipping to %f" % (np.size(idx),THRESH)
      pot[ idx ] = -THRESH

    # to large 
    idx = (np.where(pot > THRESH))[0]
    if(np.size(idx)>0):
      print "Encountered %d values with very high potential values. Clipping to %f" % (np.size(idx),THRESH)
      pot[ idx ] = THRESH

    apbs.values[:] = pot[:]

    # now interp
    mcoordinates = mesh.coordinates()
    numCoor = len( mcoordinates[:,0])

    # replace subset
    if(usesubset==1):
      # find min/max of apbs 
      apbs.min = np.min(apbs.coordinates,axis=0)
      apbs.max = np.max(apbs.coordinates,axis=0)
      print "APBS min"
      print apbs.min
      print "APBS max"
      print apbs.max
       
      #apbs.min = np.array([-9,-9,-9])  
      #apbs.max = -1 * apbs.min
     
      inside=[]
      for i in range(0,numCoor):
        p = mcoordinates[i,:]
        if((p[0] > apbs.min[0] and p[0] < apbs.max[0]) and (p[1] > apbs.min[1] and p[1] < apbs.max[1]) and ( p[2] > apbs.min[2] and p[2] < apbs.max[2])):
          #print p
          inside.append(i)

      print "Found %d points inside " % len(inside)
      mcoordsSub = mcoordinates[inside]
      mvaluesSub = griddata(apbs.coordinates, apbs.values,(mcoordsSub),method='linear')
      mvalues[ inside ] = mvaluesSub

      if(debugValues!=0):
        mvalues[ inside ] = debugValues
   
    else:
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

def do_read_write(mcsffilename,apbsfilenames,skipAPBS=0):
    #read gamer
    #mcoordinates, mcells, mmarkers,mvertmarkers= read_mcsf_file(mcsffilename)
    #mesh= read_and_mark(mcsffilename,nomark=1)
    mesh= read_and_mark(mcsffilename)
    mcoordinates = mesh.coordinates()
    mvalues = np.zeros( len(mcoordinates[:,0]) ) 

    prevRes =9999;
    for apbsfilename in apbsfilenames:
      #print "Reading %s" % apbsfilename 
     
      # read apbs finite difference mesh 
      #acoordinates,avalues = read_fd_apbs_file(apbsfilename);
      apbs = read_fd_apbs_file(apbsfilename);
      # TODO goes in read_fd_apbs
      #apbs = empty()
      #apbs.coordinates = acoordinates
      #apbs.values= avalues
      res = min(apbs.res)
      if(res > prevRes):
        msg = "Grid resolution must be HIGHER (smaller grid spacing) for each listed file"
        # TODO turn on this alarm 
        raise RuntimeError(msg)
      prevRes =res
  
  
      # have markers for each cell, now need to mark vertices accordingly
      #markedvertices = assign_vertex_markers(mcoordinates, mcells, mmarkers)
  
      # interpolate apbs values onto mcoordinates grid
      #if(skipAPBS!=1):
      nvalues = interpAPBS(mesh,apbs,usesubset=1,mvalues=mvalues)#,debugValues=res)    
      mvalues = nvalues
      #print "%d->%d" % (len(zidx),len(zidxn))
    
      # hack to visualize
      from view import plotslicegeneral
      plotslicegeneral(mesh.coordinates(),mvalues,fileName=apbsfilename+".png")
      print "Interpolated potential values [kT/e]: min (%e) max (%e) " % (min(mvalues),max(mvalues))
          #print np.isnan(mvalues)
          #print "WARNING: need to verify interpolation is correct (almost sure it is not correct)"
      #else:
      #    mvalues = np.zeros( np.size(mcoordinates[:,1]))      

    # check that i can still create a mesh/ something weird was happening in 
    # interactive session
    #print "Teston ly"
    #from dolfin import File,Function,FunctionSpace
    #File("somfile.pvd") << Function(FunctionSpace(mesh, "CG", 1))
    #r = Function(FunctionSpace(mesh, "CG", 1))
    #r.vector()[:] = mvalues
    #File("sominterp.pvd") << r  


    # create mesh, etc
    #mesh = generate_dolfin_mesh(mcoordinates, mcells)
    values = generate_dolfin_function(mesh, mvalues)


    # save
    write_smol_files(mcsffilename.replace(".m", ""), mesh, values) #,markedvertices,mmarkers)


if __name__ == "__main__":

    # ~/localTemp/NBCR/smol/apbs_fe/potential-0.dx
    apbsfilename = ["example/molecule/potential-0.dx"]
    # ~/localTemp/NBCR/smol/gamer/p.pqr.output.out.m
    mcsffilename = "example/molecule/p.pqr.output.out.m"

    import sys
    if len(sys.argv) <  3:
        raise RuntimeError("expected an 1) mcsf and 2) apbs file(s) [in order of increasing resolution]")

    mcsffilename = sys.argv[1]

    # list all possible apbs args
    #apbsfilename = sys.argv[2]
    apbsfilenames = [sys.argv[i] for i in range(2, len(sys.argv))]

    print "WARNING: .M files are repositioned. not sure if APBS is in perfect alignment"

    do_read_write(mcsffilename,apbsfilenames)
    #do_read_write(mcsffilename,apbsfilename,skipAPBS=1)

