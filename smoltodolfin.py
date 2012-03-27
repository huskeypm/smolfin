#
# Revisions
# 120312 Removing Nans in interpolated mesh 
#
import numpy
from apbstodolfin import *
from mcsftodolfin import *
from smol import * 

writePotentialOnly=0
interpMethod="linear"
#interpMethod="nearest"

def write_smol_files(filename, mesh,u,writePotentialOnly=0): # , u,vertexmarkers,facetmarkers):
    from dolfin import File

    # for dolfin
    if(writePotentialOnly==0):
      File(filename+"_mesh.xml.gz") << mesh

    File(filename+"_values.xml.gz") << u
    # for viewing in paraview 
    File(filename+"_values.pvd") << u

def ClipValues(pot, THRESH=999):
    idx = (np.where(pot < -THRESH))[0]
    if(np.size(idx)>0):
      print "Encountered %d values with very small potential valus. Clipping to %f" % (np.size(idx),-THRESH)
      pot[ idx ] = -THRESH

    # to large 
    idx = (np.where(pot > THRESH))[0]
    if(np.size(idx)>0):
      print "Encountered %d values with very high potential values. Clipping to %f" % (np.size(idx),THRESH)
      pot[ idx ] = THRESH
  


# interpolate Finite difference APBS mesh onto finite element mesh
# replaces nan with 0.0
#
# useSubset - only interpolates FE points that are within the APBS grid dimensions (e.g. will not assign 0 elsewhere) 
# mvalues - add in potential to preexisting mesh. Allows multiple FD grids to be applied
# mgridloc - describes where center of molecule is in the FE grid. Useful when the molecule is not centered
# debugValues - replaces potential with scalar to test interpolation 
# 
def interpAPBS(mesh,apbs,usesubset=0,mvalues=0,debugValues=0,mgridloc=[0,0,0]):
    from scipy.interpolate import griddata

    ## cleaning up input data 
    # I suspect that very large/small data points are affecting the interpolation   
    # I am trying to rescale the potential here to clip any points beyond a certain range 
    pot = np.array(apbs.values)
    print "Min/Max of potential are (%e,%e)" % (
        np.min(pot),
        np.max(pot))


    # get too small
    ClipValues(pot,THRESH=999)

    apbs.values[:] = pot[:]

    ## shifting apbs values to align with grid 
    apbs.meshcoordinates = apbs.coordinates + mgridloc
    

    ## interpolations 
    mcoordinates = mesh.coordinates()
    numCoor = len( mcoordinates[:,0])

    # replace subset
    if(usesubset==1):
      # find min/max of apbs 
      apbs.meshmin = np.min(apbs.meshcoordinates,axis=0)
      apbs.meshmax = np.max(apbs.meshcoordinates,axis=0)
      print "APBS(mesh) min: "
      print apbs.meshmin
      print np.min(apbs.coordinates,axis=0)
      print "APBS(mesh) max"
      print apbs.meshmax
      print np.max(apbs.coordinates,axis=0)
       
      #apbs.min = np.array([-9,-9,-9])  
      #apbs.max = -1 * apbs.min
     
      inside=[]
      for i in range(0,numCoor):
        p = mcoordinates[i,:]
        if((p[0] > apbs.meshmin[0] and p[0] < apbs.meshmax[0]) and (p[1] > apbs.meshmin[1] and p[1] < apbs.meshmax[1]) and ( p[2] > apbs.meshmin[2] and p[2] < apbs.meshmax[2])):
          #print p
          inside.append(i)

      print "Found %d points inside " % len(inside)
      mcoordsSub = mcoordinates[inside]
      mvaluesSub = griddata(apbs.meshcoordinates, apbs.values,(mcoordsSub),method=interpMethod)
      mvalues[ inside ] = mvaluesSub

      if(debugValues!=0):
        mvalues[ inside ] = debugValues
   
    else:
      mvalues = griddata(apbs.meshcoordinates, apbs.values,(mcoordinates),method=interpMethod)

    print "Interpolated potential values [kT/e]: min (%e) max (%e) " % (min(mvalues),max(mvalues))


    ## cleaning up bad pixes 
    bad = np.where(np.isnan(mvalues))        
    bad = bad[0]

    
    if(np.size(bad) > 0):
      print "Found %d bad entries, which suggest APBS grid does not overlap with FE grid." % np.size(bad)
      print "Replacing these regions with 0.0 for the time being"
      mvalues[bad] = 0.0 
      #print mvalues

    return mvalues

def do_read_write(problem,apbsfilenames,skipAPBS=0,mgridloc=[0,0,0]):
    #read gamer
    #mcoordinates, mcells, mmarkers,mvertmarkers= read_mcsf_file(mcsffilename)
    #mesh= read_and_mark(mcsffilename,nomark=1)
    mesh = problem.mesh
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
      nvalues = interpAPBS(mesh,apbs,usesubset=1,mvalues=mvalues,mgridloc=mgridloc)#,debugValues=res)    
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

    # Constrain potential to give values no higher than 55M (next to protein in stern layer) 
    # 55 = 1 exp(- q E /kT)
    maxPotent = np.log(55) * (1/parms.beta) / parms.valence
    ClipValues(mvalues,THRESH=maxPotent)
 
     


    # create mesh, etc
    #mesh = generate_dolfin_mesh(mcoordinates, mcells)
    values = generate_dolfin_function(mesh, mvalues)


    # save
    write_smol_files(mcsffilename.replace(".m", ""), mesh, values,writePotentialOnly=writePotentialOnly) #,markedvertices,mmarkers)


if __name__ == "__main__":

    # ~/localTemp/NBCR/smol/apbs_fe/potential-0.dx
    #apbsfilename = ["example/molecule/potential-0.dx"]
    apbsfilenames=[]       
    # ~/localTemp/NBCR/smol/gamer/p.pqr.output.out.m
    #mcsffilename = "example/molecule/p.pqr.output.out.m"
    mcsffilename = "none"
    mgridloc=[0,0,0] # location of molecular center within finite element mesh 

    import sys
    if len(sys.argv) <  3:
        raise RuntimeError("expected an 1) -mcsf mcsf and 2) -p apbs file(s) [in order of increasing resolution] <3) -mgridloc '0 0 0' >")

    for i in np.arange(len(sys.argv)):
      if(sys.argv[i] == '-mcsf'):
        mcsffilename = sys.argv[i+1]
        print "Here"
        problem.mesh= read_and_mark(mcsffilename)

      if(sys.argv[i] == '-mesh'):
        meshfilename = sys.argv[i+1]
        mcsffilename = meshfilename.replace("_mesh.xml.gz",".m")
        print "Mostly for debugging at this point since not fully implemented"
        from dolfin import Mesh
        problem.mesh = Mesh(meshfilename)
        writePotentialOnly=1

        
      if(sys.argv[i] == '-p'):
        apbsfilenames.append(sys.argv[i+1])

      if(sys.argv[i] == '-mgridloc'):
        spl = sys.argv[i+1].split(' ')
        mgridloc=[ float(spl[0]), float(spl[1]), float(spl[2]) ]


    ## report 
    #print "mcsfilename:"
    #print mcsffilename
    print "potentials:"
    print apbsfilenames
    print "grid loc (optional)"
    print mgridloc


    do_read_write(problem,apbsfilenames,mgridloc=mgridloc)
    #do_read_write(mcsffilename,apbsfilename,skipAPBS=1)

