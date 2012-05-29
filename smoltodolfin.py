#
# Revisions
# 120312 Removing Nans in interpolated mesh 
#
import numpy
from apbstodolfin import *
from mcsftodolfin import *
from smol import * 
from view import plotslicegeneral

writePotentialOnly=0
interpMethod="linear"
#interpMethod="nearest"

def get_range(mesh):
  return np.ceil(max(np.max(mesh.coordinates(),0) - np.min(mesh.coordinates(),0))/2)

def write_smol_files(filename, mesh,u,writePotentialOnly=0): # , u,vertexmarkers,facetmarkers):
    from dolfin import File

    # for dolfin
    if(writePotentialOnly==0):
      File(filename+"_mesh.xml.gz") << mesh

    File(filename+"_values.xml.gz") << u
    # for viewing in paraview 
    pot = u
    File(filename+"_values.pvd") << pot

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
  


# apbsfilename - dx file from apbs
# coordinates - array of 3d coordinates 
def interpolate_dx(apbsfilename,coordinates,mvalues=-1):
    print "Reading APBS file %s" % apbsfilename
    if not ".dx" in apbsfilename:
      import sys
      raise RuntimeError("File format not recognized (need dx)")

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

    nEle = np.shape(coordinates)[0]
    if(mvalues==-1):
      mvalues = np.zeros((nEle,3))

    ## do interpolation 
    for i in range(nEle):
      #print coordinates[i,:]
      #print i 
      (x,y,z) = (coordinates[i,0],coordinates[i,1],coordinates[i,2])
      #mvalues.append(vgrid.value((x,y,z)))
      mvalues[i] = vgrid.value((x,y,z))

    ### cleaning up bad pixes 
    #bad = np.where(np.isnan(mvalues))
    #bad = bad[0]
#
#
#    if(np.size(bad) > 0):
#      print "Found %d bad entries, which suggest APBS grid does not overlap with FE grid." % np.size(bad)
#      print "Replacing these regions with 0.0 for the time being"
#      mvalues[bad] = 0.0
#      #print mvalues


    return mvalues

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
    if(len(mgridloc)>1):
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

def InterpolateAPBSFiles(mesh,apbsfilenames,mgridloc=-1):
    if(np.linalg.norm(mgridloc) > 0):
      print "Need to reimplement mgridloc"
      quit()


    print "WARNING: there is still something funkly going on with the interpolation. This needs to be fixed!"
    # Need to make simple (small number of vertices) dx file and coordinate, then compare smol version versus my version, value by value 
    # It's possible that there's something about the dx reader that I don't understand or am not using correctly, 
    # so maybe even need to very that the coordinates are the same between both dx files when read in  
    quit()
    
    mcoordinates = mesh.coordinates()
    mvalues = np.zeros( len(mcoordinates[:,0]) ) 
    coverage= np.zeros( len(mcoordinates[:,0]) ) 
    prevRes =9999;
    range = get_range(mesh)
 
    print "WARNING: Be sure to read lowest resolution mesh first, second lowest next, etc"


    # starting from lowest resolution and working toward high, 
    # interpolate coordinates that exist within the apbs grid  
    # all other values are assigned nan. We'll replace all non-nan
    # values
    i=0
    for apbsfilename in apbsfilenames:
      # do interpolation 
      values = interpolate_dx(apbsfilename,mcoordinates)
  
      # only pull out the values that were interpolated 
      goodidx = np.where(np.isnan(values)==False)
      goodidx = goodidx[0]
     
      # store interpolated values 
      mvalues[goodidx] = values[goodidx]

      # here we just mark the interpolated mesh do determine where we actually 
      i+=1
      coverage[goodidx] = i * 1.0    

      plotslicegeneral(mesh.coordinates(),mvalues,fileName=apbsfilename+".png",range=range)
      print "Interpolated potential values [kT/e]: min (%e) max (%e) " % (min(mvalues),max(mvalues))

   # plot coverage
    plotslicegeneral(mesh.coordinates(),coverage,fileName="coverage.png",range=range)

    return mvalues



# use for interpolating form FD apbs files 
def InterpolateAPBSFiles_SCIPY(apbsfilenames,mgridloc=-1):
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

      return mvalues 

def readapbscsv(mesh,csvfilename):
    #print "WARNING: this does not attempt to check that your vertices match between apbs and .m!!!!!"
    print "WARNING: expects spaces"
    lines = open(csvfilename).readlines()
    all = []
    for line in lines:    
      all.append(map(float, line.split()))  

    all= np.array(all, dtype="d")
    #coordinates = all[:,
    avalues = all[:,3]
    acoords = all[:,0:3]

    from scipy.interpolate import griddata
    mvalues= griddata(acoords, avalues, (mesh.coordinates()),method="nearest")

    bad = np.where(np.isnan(mvalues))
    bad = bad[0]


    if(np.size(bad) > 0):
      print "Found %d bad entries, which suggest APBS grid does not overlap with FE grid." % np.size(bad)
      # shouldn't have this issue do definitely die
      quit()



    print "min/max %f %f " % (np.min(mvalues), np.max(mvalues))
    return mvalues

# csvfilename - provide csv file of interpolated values (see 120327_troubleshoot.tex)
def do_read_write(problem,apbsfilenames,skipAPBS=0,mgridloc=-1,csvfilename="none"):
    #read gamer
    #mcoordinates, mcells, mmarkers,mvertmarkers= read_mcsf_file(mcsffilename)
    #mesh= read_and_mark(mcsffilename,nomark=1)
    mesh = problem.mesh
    mcoordinates = mesh.coordinates()
    mvalues = np.zeros( len(mcoordinates[:,0]) ) 


    # hack
    if(csvfilename!="none"): 
      mvalues = readapbscsv(mesh,csvfilename)
    
    else:
      mvalues = InterpolateAPBSFiles(mesh,apbsfilenames,mgridloc=mgridloc)


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
 
     
    from view import plotslicegeneral
    range = get_range(mesh)
    plotslicegeneral(mesh.coordinates(),mvalues,fileName="final.png",range=range)
    print "Final: Interpolated potential values [kT/e]: min (%e) max (%e) " % (min(mvalues),max(mvalues))


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
    csvfilename = "none"
    mgridloc=[0,0,0] # location of molecular center within finite element mesh 

    import sys
    if len(sys.argv) <  3:
        raise RuntimeError("expected an 1) -mcsf/mesh mcsf/mesh file and 2) -p apbs file(s) [in order of increasing resolution] < -mgridloc '0 0 0' -csv csvapbs>")

    for i in np.arange(len(sys.argv)):
      if(sys.argv[i] == '-mcsf'):
        mcsffilename = sys.argv[i+1]
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

      if(sys.argv[i] == '-csv'):
        csvfilename = sys.argv[i+1]

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


    do_read_write(problem,apbsfilenames,mgridloc=mgridloc,csvfilename=csvfilename)
    #do_read_write(mcsffilename,apbsfilename,skipAPBS=1)

