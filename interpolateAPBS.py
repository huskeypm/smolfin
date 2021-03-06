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
#
# Revisions
# 120727 Major bug fixed with interpolating potential values 
# 120312 Removing Nans in interpolated mesh 
#
import numpy
from apbstodolfin import *
from mcsftodolfin import *
from smol import * 
from view import plotslicegeneral

writePotentialOnly=False
interpMethod="linear"
#interpMethod="nearest"

def get_range(mesh):
  return np.ceil(max(np.max(mesh.coordinates(),0) - np.min(mesh.coordinates(),0))/2)

def write_smol_files(filename, mesh,u,writePotentialOnly=False): # , u,vertexmarkers,facetmarkers):
    from dolfin import File

    # for dolfin
    if(writePotentialOnly==False):
      fileMesh = filename+"_mesh.xml.gz"
      File(fileMesh) << mesh

    filePotential = filename+"_values.xml.gz"
    File(filePotential) << u
    # for viewing in paraview 
    pot = u
    File(filename+"_values.pvd") << pot


    # Test written file 
    # PKH 120731
    #print "WARNING: there is some imprecision I think in how interplation is done in dolfin"
    #print "If I load mesh, then I get different inpolation results then If i use the current"
    #print "DEBUG"
    #mesh = Mesh(fileMesh)
    # THIS WORKS 
    #V = FunctionSpace(mesh, "CG", 1)
    #psi = Function(V,filePotential);
    #z = Function(V)
    #z.vector()[:] = psi.vector()[:]
    #File("psimcsf.pvd") << z
    #File("psimcsf_orig.pvd") << z

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

    # get mesh dimensions 
    femDim = np.shape(coordinates)[1]

    nEle = np.shape(coordinates[:,0])[0]
    if(mvalues==-1):
      mvalues = np.zeros(nEle)

    ## do interpolation 
    if(femDim==3):
      for i in range(nEle):
        (x,y,z) = (coordinates[i,0],coordinates[i,1],coordinates[i,2])
        #mvalues.append(vgrid.value((x,y,z)))
        mvalues[i] = vgrid.value((x,y,z))
    elif(femDim==2):
      for i in range(nEle):
        (x,y,z) = (coordinates[i,0],coordinates[i,1],0)
        mvalues[i] = vgrid.value((x,y,z))
    else: 
      raise RuntimeError("Unknown dimension")


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
def InterpolateAPBSdx(mesh,apbs,usesubset=0,mvalues=0,debugValues=0,mgridloc=[0,0,0]):
    from scipy.interpolate import griddata

    ## cleaning up input data 
    # I suspect that very large/small data points are affecting the interpolation   
    # I am trying to rescale the potential here to clip any points beyond a certain range 
    pot = np.array(apbs.values)
    print "Min/Max of potential are (%e,%e)" % (
        np.min(pot),
        np.max(pot))


    # get too small
    if clipValues:
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

def InterpolateAPBSdxs(mesh,apbsfilenames,mgridloc=-1):

    # CHeck on offset parameter 
    if(np.size(mgridloc) != 1):
      print "Need to reimplement mgridloc"
      print mgridloc
      quit()


    # load in coordinates 
    mcoordinates = mesh.coordinates()

    # debug 
    if(0):
      mcoordinates = np.array([[-6.1909236908e+00,    1.2481193542e+01,   -5.8969097137e+00]])

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

      #plotslicegeneral(mesh.coordinates(),mvalues,fileName=apbsfilename+".png",range=range)
      print "Interpolated potential values [kT/e]: min (%e) max (%e) " % (min(mvalues),max(mvalues))

    #print values

   # plot coverage
    #plotslicegeneral(mesh.coordinates(),coverage,fileName="coverage.png",range=range)

    return mvalues



# use for interpolating form FD apbs files 
def InterpolateAPBSdxs_SCIPY(apbsfilenames,mgridloc=-1):
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
      nvalues = InterpolateAPBSdx(mesh,apbs,usesubset=1,mvalues=mvalues,mgridloc=mgridloc)#,debugValues=res)    
      mvalues = nvalues
      #print "%d->%d" % (len(zidx),len(zidxn))
    
      # hack to visualize
      #from view import plotslicegeneral
      #plotslicegeneral(mesh.coordinates(),mvalues,fileName=apbsfilename+".png")
      print "Interpolated potential values [kT/e]: min (%e) max (%e) " % (min(mvalues),max(mvalues))
          #print np.isnan(mvalues)
          #print "WARNING: need to verify interpolation is correct (almost sure it is not correct)"
      #else:
      #    mvalues = np.zeros( np.size(mcoordinates[:,1]))      

      return mvalues 

def InterpolateAPBScsv(mesh,csvfilename):
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
def convertAllAPBS(problem,apbsfilenames,skipAPBS=0,mgridloc=-1,\
    csvfilename=None,writePotentialOnly=writePotentialOnly,clipValues=True):

    #read gamer
    #mcoordinates, mcells, mmarkers,mvertmarkers= read_mcsf_file(mcsffilename)
    #mesh= read_and_mark(mcsffilename,nomark=1)
    mesh = problem.mesh

    # bug hack - I think there is a discrepency in how interpolations are done
    # for orig. mesh. Here I will save the current mesh, then reload to ensure
    # the coords are the same
    # PKH 120731
    print "NOTE: Need to save mesh, then load 'interpolated', otherwise potential will disagree"
    fileMesh = "temp_mesh.xml.gz"
    File(fileMesh) << mesh
    mesh = Mesh(fileMesh)


    mcoordinates = mesh.coordinates()
    mvalues = np.zeros( len(mcoordinates[:,0]) ) 

    #print apbsfilenames
    #print mgridloc
    #quit()

    # hack
    if(csvfilename!=None):    
      mvalues = InterpolateAPBScsv(mesh,csvfilename)
    
    # correct way 
    else:
      mvalues = InterpolateAPBSdxs(mesh,apbsfilenames,mgridloc=mgridloc)


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
    # PKH 120731
    maxPotent = np.log(55.) * (1/parms.beta) / np.abs(parms.valence)
    if clipValues:
      ClipValues(mvalues,THRESH=maxPotent)
 
     
    from view import plotslicegeneral
    range = get_range(mesh)
    #plotslicegeneral(mesh.coordinates(),mvalues,fileName="final.png",range=range)
    print "Final: Interpolated potential values [kT/e]: min (%e) max (%e) " % (min(mvalues),max(mvalues))


    # create mesh, etc
    #mesh = generate_dolfin_mesh(mcoordinates, mcells)
    values = generate_dolfin_function(mesh, mvalues)


    # save
    write_smol_files(problem.root, mesh, values,writePotentialOnly=writePotentialOnly) #,markedvertices,mmarkers)


if __name__ == "__main__":

    clipValues=True;
    # ~/localTemp/NBCR/smol/apbs_fe/potential-0.dx
    #apbsfilename = ["example/molecule/potential-0.dx"]
    apbsfilenames=[]       
    # ~/localTemp/NBCR/smol/gamer/p.pqr.output.out.m
    #mcsffilename = "example/molecule/p.pqr.output.out.m"
    mcsffilename = None   
    csvfilename = None   
    mgridloc=[0,0,0] # location of molecular center within finite element mesh 
    mgridloc=-1      # location of molecular center within finite element mesh 

    import sys
    if len(sys.argv) <  3:
        raise RuntimeError("expected an 1) -mcsf/mesh mcsf/mesh file and 2) -p apbs file(s) [in order of increasing resolution] < -mgridloc '0 0 0' -csv csvapbs -noclipping>")

    for i in np.arange(len(sys.argv)):
      if(sys.argv[i] == '-mcsf'):
        mcsffilename = sys.argv[i+1]
        problem.root = mcsffilename.replace(".m", "")

      if(sys.argv[i] == '-mesh'):
        meshfilename = sys.argv[i+1]
        #mcsffilename = meshfilename.replace("_mesh.xml.gz",".m")
        from dolfin import Mesh
        problem.root = meshfilename.replace("_mesh.xml.gz", "")
        problem.mesh = Mesh(meshfilename)
        writePotentialOnly=True

        
      if(sys.argv[i] == '-p'):
        apbsfilenames.append(sys.argv[i+1])

      if(sys.argv[i] == '-csv'):
        csvfilename = sys.argv[i+1]

      if(sys.argv[i] == '-mgridloc'):
        spl = sys.argv[i+1].split(' ')
        mgridloc=[ float(spl[0]), float(spl[1]), float(spl[2]) ]
      if(sys.argv[i] == '-noclipping'):
        clipValues=False

      
    if(len(apbsfilenames)==0 and csvfilename=='none'):
        raise RuntimeError("Must provide apbs electrostatic potential files. Otherwise use mcsftodolfin.py")

    if(mcsffilename!=None):       
        problem.mesh= read_and_mark(mcsffilename)


    ## report 
    #print "mcsfilename:"
    #print mcsffilename
    print "Root name ", problem.root
    print "potentials:"
    print apbsfilenames
    print "grid loc (optional)"
    print mgridloc


    convertAllAPBS(problem,apbsfilenames,mgridloc=mgridloc,\
      csvfilename=csvfilename,writePotentialOnly=writePotentialOnly,clipValues=clipValues)
    #convertAllAPBS(mcsffilename,apbsfilename,skipAPBS=1)

