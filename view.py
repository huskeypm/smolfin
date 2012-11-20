import matplotlib.pyplot as plt
import numpy as np
class empty:pass

# mesh       - mesh  file (usually file_mesh.xml.gz)
# subdomains - facet file (usually file_subdomains.xml.gz)
def PrintSubdomains(mesh,subdomains,file="subdomains"):
  from dolfin import FunctionSpace, Function, File,DirichletBC,Constant
  V = FunctionSpace(mesh, "CG", 1)

  print "WARNING: still doesn't seem to mark boundaries correctly "

  subdomidx =  list( set(subdomains.array()[:]) )
  bcs = []
  for i, idx in enumerate(subdomidx):
    val = idx * 2.0
    print "Found subdomain idx %d - marking as %f" % (idx,val)
    bci = DirichletBC(V,Constant(val),subdomains,idx)
    bcs.append(bci)

  PrintBoundary(mesh,bcs,file=file)

  


def PrintBoundary(mesh, boundaries,file="marked",dbg=0):
  from dolfin import FunctionSpace, Function, File, plot
  V = FunctionSpace(mesh, "CG", 1)
 
  # make into list
  if(isinstance(boundaries,(list))):
    dalist = boundaries
  else:
    dalist = [boundaries]
  
  marked = Function(V)
  marked.vector()[:]=0
  for i,boundary in enumerate(dalist):
    print "Printing boundary %d" % i
    marked1 = Function(V)
    boundary.apply(marked1.vector())
    tot = np.size(marked1.vector().array())
    whereIdx = np.where(marked1.vector() >0)
    n = np.size(whereIdx) # assuming 0 is unmarked
    print "Marked %d/%d " % (n,tot)
    if(dbg!=0):
      print mesh.coordinates()[whereIdx,:]
  
    marked.vector()[:] +=marked1.vector()[:]

  #plot(marked, interactive=1)
  
  if("pvd" not in file):
    file = file+".pvd"

  print "Printing %s for viewing" % file
  File(file) << marked

  return (n) 

# create figures in style of jounrla formats 
def JournalFig(scale=4): # this is the scaleup needed for DPI
				 # requirements. set to 1 for viewing
				 # interactively
    # size in publication 
    # single col - 3.085 inches 
    w = 3.085

    # dpi 
    dpi = 300

    # golden ratio
    gr = 1.61803
    h = w/gr

    # make fig 
    journalfig = empty()
    journalfig.fig=plt.figure(figsize=(scale*w,scale*h),dpi=dpi)
    journalfig.fontSize=30
 

    ax = plt.gca()
    journalfig.ax = ax 
    fontsize = 18
    for tick in ax.xaxis.get_major_ticks():
      tick.label1.set_fontsize(fontsize)

    for tick in ax.yaxis.get_major_ticks():
      tick.label1.set_fontsize(fontsize)
    
    return (journalfig)
    


def plotslice(problem,result,title="no title",fileName="slice.png",show=0):
    plotslicegeneral(problem.mesh.coordinates(),result.up.vector(),title=title,fileName=fileName,show=show)

def plotslicegeneral(meshcoor,vals,title="no title",fileName="slice.png",show=0,range=50):
    import numpy as np
    #meshcoor = problem.mesh.coordinates()
    femDim = np.shape(meshcoor)[1]
    
    # assuming molecule is within 50 of middle of grid 
    # want 500 points in each dir (resolution)
    numpt = 500
    incr = numpt/ (2 * range) 
    #(grid_x,grid_y,grid_z) = np.mgrid[0:0:1j,-range:range:(incr*1j),-range:range:(incr*1j)]
    (grid_x,grid_y,grid_z) = np.mgrid[0:0:1j,-range:range:(numpt*1j),-range:range:(numpt*1j)]
    from scipy.interpolate import griddata
    #slice = griddata(meshcoor, result.up.vector(), (grid_x, grid_y,grid_z),method="linear")
    if(femDim==3):
      slice = griddata(meshcoor, vals, (grid_x, grid_y,grid_z),method="linear")
    elif(femDim==2):
      slice = griddata(meshcoor, vals, (grid_x, grid_y),method="linear")

    slice[np.isnan(slice)]=0
    x1 = np.linspace(-range,range,numpt)
    x2 = np.linspace(-range,range,numpt)
    X1,X2 = np.meshgrid(x1,x2)

    

    #fig = plt.figure(figsize=(14,10))
    journalfig = JournalFig() 
    plt.title(title,fontsize=journalfig.fontSize)
    plt.ylabel("y [$\AA$]",fontsize=journalfig.fontSize)
    plt.xlabel('z [$\AA$]',fontsize=journalfig.fontSize)
    plt.pcolormesh( X1.T,X2.T,slice.T, shading='flat' )
    plt.colorbar()
    if(show==1):
      plt.show()

    F = plt.gcf()
    F.savefig(fileName)
    print "Plotted %s" % fileName

