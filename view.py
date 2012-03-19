def PrintBoundary(mesh, boundary,file="marked"):
  from dolfin import FunctionSpace, Function, File, plot
  V = FunctionSpace(mesh, "CG", 1)
  
  marked = Function(V)
  boundary.apply(marked.vector())

  #nummrked = (subdomains.array()==1).sum()
  #print "Num marked entries:"% nummrked


  #plot(marked, interactive=1)
  
  if("pvd" not in file):
    file = file+".pvd"

  print "Printing %s for viewing" % file
  File(file) << marked


def plotslice(problem,result,title="no title",fileName="slice.png"):
    plotslicegeneral(problem.mesh.coordinates(),result.up.vector(),title=title,fileName=fileName)

def plotslicegeneral(meshcoor,vals,title="no title",fileName="slice.png"):
    import numpy as np
    #meshcoor = problem.mesh.coordinates()
    
    # assuming molecule is within 50 of middle of grid 
    # want 500 points in each dir (resolution)
    range = 50
    numpt = 500
    incr = numpt/ (2 * range) 
    #(grid_x,grid_y,grid_z) = np.mgrid[0:0:1j,-range:range:(incr*1j),-range:range:(incr*1j)]
    (grid_x,grid_y,grid_z) = np.mgrid[0:0:1j,-range:range:(numpt*1j),-range:range:(numpt*1j)]
    from scipy.interpolate import griddata
    #slice = griddata(meshcoor, result.up.vector(), (grid_x, grid_y,grid_z),method="linear")
    slice = griddata(meshcoor, vals, (grid_x, grid_y,grid_z),method="linear")
    slice[np.isnan(slice)]=0
    x1 = np.linspace(-range,range,numpt)
    x2 = np.linspace(-range,range,numpt)
    X1,X2 = np.meshgrid(x1,x2)

    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(14,10))
    gFontSize = 15
    plt.title(title)
    plt.ylabel("y [\AA]",fontsize=gFontSize)
    plt.xlabel('z [\AA]',fontsize=gFontSize)
    plt.pcolormesh( X1.T,X2.T,slice.T, shading='flat' )
    plt.colorbar()
    F = plt.gcf()
    F.savefig(fileName)
    print "Plotted %s" % fileName
