import numpy as np

def read_fe_apbs_file(filename):
    print "Reading APBS file %s" % filename
    # Grab the lines from the file
    lines = open(filename).readlines()
    coordinates = []
    cells = []
    values = []
    line = lines.pop(0)

    # Remove info
    while "data follows" not in line:
        line = lines.pop(0)
        
    # Grab vertices
    line = lines.pop(0).strip()
    while line and "Simplices" not in line:
        coordinates.append(map(float, line.split()))
        if len(coordinates[-1]) != 3:
            raise RuntimeError("expected coordinate line of length 3, "\
                               "got %s"%line)
        line = lines.pop(0).strip()

    print "Read %d vertices"%len(coordinates)

    # Remove info
    while "data follows" not in line:
        line = lines.pop(0)
    line = lines.pop(0).strip()
    while line and "attribute" not in line:
        cells.append(map(int, line.split()))
        if len(cells[-1]) != 4:
            raise RuntimeError("expected cell line of length 4, "\
                               "got %s"%cells[-1])
        
        line = lines.pop(0).strip()

    print "Read %d cells"%len(cells)

    # Remove info
    while "data follows" not in line:
        line = lines.pop(0)
    line = lines.pop(0).strip()
    while "attribute" not in line:
        values.append(float(line))
        line = lines.pop(0).strip()

    # Convert data into NumPy arrays
    coordinates = np.array(coordinates, dtype="d")
    cells = np.array(cells, dtype="I")
    values = np.array(values, dtype="d")

    # Check data
    if coordinates.shape[0] != len(values):
        raise RuntimeError("Expected equal number of values as coordinates")
    if cells.min() != 0:
        raise RuntimeError("Expected the lowest vertex number in the "\
                           "cell array to be 0.")
    if cells.max() != len(values)-1:
        raise RuntimeError("Expected the largest vertex number in the "\
                           "cell array to be %d."%(len(values)-1))
    if len(set(cells.flatten())) != len(values):
        raise RuntimeError("Expected all vertices to belong to a cell.")

    # Return data
    return coordinates, cells, values

def generate_dolfin_mesh(coordinates, cells):
    from dolfin import Mesh, MeshEditor

    # Create a mesh with correct number of vertices and cells
    mesh = Mesh()
    editor = MeshEditor()
    editor.open(mesh, 3, 3)
    editor.init_vertices(coordinates.shape[0])
    editor.init_cells(cells.shape[0])
    editor.close()

    # Assign coordinates and cells
    mesh.coordinates()[:] = coordinates
    mesh_cells = mesh.cells()
    mesh_cells.flags.writeable = True
    mesh_cells[:] = cells

    # Return mesh
    return mesh

def generate_dolfin_function(mesh, values):
    from dolfin import FunctionSpace, Function

    V = FunctionSpace(mesh, "CG", 1)
    u = Function(V)
    u.vector()[:] = values
    
    return u

def write_dolfin_files(filename, mesh, u):
    from dolfin import File
    
    File(filename+"_mesh.xml.gz") << mesh
    File(filename+"_values.xml.gz") << u

def interpolate_to_grid(coordinates, values):
    import numpy 
    from scipy.interpolate import griddata

    z = numpy.zeros((2,3))
    z[0,:]=(4,21,23)

    line = griddata(coordinates, values, (z),method='linear')


def read_fd_apbs_file(apbsfilename):
    print "Reading APBS file %s" % apbsfilename
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

    

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        raise RuntimeError("expected an apbs file as second argument")
    filename = sys.argv[1]
    gridType = "fd"

    if(gridType=="fe"):
      coordinates, cells, values = read_fe_apbs_file(filename)
    if(gridType=="fd"):
      coordinates, cells, values = read_fd_apbs_file(filename)

    mesh = generate_dolfin_mesh(coordinates, cells)
    values = generate_dolfin_function(mesh, values)
    write_dolfin_files(filename.replace(".dx", ""), mesh, values)

    # attempting interpolation here 
    # filename = "example/molecule/potential-0.dx"
 

    
