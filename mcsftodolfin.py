# Pete Kekenes-Huskey
import numpy as np
def read_mcsf_file(filename):
  lines    =open(filename).readlines()
  
  # header 
  nheader = 10
  # verts
  split   = lines[4].split()
  nverts  = int(split[2].replace(';',''))
  # simplices
  split   = lines[5].split()
  nsimps  = int(split[2].replace(';',''))
  
  coordinates = []
  
  # read coords
  startidx = nheader 
  for i in np.arange(nverts):
    line = lines[i+startidx]
    split = line.split()
  
    coordinates.append(map(float, split[2:5]))
    if len(coordinates[-1]) != 3:
      raise RuntimeError("expected coordinate line of length 3, "\
                         "got %s"%line)
  
  print "Read %d vertices"%len(coordinates)
   
  
  
  # read cells 
  cells = []
  markers = []
  startidx = startidx + nverts + 5 
  for i in np.arange(nsimps):
    line = lines[i+startidx]
    split = line.split()
  
    markers.append(int(split[6]))
    cells.append(map(int, split[7:11]))
    if len(cells[-1]) != 4:
      raise RuntimeError("expected coordinate line of length 3, "\
                         "got %s"%line)
  
  
  print "Read %d cells"%len(cells)
  print "Read %d markers"%len(markers)
  
  
  # convert into numpy arrays 
  coordinates = np.array(coordinates, dtype="d")
  cells = np.array(cells, dtype="I")
  markers = np.array(markers, dtype="I")
  
  # sanity checks
  #print coordinates
  #print cells
  
  # return data 
  return coordinates, cells,markers

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
    File(filename+"_markers.xml.gz") << u
    

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        raise RuntimeError("expected an mcsf file as second argument")
    filename = sys.argv[1]
    coordinates, cells, markers= read_mcsf_file(filename)
    print markers[:]
    #quit()
    mesh = generate_dolfin_mesh(coordinates, cells)
    # can't store this just yet 
    #values = generate_dolfin_function(mesh, markers)
    write_dolfin_files(filename.replace(".m", ""), mesh, markers)
    
#mcsffile = "p.pqr.output.all.m"
#coordinates,cells, markers = read_mcsf(mcsffile)

