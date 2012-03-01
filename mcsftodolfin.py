# Pete Kekenes-Huskey
import numpy as np
#from params import *
import params 
from params import * # must do this for class
parms = params()

def mark_neumann_facets(mesh, cellmarkers):
    #from dolfin import cells, facets,mesh
    from dolfin import cells, facets, MeshFunction, File
    md = mesh.domains()
    facet_markers = md.markers(2)
    
    # 2 for top dim 2 (facets)
    for j, cell in enumerate(cells(mesh)):
       for i, face in enumerate(facets(cell)):
           marker = int(cellmarkers[j,i])
           facet_markers.set_value(cell.index(), i, marker)

    subdomains = MeshFunction("uint", mesh, facet_markers)
    print "Num faces at active_site_marker:", \
          (subdomains.array()==parms.active_site_marker).sum()

    File("subdomains.pvd") << subdomains
    
    return subdomains

def mark_dirichlet_vertices(mesh, vertmarkers):
  from dolfin import VertexFunction

  # This function is not needed, but you can do everything by:
  boundary_markers = VertexFunction("uint", mesh)  
  boundary_markers.set_all(0)
  boundary_markers.array()[:] = vertmarkers

def read_mcsf_file(filename):
  print "Reading %s" % filename
  if not ".m" in filename:
    import sys
    raise RuntimeError("File format not recognized (need .m)")

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
  
  ## read coords
  startidx = nheader 
  for i in np.arange(nverts):
    line = lines[i+startidx]
    split = line.split()
  
    coordinates.append(map(float, split[2:5]))
    if len(coordinates[-1]) != 3:
      raise RuntimeError("expected coordinate line of length 3, "\
                         "got %s"%line)
  
  print "Read %d vertices (of %d) "% (len(coordinates),nverts)
   
  
  
  ## read cells 
  cells = []
  cellmarkers = []
  vertexmarkers = np.zeros( nverts)

  startidx = startidx + nverts + 5 
  for i in np.arange(nsimps):
    line = lines[i+startidx]
    split = line.split()
    markerlist = np.array(split[3:7],dtype="I")
    vertexlist = np.array(split[7:11],dtype="I")

    ## handle vertex markers 
    # not all markers start at '6', so we need to check all
    #markers.append(int(split[6]))
    marker = max(markerlist)
    cellmarkers.append(map(int,markerlist))
    #print markerlist
    #print map(int,markerlist)

    # (zero indexed both here and in 'm' file)
    #print vertexlist
    #print markerlist
    markedmarkerlist = np.where(markerlist!=marker) #
    #print markedmarkerlist
    markedvertices  = vertexlist[ markedmarkerlist ]
    if(marker!=0):  # do this to avoid overwriting others
      #print marker 
      #print markedvertices
      vertexmarkers[ markedvertices ] = marker
   
    

    ## handle cell markers 
    cells.append(map(int, vertexlist))  
    if len(cells[-1]) != 4:
      raise RuntimeError("expected coordinate line of length 3, "\
                         "got %s"%line)
  
  
  print "Read %d cells"%len(cells)
  #print "Read %d nonzero, marked vertices" % len( np.where(vertexmarkers!=0)[0] )
  #print vertexmarkers
  
  
  # convert into numpy arrays 
  coordinates = np.array(coordinates, dtype="d")
  cells = np.array(cells, dtype="I")
  cellmarkers = np.array(cellmarkers, dtype="I")
  vertexmarkers = np.array(vertexmarkers, dtype="I")
  
  # sanity checks
  #print coordinates
  #print cells
  
  # return data 
  return coordinates, cells, cellmarkers, vertexmarkers

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

    print "I think cells need to be marked with '1' here"

    # Return mesh
    return mesh

#def generate_dolfin_function(mesh, values):
#    from dolfin import FunctionSpace, Function
#
#    V = FunctionSpace(mesh, "CG", 1)
#    u = Function(V)
#    u.vector()[:] = values
    
#    return u

def write_dolfin_files(filename, mesh):      
    from dolfin import File,MeshFunction

    print "Filename:", filename
    print "mesh:", mesh.num_vertices(), mesh.num_cells()
    File(filename+"_mesh.xml.gz") << mesh
    #File(filename+"_cellmarkers.xml.gz") << u
    #File(filename+"_vertmarkers.xml.gz") << uv

    # write subdomains
    #sub_domains = MeshFunction("uint", mesh, mesh.topology().dim() - 1)
    sub_domains = MeshFunction("uint", mesh, mesh.domains().markers(2))
    File(filename+"_subdomains.xml.gz") << sub_domains

    # write cell marker file 
    from dolfin import CellFunction
    cells = CellFunction("uint", mesh)
    File(filename+"_cells.xml.gz") << cells       



# no mark skips the domain marking, etc (workaround for a local problem)
def read_and_mark(filename, nomark=0):
    coordinates, cells, cellmarkers, vertmarkers = read_mcsf_file(filename)

    ## 
    mesh = generate_dolfin_mesh(coordinates, cells)
    if (nomark==1):
      print "WARNING: only returning mesh for debugfgin purposes"
      return(mesh)
    
    # generate neumann/dirichlet here
    subdomains = mark_neumann_facets(mesh, cellmarkers)
    #mark_dirichlet_vertices(mesh, vertmarkers)

    # TEST!!      
    from dolfin import MeshFunction, DirichletBC, Constant, FunctionSpace
    V = FunctionSpace(mesh, "CG", 1)

    # marking is actually done inside smo.py
    #bc0 = DirichletBC(V, Constant(active_site_absorb), subdomains,active_site_marker)

    write_dolfin_files(filename.replace(".m", ""), mesh)


    #return mesh,coordinates
    return mesh
    

if __name__ == "__main__":
    ## read data 
    filename = r"xxxx/x.m"
    import sys
    if len(sys.argv) != 2:
        raise RuntimeError("expected an mcsf file as second argument")
    filename = sys.argv[1]
    mesh = read_and_mark(filename)
    
#mcsffile = "p.pqr.output.all.m"
#coordinates,cells, markers = read_mcsf(mcsffile)

