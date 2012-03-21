# Pete Kekenes-Huskey

#
# Revisions
# 120312 added function to check mesh dimensions to verify agreement with APBS 
#


import numpy as np
from params import * # must do this for class
parms = params()

class empty:pass

# 
from dolfin import SubDomain, DOLFIN_EPS
class NeumannBoundaryHack(SubDomain): 
  def inside(self, x, on_boundary):
    return x[0] < 5 + DOLFIN_EPS and on_boundary
    #return on_boundary

# bit confused her, think all markers are being stored in subdmains, but 
# we are doing somethin separete for subdomains here as well?
def check_facets(mesh, cellmarkers):
    # NOT SURE IF WORKING 

    #from dolfin import cells, facets,mesh
    from dolfin import cells, facets, MeshFunction, File
    md = mesh.domains()
    facet_markers = md.markers(2)
    
    # 2 for top dim 2 (facets)
    unmarked=0
    for j, cell in enumerate(cells(mesh)):
       for i, face in enumerate(facets(cell)):
           marker = int(cellmarkers[j,i])
           if(marker==0): 
             unmarked= unmarked+1
             facet_markers.set_value(cell.index(), i, parms.molecular_boundary_marker)

    print "Found %d unmarked facets. NOT DOING ANYTHING RIGHT NOW" % unmarked

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

    # TESTING 
    # check for errors 
    #numUnmrked = (subdomains.array()==0).sum()
    #print "Num unmarked enetries:"% numUnmrked 
    
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

    # report on mesh 
    # not importing the right things so doesnt work here (worked in smol though)  
    #from dolfin import *         
    #cell = mesh.ufl_cell()
    #c  = Constant(1, cell)
    #print "Total volume %e [A^3]" % assemble(c*dx, mesh=mesh)
    ## each cell
    #cv = []
    #for c in cells(mesh):
    #  v.append( c.volume() )
    #print "Cell vol: Min %e Max %e [A^3]" % (min(v),max(v))



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
    #print "mesh:", mesh.num_vertices(), mesh.num_cells()
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
    cells.set_all(1) # not sure if this is right 
    File(filename+"_cells.xml.gz") << cells       

# gives basic gemetric information for the marked boundaries 
def meshgeoms(mesh,idx=":"):
    coor = mesh.coordinates()
    mol= empty()
    mol.coor = coor[idx]
    mol.min = np.min(mol.coor,axis=0)
    mol.max = np.max(mol.coor,axis=0)
    mol.dim = mol.max - mol.min 
    mol.mid = 0.5 * ( mol.max + mol.min ) 

    return mol 

# checks that mesh is centered and of appropriate dimensions (to ensure 
# we are using appropriate scale and FE mesh matches with APBS FD mesh  
def do_checks(mesh,subdomains):
    from dolfin import DirichletBC, Constant, FunctionSpace, Function, File

    ## get all boundaries corresponding to molecular 
    V = FunctionSpace(mesh, "CG", 1)
    markedall = Function(V)

    # active site
    print "Checking active site..."
    bc0 = DirichletBC(V, Constant(1), subdomains,parms.active_site_marker)
    marked0 = Function(V)
    bc0.apply(marked0.vector())
    markedall.vector()[:] = markedall.vector()[:] + marked0.vector()[:]
    
    # molecular boundary 
    print "Checking molecule boundary ... (may be zero)"
    bc1 = DirichletBC(V, Constant(1), subdomains,parms.molecular_boundary_marker)
    marked1 = Function(V)
    bc1.apply(marked1.vector())
    markedall.vector()[:] = markedall.vector()[:] + marked1.vector()[:]

    # debug
    #File("test.pvd") << markedall
 
    ## get limits of molecule  
    
    # mesh coords
    outeridx = np.where(markedall.vector()[:]  < 10000) # always true..
    outer = meshgeoms(mesh,outeridx)

    # find all indices where molecule exists
    molidx = np.where(markedall.vector()[:] > 0)

    # extremes
    #mol = empty()
    #mol.coor = coor[molidx]
    #mol.min = np.min(mol.coor,axis=0)
    #mol.max = np.max(mol.coor,axis=0)
    #mol.dim = mol.max - mol.min 
    mol = meshgeoms(mesh,idx=molidx)

    actidx = np.where(marked0.vector()[:] > 0)
    act = meshgeoms(mesh,actidx) 

    # check if range contains (0,0,0) which indicates we pass through origin
    # I check for this by noting that the min must be < 0, and max > 0, therefore
    # product of min/max s.b. negative for all axes
    passing = np.size(np.where(mol.min * mol.max < 0))
    if(passing < 3):
      print "WARNING: It appears that the grid does not pass through origin. Please verify and recenter if need be"

    # check sie of grid
    print "Grid size [A]"
    print outer.dim

    # check size of molecule 
    print "Molecule size [A] is "
    print mol.dim
    print "Compare this with original molecule. If incorrect, rescale mesh using XXX and repeat conversion"

    # active site 
    print "Active site center and dims"
    print act.mid
    print act.dim

    




# no mark skips the domain marking, etc (workaround for a local problem)
# rescaleCoor requests isotropic scaling of the mesh coordinates (meshNew = meshOld * rescaleCoor)
def read_and_mark(filename, nomark=0,rescaleCoor=0):
    from dolfin import FacetFunction,MeshFunction, DirichletBC, Constant, FunctionSpace

    # read mesh 
    coordinates, cells, cellmarkers, vertmarkers = read_mcsf_file(filename)

    # rescale mesh 
    if(rescaleCoor > 0):
      print "Rescaling coordinates of molecule by %f" % rescaleCoor
      coordinates = coordinates * rescaleCoor
      
    
    # create dolfin 
    mesh = generate_dolfin_mesh(coordinates, cells)

    # debuf 
    if (nomark==1):
      print "WARNING: only returning mesh for debugfgin purposes"
      return(mesh)

    print "Need to verify meshing is correct again.."

    # generate neumann/dirichlet here
    subdomains = mark_neumann_facets(mesh, cellmarkers)
    #mark_dirichlet_vertices(mesh, vertmarkers)
    #check_facets(mesh,cellmarkers)
    #sub_domains = MeshFunction("uint", mesh, mesh.domains().markers(2))
    #numUnmrked = (subdomains.array()==0).sum()
    #print "sdfs %d" % numUnmrked

    # mark all boundaries 
    fixhack=0
    if(fixhack==1):
      from dolfin import UnitCube
      mesh = UnitCube(8,8,8)
      print "Forcing boundaries to have value at exterior [overridden by file]"
      neumann_boundary = NeumannBoundaryHack() 
      exterior_facet_domains = FacetFunction("uint", mesh) 
      exterior_facet_domains.set_all(parms.unmarked_marker) 
      neumann_boundary.mark(exterior_facet_domains, parms.molecular_boundary_marker)



    V = FunctionSpace(mesh, "CG", 1)

    # marking is actually done inside smo.py
    #bc0 = DirichletBC(V, Constant(active_site_absorb), subdomains,active_site_marker)

    test =1 
    from view import PrintBoundary
    if(test==1):
      # surface w no marker
      bc0 = DirichletBC(V, Constant(1), subdomains,parms.unmarked_marker )
      PrintBoundary(mesh, bc0,file="unmarked")
      # active site
      bc0 = DirichletBC(V, Constant(1), subdomains,parms.active_site_marker)
      PrintBoundary(mesh, bc0,file="active")
      # molecular boundary 
      bc0 = DirichletBC(V, Constant(1), subdomains,parms.molecular_boundary_marker)
      PrintBoundary(mesh, bc0,file="molecular")
      # outer_boundary
      bc0 = DirichletBC(V, Constant(1), subdomains,parms.outer_boundary_marker)
      PrintBoundary(mesh, bc0,file="outer")
      #quit()

    do_checks(mesh,subdomains)


    #write_dolfin_files(filename.replace(".m", ""), mesh, vertmarkers)
    # subdomains.set_all(3) # mark facets as sub domain 3
    write_dolfin_files(filename.replace(".m", ""), mesh)

    # test (this works) 
    #from dolfin import Mesh, FunctionSpace, DirichletBC
    #fileMesh  ="p.pqr.output.out_mesh.xml.gz"
    #mesh = Mesh(fileMesh);
    #V = FunctionSpace(mesh, "CG", 1)
    #fileSubdomains="p.pqr.output.out_subdomains.xml.gz"
    #subdomains = MeshFunction("uint", mesh, fileSubdomains)
    #bc1 = DirichletBC(V,Constant(bulk_conc),subdomains,outer_boundary_marker)
    ##bc1 = DirichletBC(V, Constant(1), subdomains,active_site_marker)
    #PrintBoundary(mesh,bc1,file="test.pvd")


    #return mesh,coordinates
    return mesh
    

if __name__ == "__main__":
    ## read data 
    filename = r"xxxx/x.m"
    import sys
    if len(sys.argv) < 2:
        raise RuntimeError("expected an mcsf file as second argument (or optional rescale value as third arg)")
    filename = sys.argv[1]

    if(sys.argv==3):
      mesh = read_and_mark(filename)
    else:
      mesh = read_and_mark(filename,rescaleCoor=sys.argv[2])
    
#mcsffile = "p.pqr.output.all.m"
#coordinates,cells, markers = read_mcsf(mcsffile)

