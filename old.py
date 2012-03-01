

# dont think this is needed before, since I should usually have a mesh from APBS 

# boundary is given as input, here we compute potential (for sphere)
# given input boundary. Note: normally youd want to use APBS to get
# the electrostatic potential
def MeshNoPotential(debug=0):

  ## load data
  # data from ~/localTemp/NBCR/smol/born_ex
  fileMesh = "example/sphere/p.pqr.output.all_mesh.xml.gz"
  mesh = Mesh(fileMesh)
  

  # Function space
  V = FunctionSpace(mesh, "CG", 1)


  ## load markers
  useMarkers=1
  # Need to pull these from Gamer output at some point 

  ## define Dirichlet boundary
  bcs=1
  if(useMarkers==1):
    bcs=1
  else:
    # PKH: (Found no facets matching domain for boundary condition.) <-- probably from Dirichlet, since Neumann BC expressed in weak form of PDE
    bc_active = DirichletBC(V, Constant(parms.active_site_absorb), Sphere.DirichletActiveSite())
    bc_bulk = DirichletBC(V, Constant(parms.bulk_conc), Sphere.DirichletBulkBoundary())
    
    bcs = [bc_active, bc_bulk]

  if(debug==0):
      return mesh

  ## Define Neumann boundary
  if(useMarkers==1):
    1
  else:      
    # PKH: Verify that Neumann BC is implemented correctly
    subdomains = MeshFunction("uint",mesh,2)
    subdomain = Sphere.NeumannMolecularBoundary()
    subdomain.mark(subdomains,molecular_boundary_marker)

    # PKH: do I need to mark BulkBoundary as well in order to use assemble call later?
    # NOTE: is this for sure where I should be computing kon? Seems like the flux along this boundary s.b. zero!
    subdomainOuter = Sphere.DirichletBulkBoundary()
    subdomainOuter.mark(subdomains,outer_boundary_marker)

  haveAPBS=0
  if(haveAPBS):
    #psi = GetPotential()
    1
  ## compute potential
  else:
    # PKH Is this the best way of computing Sphere potential over mesh coordinates? 
    x = mesh.coordinates() 
    psi = interpolate(Sphere.Potential(x),V)
    # PKH: I need to check potential, since code x-plodes 
    #psi.vector()[:]=0;

  problem.subdomains = subdomains
  problem.mesh=mesh
  problem.psi = psi
  problem.bcs = bcs
  problem.V   = V
  return problem
