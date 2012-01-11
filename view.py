def PrintBoundary(mesh,boundary):
  from dolfin import FunctionSpace,Function,File
  V = FunctionSpace(mesh, "CG", 1)
  
  marked = Function(V)
  boundary.apply(marked.vector())

  #plot(u, interactive=TRUE)
  #plot(marked, interactive=1)
  #interactive()

  File("marked.pvd") << marked

