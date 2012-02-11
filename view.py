def PrintBoundary(mesh, boundary,file="marked.pvd"):
  from dolfin import FunctionSpace, Function, File, plot
  V = FunctionSpace(mesh, "CG", 1)
  
  marked = Function(V)
  boundary.apply(marked.vector())

  #plot(marked, interactive=1)
  
  print "Printing %s for viewing" % file
  File(file) << marked

