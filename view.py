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

