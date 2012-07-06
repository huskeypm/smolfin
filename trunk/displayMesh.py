from dolfin import *
mesh = Mesh("potential-0_mesh.xml.gz")
V = FunctionSpace(mesh, "CG", 1)
u = Function(V, "potential-0_values.xml.gz")
plot(u, interactive=1)

