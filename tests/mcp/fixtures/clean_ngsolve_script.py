"""Clean NGSolve script -- should produce zero lint findings.

Used by --selftest to verify no false positives.
"""
from ngsolve import *
from netgen.occ import *

# Correct: OCCGeometry with dim=2 for 2D
shape = Rectangle(1, 1).Face()
geo = OCCGeometry(shape, dim=2)

# Correct: Preconditioner before Assemble
mesh = Mesh(geo.GenerateMesh(maxh=0.1))
fes = H1(mesh, order=2)
u, v = fes.TnT()
a = BilinearForm(fes)
a += grad(u) * grad(v) * dx
pre = Preconditioner(a, "bddc")
a.Assemble()

# Correct: using xi instead of x for loop variable
for xi in range(10):
    pass

# Correct: using .vec.data for assignment
# gfu.vec.data = f.vec

# Plain numerical computation (no domain-specific keywords)
import numpy as np
x_vals = np.linspace(0, 1, 100)
