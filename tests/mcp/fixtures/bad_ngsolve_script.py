"""Intentionally bad NGSolve script for lint testing.

This script contains violations for NGSolve-specific rules.
Used by --selftest when examples/ directory is not available.

Expected findings: 7+ (one per triggerable NGSolve rule)
"""
from ngsolve import *
from netgen.occ import *

# Rule: ngsolve-overwrite-xyz (MODERATE)
# 'x' is a built-in NGSolve coordinate variable
for x in range(10):
    pass

# Rule: ngsolve-dim2-occ (MODERATE)
shape = Rectangle(1, 1).Face()
geo = OCCGeometry(shape)

# Rule: ngsolve-vec-assign (MODERATE)
# gfu.vec = f.vec  # should be gfu.vec.data = f.vec

# Rule: ngsolve-precond-after-assemble (MODERATE)
mesh = Mesh(geo.GenerateMesh(maxh=0.1))
fes = H1(mesh, order=2)
u, v = fes.TnT()
a = BilinearForm(fes)
a += grad(u) * grad(v) * dx
a.Assemble()
pre = Preconditioner(a, "bddc")

# Rule: hcurl-missing-nograds (HIGH)
# magnetostatic solver
fes_mag = HCurl(mesh, order=2)

# Rule: eddy-current-missing-complex (HIGH)
# eddy current simulation
fes_eddy = HCurl(mesh, order=1)

# Rule: ngsolve-kelvin-missing-bonus-intorder (MODERATE)
# a += mu * curl(u) * curl(v) * dx("Kelvin")
