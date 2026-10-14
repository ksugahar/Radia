"""Intentionally bad NGSolve script for lint testing.

This script contains violations for NGSolve-specific rules.
Used by --selftest when examples/ directory is not available.

Expected findings: 18+ (one per triggerable NGSolve/BEM/HDiv rule)
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
gfu.vec = f.vec  # should be gfu.vec.data = f.vec

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
a += mu * curl(u) * curl(v) * dx("Kelvin")

# Rule: ngsolve-missing-trace-bem (CRITICAL)
fes_bem = HDivSurface(mesh, order=1)
u_b, v_b = fes_bem.TnT()
a_bem = BilinearForm(LaplaceSL(3) * u_b * v_b * ds)

# Rule: ngsolve-cg-on-saddle-point (MODERATE)
fes_saddle = fesA * fesOmega
gfu_saddle = solvers.CG(a.mat, rhs.vec)

# Rule: ngsolve-vectorh1-for-em (MODERATE)
# magnetostatic simulation
fes_vh1 = VectorH1(mesh, order=2)

# Rule: ngsolve-pinvit-no-projection (MODERATE)
fes_pinvit = HCurl(mesh, order=2)
evals, evecs = solvers.PINVIT(a.mat, m.mat, pre, num=10)

# Rule: joule-heat-missing-conj (MODERATE)
Q = sigma * InnerProduct(E, E)

# Rule: axisymmetric-naive-1-over-r (MODERATE)
# axisymmetric problem
a += nu / r * grad(u) * grad(v) * dx

# Rule: ngsbem-volume-mesh (HIGH)
from ngsolve.comp import LaplaceSL
box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
geo_bem = OCCGeometry(box)

# Rule: scattered-sibc-missing-nxH (HIGH)
# scattered field SIBC formulation
a += sibc * A_inc * v * ds("skin")

# Rule: efie-sibc-use-scalar-bie (HIGH)
fes_efie = HDivSurface(mesh, order=1)
a_efie = BilinearForm(LaplaceSL(3) * u * v * ds + Z_s * u * v * ds)

# Rule: ngsbem-missing-curvaturesafety (MODERATE)
from ngsolve.comp import LaplaceSL
geo_curved = OCCGeometry(sphere)
mesh_curved = geo_curved.GenerateMesh(maxh=0.1)

# Rule: ngsbem-taskmanager-nondeterministic (LOW)
from ngsolve.comp import LaplaceSL
with TaskManager():
    a_task = BilinearForm(LaplaceSL(3) * u * v * ds)

# Rule: hdiv-wrong-sign (HIGH)
N_matrix = -np.diag(inv_chi) - N

# Rule: hdiv-eval-face-center (MODERATE)
eval_point = face_center
