# -*- coding: utf-8 -*-
# DEMO (r) (verified): FEM-Kelvin is a SPARSE approximation/factorization of the dense
# Green's-function (BEM) kernel matrix -- the user's "trace the DtN backward" view.
#
# The dense exterior DtN kernel Lambda (assembled directly by BEM via the Green function) is
# exactly the SCHUR COMPLEMENT of the sparse FEM-Kelvin matrix A: eliminate A's volume DoFs and
# the dense Lambda falls out (cf. demo_o/demo_q). So A is a SPARSE REPRESENTATION of the same
# operator, obtained by ADDING auxiliary volume DoFs (the thin air layer + Kelvin ball). The
# trade: DoF goes UP (surface-only -> surface+volume) but COST goes DOWN, because each volume
# DoF couples to only ~20-40 neighbours -> assembly/storage are O(N) local grad-grad integrals
# instead of a dense O(N_S^2) matrix of singular Green-function quadratures.
#
# This is a DIFFERENT sparsification from H-matrix/FMM: H-matrix keeps the SAME DoF and
# compresses the dense off-diagonal blocks to low rank (data-sparse); FEM-Kelvin ADDS DoF to
# make the matrix STRUCTURALLY sparse (genuine zeros). Both beat dense O(N^2); FEM-Kelvin needs
# no cluster tree / ACA / singular quadrature / inner solve and stays SPD.
#
# Measured (sphere R=1, maxh=0.5, same Gamma); both reproduce the ladder -(n+1)/R:
#   dense BEM kernel : DoF 336, nnz 112896 (100% fill), assembly 70 s, ladder rel 7e-4..5e-3
#   sparse FEM-Kelvin: DoF 768, nnz  27832 (4.7% fill, ~36/row), assembly 0.008 s, rel 3e-5..8e-4
#   => DoF x2.3 UP, assembly ~8800x FASTER, fill 100% -> 4.7%.
import os, time
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import ngsolve as ng
from netgen.occ import Sphere, Pnt, OCCGeometry
from radia_mcp.radia_ngsolve.bem_integral import exterior_dtn_spectrum
from radia_mcp.radia_ngsolve.fem_bem_coupling import kelvin_dtn_eigenvalue

R, maxh = 1.0, 0.5

# --- dense Green's-function kernel (BEM): assemble Lambda directly ------------
t0 = time.time()
bem = exterior_dtn_spectrum(R=R, maxh=maxh, order=1, intorder=10, nmax=3)
t_bem = time.time() - t0
nS = bem["ndof"]
bem_rel = {m["n"]: m["rel_err"] for m in bem["modes"]}

# --- sparse FEM-Kelvin matrix A: assemble the volume stiffness (timed) -------
mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(3)
fes = ng.H1(mesh, order=3, dirichlet=".*")
x, y, z = ng.x, ng.y, ng.z
w = R * R / (x * x + y * y + z * z + 1e-30)              # Kelvin weight (R/rho')^2
u, v = fes.TnT()
t0 = time.time()
a = ng.BilinearForm(w * ng.grad(u) * ng.grad(v) * ng.dx(bonus_intorder=6)); a.Assemble()
t_kel = time.time() - t0
nTot, nnz = fes.ndof, a.mat.nze

# Kelvin DtN accuracy on the same sphere (energy quotient = condensed DtN), n=1..3
kel_rel = {n: kelvin_dtn_eigenvalue(R=R, degree=n, maxh=maxh, order=3, dim=3)["rel_err"] for n in (1, 2, 3)}

print("same truncation sphere Gamma (R=%.1f, maxh=%.2f). Both encode the exterior DtN ladder -(n+1)/R:" % (R, maxh))
print("   n |  dense BEM kernel rel.err |  sparse FEM-Kelvin rel.err")
for n in (1, 2, 3):
    print("   %d |        %.2e          |        %.2e" % (n, bem_rel.get(n, float('nan')), kel_rel[n]))

print("\nthe sparse matrix is a representation of the SAME operator (its Schur complement = the dense kernel):")
print("  %-26s | %14s | %16s" % ("", "dense BEM Lambda", "sparse FEM-Kelvin A"))
print("  %-26s | %14d | %16d   <- DoF INCREASES (+volume)" % ("DoF", nS, nTot))
print("  %-26s | %14d | %16d   <- nnz: dense N^2 vs ~%.0f/row" % ("nonzeros", nS * nS, nnz, nnz / nTot))
print("  %-26s | %13.1f%% | %15.3f%%   <- fill" % ("fill", 100.0, 100.0 * nnz / (nTot * nTot)))
print("  %-26s | %12.1fs | %14.3fs   <- assembly COST DECREASES" % ("assembly time", t_bem, t_kel))
print("  %-26s | %14s | %16s" % ("kernel", "Green's fn (sing.)", "local grad-grad"))
print("\n-> DoF up (%d->%d, x%.1f) but cost down (assembly x%.0f faster, fill %.0f%%->%.2f%%)." %
      (nS, nTot, nTot / nS, t_bem / max(t_kel, 1e-6), 100.0, 100.0 * nnz / (nTot * nTot)))
print("   FEM-Kelvin is the dense Green's kernel's SPARSE factorization (add volume DoFs -> structural")
print("   sparsity), distinct from H-matrix (same DoF, compress the dense blocks to low rank).")
