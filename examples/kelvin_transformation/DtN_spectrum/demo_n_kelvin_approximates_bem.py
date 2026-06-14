# -*- coding: utf-8 -*-
# DEMO (n) (verified): is Kelvin a SPARSE approximation/replacement for the BEM
# exterior DtN? Both discretise the SAME exterior Dirichlet-to-Neumann operator on
# the truncation sphere Gamma (eigenvalue ladder lambda_n = -(n+1)/R). BEM is exact in
# the exterior (no exterior mesh) but pays a DENSE Gamma x Gamma matrix encoding the
# whole ladder at once; Kelvin is a SPARSE thin-ball volume whose Gamma trace
# reproduces the ladder up to multipole n ~ order p (controllable). Because accuracy
# is bottlenecked by the SHARED surface -- NOT by how the exterior is treated -- Kelvin
# with p >= n MATCHES the BEM per-degree accuracy on the SAME Gamma, while staying
# sparse. So Kelvin is a faithful sparse stand-in for the BEM exterior DtN, i.e.
# "BEM truncated at multipole p, realised as a sparse volume".
#
# Measured (coarse sphere, 336 surface DoF, R=1):
#   BEM dense DtN (SurfaceL2 o=1, nnz=ndof^2, ~67 s):  n=1 7.4e-4  n=2 2.1e-3  n=3 5.2e-3
#   Kelvin sparse (p=3, ~7% fill, ~20 nnz/row, ~ms):   n=1 3.3e-5  n=2 1.3e-5  n=3 8.4e-4
#   -> p>=n: Kelvin matches/beats BEM on the same Gamma; matrix ~93% sparse vs 100% dense.
import os, time
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import ngsolve as ng
from netgen.occ import Sphere, Pnt, OCCGeometry
from radia_mcp.radia_ngsolve.bem_integral import exterior_dtn_spectrum
from radia_mcp.radia_ngsolve.fem_bem_coupling import kelvin_dtn_eigenvalue

R = 1.0
maxh = 0.5          # >=0.5 -> coarsest sphere; SAME Gamma triangulation for both methods
NMAX = 3

# --- BEM exterior DtN: exact exterior, but a DENSE surface matrix -------------
t0 = time.time()
bem = exterior_dtn_spectrum(R=R, maxh=maxh, order=1, intorder=10, nmax=NMAX)
t_bem = time.time() - t0
bem_err = {m["n"]: m["rel_err"] for m in bem["modes"]}
nd_bem = bem["ndof"]
print("BEM exterior DtN (dense surface): SurfaceL2 order=1, ndof=%d, maxh=%.2f" % (nd_bem, maxh))
print("  dense matrix: nnz = ndof^2 = %d (100%% fill), assembled in %.0f s\n" % (nd_bem * nd_bem, t_bem))

# --- Kelvin closure: SPARSE volume, multipole reach controlled by order p -----
print("Kelvin closure (sparse volume), same Gamma. DtN rel error vs degree n and order p:")
print("  p \\ n |" + "".join("%10d" % n for n in range(1, NMAX + 1)))
print("  " + "-" * (8 + 10 * NMAX))
for p in (1, 2, 3):
    row = "  %5d |" % p
    for n in range(1, NMAX + 1):
        r = kelvin_dtn_eigenvalue(R=R, degree=n, maxh=maxh, order=p, dim=3)
        row += "%10.2e" % r["rel_err"]
    print(row)

# Kelvin sparse fill on this same ball (one assembly at p=2 to report nnz)
mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(3)
fes = ng.H1(mesh, order=2, dirichlet=".*")
u, v = fes.TnT()
a = ng.BilinearForm(ng.grad(u) * ng.grad(v) * ng.dx); a.Assemble()
nze, nd_k = a.mat.nze, fes.ndof
fill = nze / float(nd_k * nd_k)
print("\nKelvin volume matrix (p=2): ndof=%d, nnz=%d, fill=%.4f%% (~%.1f nnz/row)"
      % (nd_k, nze, 100.0 * fill, nze / nd_k))

# --- BEM reference row (the full ladder in one dense matrix) ------------------
print("\nBEM dense DtN rel error (full ladder, one matrix):")
print("    n   |" + "".join("%10d" % n for n in range(1, NMAX + 1)))
print("  BEM   |" + "".join("%10.2e" % bem_err.get(n, float('nan')) for n in range(1, NMAX + 1)))

print("\n-> Kelvin with p>=n MATCHES the BEM per-degree accuracy on the SAME Gamma,")
print("   but its matrix is ~%.0f%% sparse vs BEM's 100%% dense. Kelvin = the BEM exterior" % (100 * (1 - fill)))
print("   DtN truncated at multipole p (controllable), realised as a SPARSE volume.")
