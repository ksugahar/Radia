# -*- coding: utf-8 -*-
# DEMO (o) (verified): is Kelvin BETTER than H-matrix/FMM BEM? Both turn the O(N^2)
# dense exterior coupling into a sub-quadratic one -- but by OPPOSITE mechanisms, and
# Kelvin's is simpler on every axis WITHIN its applicability domain.
#
#   H-matrix/FMM: the dense BEM matrix has admissible (well-separated) off-diagonal
#     blocks that are NUMERICALLY LOW RANK (the Green's function is smooth across the
#     gap) -> store each as r<<m,n columns -> O(N log N). It must BUILD this: cluster
#     tree + ACA/analytic low-rank + singular near-field quadrature + an inner V^-1
#     solve + an indefinite FEM-BEM coupling.
#   Kelvin: the SAME geometric far-DoF pair is STRUCTURALLY ZERO in the sparse volume
#     Laplacian -- nothing to compress. Sparse SPD, local element integrals, ~ms
#     assembly, composes with the interior FEM as ONE sparse solve.
#
# Measured (sphere R=1): the dense single-layer block between two well-separated polar
# caps (50x56, gap 0.91) has rel singular values 1, 9.5e-2, 8.1e-2, ... -> numerical
# rank 5 @1e-2 (90% saved), 14 @1e-4 (72% saved) = H-matrix-compressible. The Kelvin
# volume matrix on the same ball is ~21 nnz/row and the analogous far block holds
# EXACTLY 0 nonzeros (structural sparsity), reached with no compression machinery.
#
# VERDICT: for a sphere-able truncation where the field is needed INSIDE, Kelvin is the
# simpler/cheaper data-sparse route (no cluster tree, no ACA, no singular quadrature, no
# inner solve, SPD). H-matrix's edge is GENERALITY: arbitrary closed Gamma and the
# exterior field at any point -- which Kelvin (a spherical inversion) does not give.
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import ngsolve as ng
from netgen.occ import Sphere, Pnt, OCCGeometry
import ngsolve.bem as bem
from radia_mcp.radia_ngsolve.bem_integral import _dense_mat

R, maxh = 1.0, 0.4
mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(3)

# --- (1) dense single-layer V: the object H-matrix compresses -----------------
fes = ng.SurfaceL2(mesh, order=0, dual_mapping=True)
nd = fes.ndof
V = bem.SingleLayerPotentialOperator(fes, intorder=12)
Vd = np.real(_dense_mat(V.mat, nd))
print("dense single-layer V on sphere: ndof=%d, fill=100%% (every DoF couples to every DoF)" % nd)

xyz = np.array(bem.GetDofCoordinates(fes))
if xyz.shape[0] != nd:
    xyz = xyz.T
z = xyz[:, 2]
c = 0.45
A = np.where(z > c)[0]
B = np.where(z < -c)[0]
print("two well-separated caps: |A|=%d (z>%.2f), |B|=%d (z<-%.2f), vertical gap=%.2f"
      % (len(A), c, len(B), c, z[A].min() - z[B].max()))

blk = Vd[np.ix_(A, B)]
s = np.linalg.svd(blk, compute_uv=False); s0 = s[0]
print("\noff-diagonal block V[A,B]: %dx%d DENSE in plain BEM; H-matrix low-rank-compresses it" % blk.shape)
print("  rel singular values: " + "  ".join("%.1e" % (si / s0) for si in s[:8]) + " ...")
for eps in (1e-2, 1e-4, 1e-6, 1e-8):
    r = int(np.sum(s > eps * s0)); full = min(blk.shape)
    print("  eps=%.0e -> rank r=%d of %d  (%.0f%% saved)" % (eps, r, full, 100 * (1 - r / float(full))))

# --- (2) Kelvin volume: the same far block is STRUCTURALLY ZERO ---------------
fesK = ng.H1(mesh, order=1, dirichlet=".*")   # order 1 -> DoF i = vertex i (easy coords)
u, v = fesK.TnT()
a = ng.BilinearForm(ng.grad(u) * ng.grad(v) * ng.dx); a.Assemble()
ndK, nzeK = fesK.ndof, a.mat.nze
vz = np.array([p.point[2] for p in mesh.vertices])
AK = set(np.where(vz > c)[0]); BK = set(np.where(vz < -c)[0])
ri, ci, _ = a.mat.COO()
ri = np.array(ri); ci = np.array(ci)
far_nnz = int(np.sum([(i in AK and j in BK) or (i in BK and j in AK) for i, j in zip(ri, ci)]))
print("\nKelvin volume Laplacian (same ball, p=1): ndof=%d, ~%.0f nnz/row (fill->0 as N grows)"
      % (ndK, nzeK / ndK))
print("  far block Kelvin[A,B] (north-cap vs south-cap DoFs): %d nonzeros  <-- structural ZERO"
      % far_nnz)
print("  (vs V[A,B] = a dense %dx%d block; Kelvin has NOTHING for H-matrix to compress)" % blk.shape)

print("\nVERDICT: same sub-quadratic endpoint, opposite mechanisms.")
print("  H-matrix: COMPRESS the dense block to rank r (cluster tree + ACA + singular near-field")
print("            quadrature + inner V^-1 + indefinite FEM-BEM coupling).")
print("  Kelvin:   the block is ZERO by construction -> sparse SPD, local integrals, ~ms assembly,")
print("            one sparse solve with the interior FEM. Simpler/cheaper for sphere-able truncation")
print("            + interior field; H-matrix's edge is generality (arbitrary Gamma, exterior field).")
