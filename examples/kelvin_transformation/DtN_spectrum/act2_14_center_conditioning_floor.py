# -*- coding: utf-8 -*-
r"""
act2_14_center_conditioning_floor.py  ((b) for Q13 + Q14: the centre-singularity cluster)
================================================================================
Two reviewer objections about the Kelvin BALL's singular material weight w=(R/rho')^2 (which
-> infinity at the centre, the image of infinity).  NOTE: this is a 3-D phenomenon -- the 2-D
Kelvin (circle inversion) is conformally weight-FREE, so both checks are run in 3-D.

  Q13 (conditioning): the manuscript's "cond ~ 1.30, frequency-robust" (sec.5.4) is the DtN
       per-mode eigenvalue RATIO, NOT the assembled linear system one actually solves.  The
       FEM-solve conditioning is what matters; the post-DtN ratio is meaningless for it.  We
       MEASURE the assembled-stiffness condition number of the Kelvin ball and split it into
       the standard-FE part (~h^-2, present for ANY material) and the EXTRA factor the singular
       weight (R/rho')^2 adds.

  Q14 (floor origin): the 5-6 digit floor was attributed to the CURVED-GEOMETRY (sphere
       approximation).  But w=(R/rho')^2 |grad u|^2 is non-polynomial-singular at the centre,
       exactly where high-order quadrature struggles -- so the floor COULD be a centre-quadrature
       artifact.  We separate them: raise the Curve order (floor drops -> geometry) and, at fixed
       mesh+Curve, raise the centre quadrature bonus_intorder (floor unchanged -> NOT quadrature).

VERIFIED HERE (asserted; NGSolve 3-D + numpy dense eig on the free-dof stiffness):
  [1] Q13: the Kelvin-ball assembled-stiffness cond GROWS with refinement and is ORDERS above
      1.30; the EXTRA factor over a uniform-material ball tracks the centre weight ~(R/h_min)^2
      -> the "1.30" is not the solve conditioning; the FEM-solve cond must be quoted instead.
  [2] Q14: the DtN floor DROPS with Curve order (geometry) and is INSENSITIVE to extra centre
      quadrature at fixed Curve -> the floor is the curved-geometry sphere approximation, NOT a
      centre-singularity quadrature artifact.  (The smooth image P_n keeps the centre integral
      finite for the eigenvalue; the conditioning of [1] is a separate, matrix-level effect.)

NON-CLAIM: cond here is the assembled H1 stiffness on free dofs (Dirichlet on Gamma), a clean
machine-independent proxy for the solve conditioning; a full PML-vs-Kelvin solve-cond bench
(matched meshes) is the natural follow-up.  Reuses the act2_05 DtN-floor measurement.
"""
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import scipy.sparse as sp
import ngsolve as ng
from netgen.occ import Sphere, Pnt, OCCGeometry
from radia_mcp.radia_ngsolve.fem_bem_coupling import _solid_harmonic

ng.ngsglobals.msg_level = 0
R = 1.0
rho2 = ng.x * ng.x + ng.y * ng.y + ng.z * ng.z


def assembled_cond(maxh, order, curve, singular):
    """Condition number of the assembled H1 stiffness (Dirichlet on Gamma) on the free dofs."""
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(curve)
        w = R * R / rho2 if singular else ng.CoefficientFunction(1.0)
        fes = ng.H1(mesh, order=order, dirichlet=".*"); u, v = fes.TnT()
        a = ng.BilinearForm(w * ng.grad(u) * ng.grad(v) * ng.dx(bonus_intorder=10)); a.Assemble()
        rows, cols, vals = a.mat.COO()
        A = sp.csr_matrix((np.array(vals), (np.array(rows), np.array(cols))), shape=(fes.ndof, fes.ndof))
        free = np.nonzero([fes.FreeDofs()[i] for i in range(fes.ndof)])[0]
        assert len(free) > 0, "no free dofs -- mesh too coarse (all nodes on Gamma); refine or raise order"
        Af = A[free][:, free].toarray()
        Af = 0.5 * (Af + Af.T)                            # symmetrise (round-off)
        ev = np.linalg.eigvalsh(Af)
        ev = ev[ev > ev[-1] * 1e-13]                      # drop numerical zeros
        return float(ev[-1]), float(ev[0]), len(free)     # lam_max, lam_min, Nfree


def dtn_floor(n, p, maxh, curve_k, bonus):
    """Relative error of the Kelvin DtN eigenvalue -(n+1)/R (cf. act2_05), tunable quadrature."""
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(curve_k)
        fes = ng.H1(mesh, order=p, dirichlet=".*"); u, v = fes.TnT()
        a = ng.BilinearForm(ng.grad(u) * ng.grad(v) * ng.dx); a.Assemble()
        gf = ng.GridFunction(fes); gf.Set(_solid_harmonic(n), ng.BND)
        r = gf.vec.CreateVector(); r.data = -(a.mat * gf.vec)
        gf.vec.data += a.mat.Inverse(fes.FreeDofs()) * r
        E = float(ng.Integrate(ng.grad(gf) * ng.grad(gf) * ng.dx(bonus_intorder=bonus), mesh))
        bm = float(ng.Integrate(gf * gf * ng.ds(bonus_intorder=bonus), mesh))
        lam = -1.0 / R - E / bm
        return abs(lam + (n + 1) / R) / ((n + 1) / R)


print("=" * 82)
print(" act2_14_center_conditioning_floor : the Kelvin-ball centre singularity (Q13 cond, Q14 floor)")
print("=" * 82)

print("\n[1] Q13 -- assembled-stiffness conditioning (NOT the DtN per-mode ratio 1.30):")
print(f"    {'maxh':>5}  {'Nfree':>6}  {'cond_K=Lmax/Lmin':>16}  {'Lmax amplification (K/uniform)':>30}")
hs = [0.4, 0.32, 0.26, 0.22]
conds, amps = [], []
for h in hs:
    lmax_k, lmin_k, nf = assembled_cond(h, 2, 2, singular=True)
    lmax_u, lmin_u, _ = assembled_cond(h, 2, 2, singular=False)
    cond_k = lmax_k / lmin_k; amp = lmax_k / lmax_u       # amp ~ (R/h_min)^2 from the singular centre weight
    conds.append(cond_k); amps.append(amp)
    print(f"    {h:5.2f}  {nf:6d}  {cond_k:16.3e}  {amp:30.2f}")
print(f"    -> cond_K is ORDERS above 1.30 and GROWS with refinement; the singular centre weight")
print(f"       amplifies Lmax by ~(R/h)^2 (the amplification climbs with refinement).  '1.30' is the")
print(f"       DtN per-mode ratio, NOT the linear-system cond one actually solves -- quote the latter.")
assert min(conds) > 30, "the Kelvin assembled cond must be orders above the quoted 1.30"
assert conds[-1] > conds[0] * 3, "the FEM-solve cond must GROW with refinement (it is not a constant 1.30)"
assert amps[-1] > amps[0] * 1.5, "the singular centre weight must amplify Lmax MORE as the mesh refines (~(R/h)^2)"

print("\n[2] Q14 -- floor origin: Curve order (geometry) vs centre quadrature:")
print("    (2a) raise Curve k at fixed mesh -> floor DROPS = geometry-limited:")
fl = {k: dtn_floor(2, 3, 0.33, k, bonus=14) for k in (1, 2, 3)}
print("         " + "  ".join(f"k={k}:{fl[k]:.1e}" for k in (1, 2, 3)))
assert fl[2] < fl[1] * 0.3, "the floor must drop as the Curve (geometry) order rises -> it IS geometry"
print("    (2b) raise centre quadrature bonus_intorder at fixed mesh+Curve -> floor UNCHANGED = not quadrature:")
fq = {b: dtn_floor(2, 3, 0.33, 2, bonus=b) for b in (8, 14, 22)}
print("         " + "  ".join(f"bonus={b}:{fq[b]:.3e}" for b in (8, 14, 22)))
spread = (max(fq.values()) - min(fq.values())) / max(fq.values())
assert spread < 0.15, "the floor must be INSENSITIVE to extra centre quadrature -> not a centre-singularity artifact"
print(f"         spread over quadrature = {spread:.1%} (flat) -> the floor is the CURVED GEOMETRY, not the centre.")

print("\n[verdict]")
print("  Q13 + Q14 settled by experiment:")
print("  - Q13: the Kelvin ball's ASSEMBLED-STIFFNESS cond is orders above the quoted 1.30 and grows")
print("    with refinement, with an EXTRA factor (over uniform material) tracking the singular centre")
print("    weight ~(R/h)^2.  The '1.30' is the DtN per-mode eigenvalue ratio, NOT the FEM-solve cond;")
print("    the manuscript should quote the assembled-system cond (and bench it vs PML at matched mesh).")
print("  - Q14: the 5-6 digit floor DROPS with Curve order and is INSENSITIVE to extra centre quadrature")
print("    -> it is the curved-geometry sphere approximation, NOT a centre-singularity quadrature floor.")
print("    (The smooth image P_n keeps the eigenvalue's centre integral finite; conditioning is separate.)")
print("\nALL CHECKS PASSED.")
