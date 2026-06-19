# -*- coding: utf-8 -*-
r"""
demo_pml_dof_cost.py  (Track A -- kelvin branch)
================================================
The COMMITTED, asserted backing for the manuscript claim that, among the open-BC
closures, **PML reaches a given accuracy only by spending more exterior DoF**,
whereas the **Kelvin closure buys accuracy through the isoparametric geometry order
at essentially CONSTANT DoF**.  This was previously only a structural argument in
the text; here it is measured on two independent axes.

The honest framing this pins (no overclaim -- every 'ok' gated on a number):
  Kelvin's distinguishing property is CONTROLLABILITY, not unconditional accuracy.
  Its open-BC error is set by two knobs you choose -- polynomial order p (>= the
  multipole degree n present; see demo_n_kelvin_approximates_bem.py) and the
  isoparametric Curve order (the residual floor IS the curved-sphere geometry
  error).  Crucially the Curve lever costs NO extra H1 DoF: .Curve(c) changes the
  geometry map, not the FE space dimension.  PML's ONLY lever is more absorber DoF.

PROVEN HERE:
  [P] PML radial cost (1-D radial DtN, complex stretch, low frequency kR=0.2 -- the
      regime where PML still applies as an open BC; for the truly static ladder a
      wave-absorbing PML degenerates, which is itself a point for Kelvin/BEM).
      To gain two digits (rel DtN err 1e-3 -> 1e-5) the radial DoF must grow >=5x
      (n=1: 31 -> 399); and at a FIXED absorber thickness a tight target can be
      unreachable altogether (n=9 floors at ~2.4e-5 > 1e-5 even at 400 elements).
  [K] Kelvin geometry lever at CONSTANT DoF: at identical (maxh, order) -> identical
      ndof, raising Curve 3 -> 5 improves the DtN error >=10x.  So Kelvin reaches a
      tighter target WITHOUT adding DoF -- the opposite of PML.
  [G] 3-D geometry DoF: the PML absorber shell (1<=r<=2) carries ~6x the Kelvin
      ball's (r<=1) DoF at matched maxh -- its nonzero thickness alone is a DoF tax,
      before any radial-resolution argument.

PML radial part = pure numpy/scipy (fast).  [K] and [G] call the committed library
+ NGSolve (a few small meshes): under a minute total.  Run with the worktree on
PYTHONPATH so the library is in effect.
"""
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
from scipy.special import jv, yv, hankel1, hankel2
from radia_mcp.radia_ngsolve.fem_bem_coupling import kelvin_dtn_eigenvalue

np.seterr(all="ignore")
N_FAIL = 0
def check(name, cond, detail=""):
    global N_FAIL
    if not cond: N_FAIL += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}{('  -- ' + detail) if detail else ''}")

print("=" * 78)
print(" demo_pml_dof_cost: PML buys accuracy with DoF; Kelvin with geometry order")
print("=" * 78)

# ===========================================================================
# [P] PML radial DoF to reach a target DtN accuracy (low frequency kR=0.2)
#     1-D radial complex-stretch DtN, lifted from demo_lf9_spectrum_n1to9.py [C].
# ===========================================================================
_GX = np.array([-np.sqrt(3/5), 0.0, np.sqrt(3/5)]); _GW = np.array([5/9, 8/9, 5/9])
def _sph(kind, n, z):
    z = np.asarray(z, dtype=complex); pref = np.sqrt(np.pi/(2.0*z))
    return pref * {'j': jv, 'y': yv, 'h1': hankel1, 'h2': hankel2}[kind](n+0.5, z)
def _sphp(kind, n, z):
    z = np.asarray(z, dtype=complex); return _sph(kind, n-1, z) - (n+1.0)/z*_sph(kind, n, z)
def dtn_exact_helm(n, z):
    z = complex(z); return (z*_sphp('h1', n, z)/_sph('h1', n, z)).item()
def _A_pml(n, k, a, d, M, s0):
    def s(r): return 1 + 1j*s0*(r-a)**2/(k*d*d)
    def rt(r): return r + 1j*s0*(r-a)**3/(3*k*d*d)
    nod = np.linspace(a, a+d, M+1); Am = np.zeros((M+1, M+1), complex)
    for e in range(M):
        r0, r1 = nod[e], nod[e+1]; h = r1-r0; dp = np.array([-1/h, 1/h])
        K = np.zeros((2, 2), complex); Cc = np.zeros((2, 2), complex); Mm = np.zeros((2, 2), complex)
        for gp, gw in zip(_GX, _GW):
            r = .5*(r0+r1)+.5*h*gp; wq = .5*h*gw
            ph = np.array([(r1-r)/h, (r-r0)/h]); sr = s(r); rr = rt(r)
            K += np.outer(dp, dp)*(rr*rr/sr)*wq
            Cc += np.outer(ph, ph)*sr*wq; Mm += np.outer(ph, ph)*sr*rr*rr*wq
        Am[e:e+2, e:e+2] += K + n*(n+1)*Cc - k*k*Mm
    return Am
def dtn_pml(n, k, a, d, M, s0):
    Am = _A_pml(n, k, a, d, M, s0); u = np.zeros(M+1, complex); u[0] = 1.0; idx = list(range(1, M))
    u[idx] = np.linalg.solve(Am[np.ix_(idx, idx)], -Am[np.ix_(idx, [0])][:, 0])
    return -(Am[0, :] @ u)/a

print("\n[P] PML radial DoF to reach a target DtN accuracy (kR=0.2, d=1, s0=15):")
k, a, d, s0 = 0.2, 1.0, 1.0, 15.0
Ms = [2, 4, 8, 16, 32, 64, 128, 256, 400]
def min_elems(errs, tol):
    for M in Ms:
        if errs[M] <= tol:
            return M
    return None
pml_err = {}
for n in (1, 4, 9):
    ex = dtn_exact_helm(n, k)
    errs = {M: abs(dtn_pml(n, k, a, d, M, s0) - ex)/abs(ex) for M in Ms}
    pml_err[n] = errs
    m3, m5 = min_elems(errs, 1e-3), min_elems(errs, 1e-5)
    print(f"    n={n}: " + "  ".join(f"M{M:>3}:{errs[M]:.1e}" for M in Ms))
    print(f"          min elements  1e-3: {m3}   1e-5: {m5}")
# monotone refinement (error decreases as the absorber is resolved)
mono = all(pml_err[1][Ms[i+1]] < pml_err[1][Ms[i]] for i in range(len(Ms)-1))
m3_1, m5_1 = min_elems(pml_err[1], 1e-3), min_elems(pml_err[1], 1e-5)
check("PML error decreases monotonically as the absorber is radially refined (n=1)",
      bool(mono), f"{pml_err[1][Ms[0]]:.1e} -> {pml_err[1][Ms[-1]]:.1e}")
check("PML must grow DoF >=5x to gain 2 digits (n=1: 1e-3 -> 1e-5)",
      m3_1 is not None and m5_1 is not None and m5_1 / m3_1 >= 5,
      f"elements 1e-3:{m3_1} -> 1e-5:{m5_1}  ({(m5_1/m3_1):.0f}x)")
check("PML at fixed thickness cannot reach 1e-5 for n=9 (saturates ~2.4e-5 at 400 elems)",
      min_elems(pml_err[9], 1e-5) is None and pml_err[9][400] > 1e-5,
      f"floor(n=9, M=400) = {pml_err[9][400]:.1e}")

# ===========================================================================
# [K] Kelvin geometry lever at CONSTANT DoF: Curve 3 -> 5 improves error, ndof fixed
# ===========================================================================
print("\n[K] Kelvin: Curve order improves accuracy at IDENTICAL DoF (library):")
c3 = kelvin_dtn_eigenvalue(R=1.0, degree=1, maxh=0.45, order=9, intorder=16, curve=3)
c5 = kelvin_dtn_eigenvalue(R=1.0, degree=1, maxh=0.45, order=9, intorder=16, curve=5)
print(f"    Curve 3: rel_err={c3['rel_err']:.2e}  ndof={c3['ndof']}")
print(f"    Curve 5: rel_err={c5['rel_err']:.2e}  ndof={c5['ndof']}  (same FE space)")
check("Curve 3->5 leaves the H1 DoF count UNCHANGED (geometry map, not FE space)",
      c3["ndof"] == c5["ndof"], f"ndof {c3['ndof']} == {c5['ndof']}")
check("Curve 3->5 improves DtN error >=10x at that constant DoF (the no-DoF lever)",
      c5["rel_err"] * 10 <= c3["rel_err"], f"{c3['rel_err']:.1e} -> {c5['rel_err']:.1e}")

# ===========================================================================
# [G] 3-D geometry: PML absorber shell carries ~6x the Kelvin ball's DoF
# ===========================================================================
print("\n[G] 3-D exterior DoF at matched maxh (order 2, Curve 3): ball r<=1 vs shell 1..2:")
import ngsolve as ng
from netgen.occ import Sphere, Pnt, OCCGeometry
def _ndof(geo, mh):
    m = ng.Mesh(OCCGeometry(geo).GenerateMesh(maxh=mh)).Curve(3)
    return ng.H1(m, order=2, dirichlet=".*").ndof
ratios = []
for mh in (0.5, 0.4, 0.3):
    nk = _ndof(Sphere(Pnt(0, 0, 0), 1.0), mh)
    ns = _ndof(Sphere(Pnt(0, 0, 0), 2.0) - Sphere(Pnt(0, 0, 0), 1.0), mh)
    ratios.append(ns / nk)
    print(f"    maxh={mh}:  ball ndof={nk:5d}   shell ndof={ns:6d}   shell/ball = {ns/nk:.1f}x")
check("PML shell carries ~6x the Kelvin ball DoF at matched maxh (ratio in [4,7])",
      all(4.0 <= r <= 7.0 for r in ratios), f"ratios = {[round(r,1) for r in ratios]}")

print("\n" + "=" * 78)
if N_FAIL == 0:
    print(" ALL CHECKS PASSED -- claim (3) is measured: PML's accuracy lever is DoF")
    print(" (grows with the target, can saturate); Kelvin's is the geometry order at")
    print(" constant DoF.  Hence Kelvin reaches tighter targets at lower exterior DoF.")
else:
    print(f" {N_FAIL} CHECK(S) FAILED")
print("=" * 78)
assert N_FAIL == 0, f"{N_FAIL} checks failed"
