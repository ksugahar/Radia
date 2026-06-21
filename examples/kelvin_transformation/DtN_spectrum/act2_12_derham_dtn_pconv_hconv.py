# -*- coding: utf-8 -*-
"""
act2_12_derham_dtn_pconv_hconv.py  (Act 2 -- the DtN spectrum is FORM-dependent; p- and h-convergence)
=====================================================================================================
act2_11 measured the Kelvin DtN spectrum for the SCALAR (H1) field only.  But the open-boundary DtN /
Steklov operator depends on the de Rham FORM of the unknown -- H1 (0-form), H(curl) (1-form), H(div)
(2-form) do NOT share one spectrum.  The de Rham/Hodge structure gives exactly TWO distinct radial
ladders (act7_30):

    gradient / irrotational   (H1 0-form ; H(div) NORMAL trace)      DtN = -(n+1)/R   (scalar ladder)
    toroidal / transverse     (H(curl) TANGENTIAL trace)             DtN = -n/R       (vector ladder)

so H1 and the H(div) normal trace coincide at -(n+1)/R (Hodge duals), while H(curl) is the distinct
-n/R.  Kelvin closes ALL of them automatically (it is a geometric map, form-agnostic); a decay-basis
infinite element needs a separate construction per form (the manuscript's de Rham point).

This act shows, for each form, BOTH convergence axes the manuscript leans on:
  * p-convergence : raise the FE / radial order at a fixed mesh -- Kelvin is a p-method, so the DtN
                    defect drops fast (radial IE: EXACT once order >= n+1; FEM: ~(h/R)^(2p)).
  * h-convergence : refine the mesh at fixed order -- the algebraic rate ~(h/R)^(2p) (act2_11).
p-refinement reaches a target at far fewer DoF than h-refinement (the manuscript's 20-80x).

Panels:
  [A] the two ladders, ALL THREE forms, RADIAL p-convergence -- the rigorous de Rham radial IE
      (act7_30 energy matrices, pure numpy): scalar -> -(n+1), toroidal -> -n, exact for order>=n+1.
  [B] on a real CURVED MESH (the Kelvin closure, the manuscript's method):
        H1 (0-form)  via the condensed surface DtN (act2_11 machinery): p-conv (FE order) + h-conv;
        H(curl) (1-form) via the Kelvin-ball curl-curl FEM (act3_03): the -n ladder on a mesh, p-conv.
      (H(div) normal-trace FEM == the H1 problem with B=grad phi, so H1's p/h curves cover it too.)

Reuses radia.infinite_element + the act7_30 radial energies + the act3_03 vector Kelvin ball.
Needs NGSolve + numpy + scipy.  Self-asserting; writes JSON.
"""
import os
import json
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp

import ngsolve as ng
from netgen.occ import OCCGeometry, Sphere, Pnt
from radia.infinite_element import dtn_surface_matrix

N_FAIL = 0


def check(name, cond, detail=""):
    global N_FAIL
    if not cond:
        N_FAIL += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}{('  -- ' + detail) if detail else ''}")


# ===========================================================================
# [A] the de Rham radial ladders -- pure-numpy IE energy matrices (act7_30)
# ===========================================================================
def scalar_energy(n, P):            # gradient / H1 / H(div)-normal  ->  ladder -(n+1)
    k = np.arange(1, P + 1)
    return (np.outer(k, k) + n * (n + 1)) / (k[:, None] + k[None, :] - 1.0)


def toroidal_energy(n, P):          # toroidal / H(curl)-tangential  ->  ladder -n
    k = np.arange(1, P + 1)
    return (np.outer(1 - k, 1 - k) + n * (n + 1)) / (k[:, None] + k[None, :] - 1.0)


def ie_dtn(E):
    g = np.ones(E.shape[0])
    return -1.0 / (g @ np.linalg.solve(E, g))


# ===========================================================================
# [B] H1 (0-form) DtN on a curved mesh -- condensed surface operator (act2_11)
# ===========================================================================
def _to_dense(mat, ndof):
    r, c, v = mat.COO()
    A = sp.coo_matrix((np.array(v), (np.array(r), np.array(c))), shape=(ndof, ndof)).toarray()
    return 0.5 * (A + A.T)


def h1_surface_dtn_defect(R, maxh, curve, order, n, P_radial=6):
    """Relative defect of the discrete H1 (0-form) DtN energy at degree n: |eig_n(S,MS) - (n+1)R|/((n+1)R)."""
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh))
    mesh.Curve(curve)
    fes = ng.H1(mesh, order=order)
    u, v = fes.TnT()
    nrm = ng.specialcf.normal(3)
    gu = ng.grad(u).Trace(); gv = ng.grad(v).Trace()
    gut = gu - (gu * nrm) * nrm; gvt = gv - (gv * nrm) * nrm
    ds = ng.ds(bonus_intorder=2 * order + 2 * curve)
    m = ng.BilinearForm(fes, check_unused=False); m += u * v * ds; m.Assemble()
    k = ng.BilinearForm(fes, check_unused=False); k += gut * gvt * ds; k.Assemble()
    MS = _to_dense(m.mat, fes.ndof); KS = _to_dense(k.mat, fes.ndof)
    d = np.diag(MS); bnd = np.where(d > 1e-10 * d.max())[0]
    MS = MS[np.ix_(bnd, bnd)]; KS = KS[np.ix_(bnd, bnd)]
    S = dtn_surface_matrix(MS, KS, P_radial, a=R)
    lam = sla.eigh(S, MS, eigvals_only=True)
    idx = sum(2 * j + 1 for j in range(n))          # start of the degree-n bucket
    grp = lam[idx:idx + 2 * n + 1]
    return abs(float(np.mean(grp)) - (n + 1) * R) / ((n + 1) * R), int(len(bnd))


# ===========================================================================
# [B] H(curl) (1-form) DtN on a curved mesh -- Kelvin-ball curl-curl FEM (act3_03)
# ===========================================================================
def hcurl_ball_dtn(R, maxh, order):
    """Effective H(curl) DtN eigenvalue lambda_vec = int nu'|curl A|^2 / oint |A_t|^2 for the z-dipole
    (degree n=1) in the Kelvin ball (reluctivity nu'=(rho'/R)^2).  Recovers the vector ladder n/R=1/R."""
    xx, yy, zz = ng.x, ng.y, ng.z
    r2 = xx * xx + yy * yy + zz * zz; r1 = ng.sqrt(r2)
    nup = r2 / (R * R)
    Adip = ng.CoefficientFunction((-yy, xx, 0.0)) / (4.0 * np.pi * r2 * r1)
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(min(order + 1, 3))
    fes = ng.HCurl(mesh, order=order, dirichlet=".*")
    u, v = fes.TnT()
    a = ng.BilinearForm(nup * ng.curl(u) * ng.curl(v) * ng.dx(bonus_intorder=12)
                        + 1e-6 * u * v * ng.dx(bonus_intorder=8)); a.Assemble()
    gf = ng.GridFunction(fes); gf.Set(Adip, ng.BND)
    rhs = gf.vec.CreateVector(); rhs.data = -(a.mat * gf.vec)
    gf.vec.data += a.mat.Inverse(fes.FreeDofs()) * rhs
    E = float(ng.Integrate(nup * ng.curl(gf) * ng.curl(gf) * ng.dx(bonus_intorder=12), mesh))
    bm = float(ng.Integrate(gf.Trace() * gf.Trace() * ng.ds(bonus_intorder=12), mesh))
    return E / bm, fes.ndof


# ===========================================================================
R = 1.0
MODES = [1, 2, 3]
print("=" * 98)
print(" act2_12 : the DtN spectrum is FORM-dependent -- two de Rham ladders, p- & h-convergence")
print("=" * 98)

# ---- [A] the two ladders, all three forms, radial-order p-convergence (rigorous IE) ----
print("\n[A] de Rham radial ladders (pure-numpy IE; order P=6, exact for P>=n+1):")
print("    n |  H1 / Hdiv-normal (scalar)  exact |  Hcurl-tangential (toroidal)  exact")
A_ok = True
for n in MODES:
    ds_ = ie_dtn(scalar_energy(n, 6)); dt_ = ie_dtn(toroidal_energy(n, 6))
    print(f"    {n} |        {ds_:8.4f}            {-(n+1):3d}  |          {dt_:8.4f}             {-n:3d}")
    A_ok = A_ok and abs(ds_ + (n + 1)) < 1e-8 and abs(dt_ + n) < 1e-8
check("[A] H1/Hdiv-normal radial IE DtN = -(n+1) (scalar ladder), n=1..3", A_ok)
tor_dist = all(abs(ie_dtn(toroidal_energy(n, 6)) + n) < 1e-8 and (n != n + 1) for n in MODES)
check("[A] the two ladders are DISTINCT: Hcurl -n != H1/Hdiv -(n+1)", tor_dist,
      "e.g. n=1: Hcurl -1 vs H1 -2")
print("\n[A.p] RADIAL p-convergence (the IE is SPECTRAL -- exact once order >= n+1):")
pconv_ie = {}
for form, energy, lad in (("scalar(H1/Hdiv)", scalar_energy, lambda n: n + 1),
                          ("toroidal(Hcurl)", toroidal_energy, lambda n: n)):
    defs = [abs(ie_dtn(energy(2, Pp)) + lad(2)) / lad(2) for Pp in range(1, 7)]
    pconv_ie[form] = defs
    print(f"    {form:18s} n=2 reldef vs P=1..6: " + ", ".join(f"{x:.0e}" for x in defs))
    check(f"[A.p] {form} exact for P>=n+1 (n=2 -> P=3 reldef<1e-9)", defs[2] < 1e-9, f"{defs[2]:.0e}")

# ---- [B] on a curved MESH (the Kelvin closure): H1 p- and h-convergence ----
print("\n[B.H1.p] H1 (0-form) DtN on the curved Kelvin mesh -- FE-order p-convergence (maxh=0.30, Curve=p+1):")
print("    p   d_1 (dipole defect)   ndof")
h1_pconv = {}
with ng.TaskManager():
    for p in (1, 2, 3, 4):
        d1, nd = h1_surface_dtn_defect(R, 0.30, min(p + 1, 5), p, n=1)
        h1_pconv[p] = {"d1": d1, "ndof": nd}
        print(f"    {p}   {d1:.2e}            {nd}")
check("[B.H1.p] H1 p-convergence: order 4 beats order 1 by >30x", h1_pconv[4]["d1"] < h1_pconv[1]["d1"] / 30.0,
      f"d1(p1)={h1_pconv[1]['d1']:.1e} -> d1(p4)={h1_pconv[4]['d1']:.1e}")
check("[B.H1.p] H1 p>=3 reaches < 1e-3 on a coarse mesh", h1_pconv[3]["d1"] < 1e-3, f"{h1_pconv[3]['d1']:.1e}")

print("\n[B.H1.h] H1 (0-form) DtN -- h-convergence (FE order p=2, Curve=3):")
print("    maxh   d_1        ndof")
MAXH = [0.5, 0.4, 0.3, 0.25]
h1_hconv = {}
with ng.TaskManager():
    for h in MAXH:
        d1, nd = h1_surface_dtn_defect(R, h, 3, 2, n=1)
        h1_hconv[h] = {"d1": d1, "ndof": nd}
        print(f"    {h:.2f}   {d1:.2e}   {nd}")
hs = np.array(MAXH); ds_h = np.array([h1_hconv[h]["d1"] for h in MAXH])
slope_h = float(np.polyfit(np.log(hs), np.log(ds_h), 1)[0])
print(f"    h-convergence slope (dlog d / dlog h) = {slope_h:.2f}  (expect ~2p = 4)")
check("[B.H1.h] H1 h-convergence algebraic slope ~ 2p=4 (in [3.0,5.0])", 3.0 <= slope_h <= 5.0,
      f"slope={slope_h:.2f}")

# ---- [B] the VECTOR ladder on a mesh: H(curl) p-convergence -> n/R ----
print("\n[B.Hcurl.p] H(curl) (1-form) DtN on the Kelvin-ball mesh -- the -n ladder (dipole 1/R), p-conv:")
print("    order   lambda_vec   ladder n/R=1   rel        ndof")
hc_pconv = {}
with ng.TaskManager():
    for order in (2, 3):
        lam, nd = hcurl_ball_dtn(R, 0.20, order)
        rel = abs(lam - 1.0 / R) / (1.0 / R)
        hc_pconv[order] = {"lambda": lam, "rel": rel, "ndof": nd}
        print(f"    {order}       {lam:.5f}      1.00000      {rel:.2e}   {nd}")
check("[B.Hcurl] H(curl) recovers the VECTOR ladder n/R=1/R (dipole, rel<5e-3)",
      hc_pconv[2]["rel"] < 5e-3, f"order2 rel={hc_pconv[2]['rel']:.1e}")
check("[B.Hcurl.p] H(curl) p-convergence: order 3 beats order 2",
      hc_pconv[3]["rel"] < hc_pconv[2]["rel"], f"{hc_pconv[2]['rel']:.1e} -> {hc_pconv[3]['rel']:.1e}")

# ---- p beats h: DoF the h-method needs to MATCH the p-method accuracy (manuscript's 20-80x) ----
print("\n[C] p- vs h-refinement DoF efficiency at MATCHED accuracy (H1 dipole defect):")
p_pt = h1_pconv[3]                                   # p=3, maxh=0.30 -- the p-method reference point
d_target = p_pt["d1"]; ndof_p = p_pt["ndof"]
h0 = 0.25; d0 = h1_hconv[h0]["d1"]; ndof0 = h1_hconv[h0]["ndof"]
h_match = h0 * (d_target / d0) ** (1.0 / slope_h)    # mesh size the h-method needs to reach d_target
ndof_h = ndof0 * (h0 / h_match) ** 2                 # surface DoF scale ~ h^-2
ratio = ndof_h / ndof_p
print(f"    p-method (order 3, maxh 0.30): d_1={d_target:.1e} at ndof={ndof_p}")
print(f"    h-method (order 2) to MATCH d_1={d_target:.1e}: maxh~{h_match:.3f} -> ndof~{ndof_h:.0f}")
print(f"    => p-refinement is ~{ratio:.0f}x more DoF-efficient at matched accuracy")
check("[C] p-refinement is >5x more DoF-efficient than h at matched accuracy (manuscript 20-80x)",
      ratio > 5.0, f"{ratio:.0f}x (p:{ndof_p} DoF vs h:{ndof_h:.0f} DoF)")

# ---- JSON ----
RESULTS = {
    "R": R, "modes": MODES,
    "ladders": {"scalar_H1_Hdiv_normal": "-(n+1)/R", "toroidal_Hcurl_tangential": "-n/R"},
    "A_radial_ie": {
        "scalar_dtn": [ie_dtn(scalar_energy(n, 6)) for n in MODES],
        "toroidal_dtn": [ie_dtn(toroidal_energy(n, 6)) for n in MODES],
        "p_convergence_n2": pconv_ie,
    },
    "B_h1_fem": {"p_convergence": {str(p): h1_pconv[p] for p in h1_pconv},
                 "h_convergence": {str(h): h1_hconv[h] for h in MAXH},
                 "h_slope": slope_h,
                 "p_vs_h_dof_efficiency": {"d_target": d_target, "ndof_p": ndof_p,
                                           "ndof_h_matched": ndof_h, "ratio": ratio}},
    "B_hcurl_fem": {str(o): hc_pconv[o] for o in hc_pconv},
    "n_fail": N_FAIL,
}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "act2_12_derham_dtn_pconv_hconv.json")
with open(out, "w") as fh:
    json.dump(RESULTS, fh, indent=2)
print(f"\n  wrote {os.path.basename(out)}")

print("\n" + "=" * 98)
print(" ALL CHECKS PASSED" if N_FAIL == 0 else f" {N_FAIL} CHECK(S) FAILED")
print("=" * 98)
assert N_FAIL == 0, f"{N_FAIL} checks failed"
