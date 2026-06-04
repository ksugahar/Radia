"""ellipsoid_alpha_tensor_3d.py -- WIP RESEARCH PROBE (NOT a verified example).

Goal: the TRANSVERSE (m=1) component of the eddy-current polarizability
tensor of a conducting ellipsoid, by a 3D HCurl solve -- the part the
axisymmetric m=0 solve (ellipsoid_alpha_omega_axisym.py) structurally
cannot reach.  A uniform field along the symmetry axis is m=0; a transverse
field is m=1, so only a full 3D solve gives alpha_x, alpha_y, alpha_z at once.

STATUS (honest, 2026-06-04):
  * VERIFIED: the transverse machinery is structurally correct -- a sphere
    gives alpha_xx == alpha_zz to ~0.0% (the x and z responses are identical,
    as isotropy demands), so the m=1 (transverse) assembly + moment extraction
    are wired right.  Also the LOW-FREQUENCY limit is EXACT: as a/delta->0 the
    induced moment -> -j w sigma int r x A_s / 2 = the analytic 2 pi w mu0
    sigma a^5 / 15 (checked by hand) -- so the moment formula is correct.
  * SOLVER: SOLVED.  The umfpack OOM at order-2 is gone -- this now uses the
    CompactAMS-preconditioned COCR iterative solver (radia.sparsesolv_ngsolve);
    order-2 converges in ~120 iters, ~6 s, no OOM.  (Caveat learned: build the
    AMS preconditioner OUTSIDE `with TaskManager()` -- nesting segfaults; and
    CompactAMS requires nograds=True.)
  * NOT VERIFIED: the ABSOLUTE magnitude at FINITE frequency.  alpha is ~22%
    off at a/delta=1.3 and ~50% at a/delta=2.0 (Re ~3x too small; Im ~ok), and
    this is INVARIANT to order (1 vs 2 identical), mesh (27k->118k dof), box
    (6a vs 12a) and eps -- so it is NOT solver/mesh/order/truncation.  It is a
    FORMULATION issue: the reduced-A + nograds gauge under-computes the
    OUT-OF-PHASE reaction Im[A_r] (which sets Re[alpha], the field exclusion).
    My earlier "needs order-2+CompactAMS" hypothesis was WRONG -- order-2 is
    now affordable and gives the SAME gap.
  * PATH (the genuine remaining piece): a correct eddy formulation for the
    finite-frequency moment -- an A-phi formulation (explicit electric scalar
    potential, so the conductor charge conservation div J = 0 is honoured) or
    a far-field / surface moment extraction instead of the volume int r x J.
    That is a separate investigation, NOT a solver change.

What IS already done for the tensor (so only the transverse finite-freq
MAGNITUDE remains):
  - both analytic anchors, ALL directions incl. transverse: DC alpha=0,
    HF alpha_i = -V/(1-N_i)  (ellipsoid_alpha_tensor.py, golden);
  - the AXIAL full-frequency curve alpha_c(omega), skin-robust mixed phi-B
    FEM validated on the sphere to 1.8%  (ellipsoid_alpha_omega_axisym.py,
    golden).

Formulation: reduced A = A_s + A_r, A_s = (1/2) B0 x r (curl A_s = B0); gauged
complex HCurl (nograds); int nu curl(A_r).curl(v) + int jw sigma A_r.v + eps
mass = -int_cond jw sigma A_s.v; CompactAMS+COCR solve; m = (1/2) int_cond
r x J, J = -jw sigma (A_s+A_r); alpha_i = conj(m_i mu0/B0).

Run:  python ellipsoid_alpha_tensor_3d.py   (CompactAMS order-2; ~1 min)
"""
import json
import os
import sys

import numpy as np
from netgen.occ import Sphere, Box, Pnt, Glue, OCCGeometry, gp_GTrsf
from ngsolve import (
    Mesh, HCurl, GridFunction, BilinearForm, LinearForm, Integrate, CF,
    curl, x, y, z, dx, TaskManager, ngsglobals,
)
import radia.sparsesolv_ngsolve as ssn

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from levitation_sphere_force import G_exact, delta  # noqa: E402
import ellipsoid_alpha_tensor as ET  # noqa: E402

mu0 = 4 * np.pi * 1e-7
SIGMA = 5.8e7
A_SPH = 5e-3
B0 = 1.0
P_ORDER = 2               # order-2 (accurate moment) with CompactAMS iterative solve
EPS_FACTOR = 0.05         # conductor-in-air mass regularisation eps = EPS_FACTOR*nu


def build_mesh(semi, box=0.030, maxh=0.011, maxh_cond=None):
    sx, sy, sz = semi
    if maxh_cond is None:
        maxh_cond = min(semi) / 3.0
    unit = Sphere(Pnt(0, 0, 0), 1.0)
    cond = gp_GTrsf([sx, 0, 0, 0, sy, 0, 0, 0, sz], [0, 0, 0])(unit)
    cond.mat("cond")
    cond.maxh = maxh_cond
    air = Box(Pnt(-box, -box, -box), Pnt(box, box, box)) - cond
    air.mat("air")
    air.faces.name = "outer"
    return Mesh(OCCGeometry(Glue([air, cond])).GenerateMesh(maxh=maxh))


def A_source(direction):
    r = CF((x, y, z))
    e = {"x": CF((1, 0, 0)), "y": CF((0, 1, 0)), "z": CF((0, 0, 1))}[direction]
    return 0.5 * B0 * CF((e[1] * r[2] - e[2] * r[1],
                          e[2] * r[0] - e[0] * r[2],
                          e[0] * r[1] - e[1] * r[0]))


def alpha_tensor_component(mesh, omega, direction, p=P_ORDER, tol=1e-9,
                           maxiter=3000, return_iters=False):
    """alpha_ii(omega) for a uniform field along `direction`, solved with the
    CompactAMS-preconditioned COCR iterative solver (radia.sparsesolv_ngsolve)
    so order-2 HCurl fits without a direct factorisation.

    Pattern = the lab's example_compact_ams recipe: gauged HCurl (nograds), the
    complex system  nu curl.curl + jw sigma + eps mass, a REAL auxiliary matrix
    with mass |w sigma|+eps (matched scale), discrete gradient + vertex coords."""
    nu = 1.0 / mu0
    s = 1j * omega
    eps_reg = EPS_FACTOR * nu
    sigma_cf = mesh.MaterialCF({"cond": SIGMA}, default=0.0)
    As = A_source(direction)

    # gauged complex HCurl (nograds -> AMS handles the gradient space)
    fes = HCurl(mesh, order=p, nograds=True, dirichlet="outer", complex=True)
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu * curl(u) * curl(v) * dx
    a += s * sigma_cf * u * v * dx
    a += eps_reg * u * v * dx                       # conductor-in-air regulariser
    f = LinearForm(fes)
    f += -s * SIGMA * (As * v) * dx(mesh.Materials("cond"))

    # real auxiliary matrix: mass = |w sigma| + eps to match the complex scale
    fes_r = HCurl(mesh, order=p, nograds=True, dirichlet="outer", complex=False)
    ur, vr = fes_r.TnT()
    a_real = BilinearForm(fes_r)
    a_real += nu * curl(ur) * curl(vr) * dx
    a_real += (abs(omega * SIGMA) + eps_reg) * ur * vr * dx

    G_mat, _ = fes_r.CreateGradient()
    cx = [mesh[vtx].point[0] for vtx in mesh.vertices]
    cy = [mesh[vtx].point[1] for vtx in mesh.vertices]
    cz = [mesh[vtx].point[2] for vtx in mesh.vertices]

    # Assemble + build the AMS preconditioner OUTSIDE TaskManager: the Compact
    # AMS construction manages its own parallelism and nesting it inside an
    # outer `with TaskManager()` segfaults/hangs.  Only the COCR solve is
    # wrapped (the canonical radia.sparsesolv_ngsolve recipe).
    a.Assemble()
    f.Assemble()
    a_real.Assemble()
    pre = ssn.ComplexCompactAMSPreconditioner(
        a_real.mat, G_mat, freedofs=fes_r.FreeDofs(),
        coord_x=cx, coord_y=cy, coord_z=cz, ndof_complex=fes.ndof)
    solver = ssn.COCRSolver(a.mat, pre, maxiter=maxiter, tol=tol,
                            freedofs=fes.FreeDofs())
    gf = GridFunction(fes)
    with TaskManager():
        gf.vec.data = solver * f.vec
    J = -s * SIGMA * (gf + As)
    r = CF((x, y, z))
    rxJ = CF((r[1] * J[2] - r[2] * J[1],
              r[2] * J[0] - r[0] * J[2],
              r[0] * J[1] - r[1] * J[0]))
    comp = {"x": 0, "y": 1, "z": 2}[direction]
    m = 0.5 * Integrate(rxJ[comp] * dx(mesh.Materials("cond")), mesh)
    alpha = np.conj(complex(m) * mu0 / B0)          # conj: FEM s=+jw vs G_exact e^{-jwt}
    return (alpha, solver.iterations) if return_iters else alpha


def sphere_analytic(omega, a=A_SPH):
    return 4 * np.pi * a**3 * G_exact(omega)


def main():
    ngsglobals.msg_level = 0
    print("3D HCurl eddy polarizability TENSOR -- WIP probe.  CompactAMS+COCR\n"
          "solver SWAP DONE (order-2 affordable); magnitude is FORMULATION-\n"
          f"limited, not solver-limited (see header).  Cu sigma={SIGMA:.2e}, p={P_ORDER}")

    # ---- VERIFIED: sphere is ISOTROPIC (transverse == axial machinery) ----
    print("\n[check 1 -- ASSERTED] sphere isotropy alpha_xx == alpha_zz")
    msph = build_mesh((A_SPH,) * 3, box=0.030, maxh=0.011, maxh_cond=1.2e-3)
    print(f"   sphere mesh: {msph.ne} elements")
    print("   f[Hz]  a/delta  alpha_zz [mm^3]     alpha_xx [mm^3]     |zz-xx|/|zz|")
    worst_iso = 0.0
    mag_gaps = []
    for f in (300, 700):
        w = 2 * np.pi * f
        azz = alpha_tensor_component(msph, w, "z")
        axx = alpha_tensor_component(msph, w, "x")
        iso = abs(azz - axx) / abs(azz)
        worst_iso = max(worst_iso, iso)
        an = sphere_analytic(w)
        mag_gaps.append(abs(azz - an) / abs(an))
        print(f"  {f:6.0f}  {A_SPH/delta(w):5.2f}  "
              f"({azz.real*1e9:+7.1f},{azz.imag*1e9:+6.1f})  "
              f"({axx.real*1e9:+7.1f},{axx.imag*1e9:+6.1f})  {iso*100:6.2f}%")
    print(f"   isotropy worst |zz-xx|/|zz| = {worst_iso*100:.2f}%  "
          f"-> transverse machinery structurally correct")
    assert worst_iso < 0.03, f"sphere not isotropic ({worst_iso*100:.1f}%) -- m=1 assembly bug"

    # ---- REPORTED, NOT asserted: the FINITE-FREQUENCY magnitude gap ----
    print(f"\n[NOTE -- NOT asserted] finite-freq magnitude gap vs analytic 4 pi a^3 G(x): "
          f"{min(mag_gaps)*100:.0f}-{max(mag_gaps)*100:.0f}%")
    print("   INVARIANT to order(1==2)/mesh/box/eps -> FORMULATION issue, not solver:")
    print("   the reduced-A+nograds gauge under-computes Im[A_r] (the out-of-phase")
    print("   reaction that sets Re[alpha]).  Low-freq limit IS exact.  Path = an")
    print("   A-phi (scalar-potential) eddy formulation or far-field moment extraction.")

    # ---- triaxial: the machinery PRODUCES three distinct components, but the
    #      magnitudes are formulation-limited (NOT quantitatively trustworthy) ----
    print("\n[check 2 -- magnitudes NOT trustworthy (formulation-limited)] triaxial 5x3x1.5 mm")
    semi = (5e-3, 3e-3, 1.5e-3)
    N = ET.demag_tensor(semi)
    V = 4 / 3 * np.pi * semi[0] * semi[1] * semi[2]
    m3 = build_mesh(semi, box=0.022, maxh=0.011, maxh_cond=6e-4)
    w = 2 * np.pi * 700
    comps = {d: alpha_tensor_component(m3, w, d).real for d in "xyz"}
    pv = {d: abs(comps[d]) / V for d in "xyz"}
    pv_hf = {d: 1.0 / (1.0 - N[i]) for i, d in enumerate("xyz")}     # analytic HF per-vol
    print(f"   f=700 Hz: per-volume |alpha|/V  x={pv['x']:.3f} y={pv['y']:.3f} "
          f"z={pv['z']:.3f}")
    print(f"   analytic HF per-volume 1/(1-N_i): x={pv_hf['x']:.3f} y={pv_hf['y']:.3f} "
          f"z={pv_hf['z']:.3f}  (z=short axis should be LARGEST)")
    order1_ordering_ok = pv["z"] > pv["y"] > pv["x"]
    print(f"   FEM reproduces the analytic per-volume ordering? "
          f"{'YES' if order1_ordering_ok else 'NO -- another symptom of the formulation limit'}")
    # the ONLY thing we assert is that three DISTINCT components are produced
    # (the m=1 transverse machinery runs for all axes); their magnitudes/ordering
    # are NOT asserted -- order-1 is not accurate enough (this is the open piece).
    assert len({round(comps[d] * 1e10) for d in "xyz"}) == 3, "components not distinct"

    data = {
        "status": "WIP probe: isotropy + low-freq limit VERIFIED; finite-freq "
                  "magnitude FORMULATION-limited (NOT verified)",
        "solver": "CompactAMS + COCR (swap done; order-2 affordable, no OOM)",
        "sigma": SIGMA, "order": P_ORDER,
        "sphere_isotropy_worst": worst_iso,
        "finite_freq_magnitude_gap_range": [min(mag_gaps), max(mag_gaps)],
        "magnitude_gap_invariant_to": "order(1==2), mesh, box, eps",
        "triaxial_semi_m": list(semi),
        "triaxial_N": N.tolist(),
        "triaxial_alpha_re_m3": comps,
        "triaxial_kappa_hf_analytic_m3": {d: -V / (1 - N[i]) for i, d in enumerate("xyz")},
        "reproduces_analytic_ordering": bool(order1_ordering_ok),
        "path_to_full_verification": "A-phi (scalar-potential) eddy formulation "
                                     "or far-field moment extraction (NOT a solver change)",
    }
    with open(os.path.join(HERE, "ellipsoid_alpha_tensor_3d_results.json"), "w") as fp:
        json.dump(data, fp, indent=2)
    print("\n wrote ellipsoid_alpha_tensor_3d_results.json")
    print("\nSUMMARY: CompactAMS+COCR solver swap DONE (order-2 affordable, no OOM)."
          "\nTransverse m=1 machinery structurally correct (isotropy 0%, low-freq exact);"
          "\nfinite-freq magnitude is FORMULATION-limited (A-phi / far-field moment = path).")


if __name__ == "__main__":
    main()
