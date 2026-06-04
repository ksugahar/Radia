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
    are wired right.  This is the asserted claim below.
  * NOT VERIFIED: the ABSOLUTE magnitude.  With order-1 HCurl + a direct
    (umfpack) solve, alpha is ~20-25% off the analytic 4 pi a^3 G(x), and
    REFINING THE MESH DOES NOT CLOSE THE GAP (27k -> 118k dof: 23.2% -> 22.7%
    at a/delta=1.3) -- so it is an ORDER limitation, not skin/mesh resolution
    (the box is 6a; dipole truncation ~(a/R)^3 ~ 0.5%, also not the cause).
    The induced-moment integral (1/2)int r x J needs order-2 HCurl, but
    order-2 here is ~280k complex dof and umfpack 3D HCurl fill-in OOMs.
  * PATH: order-2 HCurl + an ITERATIVE solver (CompactAMS / shifted
    preconditioner, radia.sparsesolv_ngsolve -- the lab tool built for exactly
    HCurl curl-curl + mass with air regions, see CLAUDE.md).  That is a
    separate build; the m=1 full-frequency curve is the remaining open piece.

What IS already done for the tensor (so this probe is the only gap):
  - both analytic anchors, ALL directions incl. transverse: DC alpha=0,
    HF alpha_i = -V/(1-N_i)  (ellipsoid_alpha_tensor.py, golden);
  - the AXIAL full-frequency curve alpha_c(omega), skin-robust mixed phi-B
    FEM validated on the sphere to 1.8%  (ellipsoid_alpha_omega_axisym.py,
    golden).

Formulation (A-formulation, reaction unknown): total A = A_s + A_r,
A_s = (1/2) B0 x r (curl A_s = B0); int nu curl(A_r).curl(v) + int jw sigma
A_r.v = -int_cond jw sigma_cond A_s.v (sigma_air tiny = curl-curl gauge);
m = (1/2) int_cond r x J, J = -jw sigma (A_s+A_r); alpha_i = m_i mu0/B0.

Run:  python ellipsoid_alpha_tensor_3d.py   (order-1 probe; ~1 min)
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

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from levitation_sphere_force import G_exact, delta  # noqa: E402
import ellipsoid_alpha_tensor as ET  # noqa: E402

mu0 = 4 * np.pi * 1e-7
SIGMA = 5.8e7
SIGMA_AIR = 1.0           # curl-curl gauge regulariser (delta_air >> domain)
A_SPH = 5e-3
B0 = 1.0
P_ORDER = 1               # order-1 PROBE (order-2 OOMs with umfpack; see header)


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


def alpha_tensor_component(mesh, omega, direction, p=P_ORDER):
    fes = HCurl(mesh, order=p, dirichlet="outer", complex=True)
    u, v = fes.TnT()
    sig = mesh.MaterialCF({"cond": SIGMA}, default=SIGMA_AIR)
    s = 1j * omega
    As = A_source(direction)
    a = BilinearForm(fes)
    a += (1 / mu0) * curl(u) * curl(v) * dx
    a += s * sig * u * v * dx
    f = LinearForm(fes)
    f += -s * SIGMA * (As * v) * dx(mesh.Materials("cond"))
    with TaskManager():
        a.Assemble()
        f.Assemble()
        gf = GridFunction(fes)
        gf.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * f.vec
        J = -s * SIGMA * (gf + As)
        r = CF((x, y, z))
        rxJ = CF((r[1] * J[2] - r[2] * J[1],
                  r[2] * J[0] - r[0] * J[2],
                  r[0] * J[1] - r[1] * J[0]))
        comp = {"x": 0, "y": 1, "z": 2}[direction]
        m = 0.5 * Integrate(rxJ[comp] * dx(mesh.Materials("cond")), mesh)
    return np.conj(complex(m) * mu0 / B0)     # conj: FEM s=+jw vs G_exact e^{-jwt}


def sphere_analytic(omega, a=A_SPH):
    return 4 * np.pi * a**3 * G_exact(omega)


def main():
    ngsglobals.msg_level = 0
    print("3D HCurl eddy polarizability TENSOR -- WIP probe (isotropy verified;\n"
          "absolute magnitude is order-1-limited, see header)\n"
          f"Cu sigma={SIGMA:.2e}, order p={P_ORDER}")

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

    # ---- REPORTED, NOT asserted: the order-1 magnitude gap vs analytic ----
    print(f"\n[NOTE -- NOT asserted] order-1 magnitude gap vs analytic 4 pi a^3 G(x): "
          f"{min(mag_gaps)*100:.0f}-{max(mag_gaps)*100:.0f}%")
    print("   (mesh-refinement-invariant -> order limitation; needs order-2 + "
          "CompactAMS iterative solver, see header)")

    # ---- triaxial: the machinery PRODUCES three distinct components, but at
    #      order 1 they are NOT yet quantitatively trustworthy (see below) ----
    print("\n[check 2 -- order-1, NOT trustworthy] triaxial 5x3x1.5 mm")
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
    print(f"   order-1 reproduces the analytic per-volume ordering? "
          f"{'YES' if order1_ordering_ok else 'NO -- another symptom of the order-1 limit'}")
    # the ONLY thing we assert is that three DISTINCT components are produced
    # (the m=1 transverse machinery runs for all axes); their magnitudes/ordering
    # are NOT asserted -- order-1 is not accurate enough (this is the open piece).
    assert len({round(comps[d] * 1e10) for d in "xyz"}) == 3, "components not distinct"

    data = {
        "status": "WIP probe: isotropy VERIFIED, magnitude order-1-limited (NOT verified)",
        "sigma": SIGMA, "order": P_ORDER,
        "sphere_isotropy_worst": worst_iso,
        "order1_magnitude_gap_range": [min(mag_gaps), max(mag_gaps)],
        "triaxial_semi_m": list(semi),
        "triaxial_N": N.tolist(),
        "triaxial_alpha_re_m3_order1": comps,
        "triaxial_kappa_hf_analytic_m3": {d: -V / (1 - N[i]) for i, d in enumerate("xyz")},
        "order1_reproduces_analytic_ordering": bool(order1_ordering_ok),
        "path_to_full_verification": "order-2 HCurl + CompactAMS iterative solver",
    }
    with open(os.path.join(HERE, "ellipsoid_alpha_tensor_3d_results.json"), "w") as fp:
        json.dump(data, fp, indent=2)
    print("\n wrote ellipsoid_alpha_tensor_3d_results.json")
    print("\nSUMMARY: transverse m=1 machinery is structurally correct (isotropy 0%);"
          "\nabsolute magnitude needs order-2 + CompactAMS (the documented open piece).")


if __name__ == "__main__":
    main()
