"""coil_sphere_eddy_force.py -- quantify the POINT-DIPOLE approximation error
in AC magnetic levitation by computing the FULL axisymmetric eddy-current
Lorentz force on a conducting sphere in a coil's field and comparing it to
the induced-DIPOLE force at the same height, in the a/L ~ 0.5 regime.

Motivation (the honest caveat in coil_maglev_equilibrium.py / README.md):
    at a = 5 mm the on-axis field-gradient scale L ~ 10 mm gives a/L ~ 0.5,
    so the induced-point-dipole force
        F_z^dipole = (pi a^3 / mu0) Re[G(omega)] d|B0(z)|^2/dz
    is only APPROXIMATE that close to the coil.  "A full Radia + NGSolve
    eddy-current reaction-field solve would refine the force in the a/L ~ 0.5
    regime."  This script IS that solve.

Method -- reuse the lab-verified axisymmetric mixed phi-B eddy-current
machinery (team28_axisym_fem.py / team28_cln_force.py, which reproduce the
TEAM 28 ground-truth Lorentz force to 0.01%):

  * mixed space    fesPhi = H1(order=p, dirichlet="outer", complex),
                   fesB   = VectorL2(order=p-1, complex),  fes = fesPhi*fesB;
  * B operator     Bop(phi) = (-1/r dphi/dz, 1/r dphi/dr);
  * system         (K + s*N) X = F  with
                   K = nu*r*(B*dB - Bop(u)*dB + Bop(v)*B),
                   N = v*(sigma*u/r),  s = j*2*pi*f,
                   F = v*J_phi   (azimuthal coil source current density);
  * Lorentz force  the azimuthal eddy current Jt = -sigma*s*phi/r and the
                   radial field B_r give the axial force density Jt*B_r; the
                   PEAK-amplitude time average is
                       Fz = Integrate((1/2)Re[Jt conj(B_r)] * 2*pi*r,
                                      over the SPHERE material).
                   This is the SAME (1/2 peak-amplitude) convention as the
                   induced-dipole coefficient below, so the two forces are
                   directly comparable.  (The verbatim team28 integral
                   Re[B_r*Jt] omits the 1/2 / uses the non-conjugated product
                   -- a lab-.mat normalisation that is exactly 2x this here;
                   we use the proper conjugate form to share one convention.)
  * open boundary  the anisotropic-nu infinite-element semicircle shell
                   (rings u1..u4, outer edge "outer") so the coil field
                   decays correctly with no air-box truncation.

Two solves at f = 50 kHz:
  (a) sphere CONDUCTING (sigma = 5.8e7)  -> full eddy Lorentz force F_z^full;
  (b) sphere region sigma = 0 (pure coil field) -> |B0(z)|^2 on the axis,
      from which F_z^dipole and the local gradient scale L = |B0|^2 /
      |d|B0|^2/dz| (hence a/L) are computed.

The ONLY difference between the two forces is the approximation (same FEM
coil field feeds both), so F_z^full / F_z^dipole is the point-dipole error.

Verification (hard asserts, "No Fallbacks -- Fail Fast"):
  * coil-only on-axis field matches the analytic loop
    Bz = mu0 NI R^2 / (2 (R^2+z^2)^1.5) within ~5% (finite coil cross-
    section broadens the field; tolerance documented);
  * both forces are LIFT (same sign, point up) and the same order of
    magnitude (ratio in [0.5, 2.0] at a/L ~ 0.5 -- close but NOT identical:
    that nonzero gap is the result);
  * SHRINKING the sphere (a = 1 mm at the same height, a/L ~ 0.1) drives
    F_z^full / F_z^dipole -> 1: the dipole approximation becomes exact as
    a/L -> 0.  This is the key check that the gap is genuinely the dipole
    approximation and not a modelling artefact.

Run:  python coil_sphere_eddy_force.py     (~1-2 min: four axisym solves)
"""
import json
import os
import sys

import numpy as np
from numpy import pi
from netgen.occ import WorkPlane, MoveTo, Glue, OCCGeometry
from ngsolve import (
    Mesh, H1, L2, VectorL2, GridFunction, BilinearForm, LinearForm,
    Integrate, CF, x, dx, grad, TaskManager, ngsglobals,
)

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
# reuse the verified copper parameters + skin-depth helper.  NOTE: the
# imported G_exact() hardcodes the module-level a = 5 mm, so for the 1 mm
# sanity case we MUST recompute G at the correct radius (G depends on
# a/delta).  delta() depends only on (omega, sigma), not on a, so it is safe
# to reuse directly.
from maglev_sphere_force import delta, _cot, sigma as SIGMA_CU  # noqa: E402

mu0 = 4 * pi * 1e-7


def G_of_radius(a_sphere, omega):
    """Closed-form sphere magnetic response G(x), x=(1+i) a_sphere/delta.
    Landau-Lifshitz ECM sec. 59 (same as maglev_sphere_force.G_exact but
    with the radius passed explicitly so it is correct for ANY a_sphere)."""
    xv = (1 + 1j) * a_sphere / delta(omega)
    return -0.5 * (1.0 - 3.0 / xv**2 + 3.0 / xv * _cot(xv))

# ----- coil (matches coil_maglev_equilibrium.py: 30 mm, 10000 At) -----
R_COIL = 0.030          # coil mean radius [m]
NI = 100 * 100.0        # total amp-turns (peak) = N_turn * I_peak
COIL_W = 0.002          # coil cross-section radial width  [m] (thin)
COIL_H = 0.002          # coil cross-section axial height   [m] (thin)
COIL_Z = 0.0            # coil cross-section centre height  [m] (z = 0 plane)

F_HZ = 50e3             # AC frequency [Hz]
OMEGA = 2 * pi * F_HZ

# ----- sphere: 5 mm Cu, centred on axis at z0 where a/L ~ 0.5 -----
# (the on-axis a/L of a 30 mm loop peaks ~0.50 at z ~ 30 mm; see the probe in
#  the report.  z0 = 30 mm puts the 5 mm sphere squarely at a/L ~ 0.5.)
A_SPHERE_5MM = 5e-3
A_SPHERE_1MM = 1e-3
Z0 = 0.030              # sphere centre height [m]

# ----- open-boundary infinite shell (TEAM 28 anisotropic-nu values) -----
DOMAIN_A = 0.120        # outer air radius [m]
DOMAIN_Z = 0.010        # vertical centre of the air domain (between coil & sphere)
MAXH = 0.004            # background mesh size [m]
# Sphere mesh sizes chosen from the convergence study (see report): the
# full/dipole ratio creeps with refinement until the sphere surface (where
# the eddy current concentrates) and the skin layer are resolved.  250 um on
# the 5 mm sphere and 80 um on the 1 mm sphere put both ratios on their
# converged plateau (5 mm: ~+1.9%, stable to p4; 1 mm: ~-0.1%, near zero).
MAXH_SPHERE_5 = 2.5e-4  # 5 mm sphere (~20 elems across radius; converged)
MAXH_SPHERE_1 = 8e-5    # 1 mm sphere (a/delta=3.38: resolve the skin layer)
MAXH_COIL = 1.5e-3
P_ORDER = 3             # FE order (phi); B is order p-1 (p3=p4 at fixed h)


def Bz_loop_analytic(z, R, NI_amp):
    """On-axis Bz of an ideal thin circular loop (radius R, NI amp-turns)."""
    return mu0 * NI_amp * R**2 / (2.0 * (R**2 + z**2) ** 1.5)


def Bop(phi):
    """B from the axisymmetric stream function phi: B = (-1/r dphi/dz, 1/r dphi/dr)."""
    g = grad(phi)
    return CF((-1 / x * g[1], 1 / x * g[0]))


def build_mesh(a_sphere, maxh_sphere):
    """Axisymmetric (r>=0 half-plane) mesh:
       infinite-shell rings u1..u4 (open boundary, outer edge 'outer'),
       a thin coil rectangle at r ~ R_COIL ('Coil'),
       a half-disk sphere centred at (0, Z0) ('Sphere'), air elsewhere.
    """
    # --- infinite-element semicircle shell (TEAM 28 pattern) ---
    semicircle = []
    abc = np.linspace(1, 1.1, 5)
    for n in range(len(abc)):
        circle = MoveTo(0, DOMAIN_Z).Circle(abc[n] * DOMAIN_A).Face()
        rect = MoveTo(-abc[n] * DOMAIN_A, -abc[n] * DOMAIN_A + DOMAIN_Z) \
            .Rectangle(abc[n] * DOMAIN_A, 2 * abc[n] * DOMAIN_A).Face()
        semicircle.append(circle - rect)
        semicircle[n].bc('Air' if n == 0 else f'u{n}')
        if n == len(abc) - 1:
            semicircle[n].edges.name = 'outer'

    # --- thin coil rectangle centred at (R_COIL, COIL_Z) ---
    coil = MoveTo(R_COIL - COIL_W / 2, COIL_Z - COIL_H / 2) \
        .Rectangle(COIL_W, COIL_H).Face()
    coil.bc('Coil')
    coil.maxh = MAXH_COIL

    # --- sphere: full disk centred at (0, Z0), cut to the r>=0 half-plane ---
    sphere_full = WorkPlane().Circle(0.0, Z0, a_sphere).Face()
    margin = 0.05
    cutter = MoveTo(-a_sphere - margin, Z0 - a_sphere - margin) \
        .Rectangle(a_sphere + margin, 2 * a_sphere + 2 * margin).Face()
    sphere = sphere_full - cutter
    sphere.name = 'Sphere'
    sphere.maxh = maxh_sphere

    shape = Glue(semicircle + [coil, sphere])
    geo = OCCGeometry(shape, dim=2)
    return Mesh(geo.GenerateMesh(maxh=MAXH))


def assemble(mesh, sphere_sigma, p=P_ORDER):
    """Assemble the mixed phi-B (K + s*N) X = F system for the coil + sphere.
       sphere_sigma = SIGMA_CU (conducting) or 0.0 (pure coil field)."""
    fesPhi = H1(mesh, order=p, dirichlet="outer", complex=True)
    fesB = VectorL2(mesh, order=p - 1, complex=True)
    fes = fesPhi * fesB
    (u, B), (v, dB) = fes.TnT()

    nu = mesh.MaterialCF(
        {"u1": 1 / (0.285096 * mu0), "u2": 1 / (12.4933 * mu0),
         "u3": 1 / (0.051874 * mu0), "u4": 1 / (147.652 * mu0)},
        default=1 / mu0)
    sig = mesh.MaterialCF({"Sphere": sphere_sigma}, default=0.0)
    s = OMEGA * 1j

    # azimuthal coil current density (amp-turns spread over the cross-section)
    Jphi = mesh.MaterialCF({"Coil": NI / (COIL_W * COIL_H)}, default=0.0)

    K = BilinearForm(nu * x * (B * dB - Bop(u) * dB + Bop(v) * B) * dx)
    K += (v * (s * sig * (u / x))) * dx
    F = LinearForm((v * Jphi) * dx)
    K.Assemble()
    F.Assemble()

    gfuB = GridFunction(fes)
    gfuB.vec.data = K.mat.Inverse(fes.FreeDofs()) * F.vec
    return fes, fesPhi, fesB, sig, gfuB, s


def lorentz_force_on_sphere(fes, fesPhi, fesB, sig, gfuB, s, mesh, p=P_ORDER):
    """Time-averaged vertical Lorentz force on the sphere.

    The azimuthal eddy current Jt = -sigma*s*phi/r and the radial field B_r
    give the axial force density (J x B)_z = Jt * B_r.  The time-AVERAGE of a
    product of two PEAK-amplitude phasors is  <f_z> = (1/2) Re[Jt conj(B_r)].
    This is the convention that matches the induced-dipole formula
        <F_z>^dipole = (pi a^3/mu0) Re[G] d|B0|^2/dz
    of maglev_sphere_force.py, whose coefficient is likewise the (1/2)
    PEAK-amplitude time-average (B0 = peak field).  Using the same coil
    amp-turn amplitude NI for both sides makes the comparison apples-to-apples,
    so the ONLY difference is the dipole approximation.

    (NOTE on the verbatim team28 integral  Re[B_r]Re[Jt] - Im[B_r]Im[Jt] =
     Re[B_r * Jt]: that omits the 1/2 and uses the non-conjugated product;
     for the TEAM 28 problem it is calibrated to the lab .mat normalisation,
     but it is exactly 2x the proper PEAK-amplitude time-average here -- a
     pure normalisation convention, not the dipole error.  We use the proper
     (1/2) conjugate form so the full and dipole forces share one convention.)

    Returns SIGNED F_z [N] with positive = upward LIFT (sphere above coil,
    J_phi>0 -> diamagnetic lift is UP).  The raw (1/2)Re[Jt conj(B_r)]
    integral comes out negative for lift in this (r,z) sign setup, so we
    negate to express it as positive-up, matching the dipole sign.
    """
    gfu = GridFunction(fesPhi); gfu.Set(gfuB.components[0])
    gfB1 = GridFunction(fesB);  gfB1.Set(gfuB.components[1])
    fesL2 = L2(mesh, order=p, dirichlet="outer", complex=True)
    gfJt = GridFunction(fesL2); gfJt.Set(-sig * s * gfu / x)   # azimuthal eddy J
    Br = gfB1[0]
    # <f_z> = (1/2) Re[Jt conj(B_r)] = (1/2)(Re Jt Re Br + Im Jt Im Br)
    fz = 0.5 * Integrate(
        (gfJt.real * Br.real + gfJt.imag * Br.imag)
        * (2 * pi * x) * dx(mesh.Materials("Sphere")), mesh)
    return -float(fz.real)        # flip so positive = upward lift


def b2_on_axis_from_solution(fesB, gfuB, zs, r_probe=2e-4):
    """|B0(z)|^2 [T^2] on (near) the axis, from the SOLVED VectorL2 B field.

    The VectorL2 part stores B directly (no 1/r), so evaluating it at a
    tiny r_probe avoids the 1/r singularity of Bop on the axis; on the axis
    B is purely axial, so |B|^2 ~ B_z^2 and the residual r_probe bias is
    O(r_probe^2) -- negligible at r_probe = 0.2 mm for L ~ 10 mm.
    """
    gfB1 = GridFunction(fesB); gfB1.Set(gfuB.components[1])
    mesh = fesB.mesh
    b2 = np.empty(len(zs))
    for i, z in enumerate(zs):
        val = gfB1(mesh(r_probe, float(z)))   # (B_r, B_z) complex
        Br, Bz = val[0], val[1]
        # peak amplitude squared = |B_r|^2 + |B_z|^2 (real spatial field, but
        # store as complex; |.|^2 is the squared magnitude of the phasor)
        b2[i] = abs(Br) ** 2 + abs(Bz) ** 2
    return b2


def run_case(a_sphere, maxh_sphere, label):
    """Full vs dipole force for one sphere radius.  Returns a result dict."""
    print("\n" + "-" * 70)
    print(f"CASE {label}:  a = {a_sphere*1e3:.1f} mm sphere centred at "
          f"z0 = {Z0*1e3:.1f} mm")
    print("-" * 70)

    d = delta(OMEGA)
    a_over_delta = a_sphere / d
    ReG = G_of_radius(a_sphere, OMEGA).real     # G at THIS sphere's a/delta
    print(f"  skin depth delta = {d*1e6:.2f} um,  a/delta = {a_over_delta:.2f}")
    print(f"  Re[G(omega)] = {ReG:+.5f}  (require < 0 = diamagnetic lift)")
    assert ReG < 0.0, f"Re[G] = {ReG:.4f} not negative (no lift)"

    with TaskManager():
        # ---------- solve (a): conducting sphere -> full eddy force ----------
        mesh = build_mesh(a_sphere, maxh_sphere)
        mats = mesh.GetMaterials()
        assert "Sphere" in mats and "Coil" in mats, \
            f"mesh missing Sphere/Coil materials: {mats}"
        ne = mesh.ne
        fes, fesPhi, fesB, sig, gfuB, s = assemble(mesh, SIGMA_CU)
        ndof = fes.ndof
        F_full = lorentz_force_on_sphere(
            fes, fesPhi, fesB, sig, gfuB, s, mesh)
        print(f"  [full]  ne = {ne}, ndof = {ndof}")
        print(f"  [full]  F_z^full  = {F_full*1e3:+.5f} mN")

        # ---------- solve (b): sphere region air -> coil-only field ----------
        fes0, fesPhi0, fesB0, sig0, gfuB0, s0 = assemble(mesh, 0.0)

        # coil-only on-axis |B0|^2 profile (for the dipole force + a/L)
        h = max(2e-4, 0.05 * a_sphere)          # central-difference step
        zs = np.array([Z0 - 2 * h, Z0 - h, Z0, Z0 + h, Z0 + 2 * h])
        b2 = b2_on_axis_from_solution(fesB0, gfuB0, zs)
        # 4th-order central difference of |B0|^2 at z0
        dB2 = (-b2[4] + 8 * b2[3] - 8 * b2[1] + b2[0]) / (12 * h)
        b2_0 = b2[2]
        L_grad = b2_0 / abs(dB2)
        a_over_L = a_sphere / L_grad
        F_dipole = (pi * a_sphere**3 / mu0) * ReG * dB2

        # ---------- CHECK 1: coil-only on-axis B vs analytic loop ----------
        # sample on the axis AWAY from the coil (z = 8..20 mm) where the
        # finite cross-section matters least; allow 5% (documented).
        print("\n  [check 1] coil-only on-axis Bz vs analytic loop")
        print("     z[mm]   Bz_fem[T]      Bz_loop[T]     rel.err%")
        gfB0 = GridFunction(fesB0); gfB0.Set(gfuB0.components[1])
        max_err = 0.0
        for zc in (0.060, 0.080, 0.100, 0.120):
            Bz_fem = abs(gfB0(mesh(2e-4, zc))[1])
            Bz_an = Bz_loop_analytic(zc, R_COIL, NI)
            err = abs(Bz_fem - Bz_an) / abs(Bz_an) * 100.0
            max_err = max(max_err, err)
            print(f"    {zc*1e3:6.1f}  {Bz_fem: .6e}  {Bz_an: .6e}   {err:.3f}")
        print(f"     max field error = {max_err:.3f} %  (require < 5 %)")
        assert max_err < 5.0, \
            f"coil-only FEM field deviates {max_err:.2f}% from analytic loop"

    print(f"\n  [dipole] |B0(z0)|^2 = {b2_0:.4e} T^2, "
          f"d|B0|^2/dz = {dB2:+.4e} T^2/m")
    print(f"  [dipole] field-gradient scale L = {L_grad*1e3:.3f} mm, "
          f"a/L = {a_over_L:.3f}")
    print(f"  [dipole] F_z^dipole = {F_dipole*1e3:+.5f} mN")

    ratio = F_full / F_dipole
    disc = (ratio - 1.0) * 100.0
    print(f"\n  >>> F_z^full / F_z^dipole = {ratio:.4f}  "
          f"({disc:+.2f}% discrepancy)")

    # both forces must be lift (same sign, point up)
    assert F_full > 0.0, f"F_z^full = {F_full:.3e} not a lift (positive)"
    assert F_dipole > 0.0, f"F_z^dipole = {F_dipole:.3e} not a lift (positive)"

    return {
        "label": label,
        "a_sphere_m": a_sphere,
        "z0_m": Z0,
        "ne": int(ne),
        "ndof": int(ndof),
        "maxh_sphere_m": maxh_sphere,
        "a_over_delta": a_over_delta,
        "Re_G": ReG,
        "b2_axis_T2": float(b2_0),
        "dB2dz_T2_per_m": float(dB2),
        "L_grad_m": float(L_grad),
        "a_over_L": float(a_over_L),
        "F_full_N": F_full,
        "F_dipole_N": F_dipole,
        "ratio_full_over_dipole": ratio,
        "discrepancy_percent": disc,
        "coil_field_max_err_percent": float(max_err),
    }


def main():
    ngsglobals.msg_level = 0
    print("=" * 70)
    print("Point-dipole error in AC levitation: FULL eddy-current Lorentz")
    print("force on a Cu sphere vs the induced-dipole force, at a/L ~ 0.5")
    print("=" * 70)
    print(f"coil:   R = {R_COIL*1e3:.1f} mm, NI = {NI:.0f} At, "
          f"cross-section {COIL_W*1e3:.1f} x {COIL_H*1e3:.1f} mm at z = "
          f"{COIL_Z*1e3:.1f} mm")
    print(f"freq:   f = {F_HZ/1e3:.1f} kHz,  Cu sigma = {SIGMA_CU:.2e} S/m")
    print(f"FE:     order p = {P_ORDER} (phi), B order {P_ORDER-1}; "
          f"open-boundary infinite shell")

    # ---- 5 mm sphere: the a/L ~ 0.5 regime of interest ----
    res5 = run_case(A_SPHERE_5MM, MAXH_SPHERE_5, "5mm")

    # ---- 1 mm sphere: a/L ~ 0.1, dipole should become exact ----
    res1 = run_case(A_SPHERE_1MM, MAXH_SPHERE_1, "1mm")

    # ================= KEY CHECK: ratio -> 1 as a/L -> 0 =================
    print("\n" + "=" * 70)
    print("KEY CHECK: the full-vs-dipole gap is the DIPOLE APPROXIMATION")
    print("=" * 70)
    print(f"  a = 5 mm:  a/L = {res5['a_over_L']:.3f}, "
          f"ratio = {res5['ratio_full_over_dipole']:.4f} "
          f"({res5['discrepancy_percent']:+.2f}%)")
    print(f"  a = 1 mm:  a/L = {res1['a_over_L']:.3f}, "
          f"ratio = {res1['ratio_full_over_dipole']:.4f} "
          f"({res1['discrepancy_percent']:+.2f}%)")
    dev5 = abs(res5["ratio_full_over_dipole"] - 1.0)
    dev1 = abs(res1["ratio_full_over_dipole"] - 1.0)
    print(f"  |ratio-1|:  5 mm = {dev5*100:.2f}%   1 mm = {dev1*100:.2f}%")

    # ratio must be the right order of magnitude (close, not identical)
    assert 0.5 <= res5["ratio_full_over_dipole"] <= 2.0, \
        ("5 mm full/dipole ratio outside [0.5,2.0]: %.3f (open boundary or "
         "force sign wrong?)" % res5["ratio_full_over_dipole"])
    # the 1 mm (a/L ~ 0.1) ratio must be CLOSER to 1 than the 5 mm ratio:
    # this is the genuine signature that the gap is the dipole approximation
    assert dev1 < dev5, \
        ("shrinking the sphere did NOT move the ratio toward 1 "
         "(5 mm |ratio-1|=%.3f, 1 mm |ratio-1|=%.3f) -- the discrepancy is "
         "NOT the dipole approximation" % (dev5, dev1))
    print("  -> shrinking the sphere drives the ratio toward 1: the gap at "
          "a/L ~ 0.5\n     is genuinely the point-dipole approximation.  OK")

    print("\n" + "=" * 70)
    print("ALL CHECKS PASSED")
    print("=" * 70)
    print(f"\nRESULT: at a/L = {res5['a_over_L']:.2f} (a = 5 mm, z0 = "
          f"{Z0*1e3:.0f} mm, {F_HZ/1e3:.0f} kHz) the induced-point-dipole")
    print(f"levitation force is in error by {res5['discrepancy_percent']:+.1f}% "
          f"(full/dipole = {res5['ratio_full_over_dipole']:.3f}); the error")
    print(f"falls to {res1['discrepancy_percent']:+.1f}% at a/L = "
          f"{res1['a_over_L']:.2f} (a = 1 mm), confirming it is the dipole "
          f"approximation.")

    # ----- persist results JSON (committed location, next to script) -----
    results = {
        "description": ("Point-dipole approximation error in AC magnetic "
                        "levitation: FULL axisymmetric eddy-current Lorentz "
                        "force on a Cu sphere (team28 mixed phi-B solver) vs "
                        "the induced-dipole force, same FEM coil field."),
        "force_formula_full": ("Fz = Integrate((1/2)Re[Jt conj(B_r)]*2*pi*r "
                               "over Sphere), Jt=-sigma*s*phi/r (PEAK-amplitude "
                               "time average; positive = lift)"),
        "force_formula_dipole": "Fz = (pi a^3/mu0) Re[G(omega)] d|B0|^2/dz",
        "coil": {
            "R_m": R_COIL, "NI_amp_turns": NI, "cross_section_w_m": COIL_W,
            "cross_section_h_m": COIL_H, "centre_z_m": COIL_Z,
        },
        "frequency": {"f_Hz": F_HZ, "omega_rad_per_s": OMEGA,
                      "skin_depth_m": float(delta(OMEGA))},
        "sphere": {"sigma_S_per_m": SIGMA_CU, "z0_m": Z0},
        "fe": {"order_phi": P_ORDER, "order_B": P_ORDER - 1,
               "open_boundary": "anisotropic-nu infinite-element shell"},
        "case_5mm": res5,
        "case_1mm": res1,
        "key_check": {
            "ratio_5mm": res5["ratio_full_over_dipole"],
            "ratio_1mm": res1["ratio_full_over_dipole"],
            "abs_dev_5mm": dev5, "abs_dev_1mm": dev1,
            "ratio_moves_toward_one_as_aL_shrinks": bool(dev1 < dev5),
        },
    }
    out_json = os.path.join(HERE, "coil_sphere_eddy_force_results.json")
    with open(out_json, "w") as fp:
        json.dump(results, fp, indent=2)
    print(f"\nwrote {out_json}")

    try:
        plot(res5, res1)
    except Exception as e:               # plotting is optional, not a result
        print(f"(plot skipped: {e})")


def plot(res5, res1):
    """Force-ratio vs a/L bar figure.  NO in-figure title (lab figure policy:
    the title belongs in the LaTeX caption).  10 pt serif, ~8 cm."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.size": 10, "font.family": "serif", "axes.linewidth": 0.8,
        "mathtext.fontset": "dejavuserif",
    })
    fig, ax = plt.subplots(figsize=(8.0 / 2.54, 6.0 / 2.54))

    aL = [res1["a_over_L"], res5["a_over_L"]]
    ratio = [res1["ratio_full_over_dipole"], res5["ratio_full_over_dipole"]]
    labels = [f"a=1 mm\n(a/L={aL[0]:.2f})", f"a=5 mm\n(a/L={aL[1]:.2f})"]

    ax.axhline(1.0, color="0.6", ls=":", lw=0.8)
    ax.text(1.5, 1.0, " exact dipole", va="bottom", ha="right",
            fontsize=8, color="0.4")
    bars = ax.bar([0, 1], ratio, width=0.5, color=["C0", "C3"])
    for b, rt in zip(bars, ratio):
        ax.text(b.get_x() + b.get_width() / 2, rt,
                f"{rt:.3f}", ha="center",
                va="bottom" if rt >= 1 else "top", fontsize=8)
    ax.set_xticks([0, 1]); ax.set_xticklabels(labels)
    ax.set_ylabel(r"$F_z^{\mathrm{full}} / F_z^{\mathrm{dipole}}$")
    lo = min(ratio + [1.0]); hi = max(ratio + [1.0])
    pad = 0.04 * (hi - lo + 1e-9)
    ax.set_ylim(lo - 4 * pad, hi + 4 * pad)
    ax.tick_params(direction="in", right=True)
    fig.tight_layout(pad=0.3)

    out_png = os.path.join(HERE, "coil_sphere_eddy_force.png")
    fig.savefig(out_png, dpi=300)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
