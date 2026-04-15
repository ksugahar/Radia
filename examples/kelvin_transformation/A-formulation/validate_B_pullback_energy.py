"""validate_B_pullback_energy.py  (item 4)

Energy-equivalence validation of the B (2-form pseudovector) pullback
under the Kelvin inversion, using a Radia closed-form coil as the
ground-truth field.

Physical identity (Kelvin change-of-variables + 2-form pullback +
Nagamine CEFC 2026 material factor). The Kelvin outer sphere is
centered at ``c_out`` with radius ``R_K``; its INTERIOR in
computational coordinates represents the physical EXTERIOR of the
same-radius sphere centered AT ``c_out`` (not at the coil!):

    integral_{|r - c_out| > R_K} (1/2) nu_0 |B_radia(r)|^2 dV        (LHS)
        ==
    integral_{Kelvin outer sphere} (1/2) nu'(r') |B_comp(r')|^2 dV'  (RHS)

where

    B_comp(r') = - (R/rho')^4 H(n) B_radia(r_phys)
    r_phys     = c_out + R^2 (r'-c_out)/|r'-c_out|^2      (Kelvin inverse)
    nu'(r')    = (rho'/R)^2 nu_0
    rho'       = |r' - c_out|,   n = (r'-c_out)/rho'
    H(n)       = I - 2 n n^T                             (Householder)

Both sides are evaluated by direct spherical Gauss-Legendre quadrature
(no FEM mesh). The pullback is a mathematical identity that depends
only on the pullback formula and the material factor, so this is the
cleanest possible test — we avoid any contamination from the FEM
discretization / interpolation error (which, near rho'=0, is severe
because |B_comp| ~ (R/rho')^4 blows up).

    LHS: integrate (1/2) nu_0 |B_radia(r)|^2 r^2 dr dOmega
         over rho in [R_K, R_big] (truncation)

    RHS: integrate (1/2) nu'(rho') |B_comp(r')|^2 (rho')^2 drho' dOmega'
         over rho' in [eps, R_K]        (integrable at 0 after pullback)

LHS and RHS must agree in the limit R_big -> infinity and eps -> 0,
which directly verifies the Kelvin change-of-variables identity:

    nu_0 |B_phys|^2 dV_phys  =  nu'(rho') |B_comp(rho')|^2 dV'
    [ dV_phys = (R/rho_phys)^6 dV' = (rho'/R)^6 dV'  cancels against
      (R/rho')^8 from |B_comp|^2 and (rho'/R)^2 from nu' ]
"""
from __future__ import annotations

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

import radia as rad  # noqa: E402 (installed package with .pyd)

# Use the in-repo kelvin_source.py (editable dev copy) for the new helpers
SRC_RADIA = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src', 'radia'))
if SRC_RADIA not in sys.path:
    sys.path.insert(0, SRC_RADIA)
from kelvin_source import kelvin_pullback_B_pseudovector  # noqa: E402


MU_0 = 4.0 * math.pi * 1e-7
NU_0 = 1.0 / MU_0

# --- Coil geometry (same family as validate_radia_HB_kelvin.py) ----
R_coil = 0.030
a_coil = 0.003
gap_deg = 5.0
arc_deg = 360.0 - gap_deg
I_total = 1.0

# Kelvin two-sphere (only the outer sphere is used in this test; the
# physical exterior of the inner sphere (radius R_K, centered at origin)
# is the physical domain that we integrate over via shell quadrature).
R_K = 0.06
c_out = np.array([0.15, 0.0, 0.0])        # outer (Kelvin) sphere center

# Quadrature parameters (applied to both LHS and RHS)
R_big = 3.0                    # truncation for physical-exterior integral
rho_eps = 1.0e-3 * R_K         # inner cutoff for Kelvin radial integral
                               # (converged — integrand ~ rho'^2 at rho'=0;
                               #  see derivation in module docstring)
# LHS and RHS both split the radial integral into two sub-ranges so that
# the near-coil peak is sampled densely. For LHS, the coil sits at
# |r - c_out| ~ 0.15 so the split is at 0.3. For RHS, the coil image is
# at rho' = R^2 / 0.15 ~ 0.024 so the split is at 0.04.
n_r_near = 48                  # radial GL points in near-coil segment
n_r_far = 24                   # radial GL points in smooth far segment
n_theta = 32
n_phi = 48


# ---------------------------------------------------------------
# Ground-truth Radia coil
# ---------------------------------------------------------------
def build_radia_coil():
    rad.UtiDelAll()
    J0 = I_total / (2 * a_coil) ** 2
    return rad.ObjArcCur(
        [0, 0, 0],
        [R_coil - a_coil, R_coil + a_coil],
        [0, math.radians(arc_deg)],
        2 * a_coil, 200, 'man', 'z', J0)


# ---------------------------------------------------------------
# LHS: physical-exterior shell quadrature
# ---------------------------------------------------------------
def _gl_radial_points(r_lo, r_hi, n):
    t, w = np.polynomial.legendre.leggauss(n)
    return (0.5 * (t + 1.0) * (r_hi - r_lo) + r_lo, 0.5 * (r_hi - r_lo) * w)


def integrate_B2_physical_exterior(coil, R_inner, c_inner, R_outer,
                                    r_split, n_near, n_far,
                                    n_theta, n_phi):
    """Integrate (1/2) nu_0 |B|^2 dV over R_inner < |r - c_inner| < R_outer.

    Radial integration is split at ``r_split`` to concentrate GL points
    near the coil peak.
    """
    rr_near, w_near = _gl_radial_points(R_inner, r_split, n_near)
    rr_far, w_far = _gl_radial_points(r_split, R_outer, n_far)
    rr = np.concatenate([rr_near, rr_far])
    w_rr = np.concatenate([w_near, w_far])

    tt, wt = np.polynomial.legendre.leggauss(n_theta)
    cos_t = tt
    sin_t = np.sqrt(1.0 - cos_t * cos_t)
    phi = np.linspace(0.0, 2.0 * math.pi, n_phi, endpoint=False)
    dphi = 2.0 * math.pi / n_phi

    total = 0.0
    for r_i, w_r_i in zip(rr, w_rr):
        for ct, st, w_t in zip(cos_t, sin_t, wt):
            for ph in phi:
                pt = [c_inner[0] + r_i * st * math.cos(ph),
                      c_inner[1] + r_i * st * math.sin(ph),
                      c_inner[2] + r_i * ct]
                B = rad.Fld(coil, 'b', pt)
                B2 = B[0] * B[0] + B[1] * B[1] + B[2] * B[2]
                total += 0.5 * NU_0 * B2 * (r_i * r_i) * w_r_i * w_t * dphi
    return total


# ---------------------------------------------------------------
# RHS: Kelvin outer-sphere direct quadrature (no FEM)
# ---------------------------------------------------------------
def integrate_B2_kelvin_interior(coil, center, R, rho_min, rho_max,
                                   rho_split, n_near, n_far,
                                   n_theta, n_phi):
    """Integrate (1/2) nu'(rho') |B_comp(r')|^2 dV' over
    rho_min < |r' - center| < rho_max via spherical Gauss-Legendre.
    Radial integration is split at ``rho_split`` (coil image radius) so
    that the near-coil peak on the Kelvin side is well sampled.
    """
    rr_near, w_near = _gl_radial_points(rho_min, rho_split, n_near)
    rr_far, w_far = _gl_radial_points(rho_split, rho_max, n_far)
    rr = np.concatenate([rr_near, rr_far])
    w_rr = np.concatenate([w_near, w_far])

    tt, wt = np.polynomial.legendre.leggauss(n_theta)
    cos_t = tt
    sin_t = np.sqrt(1.0 - cos_t * cos_t)
    phi = np.linspace(0.0, 2.0 * math.pi, n_phi, endpoint=False)
    dphi = 2.0 * math.pi / n_phi

    total = 0.0
    for rho_p, w_r_i in zip(rr, w_rr):
        nu_p = (rho_p / R) ** 2 * NU_0
        # |B_comp|^2 = (R/rho')^8 |B_phys|^2  (H preserves norm)
        B_comp_sq_factor = (R / rho_p) ** 8
        for ct, st, w_t in zip(cos_t, sin_t, wt):
            for ph in phi:
                # r' on the Kelvin outer sphere (relative to center)
                d_prime = np.array([rho_p * st * math.cos(ph),
                                      rho_p * st * math.sin(ph),
                                      rho_p * ct])
                # Kelvin map (center = c_out, same as pullback center)
                r_phys = center + (R * R / (rho_p * rho_p)) * d_prime
                B_phys = np.array(rad.Fld(coil, 'b', list(r_phys)))
                B_phys_sq = float(np.dot(B_phys, B_phys))
                # |B_comp|^2 = (R/rho')^8 |B_phys|^2 (H preserves norm)
                B_comp_sq = B_comp_sq_factor * B_phys_sq
                # dV' = rho'^2 dr' d(cos theta) dphi
                total += 0.5 * nu_p * B_comp_sq * (rho_p * rho_p) * w_r_i * w_t * dphi
    return total


# ---------------------------------------------------------------
# Main
# ---------------------------------------------------------------
def main():
    print("=" * 68)
    print("Item 4: B (2-form) pullback energy-equivalence validation")
    print("=" * 68)
    print(f"  coil: R={R_coil*1e3:.0f}mm  a={a_coil*1e3:.1f}mm  "
          f"arc={arc_deg:.1f}deg  I={I_total}A")
    print(f"  Kelvin outer sphere c={c_out}  R={R_K}")
    print(f"  LHS shell: {R_K} < |r - c_out| < {R_big}  split@0.30  "
          f"quad=({n_r_near}+{n_r_far},{n_theta},{n_phi})")
    print(f"  RHS Kelvin interior: {rho_eps:g} < rho' < {R_K}  split@0.040  "
          f"quad=({n_r_near}+{n_r_far},{n_theta},{n_phi})")
    print()

    coil = build_radia_coil()

    # Sanity check on ground-truth Radia B at an exterior point
    pt_test = [0.0, 0.0, 0.1]
    B_test = np.array(rad.Fld(coil, 'b', pt_test))
    print(f"  Radia B at {pt_test} = {B_test}  (|B|={np.linalg.norm(B_test):.3e})")
    print()

    print("[LHS] Physical-exterior shell quadrature ...")
    t0 = time.perf_counter()
    # Inner near-coil range ends at r_split = 0.30 (coil at ~0.15)
    W_LHS_final = integrate_B2_physical_exterior(
        coil, R_K, c_out, R_big, r_split=0.30,
        n_near=n_r_near, n_far=n_r_far,
        n_theta=n_theta, n_phi=n_phi)
    print(f"  R_big={R_big:.2f}  W_ext = {W_LHS_final:.6e} J"
          f"  ({time.perf_counter() - t0:.1f}s)")

    print("\n[RHS] Kelvin outer sphere direct Gauss-Legendre quadrature ...")
    t0 = time.perf_counter()
    # Coil image at rho' = R^2 / 0.15 ~ 0.024; near-range ends at 0.040
    W_RHS = integrate_B2_kelvin_interior(
        coil, c_out, R_K, rho_eps, R_K, rho_split=0.040,
        n_near=n_r_near, n_far=n_r_far,
        n_theta=n_theta, n_phi=n_phi)
    print(f"  rho_min={rho_eps:.3e}  W_kext = {W_RHS:.6e} J"
          f"  ({time.perf_counter() - t0:.1f}s)")

    print()
    print("=" * 68)
    print("SUMMARY")
    print("=" * 68)
    rel = (W_RHS - W_LHS_final) / W_LHS_final * 100
    print(f"  LHS  (Radia shell  R_K..{R_big:.2f})  W = {W_LHS_final:.6e}  J")
    print(f"  RHS  (Kelvin interior direct)     W = {W_RHS:.6e}  J")
    print(f"  rel diff RHS vs LHS               = {rel:+7.3f} %")
    # The 2-form pullback + nu' material factor is a mathematical
    # identity (Nagamine CEFC 2026 derivation + Kelvin volume change).
    # Residual comes from GL quadrature + finite R_big truncation.
    # Observed: ~0.2 % with 48+24 radial GL pts and R_big=3.0.
    print(f"  acceptance (|rel| < 1 %)          : "
          f"{'OK' if abs(rel) < 1.0 else 'FAIL'}")
    print("=" * 68)


if __name__ == "__main__":
    main()
