"""Phase 2 -- equivalence-theorem reconstruction at 1 MHz (time-harmonic).

Hertzian-dipole source at origin, sample (E, H) on R = 1 m sphere,
reconstruct at external observation points via the time-harmonic
Stratton-Chu surface integral.

HISTORY (resolved 2026-05-26):
    EARLIER state -- this example was a "DEMONSTRATION + KNOWN
    LIMITATION" test that EXPECTEDLY failed the 2% threshold by ~3x
    in the deep near-field, because `evaluate()` used a SCALAR
    Stratton-Chu form (no `(1/k^2) grad-grad psi` dyadic correction).

    CURRENT state -- Phase B C++ kernel
    (`src/core/rad_equivalence_source.cpp::EvaluateHarmonic`) is now
    the default path for `evaluate(omega>0, use_cpp=True)`.  It uses
    the FULL dyadic Green's function `Ḡ_e = (I + ∇∇/k²) ψ`, which
    correctly reproduces the deep-near-field.  The static-limit test
    (ω = 1 Hz vs `evaluate_static_H`) agrees to 5.4e-15 relative
    diff -- machine precision -- and the harmonic dipole near-field
    is now reconstructed correctly.

    The Python scalar-form fallback is preserved behind `use_cpp=False`
    for regression cross-checking; it still has the historical
    undershoot in deep near-field and is documented as such in
    `NearFieldSource.evaluate()`.

Run:
    python phase2_wpt_harmonic.py

Output:
    results_phase2.json
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "src"))

from radia.equivalence_source import NearFieldSource, MU_0, EPS_0, C0


# ---------------------------------------------------------------------
# Hertzian dipole (point electric dipole) closed form
# ---------------------------------------------------------------------

def hertzian_dipole_EH(points, Il_phasor: complex, k: float):
    """Closed-form (E, H) of a Hertzian dipole at origin oriented along z.

    Reference: Balanis "Antenna Theory" eq. 4-10 / Stratton "EM Theory".
    Using exp(+jwt) convention so e^(-jkR) is the outgoing wave.

    Parameters:
        Il_phasor: I * l phasor (A * m).  For a Hertzian dipole the
                   "current moment" I*l is the parameter.
        k:         wave number = omega / c0
    """
    omega = k * C0
    eta = math.sqrt(MU_0 / EPS_0)  # ~ 376.73 ohm

    points = np.atleast_2d(np.asarray(points, dtype=np.float64))
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    r_safe = np.where(r > 1e-15, r, 1.0)
    theta = np.arccos(z / r_safe)
    phi = np.arctan2(y, x)

    kr = k * r_safe
    exp_jkr = np.exp(-1j * kr)
    # E_r, E_theta, H_phi  (Balanis 4-10)
    E_r = (eta * Il_phasor * np.cos(theta) / (2 * math.pi * r_safe ** 2)
           * (1.0 + 1.0 / (1j * kr)) * exp_jkr)
    E_theta = (1j * eta * k * Il_phasor * np.sin(theta) / (4 * math.pi * r_safe)
               * (1.0 + 1.0 / (1j * kr) - 1.0 / kr ** 2) * exp_jkr)
    H_phi = (1j * k * Il_phasor * np.sin(theta) / (4 * math.pi * r_safe)
             * (1.0 + 1.0 / (1j * kr)) * exp_jkr)

    # Convert spherical (E_r, E_theta, 0) and (0, 0, H_phi) to Cartesian
    sin_t, cos_t = np.sin(theta), np.cos(theta)
    sin_p, cos_p = np.sin(phi), np.cos(phi)
    # r-hat = (sin t cos p, sin t sin p, cos t)
    # theta-hat = (cos t cos p, cos t sin p, -sin t)
    # phi-hat = (-sin p, cos p, 0)
    Ex = E_r * sin_t * cos_p + E_theta * cos_t * cos_p
    Ey = E_r * sin_t * sin_p + E_theta * cos_t * sin_p
    Ez = E_r * cos_t - E_theta * sin_t
    Hx = H_phi * (-sin_p)
    Hy = H_phi * cos_p
    Hz = np.zeros_like(Hx)

    return (np.column_stack([Ex, Ey, Ez]),
            np.column_stack([Hx, Hy, Hz]))


# ---------------------------------------------------------------------
# Main verification
# ---------------------------------------------------------------------

def main():
    print("=" * 78)
    print("Phase 2: equivalence-theorem reconstruction (time-harmonic Hertzian dipole)")
    print(f"Date: {time.strftime('%Y-%m-%d %H:%M')}")
    print("=" * 78)

    # Source: small electric dipole at origin along z, oscillating
    # at 1 MHz, with current moment I*l = 1 A*m.
    f = 1.0e6
    omega = 2 * math.pi * f
    k = omega / C0
    wavelength = C0 / f
    Il_phasor = 1.0 + 0.0j  # I * l (A*m)
    print(f"Source: Hertzian dipole, f = {f/1e6} MHz, lambda = {wavelength:.1f} m")
    print(f"        I*l = {Il_phasor} A*m, k = {k:.4e} rad/m")

    # Extraction sphere: R = 1 m << lambda (deep near-field)
    R_sphere = 1.0
    n_theta = 40
    n_phi = 80
    centroids, normals, areas = NearFieldSource.spherical_surface_mesh(
        R_sphere, n_theta=n_theta, n_phi=n_phi)
    print(f"Extraction sphere: R = {R_sphere} m (R/lambda = {R_sphere/wavelength:.4f}), "
          f"{len(centroids)} faces")

    # Sample analytical (E, H) on the sphere
    E_surf, H_surf = hertzian_dipole_EH(centroids, Il_phasor, k)

    # Build NFS
    nfs = NearFieldSource.from_surface_mesh(
        centroids, normals, areas,
        E=E_surf.astype(np.complex128),
        H=H_surf.astype(np.complex128),
        omega=omega,
        source=f"hertzian_dipole_f={f}",
    )

    # Observation: 10 m on z-axis (typical EMC compliance) + a few extras
    obs = np.array([
        [0.0, 0.0, 10.0],
        [0.0, 0.0, -10.0],
        [10.0, 0.0, 0.0],   # equatorial, lambda/30
        [0.0, 0.0, 100.0],  # 100 m, lambda/3 (far-field)
    ])

    t0 = time.time()
    E_rec, H_rec = nfs.evaluate(obs, omega=omega)
    dt = time.time() - t0
    E_ana, H_ana = hertzian_dipole_EH(obs, Il_phasor, k)

    print(f"\nReconstruction (full Stratton-Chu, {dt*1000:.1f} ms):")
    print(f"{'obs (m)':>22}  {'|E_rec|':>11}  {'|E_ana|':>11}  "
          f"{'|H_rec|':>11}  {'|H_ana|':>11}  {'E_err':>7}  {'H_err':>7}")
    print("-" * 96)

    results = []
    max_E_err = 0.0
    max_H_err = 0.0
    max_H_abs_err_when_zero = 0.0
    h_zero_abs_threshold = 1e-12
    for p, er, hr, ea, ha in zip(obs, E_rec, H_rec, E_ana, H_ana):
        eE = float(np.linalg.norm(er - ea) / np.linalg.norm(ea))
        h_norm = float(np.linalg.norm(ha))
        h_abs = float(np.linalg.norm(hr - ha))
        if h_norm > h_zero_abs_threshold:
            eH = h_abs / h_norm
        else:
            eH = 0.0 if h_abs <= h_zero_abs_threshold else float("inf")
            max_H_abs_err_when_zero = max(max_H_abs_err_when_zero, h_abs)
        max_E_err = max(max_E_err, eE)
        max_H_err = max(max_H_err, eH)
        print(f"  ({p[0]:+6.1f},{p[1]:+6.1f},{p[2]:+6.1f})  "
              f"{np.linalg.norm(er):11.4e}  {np.linalg.norm(ea):11.4e}  "
              f"{np.linalg.norm(hr):11.4e}  {np.linalg.norm(ha):11.4e}  "
              f"{eE*100:6.2f}%  {eH*100:6.2f}%")
        results.append({
            "obs": p.tolist(),
            "E_rec_re": er.real.tolist(), "E_rec_im": er.imag.tolist(),
            "H_rec_re": hr.real.tolist(), "H_rec_im": hr.imag.tolist(),
            "E_ana_re": ea.real.tolist(), "E_ana_im": ea.imag.tolist(),
            "H_ana_re": ha.real.tolist(), "H_ana_im": ha.imag.tolist(),
            "E_rel_err": eE,
            "H_rel_err": eH,
            "H_abs_err": h_abs,
            "H_ana_norm": h_norm,
            "H_zero_abs_threshold": h_zero_abs_threshold,
        })

    print("-" * 96)
    accept = (max_E_err <= 0.02) and (max_H_err <= 0.02)
    status = "PASS" if accept else "FAIL"
    print(f"Max relative E error: {max_E_err*100:.2f}%, "
          f"H error: {max_H_err*100:.2f}%  --  {status}")
    print(f"Zero-H absolute-error max: {max_H_abs_err_when_zero:.3e} A/m "
          f"(threshold {h_zero_abs_threshold:.1e})")

    out = HERE / "results_phase2.json"
    out.write_text(json.dumps({
        "phase": 2,
        "description": "Equivalence-theorem time-harmonic reconstruction "
                       "vs Hertzian dipole analytical",
        "source": {
            "kind": "hertzian_dipole_z",
            "Il_phasor_re": Il_phasor.real, "Il_phasor_im": Il_phasor.imag,
            "frequency_Hz": f, "omega_rad_per_s": omega,
            "k_rad_per_m": k, "wavelength_m": wavelength,
        },
        "extraction": {
            "kind": "sphere", "R": R_sphere,
            "n_theta": n_theta, "n_phi": n_phi,
            "n_faces": len(centroids),
            "R_over_lambda": R_sphere / wavelength,
        },
        "reconstruction": {
            "method": "Stratton-Chu time-harmonic",
            "n_obs": len(obs),
            "max_E_rel_err": max_E_err,
            "max_H_rel_err": max_H_err,
            "max_H_abs_err_when_zero": max_H_abs_err_when_zero,
            "H_zero_abs_threshold": h_zero_abs_threshold,
            "t_recon_ms": dt * 1000,
            "accept_threshold": 0.02,
            "status": status,
        },
        "obs_results": results,
    }, indent=2))
    print(f"Wrote: {out}")
    return 0 if accept else 1


if __name__ == "__main__":
    sys.exit(main())
