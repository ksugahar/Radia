"""Null-field property check -- the defining test of the equivalence theorem.

Verifies the two-sided behaviour of a NearFieldSource constructed from
analytical source data on a closed surface:

  * EXTERIOR (obs outside the surface): reconstructs the TRUE source
    field, to within the surface-integral discretisation error.

  * INTERIOR (obs inside the surface, but NOT at the source itself):
    reconstructs ~zero (the "null-field" property).  This is what
    makes the equivalence theorem useful: the equivalent surface
    sources cancel the real source contribution everywhere inside,
    so the surface IS equivalent to the source as seen from outside.

Source: magnetic dipole at the origin, m = m_z * z-hat.
Extraction surface: sphere of radius R = 0.30 m around the dipole.

Acceptance:
  - Exterior obs (|r| >= 0.6 m): max rel err <= 1.5 %.
  - Interior obs (|r| <= 0.20 m, away from origin): |H_rec| / |H_true|
    <= 1 % (the null-field amplitude is dominated by the surface
    integral discretisation error, not by any real interior field).
  - The contrast between the two regimes is the WHOLE POINT.

Run:
    python null_field_property.py

Output:
    results_null_field.json
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

from radia.equivalence_source import NearFieldSource  # noqa: E402


# ---------------------------------------------------------------------
# Closed-form magnetic-dipole field
# ---------------------------------------------------------------------

def dipole_H(points: np.ndarray, m_vec: np.ndarray) -> np.ndarray:
    """H field of a point magnetic dipole at origin, in free space.

    H(r) = (1 / 4 pi) * [ 3 (m . r-hat) r-hat - m ] / r^3   [A/m]

    Args:
        points: (N, 3) observation points [m]
        m_vec:  (3,)   dipole moment [A m^2]
    Returns:
        H: (N, 3) [A/m]
    """
    r = np.linalg.norm(points, axis=1)
    rhat = points / r[:, None]
    m_dot_r = (rhat * m_vec[None, :]).sum(axis=1)
    return (3 * m_dot_r[:, None] * rhat - m_vec[None, :]) / (4 * math.pi * r[:, None] ** 3)


# ---------------------------------------------------------------------
# Build extraction surface (panel sphere) + sample dipole H
# ---------------------------------------------------------------------

def build_nfs(R_surface: float, n_theta: int, n_phi: int,
              m_vec: np.ndarray) -> NearFieldSource:
    centroids, normals, areas = NearFieldSource.spherical_surface_mesh(
        R=R_surface, n_theta=n_theta, n_phi=n_phi)
    H_surf = dipole_H(centroids, m_vec)
    return NearFieldSource.from_surface_mesh(
        centroids, normals, areas, E=None, H=H_surf, omega=0.0)


# ---------------------------------------------------------------------
# Main test: contrast exterior reconstruction vs interior null
# ---------------------------------------------------------------------

def main():
    R_surface = 0.30          # [m] extraction sphere radius
    n_theta = 40
    n_phi = 80
    m_vec = np.array([0.0, 0.0, 1.0])   # [A m^2]

    print(f"Building NearFieldSource on sphere R={R_surface} m, "
          f"{n_theta}x{n_phi}=2*{n_theta * n_phi} panels...")
    nfs = build_nfs(R_surface, n_theta, n_phi, m_vec)
    n_panels = nfs.n_faces
    print(f"  -> {n_panels} panels")

    # Exterior observation points (well outside the sphere)
    obs_ext = np.array([
        [0.60, 0.00, 0.00],
        [0.00, 0.60, 0.00],
        [0.00, 0.00, 0.60],
        [0.45, 0.45, 0.00],
        [0.60, 0.00, 0.40],
        [1.20, 0.00, 0.00],
        [0.00, 0.00, 1.50],
    ])
    # Interior observation points (well inside the sphere, NOT at the
    # singular dipole origin).  These are where the null-field property
    # is supposed to hold.
    obs_int = np.array([
        [0.10, 0.00, 0.00],
        [0.00, 0.15, 0.00],
        [0.00, 0.00, 0.20],
        [0.10, 0.05, 0.10],
        [0.18, 0.00, 0.00],
        [-0.12, 0.08, 0.05],
    ])

    print(f"\n=== EXTERIOR reconstruction (must match dipole field) ===")
    print(f"  {len(obs_ext)} obs points at |r| = "
          f"{[f'{r:.2f}' for r in np.linalg.norm(obs_ext, axis=1)]}")
    t0 = time.perf_counter()
    H_ext_rec = nfs.evaluate_static_H(obs_ext)
    t_ext = time.perf_counter() - t0
    H_ext_true = dipole_H(obs_ext, m_vec)
    rel_err_ext = np.linalg.norm(H_ext_rec - H_ext_true, axis=1) \
                  / np.linalg.norm(H_ext_true, axis=1)
    print(f"  reconstruction time: {t_ext * 1000:.1f} ms "
          f"(C++ kernel, {n_panels} panels x {len(obs_ext)} obs)")
    for i, p in enumerate(obs_ext):
        print(f"    r={tuple(round(x, 2) for x in p)} "
              f"|H_true|={np.linalg.norm(H_ext_true[i]):.4e}  "
              f"|H_rec|={np.linalg.norm(H_ext_rec[i]):.4e}  "
              f"rel_err={rel_err_ext[i]:.3%}")
    max_err_ext = float(rel_err_ext.max())
    print(f"  -> MAX exterior rel err = {max_err_ext:.3%}  "
          f"(acceptance: <= 1.5%)")

    print(f"\n=== INTERIOR null-field check (must reconstruct ~ 0) ===")
    print(f"  {len(obs_int)} obs points at |r| = "
          f"{[f'{r:.2f}' for r in np.linalg.norm(obs_int, axis=1)]}")
    H_int_rec = nfs.evaluate_static_H(obs_int)
    H_int_true = dipole_H(obs_int, m_vec)
    # Ratio: |H_reconstructed| / |H_true_dipole_at_interior|
    # Should be SMALL because the surface integral is supposed to cancel
    # the real interior field (null-field theorem).
    rel_int = np.linalg.norm(H_int_rec, axis=1) \
              / np.linalg.norm(H_int_true, axis=1)
    for i, p in enumerate(obs_int):
        print(f"    r={tuple(round(x, 2) for x in p)} "
              f"|H_true|={np.linalg.norm(H_int_true[i]):.4e}  "
              f"|H_rec|={np.linalg.norm(H_int_rec[i]):.4e}  "
              f"|H_rec|/|H_true|={rel_int[i]:.3%}")
    max_int = float(rel_int.max())
    print(f"  -> MAX interior |H_rec|/|H_true| = {max_int:.3%}  "
          f"(acceptance: <= 1.0%)")

    # ---- Pass / Fail ---------------------------------------------------
    ext_ok = max_err_ext <= 0.015
    int_ok = max_int <= 0.01
    overall_ok = ext_ok and int_ok

    print("\n" + "=" * 60)
    print(f"  EXTERIOR reconstruction: {'PASS' if ext_ok else 'FAIL'}  "
          f"({max_err_ext:.3%} <= 1.5%)")
    print(f"  INTERIOR null field:     {'PASS' if int_ok else 'FAIL'}  "
          f"({max_int:.3%} <= 1.0%)")
    print(f"  Overall:                 {'PASS' if overall_ok else 'FAIL'}")
    print("=" * 60)

    # ---- JSON report ---------------------------------------------------
    result = {
        "test": "null_field_property",
        "R_surface_m": R_surface,
        "n_panels": int(n_panels),
        "m_vec_Am2": m_vec.tolist(),
        "exterior": {
            "obs_points": obs_ext.tolist(),
            "H_reconstructed": H_ext_rec.tolist(),
            "H_true": H_ext_true.tolist(),
            "rel_err": rel_err_ext.tolist(),
            "max_rel_err": max_err_ext,
            "passed": ext_ok,
            "tolerance": 0.015,
        },
        "interior": {
            "obs_points": obs_int.tolist(),
            "H_reconstructed": H_int_rec.tolist(),
            "H_true_dipole": H_int_true.tolist(),
            "rel_ratio": rel_int.tolist(),
            "max_rel_ratio": max_int,
            "passed": int_ok,
            "tolerance": 0.01,
        },
        "overall_passed": overall_ok,
        "wall_time_ext_s": t_ext,
    }
    out = HERE / "results_null_field.json"
    out.write_text(json.dumps(result, indent=2))
    print(f"\nSaved: {out}")
    return 0 if overall_ok else 1


if __name__ == "__main__":
    sys.exit(main())
