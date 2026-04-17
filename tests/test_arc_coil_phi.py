"""test_arc_coil_phi.py — verify H = -grad(Phi) for rad.Fld('phi').

Tests the magnetic scalar potential Phi returned by rad.Fld(obj, 'phi')
against the H field from rad.Fld(obj, 'h') via finite differences:

    H_i = -d(Phi)/d(x_i)  (i = x, y, z)

Convention check: Radia may use H = +grad(Phi) or H = -grad(Phi).
The test determines the sign and checks consistency.

Test matrix:
    1. Full circle (gap=0): should PASS (uses CircularLoopSolidAngle)
    2. Arc coil (gap=5 deg): currently FAILS (buggy solid angle formula)
    3. Arc coil (gap=30 deg): currently FAILS
    4. Multiple sample points (on-axis, off-axis, above, below)
"""
import math
import numpy as np
import radia as rad

MU_0 = 4 * math.pi * 1e-7


def build_coil(gap_deg, R=0.08, a=0.008, z=0.10, I=1000.0, nsec=200):
    """Build an ObjArcCur coil with given gap angle."""
    rad.UtiDelAll()
    arc_deg = 360.0 - gap_deg
    J0 = I / (2 * a) ** 2
    return rad.ObjArcCur(
        [0, 0, z],
        [R - a, R + a],
        [0, math.radians(arc_deg)],
        2 * a, nsec, 'man', 'z', J0)


def check_grad_phi(coil, points, eps=1e-5, label=""):
    """Check H = -grad(Phi) at sample points via finite differences."""
    print(f"\n{'='*72}")
    print(f"  {label}")
    print(f"{'='*72}")
    print(f"  {'point':>28s}  {'|H_rad|':>10s}  {'|gradPhi|':>10s}  "
          f"{'cos_angle':>10s}  {'ratio':>8s}  {'status':>6s}")

    results = []
    for pt in points:
        pt = np.asarray(pt, dtype=float)
        H = np.array(rad.Fld(coil, 'h', pt.tolist()), dtype=float)

        # Finite-difference gradient of Phi
        grad_phi = np.zeros(3)
        for i in range(3):
            pt_p = pt.copy(); pt_p[i] += eps
            pt_m = pt.copy(); pt_m[i] -= eps
            phi_p = float(rad.Fld(coil, 'phi', pt_p.tolist()))
            phi_m = float(rad.Fld(coil, 'phi', pt_m.tolist()))
            grad_phi[i] = (phi_p - phi_m) / (2 * eps)

        mag_H = np.linalg.norm(H)
        mag_gP = np.linalg.norm(grad_phi)

        if mag_H < 1e-12 or mag_gP < 1e-12:
            print(f"  {str(tuple(pt)):>28s}  {mag_H:10.3e}  {mag_gP:10.3e}  "
                  f"{'N/A':>10s}  {'N/A':>8s}  {'SKIP':>6s}")
            continue

        # Check if grad(Phi) = +H or -H
        cos_angle = np.dot(H, grad_phi) / (mag_H * mag_gP)
        ratio = mag_gP / mag_H

        # Pass if |cos_angle| > 0.999 and ratio ~ 1.0
        ok = abs(abs(cos_angle) - 1.0) < 0.01 and abs(ratio - 1.0) < 0.05
        status = "OK" if ok else "FAIL"

        print(f"  {str(tuple(np.round(pt, 4))):>28s}  {mag_H:10.3e}  {mag_gP:10.3e}  "
              f"{cos_angle:+10.6f}  {ratio:8.4f}  {status:>6s}")
        results.append({
            'pt': pt, 'H': H, 'grad_phi': grad_phi,
            'cos': cos_angle, 'ratio': ratio, 'ok': ok
        })

    n_ok = sum(1 for r in results if r['ok'])
    n_total = len(results)
    sign = None
    if results:
        avg_cos = np.mean([r['cos'] for r in results])
        sign = "+grad" if avg_cos > 0 else "-grad"
    print(f"\n  Result: {n_ok}/{n_total} PASS  "
          f"(sign convention: H = {sign}(Phi))")
    return results


def main():
    # Sample points (coil at z=0.10, R=0.08)
    pts = [
        # On axis (below coil)
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.02),
        (0.0, 0.0, -0.05),
        # Off axis (in the equatorial plane of iron)
        (0.03, 0.0, 0.0),
        (0.0, 0.03, 0.0),
        (0.05, 0.0, 0.0),
        # Above coil
        (0.0, 0.0, 0.20),
        (0.05, 0.0, 0.15),
        # Far field
        (0.0, 0.0, -0.30),
        (0.20, 0.0, 0.0),
        # Near coil (but outside conductor)
        (0.0, 0.0, 0.08),
        (0.05, 0.0, 0.10),
    ]

    # Test 1: Full circle (gap=0)
    coil_full = build_coil(gap_deg=0.0)
    check_grad_phi(coil_full, pts, label="Full circle (gap=0 deg)")

    # Test 2: Small gap (5 deg)
    coil_5 = build_coil(gap_deg=5.0)
    check_grad_phi(coil_5, pts, label="Arc coil (gap=5 deg)")

    # Test 3: Large gap (30 deg)
    coil_30 = build_coil(gap_deg=30.0)
    check_grad_phi(coil_30, pts, label="Arc coil (gap=30 deg)")

    # Test 4: Half ring (180 deg gap)
    coil_180 = build_coil(gap_deg=180.0)
    check_grad_phi(coil_180, pts, label="Half ring (gap=180 deg)")


if __name__ == "__main__":
    main()
