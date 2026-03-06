"""Debug script to check Delta_Z sign from ShieldBEMSIBC."""

import sys, os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad
from fasthenry_parser import FastHenryParser

MU_0 = 4e-7 * np.pi


def check_delta_z_sign():
    """Check Delta_Z sign for a simple single-wire + shield case."""

    rad.UtiDelAll()
    rad.FldUnits('m')

    # Simple straight wire with close shield
    inp = """\
.Units mm
.default sigma=5.8e7

N1 x=-50 y=0 z=0
N2 x=50 y=0 z=0

E1 N1 N2 w=1 h=1

.external N1 N2

.shield
  type=box
  center=0,0,0
  size=120,60,20
  sigma=3.7e7
.endshield

.freq fmin=100e3 fmax=100e3 ndec=1
.end
"""
    parser = FastHenryParser()
    parser.parse_string(inp)

    # Build topology
    builder = parser.to_peec_builder()
    topo = builder.build_topology()

    # Build shield solver
    shield_solvers = parser.build_shield_solvers()
    shield = shield_solvers[0]

    freq = 100e3
    omega = 2 * np.pi * freq

    # Compute Delta_Z at 100 kHz
    Delta_Z = shield.compute_impedance_matrix(freq, topo)

    n_seg = Delta_Z.shape[0]
    print(f"\nDelta_Z shape: {Delta_Z.shape}")
    print(f"Frequency: {freq:.0f} Hz")
    print(f"omega: {omega:.2e}")

    # Print diagonal elements
    print(f"\nDiagonal elements of Delta_Z:")
    for i in range(n_seg):
        dz = Delta_Z[i, i]
        dL = dz.imag / omega
        dR = dz.real
        print(f"  Delta_Z[{i},{i}] = {dz.real:+.6e} + j*{dz.imag:+.6e}")
        print(f"    -> Delta_R = {dR*1e3:+.6f} mOhm, Delta_L = {dL*1e9:+.4f} nH")

    # Sum of diagonal: total self-coupling effect
    dz_sum = np.sum(np.diag(Delta_Z))
    print(f"\nSum of diagonal:")
    print(f"  Delta_R_total = {dz_sum.real*1e3:+.6f} mOhm")
    print(f"  Delta_L_total = {(dz_sum.imag/omega)*1e9:+.4f} nH")
    print()

    # Physical expectation:
    print("Expected (Lenz's law):")
    print("  Delta_R > 0 (shield absorbs energy)")
    print("  Delta_L < 0 (shield opposes flux)")
    print()

    if dz_sum.real > 0:
        print("  Delta_R: CORRECT (positive)")
    else:
        print("  Delta_R: WRONG (negative) -- need to negate Delta_Z")

    if dz_sum.imag < 0:
        print("  Delta_L: CORRECT (negative)")
    else:
        print("  Delta_L: WRONG (positive) -- need to negate Delta_Z")

    # Also check A_scat direction
    print("\n--- Checking A_scat direction for segment 0 ---")
    centers = np.asarray(topo['segment_centers'])
    directions = np.asarray(topo['segment_directions'])
    lengths = np.asarray(topo['segment_lengths'])

    from ngbem_eddy import _biot_savart_A

    # Compute A_inc from segment 0 at its own center
    A_inc_0 = _biot_savart_A(
        centers[0:1], centers[0], directions[0], lengths[0])
    print(f"  A_inc[0] at center[0]: {A_inc_0[0]}")
    print(f"  A_inc[0] . dir[0] = {np.dot(A_inc_0[0], directions[0]):.6e}")

    # Solve for shield current from segment 0
    def A_inc_func(pts):
        return _biot_savart_A(pts, centers[0], directions[0], lengths[0])

    shield.solve(freq, A_inc_func=A_inc_func)
    A_scat = shield.compute_A_scattered(centers[0:1])
    print(f"  A_scat[0] at center[0]: {A_scat[0]}")
    print(f"  A_scat[0] . dir[0] = {np.dot(A_scat[0], directions[0]):.6e}")
    print(f"  |A_scat / A_inc| . dir = "
          f"{abs(np.dot(A_scat[0], directions[0]) / np.dot(A_inc_0[0], directions[0])):.4f}")

    # Lenz check: A_scat . dir should be opposite to A_inc . dir
    a_inc_proj = np.dot(A_inc_0[0].real, directions[0])
    a_scat_proj = np.dot(A_scat[0].real, directions[0])
    print(f"\n  Real parts:")
    print(f"    Re(A_inc . dir)  = {a_inc_proj:.6e}")
    print(f"    Re(A_scat . dir) = {a_scat_proj:.6e}")
    if a_inc_proj * a_scat_proj < 0:
        print("    Lenz check: PASSED (opposite signs)")
    else:
        print("    Lenz check: FAILED (same signs!)")

    rad.UtiDelAll()


if __name__ == '__main__':
    check_delta_z_sign()
