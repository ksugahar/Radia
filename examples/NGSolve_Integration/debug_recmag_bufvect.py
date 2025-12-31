#!/usr/bin/env python
"""
Debug script: Reproduce radTRecMag::B_comp BufVect calculation
to understand the non-zero A field on z-axis.
"""
import numpy as np

MU_0 = 4 * np.pi * 1e-7
dx, dy, dz = 0.02, 0.02, 0.03  # half dimensions

def TransAtans(a, b, PiMult):
    """Reproduce the TransAtans function from Radia."""
    # This function handles arctangent branch cuts
    PI = np.pi
    result = a + b
    PiMult[0] = 0.0

    # Handle special cases
    if abs(a) > 1e10 or abs(b) > 1e10:
        if a > 1e10:
            result = 1.0
            PiMult[0] = 0.5
        elif a < -1e10:
            result = 1.0
            PiMult[0] = -0.5
        if b > 1e10:
            result = 1.0 if result > 0 else -1.0
            PiMult[0] += 0.5
        elif b < -1e10:
            result = 1.0 if result > 0 else -1.0
            PiMult[0] -= 0.5

    # Standard case
    if abs(result) < 1e-15:
        result = 1e-15

    return result

def compute_recmag_bufvect(obs_point, center, half_dim):
    """
    Reproduce radTRecMag::B_comp BufVect calculation.
    """
    PI = np.pi
    P_min_CenPo = np.array(obs_point) - np.array(center)

    # Compute BfSt (with potential perturbation for singularity avoidance)
    x0 = -P_min_CenPo[0] - half_dim[0]
    x1 = -P_min_CenPo[0] + half_dim[0]
    y0 = -P_min_CenPo[1] - half_dim[1]
    y1 = -P_min_CenPo[1] + half_dim[1]
    z0 = -P_min_CenPo[2] - half_dim[2]
    z1 = -P_min_CenPo[2] + half_dim[2]

    print(f"x0={x0:.6f}, x1={x1:.6f}")
    print(f"y0={y0:.6f}, y1={y1:.6f}")
    print(f"z0={z0:.6f}, z1={z1:.6f}")

    # Check for zeros and apply perturbation (Radia uses AbsRandMagnitude)
    # For simplicity, we use a fixed small perturbation
    AbsRand = 1e-12  # Small perturbation
    if abs(x0) < 1e-15:
        x0 = AbsRand
        print(f"Applied perturbation to x0")
    if abs(x1) < 1e-15:
        x1 = AbsRand
        print(f"Applied perturbation to x1")
    if abs(y0) < 1e-15:
        y0 = AbsRand
        print(f"Applied perturbation to y0")
    if abs(y1) < 1e-15:
        y1 = AbsRand
        print(f"Applied perturbation to y1")
    if abs(z0) < 1e-15:
        z0 = AbsRand
        print(f"Applied perturbation to z0")
    if abs(z1) < 1e-15:
        z1 = AbsRand
        print(f"Applied perturbation to z1")

    x0e2, x1e2 = x0*x0, x1*x1
    y0e2, y1e2 = y0*y0, y1*y1
    z0e2, z1e2 = z0*z0, z1*z1

    # Distances to 8 corners
    D000 = np.sqrt(x0e2 + y0e2 + z0e2)
    D100 = np.sqrt(x1e2 + y0e2 + z0e2)
    D010 = np.sqrt(x0e2 + y1e2 + z0e2)
    D110 = np.sqrt(x1e2 + y1e2 + z0e2)
    D001 = np.sqrt(x0e2 + y0e2 + z1e2)
    D101 = np.sqrt(x1e2 + y0e2 + z1e2)
    D011 = np.sqrt(x0e2 + y1e2 + z1e2)
    D111 = np.sqrt(x1e2 + y1e2 + z1e2)

    # Solid angle terms (T0 and T1)
    # These use the TransAtans function for branch cut handling
    # Simplified version - just using atan directly
    def safe_atan(y, x):
        if abs(x) < 1e-15:
            return np.sign(y) * np.pi / 2 if abs(y) > 1e-15 else 0.0
        return np.arctan(y / x)

    # T0.x computation (from radTRecMag::B_comp)
    T0_x = safe_atan(y0*z0/(x0*D000), 1.0) - safe_atan(y0*z1/(x0*D001), 1.0) \
         - safe_atan(y1*z0/(x0*D010), 1.0) + safe_atan(y1*z1/(x0*D011), 1.0)
    T1_x = -safe_atan(y0*z0/(x1*D100), 1.0) + safe_atan(y0*z1/(x1*D101), 1.0) \
         + safe_atan(y1*z0/(x1*D110), 1.0) - safe_atan(y1*z1/(x1*D111), 1.0)

    T0_y = safe_atan(x0*z0/(y0*D000), 1.0) - safe_atan(x0*z1/(y0*D001), 1.0) \
         - safe_atan(x1*z0/(y0*D100), 1.0) + safe_atan(x1*z1/(y0*D101), 1.0)
    T1_y = -safe_atan(x0*z0/(y1*D010), 1.0) + safe_atan(x0*z1/(y1*D011), 1.0) \
         + safe_atan(x1*z0/(y1*D110), 1.0) - safe_atan(x1*z1/(y1*D111), 1.0)

    T0_z = safe_atan(x0*y0/(z0*D000), 1.0) - safe_atan(x1*y0/(z0*D100), 1.0) \
         - safe_atan(x0*y1/(z0*D010), 1.0) + safe_atan(x1*y1/(z0*D110), 1.0)
    T1_z = -safe_atan(x0*y0/(z1*D001), 1.0) + safe_atan(x1*y0/(z1*D101), 1.0) \
         + safe_atan(x0*y1/(z1*D011), 1.0) - safe_atan(x1*y1/(z1*D111), 1.0)

    print(f"\nT0 = [{T0_x:.6e}, {T0_y:.6e}, {T0_z:.6e}]")
    print(f"T1 = [{T1_x:.6e}, {T1_y:.6e}, {T1_z:.6e}]")

    # Log terms (with singularity handling)
    def safe_log_ratio(num, denom):
        if abs(denom) < 1e-15:
            denom = 1e-15
        return np.log(num / denom)

    # z + D combinations
    z0plD100 = z0 + D100
    z1plD101 = z1 + D101
    z1plD001 = z1 + D001
    z0plD000 = z0 + D000
    z0plD010 = z0 + D010
    z1plD011 = z1 + D011
    z1plD111 = z1 + D111
    z0plD110 = z0 + D110

    y0plD100 = y0 + D100
    y1plD110 = y1 + D110
    y1plD010 = y1 + D010
    y0plD000 = y0 + D000
    y0plD001 = y0 + D001
    y1plD011 = y1 + D011
    y1plD111 = y1 + D111
    y0plD101 = y0 + D101

    x0plD010 = x0 + D010
    x1plD110 = x1 + D110
    x1plD100 = x1 + D100
    x0plD000 = x0 + D000
    x0plD001 = x0 + D001
    x1plD101 = x1 + D101
    x1plD111 = x1 + D111
    x0plD011 = x0 + D011

    ln_z0plD100_di_z1plD101 = safe_log_ratio(z0plD100, z1plD101)
    ln_z1plD001_di_z0plD000 = safe_log_ratio(z1plD001, z0plD000)
    ln_z0plD010_di_z1plD011 = safe_log_ratio(z0plD010, z1plD011)
    ln_z1plD111_di_z0plD110 = safe_log_ratio(z1plD111, z0plD110)
    ln_y0plD100_di_y1plD110 = safe_log_ratio(y0plD100, y1plD110)
    ln_y1plD010_di_y0plD000 = safe_log_ratio(y1plD010, y0plD000)
    ln_y0plD001_di_y1plD011 = safe_log_ratio(y0plD001, y1plD011)
    ln_y1plD111_di_y0plD101 = safe_log_ratio(y1plD111, y0plD101)
    ln_x0plD010_di_x1plD110 = safe_log_ratio(x0plD010, x1plD110)
    ln_x1plD100_di_x0plD000 = safe_log_ratio(x1plD100, x0plD000)
    ln_x0plD001_di_x1plD101 = safe_log_ratio(x0plD001, x1plD101)
    ln_x1plD111_di_x0plD011 = safe_log_ratio(x1plD111, x0plD011)

    # BufVect computation
    BufVect_x = (x0*T0_x + x1*T1_x
                + y0*(ln_z0plD100_di_z1plD101 + ln_z1plD001_di_z0plD000)
                + y1*(ln_z0plD010_di_z1plD011 + ln_z1plD111_di_z0plD110)
                + z0*(ln_y0plD100_di_y1plD110 + ln_y1plD010_di_y0plD000)
                + z1*(ln_y0plD001_di_y1plD011 + ln_y1plD111_di_y0plD101))

    BufVect_y = (y0*T0_y + y1*T1_y
                + x0*(ln_z0plD010_di_z1plD011 + ln_z1plD001_di_z0plD000)
                + x1*(ln_z0plD100_di_z1plD101 + ln_z1plD111_di_z0plD110)
                + z0*(ln_x0plD010_di_x1plD110 + ln_x1plD100_di_x0plD000)
                + z1*(ln_x0plD001_di_x1plD101 + ln_x1plD111_di_x0plD011))

    BufVect_z = (z0*T0_z + z1*T1_z
                + y0*(ln_x0plD001_di_x1plD101 + ln_x1plD100_di_x0plD000)
                + y1*(ln_x0plD010_di_x1plD110 + ln_x1plD111_di_x0plD011)
                + x0*(ln_y0plD001_di_y1plD011 + ln_y1plD010_di_y0plD000)
                + x1*(ln_y0plD100_di_y1plD110 + ln_y1plD111_di_y0plD101))

    print(f"\nBufVect = [{BufVect_x:.6e}, {BufVect_y:.6e}, {BufVect_z:.6e}]")

    return np.array([BufVect_x, BufVect_y, BufVect_z])

# Test
print("=" * 70)
print("Debug: radTRecMag BufVect calculation")
print("=" * 70)

Br = 1.2  # T
Mr = Br / MU_0  # A/m
M = np.array([0, 0, Mr])

center = [0, 0, 0]
half_dim = [dx, dy, dz]

test_points = [
    [0.05, 0.0, 0.0],   # Off-axis
    [0.0, 0.0, 0.05],   # On z-axis
]

for pt in test_points:
    print(f"\n{'='*70}")
    print(f"Observation point: {pt}")
    print(f"{'='*70}")

    BufVect = compute_recmag_bufvect(pt, center, half_dim)

    # A = (1/4pi) * M x BufVect
    INV_FOUR_PI = 1.0 / (4.0 * np.pi)
    MxBuf = np.cross(M, BufVect)
    A = INV_FOUR_PI * MxBuf

    print(f"\nM = [{M[0]:.0f}, {M[1]:.0f}, {M[2]:.0f}]")
    print(f"M x BufVect = [{MxBuf[0]:.6e}, {MxBuf[1]:.6e}, {MxBuf[2]:.6e}]")
    print(f"A = (1/4pi) * (M x BufVect) = [{A[0]:.6e}, {A[1]:.6e}, {A[2]:.6e}]")
    print(f"|A| = {np.linalg.norm(A):.6e}")
