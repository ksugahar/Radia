"""
evaluate_peec_formula_accuracy.py

Quantitative evaluation of current PEEC formula accuracy vs reference solutions.
Tests three areas:
  1. Self-inductance: 3-term Grover vs full Rosa/Grover (14+ terms)
  2. Mutual inductance: 8-pt Gauss vs analytical (parallel) and high-order quad (non-parallel)
  3. Panel potential: current implementation vs analytical (parallel plates)

Part of Radia project - FastMaxwell incorporation analysis
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

MU_0 = 4.0 * np.pi * 1e-7
EPS_0 = 8.854187817e-12
INV_FOUR_PI = 1.0 / (4.0 * np.pi)

# ==============================================================================
# 1. SELF-INDUCTANCE: 3-term vs Rosa/Grover full formula
# ==============================================================================

def self_inductance_3term(l, w, h):
    """Current Radia implementation: 3-term Grover formula."""
    d_rect = np.sqrt(w*w + h*h)
    if d_rect < 1e-15:
        d_rect = 1e-6
    if l > d_rect:
        term1 = np.log(2.0 * l / d_rect)
        term2 = 0.25
        term3 = (w*w + h*h) / (12.0 * l*l)
        return (MU_0 / (2.0 * np.pi)) * l * (term1 + term2 + term3)
    else:
        return (MU_0 / (2.0 * np.pi)) * l * 0.5


def self_inductance_rosa_grover(l, w, h):
    """Full Rosa/Grover formula with all correction terms.

    Reference: Rosa, NBS Scientific Paper 169 (1908)
               Grover, "Inductance Calculations", Dover, 1946
               Ruehli, "Inductance Calculations in a Complex Integrated Circuit
                        Environment", IBM J. Res. Dev., 16(5), 470-481, 1972

    Uses the exact formula for self-inductance of a rectangular conductor:
    L = (mu_0 / (2*pi)) * l * Y(l, w, h)

    where Y is the Newell formula (exact double volume integral of 1/R).
    """
    # Normalize by length
    wn = w / l  # normalized width
    tn = h / l  # normalized height (thickness)

    # Derived quantities
    aw = np.sqrt(wn*wn + 1.0)
    at = np.sqrt(tn*tn + 1.0)
    r = np.sqrt(wn*wn + tn*tn)
    ar = np.sqrt(wn*wn + tn*tn + 1.0)

    # Full Rosa/Grover formula (from FastMaxwell calcaoneoverr.h)
    z = 0.0

    # Base terms (1/4 coefficient)
    z += 0.25 * (
        (1.0/wn) * np.arcsinh(wn / at) +
        (1.0/tn) * np.arcsinh(tn / aw) +
        np.arcsinh(1.0 / r)
    )

    # Correction terms (1/24 coefficient) - 6 additional terms
    z += (1.0/24.0) * (
        (tn*tn/wn) * np.arcsinh(wn / (tn * at * (r + ar))) +
        (wn*wn/tn) * np.arcsinh(tn / (wn * aw * (r + ar))) +
        (tn*tn/(wn*wn)) * np.arcsinh(wn*wn / (tn * r * (at + ar))) +
        (wn*wn/(tn*tn)) * np.arcsinh(tn*tn / (wn * r * (aw + ar))) +
        (1.0/(wn * tn*tn)) * np.arcsinh(wn * tn*tn / (at * (aw + ar))) +
        (1.0/(tn * wn*wn)) * np.arcsinh(tn * wn*wn / (aw * (at + ar)))
    )

    # Arctangent corrections (1/6 coefficient)
    z -= (1.0/6.0) * (
        (1.0/(wn*tn)) * np.arctan(wn*tn / ar) +
        (tn/wn) * np.arctan(wn / (tn * ar)) +
        (wn/tn) * np.arctan(tn / (wn * ar))
    )

    # Rational corrections (1/60 coefficient)
    z -= (1.0/60.0) * (
        ((ar + r + tn + at) * tn*tn) /
            ((ar + r) * (r + tn) * (tn + at) * (at + ar)) +
        ((ar + r + wn + aw) * wn*wn) /
            ((ar + r) * (r + wn) * (wn + aw) * (aw + ar)) +
        (ar + aw + 1.0 + at) /
            ((ar + aw) * (aw + 1.0) * (1.0 + at) * (at + ar))
    )

    # Inverse distance corrections (1/20 coefficient)
    z -= (1.0/20.0) * (
        1.0/(r + ar) + 1.0/(aw + ar) + 1.0/(at + ar)
    )

    z *= (2.0 / np.pi)
    z *= l

    return z * MU_0


def evaluate_self_inductance():
    """Compare 3-term vs full Rosa/Grover for various aspect ratios."""
    print("=" * 70)
    print("1. SELF-INDUCTANCE: 3-term Grover vs Full Rosa/Grover")
    print("=" * 70)
    print()

    # Test parameters: 1mm x 1mm cross-section, varying length
    w = 1e-3  # 1 mm width
    h = 1e-3  # 1 mm height

    print(f"  Cross-section: {w*1e3:.1f} mm x {h*1e3:.1f} mm")
    print()
    print(f"  {'l/w':>6s} | {'l (mm)':>8s} | {'L_3term (nH)':>12s} | "
          f"{'L_Rosa (nH)':>12s} | {'Error (%)':>10s} | {'Impact':>10s}")
    print("  " + "-" * 68)

    l_over_w_values = [0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0]

    for l_over_w in l_over_w_values:
        l = l_over_w * w

        L_3term = self_inductance_3term(l, w, h)
        L_rosa = self_inductance_rosa_grover(l, w, h)

        error = abs(L_3term - L_rosa) / abs(L_rosa) * 100.0

        if error > 10:
            impact = "CRITICAL"
        elif error > 5:
            impact = "HIGH"
        elif error > 1:
            impact = "MODERATE"
        elif error > 0.1:
            impact = "LOW"
        else:
            impact = "NEGLIGIBLE"

        print(f"  {l_over_w:6.1f} | {l*1e3:8.2f} | {L_3term*1e9:12.4f} | "
              f"{L_rosa*1e9:12.4f} | {error:10.3f} | {impact:>10s}")

    # Also test non-square cross-sections (bus bar)
    print()
    print(f"  Non-square cross-section (bus bar): w=10mm, h=1mm")
    w_bus = 10e-3
    h_bus = 1e-3
    print(f"  {'l/w':>6s} | {'l (mm)':>8s} | {'L_3term (nH)':>12s} | "
          f"{'L_Rosa (nH)':>12s} | {'Error (%)':>10s} | {'Impact':>10s}")
    print("  " + "-" * 68)

    for l_over_w in [0.5, 1.0, 2.0, 5.0, 10.0, 20.0]:
        l = l_over_w * w_bus

        L_3term = self_inductance_3term(l, w_bus, h_bus)
        L_rosa = self_inductance_rosa_grover(l, w_bus, h_bus)

        error = abs(L_3term - L_rosa) / abs(L_rosa) * 100.0

        if error > 10:
            impact = "CRITICAL"
        elif error > 5:
            impact = "HIGH"
        elif error > 1:
            impact = "MODERATE"
        elif error > 0.1:
            impact = "LOW"
        else:
            impact = "NEGLIGIBLE"

        print(f"  {l_over_w:6.1f} | {l*1e3:8.2f} | {L_3term*1e9:12.4f} | "
              f"{L_rosa*1e9:12.4f} | {error:10.3f} | {impact:>10s}")

    print()


# ==============================================================================
# 2. MUTUAL INDUCTANCE: 8-pt Gauss vs Analytical (parallel)
# ==============================================================================

def F_neumann(x, d):
    """Neumann auxiliary function F(x, d) = x * arsinh(x/d) - sqrt(x^2 + d^2)"""
    x2d2 = x*x + d*d
    return x * np.log((x + np.sqrt(x2d2)) / d) - np.sqrt(x2d2)


def mutual_inductance_analytical(l_i, l_j, d_perp, p_axial):
    """Exact Rosa/Grover formula for parallel filaments.

    l_i, l_j: lengths
    d_perp: perpendicular separation
    p_axial: axial offset (center-to-center along direction)
    """
    b1 = l_i / 2.0
    a1 = -l_i / 2.0
    b2 = p_axial + l_j / 2.0
    a2 = p_axial - l_j / 2.0

    M = (MU_0 * INV_FOUR_PI) * (
        F_neumann(b1 - a2, d_perp) + F_neumann(a1 - b2, d_perp) -
        F_neumann(b1 - b2, d_perp) - F_neumann(a1 - a2, d_perp)
    )
    return M


def mutual_inductance_gauss(l_i, l_j, d_perp, p_axial, angle_deg, n_gauss=8):
    """N-point Gauss-Legendre quadrature for mutual inductance.

    Handles non-parallel filaments at given angle.
    """
    # Gauss-Legendre points and weights
    gp, gw = np.polynomial.legendre.leggauss(n_gauss)

    # Filament i along x-axis, centered at origin
    dir_i = np.array([1.0, 0.0, 0.0])

    # Filament j at angle, offset by d_perp in y, p_axial in x
    angle_rad = np.radians(angle_deg)
    dir_j = np.array([np.cos(angle_rad), np.sin(angle_rad), 0.0])
    center_j = np.array([p_axial, d_perp, 0.0])

    dot_ij = np.dot(dir_i, dir_j)
    if abs(dot_ij) < 1e-10:
        return 0.0

    total = 0.0
    for ki in range(n_gauss):
        ti = gp[ki] * (l_i / 2.0)
        pi_pt = np.array([ti, 0.0, 0.0])

        for kj in range(n_gauss):
            tj = gp[kj] * (l_j / 2.0)
            pj_pt = center_j + tj * dir_j

            dist = np.linalg.norm(pi_pt - pj_pt)
            if dist > 1e-15:
                total += gw[ki] * gw[kj] / dist

    M = (MU_0 * INV_FOUR_PI) * dot_ij * (l_i / 2.0) * (l_j / 2.0) * total
    return M


def evaluate_mutual_inductance():
    """Compare Gauss quadrature vs analytical for various configurations."""
    print("=" * 70)
    print("2. MUTUAL INDUCTANCE: 8-pt Gauss vs Analytical")
    print("=" * 70)
    print()

    l = 10e-3  # 10 mm segment
    print(f"  Parallel filaments, l = {l*1e3:.0f} mm, p_axial = 0")
    print()
    print(f"  {'d/l':>8s} | {'d (mm)':>8s} | {'M_exact (nH)':>12s} | "
          f"{'M_8pt (nH)':>12s} | {'Error (%)':>10s} | {'Impact':>10s}")
    print("  " + "-" * 68)

    d_over_l_values = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]

    for d_over_l in d_over_l_values:
        d = d_over_l * l

        M_exact = mutual_inductance_analytical(l, l, d, 0.0)
        M_gauss = mutual_inductance_gauss(l, l, d, 0.0, 0.0, n_gauss=8)

        if abs(M_exact) > 1e-30:
            error = abs(M_gauss - M_exact) / abs(M_exact) * 100.0
        else:
            error = 0.0

        if error > 10:
            impact = "CRITICAL"
        elif error > 5:
            impact = "HIGH"
        elif error > 1:
            impact = "MODERATE"
        elif error > 0.1:
            impact = "LOW"
        else:
            impact = "NEGLIGIBLE"

        print(f"  {d_over_l:8.3f} | {d*1e3:8.3f} | {M_exact*1e9:12.6f} | "
              f"{M_gauss*1e9:12.6f} | {error:10.4f} | {impact:>10s}")

    # Non-parallel filaments (45 degrees)
    print()
    print(f"  Non-parallel filaments (45 deg), l = {l*1e3:.0f} mm, p_axial = 0")
    print()
    print(f"  {'d/l':>8s} | {'d (mm)':>8s} | "
          f"{'M_4pt (nH)':>12s} | {'M_8pt (nH)':>12s} | {'M_16pt (nH)':>12s} | "
          f"{'Err 8vs16':>10s}")
    print("  " + "-" * 80)

    for d_over_l in [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0]:
        d = d_over_l * l

        M_4pt = mutual_inductance_gauss(l, l, d, 0.0, 45.0, n_gauss=4)
        M_8pt = mutual_inductance_gauss(l, l, d, 0.0, 45.0, n_gauss=8)
        M_16pt = mutual_inductance_gauss(l, l, d, 0.0, 45.0, n_gauss=16)

        if abs(M_16pt) > 1e-30:
            error_8vs16 = abs(M_8pt - M_16pt) / abs(M_16pt) * 100.0
        else:
            error_8vs16 = 0.0

        if error_8vs16 > 10:
            impact = "CRITICAL"
        elif error_8vs16 > 5:
            impact = "HIGH"
        elif error_8vs16 > 1:
            impact = "MODERATE"
        elif error_8vs16 > 0.1:
            impact = "LOW"
        else:
            impact = "NEGLIGIBLE"

        print(f"  {d_over_l:8.3f} | {d*1e3:8.3f} | {M_4pt*1e9:12.6f} | "
              f"{M_8pt*1e9:12.6f} | {M_16pt*1e9:12.6f} | {error_8vs16:10.4f}")

    print()


# ==============================================================================
# 3. PANEL POTENTIAL: Current implementation vs Analytical
# ==============================================================================

def panel_potential_centroid_approx(area_i, area_j, dist):
    """Current Radia far-field approximation: P = (A_i * A_j) / (4*pi*eps_0 * r)"""
    return (area_i * area_j) / (4.0 * np.pi * EPS_0 * dist)


def panel_potential_single_area(area_j, dist):
    """Standard far-field: P = A_j / (4*pi*eps_0 * r)"""
    return area_j / (4.0 * np.pi * EPS_0 * dist)


def panel_potential_gauss_3pt(vertices_i, vertices_j):
    """3-point Gauss quadrature on triangles (current Radia implementation)."""
    # Barycentric coordinates for 3-point rule
    xi = np.array([
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5]
    ])
    w = 1.0 / 6.0

    # Compute areas
    e1_i = vertices_i[1] - vertices_i[0]
    e2_i = vertices_i[2] - vertices_i[0]
    area_i = 0.5 * np.linalg.norm(np.cross(e1_i, e2_i))

    e1_j = vertices_j[1] - vertices_j[0]
    e2_j = vertices_j[2] - vertices_j[0]
    area_j = 0.5 * np.linalg.norm(np.cross(e1_j, e2_j))

    total = 0.0
    for qi in range(3):
        pt_i = xi[qi, 0]*vertices_i[0] + xi[qi, 1]*vertices_i[1] + xi[qi, 2]*vertices_i[2]
        for qj in range(3):
            pt_j = xi[qj, 0]*vertices_j[0] + xi[qj, 1]*vertices_j[1] + xi[qj, 2]*vertices_j[2]
            R = np.linalg.norm(pt_i - pt_j)
            if R > 1e-15:
                total += w * w / R

    return (area_i * area_j * total) / (4.0 * np.pi * EPS_0)


def panel_potential_gauss_high(vertices_i, vertices_j, n_pts=7):
    """High-order Gauss quadrature on triangles (reference solution).

    Uses 7-point rule (degree 5) for triangles.
    """
    # 7-point Gauss rule for triangles (Strang & Fix)
    if n_pts == 7:
        bary = np.array([
            [1/3, 1/3, 1/3],
            [0.059715871789770, 0.470142064105115, 0.470142064105115],
            [0.470142064105115, 0.059715871789770, 0.470142064105115],
            [0.470142064105115, 0.470142064105115, 0.059715871789770],
            [0.797426985353087, 0.101286507323456, 0.101286507323456],
            [0.101286507323456, 0.797426985353087, 0.101286507323456],
            [0.101286507323456, 0.101286507323456, 0.797426985353087],
        ])
        weights = np.array([
            0.225,
            0.132394152788506, 0.132394152788506, 0.132394152788506,
            0.125939180544827, 0.125939180544827, 0.125939180544827,
        ]) * 0.5  # multiply by 0.5 for triangle area factor
    else:
        # Fallback to uniform sampling
        bary = np.array([[1/3, 1/3, 1/3]])
        weights = np.array([0.5])

    e1_i = vertices_i[1] - vertices_i[0]
    e2_i = vertices_i[2] - vertices_i[0]
    area_i = 0.5 * np.linalg.norm(np.cross(e1_i, e2_i))

    e1_j = vertices_j[1] - vertices_j[0]
    e2_j = vertices_j[2] - vertices_j[0]
    area_j = 0.5 * np.linalg.norm(np.cross(e1_j, e2_j))

    total = 0.0
    for qi in range(len(weights)):
        pt_i = bary[qi, 0]*vertices_i[0] + bary[qi, 1]*vertices_i[1] + bary[qi, 2]*vertices_i[2]
        for qj in range(len(weights)):
            pt_j = bary[qj, 0]*vertices_j[0] + bary[qj, 1]*vertices_j[1] + bary[qj, 2]*vertices_j[2]
            R = np.linalg.norm(pt_i - pt_j)
            if R > 1e-15:
                total += weights[qi] * weights[qj] / R

    return (area_i * area_j * total) / (np.pi * EPS_0)  # 4*area_i*area_j factor from Jacobian


def parallel_plate_capacitance_analytical(area, distance):
    """Analytical capacitance of parallel plates: C = eps_0 * A / d."""
    return EPS_0 * area / distance


def evaluate_panel_potential():
    """Compare panel potential implementations."""
    print("=" * 70)
    print("3. PANEL POTENTIAL: Current vs High-Order Quadrature")
    print("=" * 70)
    print()

    # Two parallel square panels (10mm x 10mm each)
    side = 10e-3
    area = side * side

    # Create triangular panels (split square into 2 triangles)
    # Panel 1 at z=0
    v1_tri1 = np.array([
        [0, 0, 0],
        [side, 0, 0],
        [side, side, 0]
    ], dtype=float)

    print(f"  Parallel triangular panels, side = {side*1e3:.0f} mm")
    print(f"  Panel area = {area*1e4:.2f} cm^2 (half-square triangle)")
    print()
    print(f"  {'d/side':>8s} | {'d (mm)':>8s} | {'P_3pt':>14s} | "
          f"{'P_7pt':>14s} | {'Err 3vs7':>10s} | {'Impact':>10s}")
    print("  " + "-" * 72)

    d_over_side_values = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]

    for d_over_side in d_over_side_values:
        d = d_over_side * side

        # Panel 2 at z=d
        v2_tri1 = np.array([
            [0, 0, d],
            [side, 0, d],
            [side, side, d]
        ], dtype=float)

        P_3pt = panel_potential_gauss_3pt(v1_tri1, v2_tri1)
        P_7pt = panel_potential_gauss_high(v1_tri1, v2_tri1, n_pts=7)

        if abs(P_7pt) > 1e-30:
            error = abs(P_3pt - P_7pt) / abs(P_7pt) * 100.0
        else:
            error = 0.0

        if error > 10:
            impact = "CRITICAL"
        elif error > 5:
            impact = "HIGH"
        elif error > 1:
            impact = "MODERATE"
        elif error > 0.1:
            impact = "LOW"
        else:
            impact = "NEGLIGIBLE"

        print(f"  {d_over_side:8.3f} | {d*1e3:8.2f} | {P_3pt:14.6e} | "
              f"{P_7pt:14.6e} | {error:10.4f} | {impact:>10s}")

    # Also check: self-potential comparison
    print()
    print(f"  Self-potential: observation on same panel (singularity test)")
    print(f"  Note: 3-pt Gauss avoids singularity by skipping R<1e-10")
    print()

    # Compute self-potential with Gauss
    P_self_3pt = panel_potential_gauss_3pt(v1_tri1, v1_tri1)
    P_self_7pt = panel_potential_gauss_high(v1_tri1, v1_tri1, n_pts=7)

    print(f"  P_self (3pt):  {P_self_3pt:14.6e}")
    print(f"  P_self (7pt):  {P_self_7pt:14.6e}")
    if abs(P_self_7pt) > 1e-30:
        print(f"  Error 3 vs 7:  {abs(P_self_3pt - P_self_7pt) / abs(P_self_7pt) * 100.0:.2f}%")
    print(f"  NOTE: Both are inaccurate for self-potential (singularity not handled)")
    print(f"        -> Analytical edge integration (Hess-Smith/Wilton) is needed")

    print()


# ==============================================================================
# SUMMARY
# ==============================================================================

def print_summary():
    print("=" * 70)
    print("SUMMARY: Impact Assessment")
    print("=" * 70)
    print()
    print("  Priority | Improvement             | Error Range  | When Matters")
    print("  " + "-" * 66)
    print("  P2       | Self-L Rosa/Grover 14T   | 0-50%+       | l/w < 5 (bus bar)")
    print("  P1       | P matrix Hess-Smith      | 10-30%+      | Near-field panels")
    print("  P3       | fourfil() mutual-L       | 0-5%         | Close non-parallel")
    print()


if __name__ == '__main__':
    print()
    evaluate_self_inductance()
    evaluate_mutual_inductance()
    evaluate_panel_potential()
    print_summary()
