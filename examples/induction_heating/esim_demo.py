"""
ESIM Induction Heating Demo

This script demonstrates the ESIM (Effective Surface Impedance Method) for
analyzing induction heating of a steel workpiece with nonlinear magnetic material.

The demo:
1. Creates an ESI table from a steel BH-curve
2. Sets up a block workpiece with ESIM
3. Simulates a simple field distribution and computes power losses
4. Compares linear vs nonlinear impedance behavior

Reference:
    K. Hollaus, M. Kaltenbacher, J. Schoberl, "A Nonlinear Effective Surface
    Impedance in a Magnetic Scalar Potential Formulation," IEEE Trans. Magnetics,
    2025, DOI: 10.1109/TMAG.2025.3613932

Author: Radia Development Team
Date: 2026-01-08
"""

import sys
import os
import numpy as np

# Add radia module path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

from esim_cell_problem import (
    ESIMCellProblemSolver,
    BHCurveInterpolator,
    generate_esi_table_from_bh_curve
)
from esim_workpiece import create_esim_block, create_esim_cylinder


def print_header(title):
    """Print a formatted header."""
    print()
    print("=" * 70)
    print(title)
    print("=" * 70)


def demo_bh_curve():
    """Demonstrate BH curve interpolation."""
    print_header("1. BH Curve Interpolation")

    # Steel BH curve (typical electrical steel at operating temperature)
    bh_curve = [
        [0, 0],
        [100, 0.2],
        [250, 0.5],
        [500, 0.9],
        [1000, 1.3],
        [2500, 1.6],
        [5000, 1.8],
        [10000, 1.95],
        [50000, 2.1],
    ]

    interp = BHCurveInterpolator(bh_curve)

    print("Steel BH Curve Data:")
    print("-" * 50)
    print(f"{'H [A/m]':>12} {'B [T]':>10} {'mu_r':>10} {'Saturation':>12}")
    print("-" * 50)

    H_values = [0, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 50000]
    for H in H_values:
        B = interp.B(H)
        mu_r = interp.mu_r(H)
        sat = B / 2.1 * 100 if B > 0 else 0  # Saturation percentage
        print(f"{H:>12.0f} {B:>10.3f} {mu_r:>10.0f} {sat:>11.1f}%")

    return bh_curve


def demo_cell_problem(bh_curve):
    """Demonstrate Cell Problem solver."""
    print_header("2. Cell Problem Solver")

    sigma = 2e6  # S/m (hot steel)
    freq = 50000  # 50 kHz

    print(f"Material: Steel")
    print(f"Conductivity: {sigma/1e6:.1f} MS/m")
    print(f"Frequency: {freq/1000:.0f} kHz")
    print()

    solver = ESIMCellProblemSolver(bh_curve, sigma, freq)

    print(f"Initial skin depth: {solver.delta_initial*1e3:.3f} mm")
    print(f"Domain length: {solver.L*1e3:.2f} mm (10 skin depths)")
    print()

    # Solve for various H0 values
    print("Cell Problem Solutions:")
    print("-" * 75)
    print(f"{'H0 [A/m]':>10} {'Re(Z) [mOhm]':>14} {'Im(Z) [mOhm]':>14} "
          f"{'P [kW/m^2]':>12} {'Iter':>6} {'mu_r':>8}")
    print("-" * 75)

    H0_values = [10, 50, 100, 500, 1000, 2000, 5000, 10000, 20000]
    results = []

    for H0 in H0_values:
        result = solver.solve(H0)
        Z = result['Z']
        P = result['P_prime']
        iters = result['iterations']
        mu_r = result['mu_final'] / (4 * np.pi * 1e-7)
        conv = "Y" if result['converged'] else "N"

        results.append({
            'H0': H0,
            'Z': Z,
            'P_prime': P,
            'mu_r': mu_r
        })

        print(f"{H0:>10.0f} {Z.real*1e3:>14.4f} {Z.imag*1e3:>14.4f} "
              f"{P/1e3:>12.2f} {iters:>5}({conv}) {mu_r:>8.0f}")

    # Compare with linear SIBC
    print()
    Z_linear = solver.get_linear_sibc()
    print(f"Linear SIBC (initial mu_r): Z = {Z_linear.real*1e3:.4f} + j{Z_linear.imag*1e3:.4f} mOhm")

    return solver, results


def demo_esi_table(bh_curve):
    """Demonstrate ESI table generation."""
    print_header("3. ESI Table Generation")

    sigma = 2e6  # S/m
    freq = 50000  # 50 kHz

    print("Generating ESI table (30 points, log-spaced from 1 to 100,000 A/m)...")
    esi_table = generate_esi_table_from_bh_curve(bh_curve, sigma, freq, n_points=30)

    print("ESI table generated successfully!")
    print()

    # Test interpolation at intermediate points
    print("ESI Table Interpolation Test:")
    print("-" * 60)
    print(f"{'H0 [A/m]':>12} {'Re(Z) [mOhm]':>14} {'Im(Z) [mOhm]':>14} {'|Z| [mOhm]':>12}")
    print("-" * 60)

    test_H0 = [25, 75, 300, 750, 3000, 7500, 15000]
    for H0 in test_H0:
        Z = esi_table.get_impedance(H0)
        print(f"{H0:>12.0f} {Z.real*1e3:>14.4f} {Z.imag*1e3:>14.4f} {abs(Z)*1e3:>12.4f}")

    return esi_table


def demo_workpiece(bh_curve):
    """Demonstrate ESIM workpiece creation and analysis."""
    print_header("4. ESIM Workpiece Analysis")

    sigma = 2e6  # S/m
    freq = 50000  # 50 kHz

    # Create a steel block workpiece (simulating a billet)
    print("Creating steel billet workpiece:")
    print("  Dimensions: 100mm x 100mm x 30mm")
    print("  Panels per side: 8")
    print()

    workpiece = create_esim_block(
        center=[0, 0, -0.015],  # Center at z=-15mm
        dimensions=[0.1, 0.1, 0.03],  # 100x100x30 mm
        bh_curve=bh_curve,
        sigma=sigma,
        frequency=freq,
        panels_per_side=8
    )

    print(f"Workpiece created:")
    print(f"  Number of panels: {workpiece.num_panels}")
    print(f"  Total surface area: {workpiece.total_surface_area*1e4:.1f} cm^2")
    print()

    # Simulate different field distributions
    print("Power Analysis for Different Field Distributions:")
    print("-" * 70)

    # Case 1: Uniform low field
    H_low = 500  # A/m
    for panel in workpiece.panels:
        workpiece.set_tangential_field(panel.panel_id, H_low)
    P1, Q1 = workpiece.compute_power_losses()
    summary1 = workpiece.get_summary()

    print(f"Case 1: Uniform H_t = {H_low} A/m")
    print(f"  Total power: P = {P1:.1f} W, Q = {Q1:.1f} var")
    print(f"  Max power density: {summary1['max_P_density']/1e3:.2f} kW/m^2")
    print()

    # Case 2: Uniform high field (saturation)
    H_high = 5000  # A/m
    for panel in workpiece.panels:
        workpiece.set_tangential_field(panel.panel_id, H_high)
    P2, Q2 = workpiece.compute_power_losses()
    summary2 = workpiece.get_summary()

    print(f"Case 2: Uniform H_t = {H_high} A/m (near saturation)")
    print(f"  Total power: P = {P2:.1f} W, Q = {Q2:.1f} var")
    print(f"  Max power density: {summary2['max_P_density']/1e3:.2f} kW/m^2")
    print()

    # Case 3: Top surface only (typical induction heating)
    # Top face panels are the first n^2 panels
    n = 8  # panels_per_side
    H_top = 3000  # A/m
    for panel in workpiece.panels:
        if panel.panel_id < n * n:  # Top face
            workpiece.set_tangential_field(panel.panel_id, H_top)
        else:
            workpiece.set_tangential_field(panel.panel_id, H_top * 0.1)  # Small fringing
    P3, Q3 = workpiece.compute_power_losses()
    summary3 = workpiece.get_summary()

    print(f"Case 3: H_t = {H_top} A/m on top face only")
    print(f"  Total power: P = {P3:.1f} W, Q = {Q3:.1f} var")
    print(f"  Max power density: {summary3['max_P_density']/1e3:.2f} kW/m^2")
    print()

    # Power ratio analysis (nonlinear effect)
    print("Nonlinear Effect Analysis:")
    print(f"  Power ratio (H_high/H_low): {P2/P1:.2f}")
    print(f"  Linear expectation (H^2 ratio): {(H_high/H_low)**2:.2f}")
    print(f"  Nonlinear reduction factor: {P2/P1 / (H_high/H_low)**2:.2f}")
    print("  (Values < 1 indicate saturation reducing losses)")

    return workpiece


def demo_frequency_sweep(bh_curve):
    """Demonstrate frequency dependence of ESIM."""
    print_header("5. Frequency Sweep Analysis")

    sigma = 2e6  # S/m
    H0 = 2000  # A/m (moderate field)

    print(f"Analyzing Z(f) at H_t = {H0} A/m")
    print("-" * 70)
    print(f"{'Freq [kHz]':>12} {'delta [mm]':>12} {'Re(Z) [mOhm]':>14} "
          f"{'Im(Z) [mOhm]':>14} {'|Z| [mOhm]':>12}")
    print("-" * 70)

    frequencies = [1000, 5000, 10000, 25000, 50000, 100000, 200000]

    for freq in frequencies:
        solver = ESIMCellProblemSolver(bh_curve, sigma, freq)
        result = solver.solve(H0)
        Z = result['Z']
        delta = solver.delta_initial * 1e3  # mm

        print(f"{freq/1000:>12.0f} {delta:>12.3f} {Z.real*1e3:>14.4f} "
              f"{Z.imag*1e3:>14.4f} {abs(Z)*1e3:>12.4f}")


def main():
    """Run all ESIM demos."""
    print()
    print("*" * 70)
    print("*  ESIM (Effective Surface Impedance Method) Demonstration")
    print("*  For Induction Heating Analysis with Nonlinear Materials")
    print("*" * 70)

    # Run demos
    bh_curve = demo_bh_curve()
    demo_cell_problem(bh_curve)
    demo_esi_table(bh_curve)
    demo_workpiece(bh_curve)
    demo_frequency_sweep(bh_curve)

    print_header("Demo Complete")
    print("The ESIM implementation is ready for induction heating analysis.")
    print()
    print("Key features demonstrated:")
    print("  1. BH curve interpolation for nonlinear permeability")
    print("  2. Cell Problem solver for effective surface impedance")
    print("  3. ESI table generation and interpolation")
    print("  4. Workpiece modeling with surface panels")
    print("  5. Power loss computation with nonlinear materials")
    print()
    print("Next steps:")
    print("  - Integrate with FastImp coil for coupled analysis")
    print("  - Implement fixed-point iteration for 3D solve")
    print("  - Add VTK export for visualization")
    print()


if __name__ == "__main__":
    main()
