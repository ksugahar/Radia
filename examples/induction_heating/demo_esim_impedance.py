"""
ESIM Coupled Solver Demo: Coil-Workpiece Impedance Calculation

This example demonstrates:
1. Complex permeability support (mu' - j*mu")
2. Coil-workpiece coupled impedance calculation
3. Efficiency and power analysis for induction heating

Author: Radia Development Team
Date: 2026-01-08
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src'))

import numpy as np

# Import ESIM modules from radia package
from radia import (
    ESIMCellProblemSolver,
    BHCurveInterpolator,
    ComplexPermeabilityInterpolator,
    InductionHeatingCoil,
    ESIMCoupledSolver,
    create_esim_block,
    create_esim_cylinder,
)
import radia as rad


def demo_complex_permeability():
    """Demonstrate complex permeability support."""
    print("=" * 60)
    print("Demo 1: Complex Permeability (mu' - j*mu\")")
    print("=" * 60)
    print()

    # Ferrite material with complex permeability
    # mu'_r = 1000 (energy storage), mu"_r = 100 (loss)
    # Loss tangent = mu"/mu' = 0.1

    mu_prime_r = 1000
    mu_double_prime_r = 100
    sigma = 1e6  # 1 MS/m (steel)
    freq = 50000  # 50 kHz

    print(f"Material properties:")
    print(f"  mu'_r = {mu_prime_r}")
    print(f"  mu\"_r = {mu_double_prime_r}")
    print(f"  Loss tangent = {mu_double_prime_r/mu_prime_r:.3f}")
    print(f"  sigma = {sigma/1e6:.1f} MS/m")
    print(f"  freq = {freq/1000:.1f} kHz")
    print()

    # Create solver with complex permeability
    solver = ESIMCellProblemSolver(
        sigma=sigma,
        frequency=freq,
        complex_mu=(mu_prime_r, mu_double_prime_r)
    )

    # Solve for different surface field values
    print("Cell Problem Solutions:")
    print("-" * 70)
    print(f"{'H0 [A/m]':>12} {'Re(Z) [mOhm]':>14} {'Im(Z) [mOhm]':>14} {'P_mag [W/m^2]':>14}")
    print("-" * 70)

    for H0 in [100, 500, 1000, 5000]:
        result = solver.solve(H0)
        Z = result['Z']
        P_mag = result['P_magnetic']
        print(f"{H0:>12.0f} {Z.real*1e3:>14.4f} {Z.imag*1e3:>14.4f} {P_mag:>14.2f}")

    print()


def demo_h_dependent_complex_mu():
    """Demonstrate H-dependent complex permeability."""
    print("=" * 60)
    print("Demo 2: H-Dependent Complex Permeability")
    print("=" * 60)
    print()

    # H-dependent permeability data [H, mu'_r, mu"_r]
    # Simulates saturation with decreasing mu' and mu"
    complex_mu_data = [
        [0, 2000, 200],        # Low field: high permeability
        [100, 1800, 180],
        [500, 1500, 150],
        [1000, 1000, 100],
        [5000, 500, 50],
        [10000, 200, 20],      # High field: saturated
    ]

    sigma = 1e6  # S/m
    freq = 50000  # Hz

    print("H-dependent mu' and mu\":")
    print(f"{'H [A/m]':>10} {'mu_r':>8} {'mu\"_r':>8}")
    for row in complex_mu_data:
        print(f"{row[0]:>10.0f} {row[1]:>8.0f} {row[2]:>8.0f}")
    print()

    # Create solver
    solver = ESIMCellProblemSolver(
        sigma=sigma,
        frequency=freq,
        complex_mu=complex_mu_data
    )

    print("Cell Problem Solutions:")
    print("-" * 80)
    print(f"{'H0 [A/m]':>12} {'Re(Z) [mOhm]':>14} {'Im(Z) [mOhm]':>14} {'P_total [kW/m^2]':>18} {'Iter':>6}")
    print("-" * 80)

    for H0 in [100, 500, 1000, 5000, 10000]:
        result = solver.solve(H0)
        Z = result['Z']
        P = result['P_prime']
        iters = result['iterations']
        print(f"{H0:>12.0f} {Z.real*1e3:>14.4f} {Z.imag*1e3:>14.4f} {P/1e3:>18.2f} {iters:>6}")

    print()


def demo_impedance_calculation():
    """Demonstrate coil-workpiece coupled impedance calculation."""
    print("=" * 60)
    print("Demo 3: Coil-Workpiece Coupled Impedance")
    print("=" * 60)
    print()

    # Steel BH curve
    bh_curve = [
        [0, 0], [100, 0.2], [250, 0.5], [500, 0.9],
        [1000, 1.3], [2500, 1.6], [5000, 1.8], [10000, 1.95],
    ]
    sigma = 2e6  # 2 MS/m (hot steel)
    freq = 50000  # 50 kHz

    # Create spiral coil
    coil = InductionHeatingCoil(
        coil_type='spiral',
        center=[0, 0, 0.02],      # 20mm above workpiece
        inner_radius=0.03,        # 30mm inner
        outer_radius=0.06,        # 60mm outer
        pitch=0.005,              # 5mm pitch
        num_turns=4,
        axis=[0, 0, 1],
        wire_width=0.004,
        wire_height=0.002,
        conductivity=5.8e7,       # Copper
    )
    coil.set_current(150)  # 150 A

    # Create steel slab workpiece
    workpiece = create_esim_block(
        center=[0, 0, -0.01],          # 10mm below origin
        dimensions=[0.12, 0.12, 0.02],  # 120mm x 120mm x 20mm
        bh_curve=bh_curve,
        sigma=sigma,
        frequency=freq,
        panels_per_side=5
    )

    print("Configuration:")
    print(f"  Coil: 4-turn spiral, R_in=30mm, R_out=60mm")
    print(f"  Coil current: {coil.current} A")
    print(f"  Workpiece: 120mm x 120mm x 20mm steel slab")
    print(f"  Gap: 20mm")
    print(f"  Frequency: {freq/1000} kHz")
    print()

    # Create and solve
    solver = ESIMCoupledSolver(coil, workpiece, freq)
    result = solver.solve(tol=1e-4, max_iter=30, verbose=True)

    # Get impedance summary
    imp = result['impedance']

    print()
    print("=" * 60)
    print("Impedance Analysis Summary")
    print("=" * 60)
    print()
    print("Coil Self-Impedance:")
    print(f"  Inductance: L = {imp['L_coil_uH']:.3f} uH")
    print(f"  AC Resistance: R = {imp['R_coil_mOhm']:.3f} mOhm")
    print(f"  Z_coil = {solver.Z_coil_self.real*1e3:.3f} + j{solver.Z_coil_self.imag*1e3:.3f} mOhm")
    print()
    print("Reflected Impedance (from workpiece):")
    print(f"  R_reflected = {imp['R_reflected_mOhm']:.3f} mOhm (heating power)")
    print(f"  X_reflected = {imp['X_reflected_mOhm']:.3f} mOhm (reactive)")
    print(f"  Z_reflected = {solver.Z_reflected.real*1e3:.3f} + j{solver.Z_reflected.imag*1e3:.3f} mOhm")
    print()
    print("Total System Impedance:")
    print(f"  |Z_total| = {imp['Z_total_magnitude_mOhm']:.3f} mOhm")
    print(f"  Phase = {imp['phase_deg']:.1f} deg")
    print()
    print("Power Analysis:")
    print(f"  Total power to workpiece: P = {result['P_total']:.0f} W")
    print(f"  Reactive power: Q = {result['Q_total']:.0f} var")
    print(f"  Power factor: {result['power_factor']:.3f}")
    print(f"  Heating efficiency: {imp['efficiency']*100:.1f}%")
    print()

    # Resonance capacitor calculation
    C_res = 1 / (solver.omega**2 * solver.L_coil)
    print("Resonance Matching:")
    print(f"  Required capacitance for resonance: C = {C_res*1e6:.2f} uF")
    print()


def demo_frequency_sweep():
    """Demonstrate frequency sweep analysis."""
    print("=" * 60)
    print("Demo 4: Frequency Sweep Analysis")
    print("=" * 60)
    print()

    bh_curve = [
        [0, 0], [100, 0.2], [500, 0.9], [1000, 1.3], [5000, 1.8],
    ]
    sigma = 2e6

    # Single-turn loop coil
    coil = InductionHeatingCoil(
        coil_type='loop',
        center=[0, 0, 0.015],
        radius=0.04,
        normal=[0, 0, 1],
        wire_width=0.005,
        wire_height=0.003,
    )
    coil.set_current(300)

    print("Configuration:")
    print(f"  Coil: Single-turn loop, R=40mm, I=300A")
    print(f"  Workpiece: 80mm x 80mm steel slab")
    print()

    frequencies = [10000, 25000, 50000, 100000]  # 10, 25, 50, 100 kHz

    print(f"{'Freq [kHz]':>12} | {'P [W]':>10} | {'Eff [%]':>8} | {'|Z| [mOhm]':>12} | {'C_res [uF]':>10}")
    print("-" * 65)

    for freq in frequencies:
        workpiece = create_esim_block(
            center=[0, 0, -0.01],
            dimensions=[0.08, 0.08, 0.015],
            bh_curve=bh_curve,
            sigma=sigma,
            frequency=freq,
            panels_per_side=4
        )

        solver = ESIMCoupledSolver(coil, workpiece, freq)
        result = solver.solve(tol=1e-3, max_iter=20, verbose=False)

        imp = result['impedance']
        C_res = 1 / (solver.omega**2 * solver.L_coil) if solver.L_coil else 0

        print(f"{freq/1000:>12.0f} | {result['P_total']:>10.0f} | "
              f"{imp['efficiency']*100:>8.1f} | {imp['Z_total_magnitude_mOhm']:>12.2f} | "
              f"{C_res*1e6:>10.2f}")

    print()


if __name__ == "__main__":
    print()
    print("ESIM Induction Heating Solver - Impedance Calculation Demo")
    print("=" * 60)
    print()

    demo_complex_permeability()
    demo_h_dependent_complex_mu()
    demo_impedance_calculation()
    demo_frequency_sweep()

    print("All demos completed!")
