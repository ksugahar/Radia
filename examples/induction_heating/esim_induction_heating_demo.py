"""
ESIM Induction Heating Coupled Solver Demo

This script demonstrates the complete induction heating analysis workflow using:
1. Spiral coil model (FastImp or analytical)
2. ESIM workpiece with nonlinear ferromagnetic material
3. Fixed-point coupled solver

The demo simulates a typical induction heating setup for surface hardening
of a steel billet.

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

from esim_coupled_solver import (
    InductionHeatingCoil,
    ESIMCoupledSolver,
    solve_induction_heating
)
from esim_workpiece import create_esim_block, create_esim_cylinder


def print_header(title):
    """Print a formatted header."""
    print()
    print("=" * 70)
    print(title)
    print("=" * 70)


def demo_spiral_coil_block():
    """Demo: Spiral coil heating a rectangular steel block."""
    print_header("Demo 1: Spiral Coil + Steel Block")

    # Steel BH curve (typical hot working steel at ~800C)
    bh_curve_steel = [
        [0, 0],
        [100, 0.15],
        [250, 0.4],
        [500, 0.7],
        [1000, 1.0],
        [2500, 1.35],
        [5000, 1.55],
        [10000, 1.7],
        [50000, 1.9],
    ]

    sigma_steel = 1.5e6  # S/m (hot steel at 800C)
    freq = 30000  # 30 kHz (typical for surface hardening)

    print("Material Properties:")
    print(f"  BH curve: Typical hot steel at 800C")
    print(f"  Conductivity: {sigma_steel/1e6:.2f} MS/m")
    print(f"  Frequency: {freq/1000:.0f} kHz")
    print()

    # Coil parameters
    coil_params = {
        'type': 'spiral',
        'center': [0, 0, 0.015],  # 15mm gap between coil and workpiece
        'inner_radius': 0.025,    # 25mm inner radius
        'outer_radius': 0.055,    # 55mm outer radius
        'pitch': 0.006,           # 6mm pitch
        'num_turns': 4,
        'axis': [0, 0, 1],
        'wire_width': 0.004,      # 4mm wire width
        'wire_height': 0.003,     # 3mm wire height
        'cross_section': 'r',
        'conductivity': 5.8e7,    # Copper
    }

    # Workpiece parameters (steel billet)
    workpiece_params = {
        'geometry': 'block',
        'center': [0, 0, -0.0125],  # 12.5mm below origin
        'dimensions': [0.08, 0.08, 0.025],  # 80mm x 80mm x 25mm
        'panels_per_side': 6,
    }

    print("Coil Geometry:")
    print(f"  Type: {coil_params['type']}")
    print(f"  Turns: {coil_params['num_turns']}")
    print(f"  R_inner: {coil_params['inner_radius']*1000:.0f} mm")
    print(f"  R_outer: {coil_params['outer_radius']*1000:.0f} mm")
    print(f"  Gap to workpiece: {coil_params['center'][2]*1000:.0f} mm")
    print()

    print("Workpiece Geometry:")
    print(f"  Dimensions: {workpiece_params['dimensions'][0]*1000:.0f} x "
          f"{workpiece_params['dimensions'][1]*1000:.0f} x "
          f"{workpiece_params['dimensions'][2]*1000:.0f} mm")
    print()

    # Solve
    print("Solving coupled problem...")
    result = solve_induction_heating(
        coil_params=coil_params,
        workpiece_params=workpiece_params,
        frequency=freq,
        bh_curve=bh_curve_steel,
        sigma=sigma_steel,
        tol=1e-4,
        max_iter=30,
        verbose=True
    )

    # Summary
    print()
    print("Results Summary:")
    print("-" * 40)
    print(f"  Converged: {result['converged']}")
    print(f"  Iterations: {result['iterations']}")
    print(f"  Total active power: {result['P_total']:.1f} W")
    print(f"  Total reactive power: {result['Q_total']:.1f} var")
    print(f"  Apparent power: {result['S_total']:.1f} VA")
    print(f"  Power factor: {result['power_factor']:.3f}")
    print(f"  Max power density: {result['max_P_density']/1e3:.2f} kW/m^2")

    return result


def demo_loop_coil_cylinder():
    """Demo: Loop coil heating a cylindrical steel billet."""
    print_header("Demo 2: Loop Coil + Steel Cylinder")

    # Steel BH curve
    bh_curve_steel = [
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

    sigma_steel = 2e6  # S/m (hot steel)
    freq = 50000  # 50 kHz

    print("Material Properties:")
    print(f"  BH curve: Typical steel")
    print(f"  Conductivity: {sigma_steel/1e6:.1f} MS/m")
    print(f"  Frequency: {freq/1000:.0f} kHz")
    print()

    # Create single loop coil
    coil = InductionHeatingCoil(
        coil_type='loop',
        center=[0, 0, 0.01],  # 10mm above
        radius=0.04,          # 40mm radius
        normal=[0, 0, 1],
        wire_width=0.003,
        wire_height=0.002,
        cross_section='r',
        conductivity=5.8e7
    )
    coil.set_current(200)  # 200 A

    print("Coil Geometry:")
    print(f"  Type: Single loop")
    print(f"  Radius: {coil.radius*1000:.0f} mm")
    print(f"  Current: {coil.current} A")
    print()

    # Create cylindrical workpiece
    workpiece = create_esim_cylinder(
        center=[0, 0, 0],  # Top at origin
        radius=0.03,       # 30mm radius
        height=0.04,       # 40mm height
        bh_curve=bh_curve_steel,
        sigma=sigma_steel,
        frequency=freq,
        panels_radial=8,
        panels_axial=5
    )

    print("Workpiece Geometry:")
    print(f"  Radius: {workpiece.radius*1000:.0f} mm")
    print(f"  Height: {workpiece.height*1000:.0f} mm")
    print(f"  Panels: {workpiece.num_panels}")
    print()

    # Create and run solver
    solver = ESIMCoupledSolver(coil, workpiece, freq)
    print("Solving coupled problem...")
    result = solver.solve(tol=1e-4, max_iter=20, verbose=True)

    # Summary
    print()
    print("Results Summary:")
    print("-" * 40)
    print(f"  Converged: {result['converged']}")
    print(f"  Iterations: {result['iterations']}")
    print(f"  Total active power: {result['P_total']:.1f} W")
    print(f"  Total reactive power: {result['Q_total']:.1f} var")
    print(f"  Max power density: {result['max_P_density']/1e3:.2f} kW/m^2")

    return result


def demo_frequency_sweep():
    """Demo: Frequency sweep to find optimal heating frequency."""
    print_header("Demo 3: Frequency Sweep Analysis")

    # Steel BH curve
    bh_curve_steel = [
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

    sigma_steel = 2e6  # S/m

    # Frequencies to analyze
    frequencies = [10000, 20000, 30000, 50000, 80000, 100000]  # 10 kHz to 100 kHz

    print("Analyzing power vs frequency...")
    print()
    print(f"{'Frequency [kHz]':>15} {'P [W]':>12} {'Q [var]':>12} {'PF':>8} {'P_max [kW/m^2]':>16}")
    print("-" * 70)

    results = []

    for freq in frequencies:
        # Create coil
        coil = InductionHeatingCoil(
            coil_type='spiral',
            center=[0, 0, 0.015],
            inner_radius=0.03,
            outer_radius=0.05,
            pitch=0.005,
            num_turns=3,
            axis=[0, 0, 1],
            wire_width=0.003,
            wire_height=0.002,
            cross_section='r'
        )
        coil.set_current(100)

        # Create workpiece
        workpiece = create_esim_block(
            center=[0, 0, -0.01],
            dimensions=[0.08, 0.08, 0.02],
            bh_curve=bh_curve_steel,
            sigma=sigma_steel,
            frequency=freq,
            panels_per_side=5
        )

        # Solve
        solver = ESIMCoupledSolver(coil, workpiece, freq)
        result = solver.solve(tol=1e-3, max_iter=15, verbose=False)

        results.append({
            'frequency': freq,
            'P_total': result['P_total'],
            'Q_total': result['Q_total'],
            'power_factor': result['power_factor'],
            'max_P_density': result['max_P_density'],
        })

        print(f"{freq/1000:>15.0f} {result['P_total']:>12.1f} {result['Q_total']:>12.1f} "
              f"{result['power_factor']:>8.3f} {result['max_P_density']/1e3:>16.2f}")

    print()
    print("Analysis Notes:")
    print("  - Higher frequency -> smaller skin depth -> more surface heating")
    print("  - Power density increases with frequency (skin effect)")
    print("  - Power factor varies due to nonlinear material effects")

    # Find optimal frequency (max power for given current)
    P_values = [r['P_total'] for r in results]
    max_P_idx = np.argmax(P_values)
    optimal_freq = frequencies[max_P_idx]

    print()
    print(f"Optimal frequency for maximum power: {optimal_freq/1000:.0f} kHz")
    print(f"  P = {results[max_P_idx]['P_total']:.1f} W")

    return results


def demo_current_scaling():
    """Demo: Power scaling with current."""
    print_header("Demo 4: Power vs Current Analysis")

    # Steel BH curve
    bh_curve_steel = [
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

    sigma_steel = 2e6  # S/m
    freq = 50000  # 50 kHz

    # Currents to analyze
    currents = [10, 25, 50, 100, 200, 300]

    print(f"Frequency: {freq/1000:.0f} kHz")
    print()
    print(f"{'Current [A]':>12} {'P [W]':>12} {'P/I^2':>12} {'Comment':>20}")
    print("-" * 60)

    results = []

    for current in currents:
        # Create coil
        coil = InductionHeatingCoil(
            coil_type='spiral',
            center=[0, 0, 0.015],
            inner_radius=0.03,
            outer_radius=0.05,
            pitch=0.005,
            num_turns=3,
            axis=[0, 0, 1],
            wire_width=0.003,
            wire_height=0.002,
        )
        coil.set_current(current)

        # Create workpiece
        workpiece = create_esim_block(
            center=[0, 0, -0.01],
            dimensions=[0.08, 0.08, 0.02],
            bh_curve=bh_curve_steel,
            sigma=sigma_steel,
            frequency=freq,
            panels_per_side=5
        )

        # Solve
        solver = ESIMCoupledSolver(coil, workpiece, freq)
        result = solver.solve(tol=1e-3, max_iter=15, verbose=False)

        P_over_I2 = result['P_total'] / (current ** 2)

        # Comment on nonlinear effects
        if current <= 25:
            comment = "Linear region"
        elif current <= 100:
            comment = "Transition"
        else:
            comment = "Saturation region"

        results.append({
            'current': current,
            'P_total': result['P_total'],
            'P_over_I2': P_over_I2,
        })

        print(f"{current:>12.0f} {result['P_total']:>12.1f} {P_over_I2:>12.6f} {comment:>20}")

    print()
    print("Analysis Notes:")
    print("  - For linear materials, P/I^2 would be constant")
    print("  - Decreasing P/I^2 at high currents indicates saturation")
    print("  - Saturation reduces heating efficiency at high field levels")

    return results


def main():
    """Run all ESIM coupled solver demos."""
    print()
    print("*" * 70)
    print("*  ESIM Induction Heating Coupled Solver Demonstration")
    print("*  FastImp Coil + ESIM Workpiece Analysis")
    print("*" * 70)

    # Run demos
    demo_spiral_coil_block()
    demo_loop_coil_cylinder()
    demo_frequency_sweep()
    demo_current_scaling()

    print_header("Demo Complete")
    print("The ESIM coupled solver is ready for induction heating analysis.")
    print()
    print("Key features demonstrated:")
    print("  1. Spiral coil field computation (Biot-Savart)")
    print("  2. Loop coil for cylindrical workpiece")
    print("  3. Frequency sweep analysis")
    print("  4. Current scaling and saturation effects")
    print()
    print("Next steps for production use:")
    print("  - Integrate with Radia's FastImp API for accurate coil impedance")
    print("  - Add VTK export for visualization")
    print("  - Validate against experimental/FEM data")
    print()


if __name__ == "__main__":
    main()
