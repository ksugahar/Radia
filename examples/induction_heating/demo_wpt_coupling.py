"""
WPT (Wireless Power Transfer) Coil Coupling Analysis Demo

This example demonstrates:
1. Mutual inductance calculation using Neumann integral
2. Coupling coefficient k = M / sqrt(L1*L2)
3. Mutual resistance Rm (proximity effect between coils)
4. WPT system efficiency analysis
5. Resonant capacitor design for S-S topology

Mutual Resistance (Rm):
    The mutual resistance accounts for eddy current losses induced in one coil
    by the alternating magnetic field of the other coil (proximity effect).
    The impedance matrix off-diagonal elements are: Z12 = Z21 = Rm + j*omega*M

Author: Radia Development Team
Date: 2026-01-08
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src'))

import numpy as np

# Import WPT modules from radia package
from radia import (
    InductionHeatingCoil,
    WPTCoupledSolver,
    compute_mutual_inductance,
    compute_coupling_coefficient,
    analyze_coil_coupling,
)
import radia as rad
rad.FldUnits('m')


def demo_coupling_vs_distance():
    """Demonstrate coupling coefficient vs. coil separation distance."""
    print("=" * 60)
    print("Demo 1: Coupling Coefficient vs. Distance")
    print("=" * 60)
    print()

    # Create primary coil (transmitter) - 5-turn spiral
    coil_tx = InductionHeatingCoil(
        coil_type='spiral',
        center=[0, 0, 0],
        inner_radius=0.02,  # 20mm
        outer_radius=0.05,  # 50mm
        pitch=0.002,        # 2mm
        num_turns=5,
        axis=[0, 0, 1],
        wire_width=0.003,
        wire_height=0.001,
        conductivity=5.8e7,
    )

    frequency = 85000  # 85 kHz (Qi standard)

    print(f"Primary coil: 5-turn spiral, R_in=20mm, R_out=50mm")
    print(f"Frequency: {frequency/1000} kHz")
    print()

    # Vary distance
    distances = [5, 10, 15, 20, 30, 40, 50]  # mm

    print(f"{'Distance [mm]':>14} | {'k':>8} | {'M [uH]':>10} | {'L1 [uH]':>10} | {'L2 [uH]':>10}")
    print("-" * 70)

    for d_mm in distances:
        d = d_mm / 1000  # Convert to meters

        # Create secondary coil (receiver) at distance d
        coil_rx = InductionHeatingCoil(
            coil_type='spiral',
            center=[0, 0, d],  # Distance along z-axis
            inner_radius=0.02,
            outer_radius=0.05,
            pitch=0.002,
            num_turns=5,
            axis=[0, 0, 1],
            wire_width=0.003,
            wire_height=0.001,
            conductivity=5.8e7,
        )

        # Compute coupling
        k, L1, L2, M = compute_coupling_coefficient(coil_tx, coil_rx, n_segments=50)

        print(f"{d_mm:>14} | {k:>8.4f} | {M*1e6:>10.3f} | {L1*1e6:>10.3f} | {L2*1e6:>10.3f}")

    print()


def demo_coaxial_loops():
    """Demonstrate coupling between two coaxial circular loops."""
    print("=" * 60)
    print("Demo 2: Coaxial Circular Loops")
    print("=" * 60)
    print()

    # Two identical circular loops
    R1 = 0.05  # 50mm radius
    R2 = 0.05  # 50mm radius
    d = 0.03   # 30mm separation

    coil1 = InductionHeatingCoil(
        coil_type='loop',
        center=[0, 0, 0],
        radius=R1,
        normal=[0, 0, 1],
        wire_width=0.002,
        wire_height=0.002,
        conductivity=5.8e7,
    )

    coil2 = InductionHeatingCoil(
        coil_type='loop',
        center=[0, 0, d],
        radius=R2,
        normal=[0, 0, 1],
        wire_width=0.002,
        wire_height=0.002,
        conductivity=5.8e7,
    )

    print(f"Loop 1: R = {R1*1000:.0f}mm at z = 0")
    print(f"Loop 2: R = {R2*1000:.0f}mm at z = {d*1000:.0f}mm")
    print()

    # Compute coupling
    k, L1, L2, M = compute_coupling_coefficient(coil1, coil2, n_segments=100)

    print("Neumann Integral Results:")
    print(f"  L1 = {L1*1e9:.2f} nH")
    print(f"  L2 = {L2*1e9:.2f} nH")
    print(f"  M  = {M*1e9:.2f} nH")
    print(f"  k  = {k:.4f}")
    print()

    # Analytical formula for coaxial loops (Maxwell's formula)
    # M = mu_0 * sqrt(R1*R2) * ((2/k_m - k_m) * K(k_m) - 2/k_m * E(k_m))
    # where k_m^2 = 4*R1*R2 / ((R1+R2)^2 + d^2)
    # K, E are complete elliptic integrals

    # Simplified formula for R1 = R2 = R:
    # M = mu_0 * R * ((2/k_m - k_m) * K(k_m) - 2/k_m * E(k_m))
    from scipy.special import ellipk, ellipe
    from scipy.constants import mu_0

    k_m_sq = 4 * R1 * R2 / ((R1 + R2)**2 + d**2)
    k_m = np.sqrt(k_m_sq)

    K_k = ellipk(k_m_sq)
    E_k = ellipe(k_m_sq)

    M_analytical = mu_0 * np.sqrt(R1 * R2) * ((2/k_m - k_m) * K_k - 2/k_m * E_k)

    print("Analytical (Maxwell's formula):")
    print(f"  M_analytical = {M_analytical*1e9:.2f} nH")
    print(f"  Relative error = {abs(M - M_analytical)/M_analytical * 100:.2f}%")
    print()


def demo_wpt_system_analysis():
    """Demonstrate complete WPT system analysis."""
    print("=" * 60)
    print("Demo 3: Complete WPT System Analysis")
    print("=" * 60)
    print()

    # EV charging coil geometry (scaled down for demo)
    # Primary: Ground Assembly (GA)
    coil_ga = InductionHeatingCoil(
        coil_type='spiral',
        center=[0, 0, 0],
        inner_radius=0.08,   # 80mm inner
        outer_radius=0.15,   # 150mm outer
        pitch=0.005,         # 5mm pitch
        num_turns=8,
        axis=[0, 0, 1],
        wire_width=0.006,    # Litz wire equivalent
        wire_height=0.003,
        conductivity=5.8e7,
    )

    # Secondary: Vehicle Assembly (VA) - at 100mm gap
    gap = 0.10  # 100mm air gap
    coil_va = InductionHeatingCoil(
        coil_type='spiral',
        center=[0, 0, gap],
        inner_radius=0.08,
        outer_radius=0.15,
        pitch=0.005,
        num_turns=8,
        axis=[0, 0, 1],
        wire_width=0.006,
        wire_height=0.003,
        conductivity=5.8e7,
    )

    frequency = 85000  # 85 kHz (SAE J2954)

    print("EV Wireless Charging System (Scaled)")
    print(f"  Ground Assembly (GA): 8-turn spiral, R_in=80mm, R_out=150mm")
    print(f"  Vehicle Assembly (VA): Same geometry, {gap*1000:.0f}mm above GA")
    print(f"  Frequency: {frequency/1000} kHz (SAE J2954)")
    print()

    # Create WPT solver and analyze
    solver = WPTCoupledSolver(coil_ga, coil_va, frequency)
    result = solver.analyze(n_segments=80, topology='SS', verbose=True)

    # Additional analysis: efficiency vs load
    print("Efficiency vs. Load Resistance:")
    print(f"{'R_load [Ohm]':>14} | {'eta [%]':>10} | {'P_load [W]':>12}")
    print("-" * 45)

    for R_load in [1, 2, 5, 10, 20, 50, 100]:
        eta, P_load, _ = solver.compute_transfer_efficiency(R_load)
        print(f"{R_load:>14} | {eta*100:>10.1f} | {P_load:>12.2f}")

    print()


def demo_misalignment():
    """Demonstrate coupling degradation with lateral misalignment."""
    print("=" * 60)
    print("Demo 4: Coupling vs. Lateral Misalignment")
    print("=" * 60)
    print()

    # Fixed vertical gap
    gap = 0.03  # 30mm

    # Primary coil (fixed)
    coil_tx = InductionHeatingCoil(
        coil_type='loop',
        center=[0, 0, 0],
        radius=0.05,  # 50mm
        normal=[0, 0, 1],
        wire_width=0.003,
        wire_height=0.002,
        conductivity=5.8e7,
    )

    print(f"Primary coil: R = 50mm at origin")
    print(f"Secondary coil: R = 50mm, {gap*1000:.0f}mm above")
    print()

    # Vary lateral offset
    offsets = [0, 10, 20, 30, 40, 50]  # mm

    print(f"{'Offset [mm]':>12} | {'k':>8} | {'k/k0':>8} | {'M [nH]':>10}")
    print("-" * 50)

    k0 = None
    for offset_mm in offsets:
        offset = offset_mm / 1000

        # Secondary coil with lateral offset
        coil_rx = InductionHeatingCoil(
            coil_type='loop',
            center=[offset, 0, gap],  # Offset in x-direction
            radius=0.05,
            normal=[0, 0, 1],
            wire_width=0.003,
            wire_height=0.002,
            conductivity=5.8e7,
        )

        k, L1, L2, M = compute_coupling_coefficient(coil_tx, coil_rx, n_segments=80)

        if k0 is None:
            k0 = k

        print(f"{offset_mm:>12} | {k:>8.4f} | {k/k0:>8.3f} | {M*1e9:>10.2f}")

    print()
    print("Note: k/k0 shows coupling relative to aligned case")
    print()


def demo_different_coil_sizes():
    """Demonstrate coupling between different sized coils."""
    print("=" * 60)
    print("Demo 5: Coupling Between Different Coil Sizes")
    print("=" * 60)
    print()

    # Large primary coil (transmitter pad)
    coil_tx = InductionHeatingCoil(
        coil_type='spiral',
        center=[0, 0, 0],
        inner_radius=0.05,   # 50mm inner
        outer_radius=0.10,   # 100mm outer
        pitch=0.003,
        num_turns=6,
        axis=[0, 0, 1],
        wire_width=0.003,
        wire_height=0.002,
        conductivity=5.8e7,
    )

    gap = 0.02  # 20mm gap

    print("Primary: 6-turn spiral, R_in=50mm, R_out=100mm")
    print(f"Gap: {gap*1000:.0f}mm")
    print()

    # Vary secondary coil size
    secondary_sizes = [
        (0.03, 0.06),   # 30-60mm (small)
        (0.04, 0.08),   # 40-80mm (medium)
        (0.05, 0.10),   # 50-100mm (matched)
        (0.06, 0.12),   # 60-120mm (large)
    ]

    print(f"{'Secondary [mm]':>18} | {'k':>8} | {'M [uH]':>10} | {'L2 [uH]':>10}")
    print("-" * 60)

    for R_in, R_out in secondary_sizes:
        coil_rx = InductionHeatingCoil(
            coil_type='spiral',
            center=[0, 0, gap],
            inner_radius=R_in,
            outer_radius=R_out,
            pitch=0.003,
            num_turns=6,
            axis=[0, 0, 1],
            wire_width=0.003,
            wire_height=0.002,
            conductivity=5.8e7,
        )

        k, L1, L2, M = compute_coupling_coefficient(coil_tx, coil_rx, n_segments=60)

        size_str = f"{R_in*1000:.0f}-{R_out*1000:.0f}"
        print(f"{size_str:>18} | {k:>8.4f} | {M*1e6:>10.3f} | {L2*1e6:>10.3f}")

    print()


def demo_mutual_resistance():
    """Demonstrate mutual resistance (proximity effect) calculation."""
    print("=" * 60)
    print("Demo 6: Mutual Resistance (Proximity Effect)")
    print("=" * 60)
    print()

    print("Mutual resistance (Rm) represents eddy current losses induced")
    print("in one coil by the magnetic field of the other coil.")
    print()
    print("Impedance matrix: Z12 = Z21 = Rm + j*omega*M")
    print()

    # Create two closely coupled coils
    gap = 0.01  # 10mm - close coupling

    coil_tx = InductionHeatingCoil(
        coil_type='spiral',
        center=[0, 0, 0],
        inner_radius=0.03,
        outer_radius=0.06,
        pitch=0.003,
        num_turns=5,
        axis=[0, 0, 1],
        wire_width=0.004,    # 4mm wire width
        wire_height=0.002,   # 2mm wire height
        conductivity=5.8e7,  # Copper
    )

    coil_rx = InductionHeatingCoil(
        coil_type='spiral',
        center=[0, 0, gap],
        inner_radius=0.03,
        outer_radius=0.06,
        pitch=0.003,
        num_turns=5,
        axis=[0, 0, 1],
        wire_width=0.004,
        wire_height=0.002,
        conductivity=5.8e7,
    )

    print(f"Coil geometry: 5-turn spiral, R_in=30mm, R_out=60mm")
    print(f"Wire: 4mm x 2mm rectangular, copper (sigma=5.8e7 S/m)")
    print(f"Gap: {gap*1000:.0f}mm")
    print()

    # Analyze at different frequencies
    frequencies = [10000, 50000, 100000, 200000]  # 10, 50, 100, 200 kHz

    print(f"{'Freq [kHz]':>12} | {'k':>8} | {'Rm [mOhm]':>12} | {'omega*M [Ohm]':>14} | {'Rm/omega*M':>12}")
    print("-" * 75)

    for freq in frequencies:
        solver = WPTCoupledSolver(coil_tx, coil_rx, freq)
        solver.compute_inductances(n_segments=60)
        solver.compute_impedance_matrix()

        omega_M = solver.omega * abs(solver.M)
        Rm = solver.Rm if solver.Rm else 0

        ratio = Rm / omega_M if omega_M > 0 else 0

        print(f"{freq/1000:>12.0f} | {solver.k:>8.4f} | {Rm*1e3:>12.6f} | {omega_M:>14.4f} | {ratio:>12.6f}")

    print()
    print("Note: Rm/omega*M ratio indicates the relative importance of proximity losses.")
    print("      Higher frequencies increase omega*M but also increase Rm (skin effect).")
    print()


def demo_impedance_matrix_components():
    """Demonstrate full impedance matrix with mutual resistance."""
    print("=" * 60)
    print("Demo 7: Full Impedance Matrix Analysis")
    print("=" * 60)
    print()

    # EV charging scenario
    frequency = 85000  # 85 kHz (SAE J2954)
    gap = 0.15  # 150mm air gap

    coil_ga = InductionHeatingCoil(
        coil_type='spiral',
        center=[0, 0, 0],
        inner_radius=0.10,
        outer_radius=0.20,
        pitch=0.008,
        num_turns=8,
        axis=[0, 0, 1],
        wire_width=0.008,
        wire_height=0.004,
        conductivity=5.8e7,
    )

    coil_va = InductionHeatingCoil(
        coil_type='spiral',
        center=[0, 0, gap],
        inner_radius=0.10,
        outer_radius=0.20,
        pitch=0.008,
        num_turns=8,
        axis=[0, 0, 1],
        wire_width=0.008,
        wire_height=0.004,
        conductivity=5.8e7,
    )

    print("EV Wireless Charging Scenario (SAE J2954)")
    print(f"  Frequency: {frequency/1000} kHz")
    print(f"  Air gap: {gap*1000:.0f}mm")
    print(f"  Coil: 8-turn spiral, R_in=100mm, R_out=200mm")
    print()

    # Run full analysis (verbose output includes Rm)
    solver = WPTCoupledSolver(coil_ga, coil_va, frequency)
    result = solver.analyze(n_segments=80, topology='SS', verbose=True)

    # Additional analysis: show breakdown of impedance components
    print("Impedance Component Breakdown:")
    print("-" * 50)
    print(f"  Self-resistance R1:     {result['R1_mOhm']:.4f} mOhm")
    print(f"  Self-resistance R2:     {result['R2_mOhm']:.4f} mOhm")
    print(f"  Mutual resistance Rm:   {result['Rm_mOhm']:.4f} mOhm")
    print()
    print(f"  Self-reactance X1 = omega*L1: {solver.omega * solver.L1:.4f} Ohm")
    print(f"  Self-reactance X2 = omega*L2: {solver.omega * solver.L2:.4f} Ohm")
    print(f"  Mutual reactance Xm = omega*M: {solver.omega * solver.M:.4f} Ohm")
    print()

    # Power loss analysis
    print("Power Loss Sources (for I1 = I2 = 1A):")
    P_R1 = result['R1_mOhm'] / 1000  # Watts for 1A
    P_R2 = result['R2_mOhm'] / 1000
    P_Rm = 2 * result['Rm_mOhm'] / 1000  # Factor of 2 for I1*I2 term
    P_total_loss = P_R1 + P_R2 + P_Rm

    print(f"  Loss in R1:  {P_R1*1000:.4f} mW")
    print(f"  Loss in R2:  {P_R2*1000:.4f} mW")
    print(f"  Loss in Rm:  {P_Rm*1000:.4f} mW  (proximity effect)")
    print(f"  Total loss:  {P_total_loss*1000:.4f} mW")
    print(f"  Rm contribution: {P_Rm/P_total_loss*100:.1f}%")
    print()


if __name__ == "__main__":
    print()
    print("WPT Coil Coupling Analysis Demo")
    print("=" * 60)
    print()

    demo_coupling_vs_distance()
    demo_coaxial_loops()
    demo_wpt_system_analysis()
    demo_misalignment()
    demo_different_coil_sizes()
    demo_mutual_resistance()
    demo_impedance_matrix_components()

    print("All demos completed!")
