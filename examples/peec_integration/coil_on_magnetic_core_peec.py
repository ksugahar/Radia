"""
Coil on Magnetic Core - CplMag Solver Analysis

This example demonstrates the CplMag solver (Coupled Conductor + Magnetic Material)
using PEEC + Magnetic Moment Method for analyzing a coil wound on a magnetic core.

Features demonstrated:
- Creating a circular loop coil conductor
- Creating a magnetic core (hexahedral magnet with soft iron material)
- Setting complex permeability for magnetic losses
- Coupled PEEC+MMM solution
- Frequency sweep for impedance analysis
- Power loss computation (conductor + magnetic core)

Physical model:
- Circular loop coil (copper conductor)
- Cylindrical/rectangular magnetic core with complex permeability
- mu = mu' - j*mu" for magnetic losses in the core
- ESIM surface impedance for conductor skin effect

Output:
- Impedance vs frequency (comparison with/without core)
- Power loss breakdown (conductor vs core losses)
- Inductance enhancement due to magnetic core

Physics background:
- Coil field induces magnetization in the core
- Core magnetization enhances flux linkage -> increased inductance
- Core magnetic loss (mu") -> increased resistance
- Frequency-dependent skin effect in conductor

Part of Radia CplMag solver examples.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Add Radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

# Physical constants
MU_0 = 4 * np.pi * 1e-7  # H/m
SIGMA_COPPER = 5.8e7     # S/m


def skin_depth(freq, sigma, mu_r=1.0):
    """Calculate skin depth [m]."""
    if freq <= 0 or sigma <= 0:
        return float('inf')
    omega = 2 * np.pi * freq
    return np.sqrt(2 / (omega * MU_0 * mu_r * sigma))


def analytical_loop_inductance(radius, wire_radius):
    """
    Analytical inductance of a circular loop in air.
    L = mu_0 * R * (ln(8*R/a) - 2)
    """
    return MU_0 * radius * (np.log(8 * radius / wire_radius) - 2)


def run_peec_only_analysis():
    """
    Run PEEC analysis for coil without magnetic core (baseline).
    """
    import radia as rad

    rad.UtiDelAll()
    rad.FldUnits('m')

    print("=" * 60)
    print("PEEC Analysis: Coil WITHOUT Magnetic Core (Baseline)")
    print("=" * 60)

    # Coil parameters
    loop_radius = 0.05     # 50 mm
    wire_width = 2e-3      # 2 mm (rectangular wire)
    wire_height = 2e-3     # 2 mm
    sigma = SIGMA_COPPER

    # Equivalent wire radius for analytical comparison
    wire_radius = np.sqrt(wire_width * wire_height / np.pi)

    # Create circular loop conductor
    coil = rad.CndLoop(
        [0, 0, 0],          # center
        loop_radius,         # radius
        [0, 0, 1],          # normal (z-axis)
        'r',                # rectangular cross-section
        wire_width,
        wire_height,
        sigma,
        8,                  # num_panels_around
        36                  # num_panels_loop
    )

    # Frequency range
    frequencies = np.logspace(2, 6, 20)  # 100 Hz to 1 MHz

    # Results storage
    results = {
        'freq': frequencies,
        'R': [],
        'L': [],
        'Q': [],
        'delta': []
    }

    print(f"\nCoil: R={loop_radius*1000:.1f}mm, wire={wire_width*1000:.1f}x{wire_height*1000:.1f}mm")
    print(f"Analytical L (air): {analytical_loop_inductance(loop_radius, wire_radius)*1e9:.2f} nH")
    print()

    for freq in frequencies:
        rad.CndSetFrequency(coil, freq)
        rad.CndSetVoltage(coil, 1.0, 0.0)  # 1V excitation
        rad.CndSolve(coil)

        Z = rad.CndGetImpedance(coil)
        R = Z.real
        X = Z.imag
        omega = 2 * np.pi * freq

        L = X / omega if omega > 0 else 0
        Q = X / R if R > 0 else 0
        delta = skin_depth(freq, sigma)

        results['R'].append(R)
        results['L'].append(L)
        results['Q'].append(Q)
        results['delta'].append(delta)

    # Convert to arrays
    for key in ['R', 'L', 'Q', 'delta']:
        results[key] = np.array(results[key])

    # Print summary
    print(f"{'Freq [Hz]':>12} {'R [mOhm]':>12} {'L [nH]':>12} {'Q':>10} {'delta [mm]':>12}")
    print("-" * 60)
    for i in [0, len(frequencies)//4, len(frequencies)//2, 3*len(frequencies)//4, -1]:
        f = results['freq'][i]
        R = results['R'][i] * 1000
        L = results['L'][i] * 1e9
        Q = results['Q'][i]
        d = results['delta'][i] * 1000
        print(f"{f:>12.0f} {R:>12.4f} {L:>12.3f} {Q:>10.1f} {d:>12.4f}")

    return results


def run_cplmag_analysis():
    """
    Run CplMag analysis for coil with magnetic core.

    NOTE: The CplMag solver is currently a stub implementation.
    The Python API bindings are in place, but the internal coupling matrix
    assembly and solve logic is incomplete. This example demonstrates the
    API usage pattern for when the full implementation is available.

    For now, the example falls back to demonstrating the expected workflow.
    """
    import radia as rad

    rad.UtiDelAll()
    rad.FldUnits('m')

    print("\n" + "=" * 60)
    print("CplMag Analysis: Coil WITH Magnetic Core")
    print("(NOTE: Currently using stub implementation - API demonstration)")
    print("=" * 60)

    # Coil parameters
    loop_radius = 0.05     # 50 mm
    wire_width = 2e-3      # 2 mm
    wire_height = 2e-3     # 2 mm
    sigma = SIGMA_COPPER

    # Create circular loop conductor
    coil = rad.CndLoop(
        [0, 0, 0],
        loop_radius,
        [0, 0, 1],
        'r',
        wire_width,
        wire_height,
        sigma,
        8,
        36
    )

    # Magnetic core parameters (cylindrical core inside the coil)
    core_height = 0.02      # 20 mm
    core_radius = 0.03      # 30 mm (smaller than coil radius)
    mu_r_real = 1000        # Real part of relative permeability
    mu_r_imag = 50          # Imaginary part (magnetic loss)

    # Create magnetic core as a hexahedral magnet (approximate cylinder with hex)
    # Using a rectangular block for simplicity
    core_Lx = 2 * core_radius * 0.8  # Approximate cylinder with square
    core_Ly = 2 * core_radius * 0.8
    core_Lz = core_height

    vertices = [
        [-core_Lx/2, -core_Ly/2, -core_Lz/2],
        [core_Lx/2, -core_Ly/2, -core_Lz/2],
        [core_Lx/2, core_Ly/2, -core_Lz/2],
        [-core_Lx/2, core_Ly/2, -core_Lz/2],
        [-core_Lx/2, -core_Ly/2, core_Lz/2],
        [core_Lx/2, -core_Ly/2, core_Lz/2],
        [core_Lx/2, core_Ly/2, core_Lz/2],
        [-core_Lx/2, core_Ly/2, core_Lz/2],
    ]

    # Create core as soft iron (zero initial magnetization)
    core = rad.ObjHexahedron(vertices, [0, 0, 0])

    # Apply linear material
    mat_core = rad.MatLin(mu_r_real)  # Soft iron with mu_r=1000
    rad.MatApl(core, mat_core)

    print(f"\nCoil: R={loop_radius*1000:.1f}mm, wire={wire_width*1000:.1f}x{wire_height*1000:.1f}mm")
    print(f"Core: {core_Lx*1000:.1f}x{core_Ly*1000:.1f}x{core_Lz*1000:.1f}mm, mu_r={mu_r_real}")
    print(f"Complex mu: mu' = {mu_r_real}, mu'' = {mu_r_imag}")
    print()

    # Create CplMag solver
    solver = rad.CplMagCreate(coil, core)

    # Set complex permeability for magnetic loss
    rad.CplMagSetMu(solver, mu_r_real, mu_r_imag)

    # Frequency range
    frequencies = np.logspace(2, 6, 20)  # 100 Hz to 1 MHz

    # Results storage
    results = {
        'freq': frequencies,
        'R': [],
        'L': [],
        'Q': [],
        'P_cond': [],
        'P_mag': []
    }

    print(f"{'Freq [Hz]':>12} {'R [mOhm]':>12} {'L [nH]':>12} {'P_cond [W]':>12} {'P_mag [W]':>12}")
    print("-" * 72)

    for freq in frequencies:
        rad.CplMagSetFrequency(solver, freq)
        rad.CplMagSetVoltage(solver, 1.0, 0.0)  # 1V excitation

        # Solve coupled system
        result = rad.CplMagSolve(solver)

        Z = result['Z']
        R = Z.real
        X = Z.imag
        omega = 2 * np.pi * freq

        L = X / omega if omega > 0 else 0
        Q = X / R if R > 0 else 0
        P_cond = result['P_conductor']
        P_mag = result['P_magnet']

        results['R'].append(R)
        results['L'].append(L)
        results['Q'].append(Q)
        results['P_cond'].append(P_cond)
        results['P_mag'].append(P_mag)

        if freq in [100, 1000, 10000, 100000, 1000000]:
            print(f"{freq:>12.0f} {R*1000:>12.4f} {L*1e9:>12.3f} {P_cond:>12.6f} {P_mag:>12.6f}")

    # Convert to arrays
    for key in ['R', 'L', 'Q', 'P_cond', 'P_mag']:
        results[key] = np.array(results[key])

    # Clean up
    rad.CplMagDelete(solver)

    return results


def compare_with_without_core():
    """
    Compare coil characteristics with and without magnetic core.
    """
    # Run both analyses
    results_no_core = run_peec_only_analysis()
    results_with_core = run_cplmag_analysis()

    # Create comparison plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    freq = results_no_core['freq']

    # Plot 1: Inductance vs frequency
    ax1 = axes[0, 0]
    ax1.semilogx(freq, results_no_core['L'] * 1e9, 'b-', linewidth=2, label='Without core')
    ax1.semilogx(freq, results_with_core['L'] * 1e9, 'r-', linewidth=2, label='With core')
    ax1.set_xlabel('Frequency [Hz]')
    ax1.set_ylabel('Inductance [nH]')
    ax1.set_title('Inductance vs Frequency')
    ax1.legend()
    ax1.grid(True, which='both', alpha=0.3)

    # Plot 2: Resistance vs frequency
    ax2 = axes[0, 1]
    ax2.loglog(freq, results_no_core['R'] * 1000, 'b-', linewidth=2, label='Without core')
    ax2.loglog(freq, results_with_core['R'] * 1000, 'r-', linewidth=2, label='With core')
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('Resistance [mOhm]')
    ax2.set_title('Resistance vs Frequency')
    ax2.legend()
    ax2.grid(True, which='both', alpha=0.3)

    # Plot 3: Q factor vs frequency
    ax3 = axes[1, 0]
    ax3.semilogx(freq, results_no_core['Q'], 'b-', linewidth=2, label='Without core')
    ax3.semilogx(freq, results_with_core['Q'], 'r-', linewidth=2, label='With core')
    ax3.set_xlabel('Frequency [Hz]')
    ax3.set_ylabel('Q Factor')
    ax3.set_title('Quality Factor vs Frequency')
    ax3.legend()
    ax3.grid(True, which='both', alpha=0.3)

    # Plot 4: Inductance ratio (enhancement due to core)
    ax4 = axes[1, 1]
    L_ratio = results_with_core['L'] / results_no_core['L']
    ax4.semilogx(freq, L_ratio, 'g-', linewidth=2)
    ax4.set_xlabel('Frequency [Hz]')
    ax4.set_ylabel('L(with core) / L(air)')
    ax4.set_title('Inductance Enhancement Ratio')
    ax4.axhline(y=1, color='k', linestyle='--', alpha=0.5, label='Air reference')
    ax4.grid(True, which='both', alpha=0.3)
    ax4.legend()

    plt.tight_layout()
    plt.savefig('coil_on_magnetic_core_peec.png', dpi=150)
    print("\nSaved: coil_on_magnetic_core_peec.png")
    plt.close()

    # Print summary
    print("\n" + "=" * 60)
    print("Summary Comparison")
    print("=" * 60)

    idx_1k = np.argmin(np.abs(freq - 1000))
    idx_100k = np.argmin(np.abs(freq - 100000))

    print(f"\n{'Parameter':30} {'@1kHz (no core)':>18} {'@1kHz (with core)':>18}")
    print("-" * 70)
    print(f"{'Inductance [nH]':30} {results_no_core['L'][idx_1k]*1e9:>18.2f} {results_with_core['L'][idx_1k]*1e9:>18.2f}")
    print(f"{'Resistance [mOhm]':30} {results_no_core['R'][idx_1k]*1000:>18.4f} {results_with_core['R'][idx_1k]*1000:>18.4f}")
    print(f"{'Q Factor':30} {results_no_core['Q'][idx_1k]:>18.1f} {results_with_core['Q'][idx_1k]:>18.1f}")

    print(f"\n{'Parameter':30} {'@100kHz (no core)':>18} {'@100kHz (with core)':>18}")
    print("-" * 70)
    print(f"{'Inductance [nH]':30} {results_no_core['L'][idx_100k]*1e9:>18.2f} {results_with_core['L'][idx_100k]*1e9:>18.2f}")
    print(f"{'Resistance [mOhm]':30} {results_no_core['R'][idx_100k]*1000:>18.4f} {results_with_core['R'][idx_100k]*1000:>18.4f}")
    print(f"{'Q Factor':30} {results_no_core['Q'][idx_100k]:>18.1f} {results_with_core['Q'][idx_100k]:>18.1f}")

    return results_no_core, results_with_core


def main():
    """Main function."""
    print("=" * 60)
    print("Coil on Magnetic Core - CplMag Solver Analysis")
    print("=" * 60)
    print()

    try:
        # Run comparison
        compare_with_without_core()
        print("\n" + "=" * 60)
        print("Analysis complete!")
        print("=" * 60)
    except AttributeError as e:
        if 'CplMagCreate' in str(e):
            print("\nNote: CplMag solver API not available.")
            print("Running PEEC-only analysis...")
            run_peec_only_analysis()
        else:
            raise


if __name__ == '__main__':
    main()
