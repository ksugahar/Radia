"""
Demonstration of the Unified Analysis Framework

This script demonstrates the three analysis types:
1. Static (DC) Analysis
2. Frequency Response Analysis
3. Transient Analysis (PRIMA method)

The example uses a simple rectangular coil PEEC model.

Author: Radia Development Team
Date: 2026-01-15
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add build path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build-msvc'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

# Import analysis framework
from analysis import (
    UnifiedAnalysis,
    step_voltage,
    pulse_voltage,
    sinusoidal_voltage,
)


def create_simple_coil_peec():
    """
    Create a simple PEEC model for a rectangular coil.

    Returns:
        L: Inductance matrix [H]
        R: Resistance vector [Ohm]
    """
    # Coil parameters
    n_turns = 10
    r_coil = 0.05       # Coil radius [m]
    wire_length = 2.0 * np.pi * r_coil * n_turns  # Total wire length [m]
    wire_radius = 0.5e-3  # Wire radius [m]
    sigma = 5.8e7       # Copper conductivity [S/m]

    # Discretize into segments
    n_segments = 20

    # Resistance per segment
    segment_length = wire_length / n_segments
    segment_area = np.pi * wire_radius**2
    R_segment = segment_length / (sigma * segment_area)
    R = np.ones(n_segments) * R_segment

    # Self-inductance per segment (simplified Neumann formula)
    mu_0 = 4.0 * np.pi * 1e-7
    L_self = (mu_0 * segment_length / (2.0 * np.pi)) * (
        np.log(2.0 * segment_length / wire_radius) - 1.0
    )

    # Build inductance matrix
    L = np.zeros((n_segments, n_segments))

    for i in range(n_segments):
        L[i, i] = L_self

        # Mutual inductance (simplified - proportional to 1/distance)
        for j in range(i + 1, n_segments):
            # Approximate distance based on coil geometry
            angle_diff = 2.0 * np.pi * abs(i - j) / n_segments
            distance = 2.0 * r_coil * np.sin(angle_diff / 2.0)
            distance = max(distance, wire_radius)  # Prevent singularity

            # Neumann formula (simplified)
            M_ij = (mu_0 * segment_length**2) / (4.0 * np.pi * distance)
            L[i, j] = M_ij
            L[j, i] = M_ij

    return L, R


def demo_static_analysis(analysis):
    """Demonstrate static (DC) analysis."""
    print("\n" + "="*60)
    print("1. STATIC (DC) ANALYSIS")
    print("="*60)

    result = analysis.static()

    print(f"Success: {result.success}")
    print(f"DC Resistance: {result.resistance:.6f} Ohm")
    print(f"DC Inductance: {result.inductance*1e6:.3f} uH")
    print(f"Computation time: {result.computation_time*1000:.2f} ms")

    return result


def demo_frequency_analysis(analysis):
    """Demonstrate frequency response analysis."""
    print("\n" + "="*60)
    print("2. FREQUENCY RESPONSE ANALYSIS")
    print("="*60)

    # Frequency sweep: 1 Hz to 10 MHz
    frequencies = np.logspace(0, 7, 200)

    result = analysis.frequency_sweep(frequencies)

    print(f"Success: {result.success}")
    print(f"Frequency range: {frequencies[0]:.0f} Hz to {frequencies[-1]/1e6:.1f} MHz")
    print(f"Computation time: {result.computation_time*1000:.2f} ms")

    # Print some key values
    idx_1k = np.argmin(np.abs(frequencies - 1e3))
    idx_100k = np.argmin(np.abs(frequencies - 1e5))
    idx_1M = np.argmin(np.abs(frequencies - 1e6))

    print(f"\nImpedance at selected frequencies:")
    print(f"  @ 1 kHz:   R = {result.resistance[idx_1k]:.4f} Ohm, "
          f"L = {result.inductance[idx_1k]*1e6:.3f} uH")
    print(f"  @ 100 kHz: R = {result.resistance[idx_100k]:.4f} Ohm, "
          f"L = {result.inductance[idx_100k]*1e6:.3f} uH")
    print(f"  @ 1 MHz:   R = {result.resistance[idx_1M]:.4f} Ohm, "
          f"L = {result.inductance[idx_1M]*1e6:.3f} uH")

    # Plot frequency response
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # Impedance magnitude
    ax = axes[0, 0]
    Z_mag = np.abs(result.impedance)
    ax.loglog(frequencies, Z_mag)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('|Z| [Ohm]')
    ax.set_title('Impedance Magnitude')
    ax.grid(True, which='both', alpha=0.3)

    # Resistance
    ax = axes[0, 1]
    ax.semilogx(frequencies, result.resistance)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('R [Ohm]')
    ax.set_title('Resistance (AC)')
    ax.grid(True, which='both', alpha=0.3)

    # Inductance
    ax = axes[1, 0]
    ax.semilogx(frequencies[1:], result.inductance[1:] * 1e6)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('L [uH]')
    ax.set_title('Inductance')
    ax.grid(True, which='both', alpha=0.3)

    # Q factor
    ax = axes[1, 1]
    Q = result.quality_factor
    ax.semilogx(frequencies, Q)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Q')
    ax.set_title('Quality Factor')
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('demo_analysis_frequency_response.png', dpi=150)
    print("\nSaved: demo_analysis_frequency_response.png")

    return result


def demo_transient_analysis(analysis):
    """Demonstrate transient analysis using PRIMA method."""
    print("\n" + "="*60)
    print("3. TRANSIENT ANALYSIS (PRIMA Method)")
    print("="*60)

    # Time parameters
    t_end = 1e-4  # 100 us
    n_points = 500
    time = np.linspace(0, t_end, n_points)

    # Step voltage excitation (1V step at t=0)
    step = step_voltage(amplitude=1.0, t_start=0.0)

    result = analysis.transient(time, step, excitation_type='voltage')

    print(f"Success: {result.success}")
    print(f"Time range: 0 to {t_end*1e6:.0f} us")
    print(f"Computation time: {result.computation_time*1000:.2f} ms")

    # Print some key values
    idx_10 = np.argmin(np.abs(time - t_end*0.1))
    idx_50 = np.argmin(np.abs(time - t_end*0.5))
    idx_90 = np.argmin(np.abs(time - t_end*0.9))

    print(f"\nCurrent at selected times:")
    print(f"  @ 10%% ({time[idx_10]*1e6:.1f} us): i = {result.current[idx_10]*1e3:.3f} mA")
    print(f"  @ 50%% ({time[idx_50]*1e6:.1f} us): i = {result.current[idx_50]*1e3:.3f} mA")
    print(f"  @ 90%% ({time[idx_90]*1e6:.1f} us): i = {result.current[idx_90]*1e3:.3f} mA")

    # Plot transient response
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # Voltage
    ax = axes[0, 0]
    ax.plot(time * 1e6, result.voltage, 'b-', linewidth=1.5)
    ax.set_xlabel('Time [us]')
    ax.set_ylabel('Voltage [V]')
    ax.set_title('Applied Voltage (Step)')
    ax.grid(True, alpha=0.3)

    # Current
    ax = axes[0, 1]
    ax.plot(time * 1e6, result.current * 1e3, 'r-', linewidth=1.5)
    ax.set_xlabel('Time [us]')
    ax.set_ylabel('Current [mA]')
    ax.set_title('Current Response')
    ax.grid(True, alpha=0.3)

    # Power
    ax = axes[1, 0]
    ax.plot(time * 1e6, result.power * 1e3, 'g-', linewidth=1.5)
    ax.set_xlabel('Time [us]')
    ax.set_ylabel('Power [mW]')
    ax.set_title('Instantaneous Power')
    ax.grid(True, alpha=0.3)

    # Energy
    ax = axes[1, 1]
    energy = result.energy
    ax.plot(time * 1e6, energy * 1e6, 'm-', linewidth=1.5)
    ax.set_xlabel('Time [us]')
    ax.set_ylabel('Energy [uJ]')
    ax.set_title('Cumulative Energy')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('demo_analysis_transient_step.png', dpi=150)
    print("\nSaved: demo_analysis_transient_step.png")

    # Additional: sinusoidal excitation
    print("\n--- Sinusoidal Excitation ---")

    freq_exc = 10e3  # 10 kHz
    n_cycles = 5
    t_end_sin = n_cycles / freq_exc
    time_sin = np.linspace(0, t_end_sin, 1000)

    sine = sinusoidal_voltage(amplitude=1.0, frequency=freq_exc)
    result_sin = analysis.transient(time_sin, sine, excitation_type='voltage')

    print(f"Frequency: {freq_exc/1e3:.0f} kHz")
    print(f"Computation time: {result_sin.computation_time*1000:.2f} ms")

    # Plot sinusoidal response
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(time_sin * 1e6, result_sin.voltage, 'b-', label='Voltage [V]', alpha=0.7)
    ax.plot(time_sin * 1e6, result_sin.current * 100, 'r-', label='Current x100 [A]', alpha=0.7)
    ax.set_xlabel('Time [us]')
    ax.set_ylabel('Amplitude')
    ax.set_title(f'Sinusoidal Response @ {freq_exc/1e3:.0f} kHz')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('demo_analysis_transient_sine.png', dpi=150)
    print("Saved: demo_analysis_transient_sine.png")

    return result


def main():
    """Main demonstration function."""
    print("="*60)
    print("Unified Analysis Framework Demo")
    print("="*60)

    # Create PEEC model
    print("\nCreating PEEC model for rectangular coil...")
    L, R = create_simple_coil_peec()
    print(f"Matrix size: {L.shape[0]} x {L.shape[1]}")
    print(f"Total DC resistance: {np.sum(R):.4f} Ohm")
    print(f"Total DC inductance: {np.sum(L)*1e6:.3f} uH")

    # Initialize analysis framework
    analysis = UnifiedAnalysis()
    analysis.set_peec_model(L, R, prima_order=5)

    # Configure skin effect (optional)
    analysis.set_skin_effect(
        sigma=5.8e7,  # Copper
        use_dowell=True,
        conductor_height=1e-3  # 1 mm wire diameter
    )

    # Run demonstrations
    demo_static_analysis(analysis)
    demo_frequency_analysis(analysis)
    demo_transient_analysis(analysis)

    print("\n" + "="*60)
    print("Demo completed successfully!")
    print("="*60)

    plt.show()


if __name__ == "__main__":
    main()
