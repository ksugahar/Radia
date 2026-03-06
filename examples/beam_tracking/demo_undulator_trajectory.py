"""
demo_undulator_trajectory.py

Demonstrates electron trajectory in a planar undulator.

This example shows how to:
1. Create a simple planar undulator (periodic magnet array)
2. Compute electron trajectory through the undulator
3. Analyze trajectory oscillations and radiation parameters

Uses the beam_tracking module for NumPy array output.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src'))

import numpy as np
import radia as rad
from radia import beam_tracking as bt

# Optional: matplotlib for visualization
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("matplotlib not available, skipping plots")


def create_planar_undulator(period=0.050, n_periods=5, gap=0.020, B0=0.8):
    """
    Create a simple planar undulator.

    Args:
        period: Undulator period [m]
        n_periods: Number of periods
        gap: Magnetic gap [m]
        B0: Peak field [T]

    Returns:
        Radia object key for the undulator
    """
    # Estimate magnetization from desired peak field
    # For Halbach-type: B0 ~ 1.7 * Br * exp(-pi*g/lambda)
    # Simplified model: use permanent magnet blocks
    M0 = B0 / (4 * np.pi * 1e-7) * 2e-3  # Rough scaling

    magnets = []

    # Magnet block dimensions (convert from mm ratios)
    block_width = 0.040   # 40 mm
    block_height = 0.030  # 30 mm
    block_length = period / 4  # Quarter period

    total_length = n_periods * period
    y_start = -total_length / 2

    for i in range(n_periods * 4):
        y_center = y_start + block_length / 2 + i * block_length

        # Alternating vertical magnetization
        # Pattern: +z, 0, -z, 0, +z, ...  (simplified Halbach)
        phase = i % 4
        if phase == 0:
            Mz = M0
        elif phase == 2:
            Mz = -M0
        else:
            continue  # Skip horizontal blocks for simplicity

        # Upper array
        upper = rad.ObjRecMag(
            [0, y_center, gap/2 + block_height/2],
            [block_width, block_length, block_height],
            [0, 0, -Mz]  # Field pointing into gap
        )
        magnets.append(upper)

        # Lower array (opposite polarity)
        lower = rad.ObjRecMag(
            [0, y_center, -(gap/2 + block_height/2)],
            [block_width, block_length, block_height],
            [0, 0, Mz]  # Field pointing into gap
        )
        magnets.append(lower)

    return rad.ObjCnt(magnets)


def analyze_trajectory(trajectory, period):
    """
    Analyze undulator trajectory parameters.

    Args:
        trajectory: NumPy array from particle_trajectory()
        period: Undulator period [m]

    Returns:
        dict with analysis results
    """
    s = trajectory[:, 0]
    x = trajectory[:, 1]
    xp = trajectory[:, 2]

    # Find oscillation amplitude
    x_amplitude = (np.max(x) - np.min(x)) / 2

    # Estimate K parameter from trajectory amplitude
    # K = gamma * theta_max, where theta_max ~ max(x')
    # Also K ~ 2*pi*x_max / lambda_u
    K_from_angle = np.max(np.abs(xp)) * 1957  # gamma for 1 GeV ~ 1957
    K_from_amplitude = 2 * np.pi * x_amplitude / period

    return {
        'x_amplitude_um': x_amplitude * 1e6,  # um
        'xp_max_urad': np.max(np.abs(xp)) * 1e6,  # urad
        'K_from_angle': K_from_angle,
        'K_from_amplitude': K_from_amplitude
    }


def main():
    print("=" * 60)
    print("Undulator Trajectory Demo")
    print("=" * 60)

    # Set units to meters (SI) - Policy compliance
    rad.FldUnits('m')

    # Undulator parameters
    period = 0.050      # 50 mm period
    n_periods = 5       # 5 periods
    gap = 0.020         # 20 mm gap
    B0 = 0.5            # 0.5 T peak field

    print(f"\nUndulator parameters:")
    print(f"  Period: {period*1000} mm")
    print(f"  Number of periods: {n_periods}")
    print(f"  Gap: {gap*1000} mm")
    print(f"  Target peak field: {B0} T")

    # Create undulator
    print("\nCreating undulator...")
    undulator = create_planar_undulator(period, n_periods, gap, B0)

    # Check field along axis
    y_test = np.linspace(-n_periods*period/2, n_periods*period/2, 101)
    Bz_axis = [rad.Fld(undulator, 'bz', [0, y, 0]) for y in y_test]
    Bz_peak = np.max(np.abs(Bz_axis))
    print(f"Measured peak field on axis: {Bz_peak:.4f} T")

    # Compute trajectory
    print("\nComputing electron trajectory...")
    energy_gev = 3.0  # 3 GeV electron

    # Start well before undulator, end well after
    margin = 2 * period
    s_start = -n_periods * period / 2 - margin
    s_end = n_periods * period / 2 + margin
    n_points = int((s_end - s_start) / 0.0005) + 1  # 0.5 mm step -> 0.0005 m

    trajectory = bt.particle_trajectory(
        undulator,
        energy_gev,
        [0, 0, 0, 0],  # On axis, no angle
        [s_start, s_end],
        n_points
    )

    print(f"Trajectory: {trajectory.shape[0]} points from s={s_start*1000:.1f} to s={s_end*1000:.1f} mm")

    # Analyze trajectory
    # Note: period is passed in meters, so analysis uses meters
    analysis = analyze_trajectory(trajectory, period)
    print(f"\nTrajectory analysis:")
    print(f"  Oscillation amplitude: {analysis['x_amplitude_um']:.2f} um")
    print(f"  Max angle: {analysis['xp_max_urad']:.2f} urad")
    print(f"  K parameter (from angle): {analysis['K_from_angle']:.3f}")
    print(f"  K parameter (from amplitude): {analysis['K_from_amplitude']:.3f}")

    # Theoretical K parameter
    # Formula: K = 0.934 * lambda_cm * B_tesla
    period_cm = period * 100
    K_theory = 0.934 * period_cm * Bz_peak 
    print(f"  K parameter (theoretical): {K_theory:.3f}")

    # Plot if matplotlib available
    if HAS_MATPLOTLIB:
        fig, axes = plt.subplots(3, 1, figsize=(12, 10))

        s = trajectory[:, 0]
        x = trajectory[:, 1]
        xp = trajectory[:, 2]

        # Convert to mm/um/urad for plotting
        s_mm = s * 1000
        x_um = x * 1e6
        xp_urad = xp * 1e6
        y_test_mm = y_test * 1000

        # Plot 1: Trajectory x vs s
        ax1 = axes[0]
        ax1.plot(s_mm, x_um, 'b-', linewidth=1)
        ax1.axvspan(-n_periods*period*1000/2, n_periods*period*1000/2, alpha=0.2, color='orange', label='Undulator')
        ax1.set_xlabel('s [mm]')
        ax1.set_ylabel('x [um]')
        ax1.set_title(f'Horizontal Trajectory (E={energy_gev} GeV, K~{K_theory:.2f})')
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        # Plot 2: Angle x' vs s
        ax2 = axes[1]
        ax2.plot(s_mm, xp_urad, 'r-', linewidth=1)
        ax2.axvspan(-n_periods*period*1000/2, n_periods*period*1000/2, alpha=0.2, color='orange')
        ax2.set_xlabel('s [mm]')
        ax2.set_ylabel("x' [urad]")
        ax2.set_title('Horizontal Angle')
        ax2.grid(True, alpha=0.3)

        # Plot 3: Field profile
        ax3 = axes[2]
        ax3.plot(y_test_mm, Bz_axis, 'g-', linewidth=2)
        ax3.set_xlabel('y [mm]')
        ax3.set_ylabel('Bz [T]')
        ax3.set_title('Vertical Field on Axis')
        ax3.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('undulator_trajectory_demo.png', dpi=150)
        print("\nPlot saved to: undulator_trajectory_demo.png")
        plt.show()

    print("\nDemo completed successfully!")


if __name__ == '__main__':
    main()
