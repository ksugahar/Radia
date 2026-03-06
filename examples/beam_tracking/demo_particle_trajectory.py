"""
demo_particle_trajectory.py

Demonstrates particle trajectory computation in magnetic fields.

This example shows how to:
1. Create a simple dipole magnet
2. Compute electron trajectory through the magnet
3. Visualize trajectory using matplotlib

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


def create_dipole_magnet(gap=0.020, width=0.100, height=0.050, length=0.200, B0=0.01):
    """
    Create a simple dipole magnet (H-type).

    The beam travels along y-axis. Bz field causes deflection in x.

    Args:
        gap: Gap between poles [m]
        width: Pole width (x) [m]
        height: Pole height (z) [m]
        length: Magnet length along beam (y) [m]
        B0: Approximate field in gap [T]

    Returns:
        Radia object key for the magnet assembly
    """
    # For a dipole with vertical field (Bz), both magnets magnetized in same direction
    # Br = 1.2 T -> M = Br / mu_0 = 954930 A/m
    # Empirical scaling: B0 ~ 0.005T when M = 954930*0.01
    M_approx = 954930 * (B0 / 0.5)

    # Upper pole - magnetized downward (-z)
    upper = rad.ObjRecMag(
        [0, 0, gap/2 + height/2],
        [width, length, height],
        [0, 0, -M_approx]  # Magnetization in -z
    )

    # Lower pole - also magnetized downward (-z) to create field in gap
    lower = rad.ObjRecMag(
        [0, 0, -(gap/2 + height/2)],
        [width, length, height],
        [0, 0, -M_approx]  # Same direction as upper
    )

    # Create container
    magnet = rad.ObjCnt([upper, lower])

    return magnet


def main():
    print("=" * 60)
    print("Particle Trajectory Demo")
    print("=" * 60)

    # Set units to meters (SI) - Policy compliance

    # Create dipole magnet
    print("\nCreating dipole magnet...")
    magnet = create_dipole_magnet(
        gap=0.020,      # 20 mm gap
        width=0.100,    # 100 mm wide poles
        height=0.050,   # 50 mm thick poles
        length=0.200,   # 200 mm long
        B0=0.01         # ~0.01 T field (weak for demo)
    )

    # Check field at center
    B_center = rad.Fld(magnet, 'b', [0, 0, 0])
    print(f"Field at center: Bx={B_center[0]:.4f}, By={B_center[1]:.4f}, Bz={B_center[2]:.4f} T")

    # Compute particle trajectory
    print("\nComputing electron trajectory...")
    energy_gev = 6.0  # 6 GeV electron (typical synchrotron)

    # Initial conditions: [x0, x'0, z0, z'0]
    # Start at y=-0.3m (-300mm), on axis, no angle
    # Units: [m, rad, m, rad]
    initial_cond = [0.0, 0.0, 0.0, 0.0]
    s_range = [-0.300, 0.300]  # y from -0.3 to +0.3 m
    n_points = 301

    trajectory = bt.particle_trajectory(
        magnet,
        energy_gev,
        initial_cond,
        s_range,
        n_points
    )

    print(f"Trajectory shape: {trajectory.shape}")
    print(f"Trajectory columns: [s, x, x', z, z']")

    # Extract trajectory data
    s = trajectory[:, 0]   # Longitudinal position [m]
    x = trajectory[:, 1]   # Horizontal position [m]
    xp = trajectory[:, 2]  # Horizontal angle [rad]
    z = trajectory[:, 3]   # Vertical position [m]
    zp = trajectory[:, 4]  # Vertical angle [rad]

    # Print summary (convert to mm for readability)
    print(f"\nTrajectory summary:")
    print(f"  Entry (s={s[0]*1000:.1f} mm): x={x[0]*1000:.4f} mm, z={z[0]*1000:.4f} mm")
    print(f"  Center (s={s[n_points//2]*1000:.1f} mm): x={x[n_points//2]*1000:.4f} mm, z={z[n_points//2]*1000:.4f} mm")
    print(f"  Exit (s={s[-1]*1000:.1f} mm): x={x[-1]*1000:.4f} mm, z={z[-1]*1000:.4f} mm")
    print(f"  Max horizontal deflection: {np.max(np.abs(x))*1000:.4f} mm")
    print(f"  Exit angle: x'={xp[-1]*1e3:.4f} mrad")

    # Batch trajectory: scan initial x positions
    print("\nComputing batch trajectories (phase space scan)...")
    x_offsets = np.linspace(-0.005, 0.005, 5)  # -5 to +5 mm
    init_conds = [[x_off, 0, 0, 0] for x_off in x_offsets]

    trajectories = bt.batch_trajectory(
        magnet,
        energy_gev,
        init_conds,
        s_range,
        n_points
    )

    print(f"Computed {len(trajectories)} trajectories")

    # Plot if matplotlib available
    if HAS_MATPLOTLIB:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # Convert to mm for plotting
        s_mm = s * 1000
        x_mm = x * 1000
        z_um = z * 1e6
        xp_mrad = xp * 1e3

        # Plot 1: x vs s (horizontal trajectory)
        ax1 = axes[0, 0]
        ax1.plot(s_mm, x_mm, 'b-', linewidth=2)
        ax1.axvspan(-100, 100, alpha=0.2, color='gray', label='Magnet')
        ax1.set_xlabel('s [mm]')
        ax1.set_ylabel('x [mm]')
        ax1.set_title('Horizontal Trajectory')
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        # Plot 2: z vs s (vertical trajectory)
        ax2 = axes[0, 1]
        ax2.plot(s_mm, z_um, 'r-', linewidth=2)  # Convert to um for visibility
        ax2.axvspan(-100, 100, alpha=0.2, color='gray', label='Magnet')
        ax2.set_xlabel('s [mm]')
        ax2.set_ylabel('z [um]')
        ax2.set_title('Vertical Trajectory (should be ~0)')
        ax2.grid(True, alpha=0.3)

        # Plot 3: x' vs s (horizontal angle)
        ax3 = axes[1, 0]
        ax3.plot(s_mm, xp_mrad, 'b-', linewidth=2)
        ax3.axvspan(-100, 100, alpha=0.2, color='gray', label='Magnet')
        ax3.set_xlabel('s [mm]')
        ax3.set_ylabel("x' [mrad]")
        ax3.set_title('Horizontal Angle')
        ax3.grid(True, alpha=0.3)

        # Plot 4: Batch trajectories
        ax4 = axes[1, 1]
        for i, (traj, x_off) in enumerate(zip(trajectories, x_offsets)):
            # traj is in meters
            ax4.plot(traj[:, 0]*1000, traj[:, 1]*1000, label=f'x0={x_off*1000:.1f}mm')
        ax4.axvspan(-100, 100, alpha=0.2, color='gray')
        ax4.set_xlabel('s [mm]')
        ax4.set_ylabel('x [mm]')
        ax4.set_title('Batch Trajectories (Initial x scan)')
        ax4.legend(loc='upper left', fontsize=8)
        ax4.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('particle_trajectory_demo.png', dpi=150)
        print("\nPlot saved to: particle_trajectory_demo.png")
        plt.show()

    print("\nDemo completed successfully!")


if __name__ == '__main__':
    main()
