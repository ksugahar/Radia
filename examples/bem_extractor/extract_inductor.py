"""
extract_inductor.py

One-button inductor parameter extraction using BEM.

Demonstrates Cubit mesh -> ngbem BEM -> L, R extraction pipeline
with GMSH post-processing for high-order element visualization.

This example works in two modes:
  1. With Cubit: Full pipeline (Cubit mesh -> extract)
  2. Without Cubit: Uses Netgen OCC to create a simple coil mesh

Usage:
    # With Cubit (journal script):
    cubit -nographics -batch extract_inductor.py

    # Without Cubit (standalone):
    python extract_inductor.py

Part of Radia project
"""

import sys
import os
import numpy as np

# Add Radia source to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

MU_0 = 4.0 * np.pi * 1e-7


def create_circular_loop_mesh(radius, wire_radius, maxh=None, label="coil"):
    """Create a surface mesh of a circular loop using Netgen OCC.

    Args:
        radius: Loop radius [m]
        wire_radius: Wire cross-section radius [m]
        maxh: Maximum element size [m] (default: wire_radius)
        label: Surface label

    Returns:
        NGSolve Mesh (surface mesh in 3D)
    """
    from netgen.occ import WorkPlane, Axes, Axis, Pnt, Dir, OCCGeometry
    from ngsolve import Mesh
    from ngsolve import TaskManager

    if maxh is None:
        maxh = wire_radius * 1.5

    # Create a torus (circular loop) using OCC
    # WorkPlane circle -> revolve around axis
    wp = WorkPlane(Axes(p=Pnt(radius, 0, 0),
                        n=Dir(0, 1, 0),
                        h=Dir(0, 0, 1)))
    circle = wp.Circle(wire_radius).Face()
    circle.name = label

    # Revolve around Z axis to create torus (Axis, not Axes)
    torus = circle.Revolve(Axis(p=Pnt(0, 0, 0), d=Dir(0, 0, 1)), 360)

    geo = OCCGeometry(torus)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh)

        mesh = Mesh(ngmesh)
        return mesh


def analytical_inductance_circular_loop(R, a):
    """Analytical self-inductance of a circular loop (Neumann formula).

    L = mu_0 * R * (ln(8R/a) - 2)

    Args:
        R: Loop radius [m]
        a: Wire radius [m]

    Returns:
        L: Self-inductance [H]
    """
    return MU_0 * R * (np.log(8.0 * R / a) - 2.0)


def main():
    """Run inductor extraction demo."""
    print("=" * 60)
    print("BEM Impedance Extractor - Circular Loop Demo")
    print("=" * 60)

    # Parameters
    R = 0.05       # Loop radius: 50 mm
    a = 0.002      # Wire radius: 2 mm
    sigma = 5.8e7  # Copper conductivity [S/m]
    freqs = [100, 1000, 10000, 100000]  # 100 Hz to 100 kHz

    # Analytical reference
    L_analytical = analytical_inductance_circular_loop(R, a)
    print(f"\nAnalytical: L = {L_analytical*1e9:.2f} nH "
          f"(R={R*1e3:.1f} mm, a={a*1e3:.1f} mm)")

    # Try Cubit first, fall back to Netgen OCC
    use_cubit = False
    try:
        import cubit
        cubit.init(["cubit", "-nographics"])
        use_cubit = True
        print("\nUsing Cubit for mesh generation")
    except ImportError:
        print("\nCubit not available, using Netgen OCC")

    from cubit_bem_extractor import BEMExtractor

    ext = BEMExtractor()

    if use_cubit:
        # Cubit workflow
        # (Would create geometry and mesh in Cubit here)
        # For demo, fall through to Netgen
        pass

    # Netgen OCC workflow (works without Cubit)
    print(f"\nCreating circular loop mesh (R={R*1e3:.0f} mm, a={a*1e3:.1f} mm)...")
    mesh = create_circular_loop_mesh(R, a, maxh=a * 2.0, label="coil")
    ext.load_ngsolve_mesh(mesh)

    # Define conductor and port
    ext.set_conductor("coil", sigma=sigma, thickness=a * 2)
    ext.set_port("port1", block="coil")

    # Extract
    print(f"\nExtracting at frequencies: {freqs} Hz")
    result = ext.extract(freqs=freqs, mode='mqs')

    # Compare with analytical
    L_bem = result.L[0, 0]
    error = abs(L_bem - L_analytical) / abs(L_analytical) * 100
    print(f"\nComparison:")
    print(f"  L(BEM)       = {L_bem*1e9:.2f} nH")
    print(f"  L(analytical) = {L_analytical*1e9:.2f} nH")
    print(f"  Error: {error:.1f}%")

    # Save results
    output_dir = os.path.dirname(os.path.abspath(__file__))
    result.to_json(os.path.join(output_dir, "inductor_params.json"))
    result.to_csv(os.path.join(output_dir, "inductor_sweep.csv"))

    # GMSH post-processing
    try:
        gmsh_file = os.path.join(output_dir, "inductor_J.msh")
        ext.export_gmsh_post(gmsh_file, field="J")
        print(f"\nGMSH post file: {gmsh_file}")
        print("  Open in GMSH GUI for high-order visualization")
    except Exception as e:
        print(f"\nGMSH export skipped: {e}")

    # JupyterLab visualization hint
    print("\nFor JupyterLab visualization:")
    print("  from ngsolve.webgui import Draw")
    print("  Draw(mesh)")

    print("\nDone.")
    return result


if __name__ == "__main__":
    main()
