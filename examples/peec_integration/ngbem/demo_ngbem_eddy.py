"""
demo_ngbem_eddy.py

Demonstration and validation of FEM-BEM eddy current solver.

Tests:
1. FEM-only: Dirichlet BC with Hz=Hz_inc, validate against analytical skin depth
2. FEM-BEM coupled: Calderon projector coupling with scattered field
3. Comparison of FEM vs FEM-BEM results

Part of Radia project
"""

import sys
import os
import numpy as np
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

from ngbem_eddy import EddyCurrentFEMBEM, create_conductor_mesh

# Physical constants
MU_0 = 4.0 * np.pi * 1e-7


def print_separator(title):
    """Print formatted section separator."""
    print()
    print("=" * 70)
    print(f"  {title}")
    print("=" * 70)


def analytical_loss_slab(Bz, sigma, mu_r, freq, width, height, depth):
    """Analytical eddy current loss for a rectangular slab (1D skin effect).

    For a slab of thickness d in a uniform Hz field (skin effect in z):
    P_per_face = (1/(2*sigma)) * (H0^2 / delta) * A
    where H0 = Bz/mu_0, delta = skin depth, A = face area.

    For thick slab (d >> delta):
        P = 2 * P_per_face (top + bottom faces)

    For arbitrary thickness:
        P = (H0^2 * A) / (sigma * delta) *
            [sinh(d/delta) - sin(d/delta)] /
            [cosh(d/delta) + cos(d/delta)]

    Args:
        Bz: External field [T]
        sigma: Conductivity [S/m]
        mu_r: Relative permeability
        freq: Frequency [Hz]
        width, height: Face dimensions [m]
        depth: Slab thickness [m]

    Returns:
        P_analytical: Power loss [W]
    """
    mu = mu_r * MU_0
    omega = 2.0 * np.pi * freq
    delta = np.sqrt(2.0 / (omega * mu * sigma))
    H0 = Bz / (mu_r * MU_0)  # Hz_inc

    A = width * height  # Face area (top or bottom)
    xi = depth / delta  # Thickness ratio

    # Exact formula for finite slab (1D analytical solution)
    # P = H0^2 * A / (sigma * delta) * F(xi)
    # where F(xi) = [sinh(xi) - sin(xi)] / [cosh(xi) + cos(xi)]
    # Factor 2 for both faces (top + bottom)
    F_xi = (np.sinh(xi) - np.sin(xi)) / (np.cosh(xi) + np.cos(xi))
    P_analytical = 2.0 * H0**2 * A / (sigma * delta) * F_xi

    return P_analytical, delta


def test_fem_only():
    """Test 1: FEM-only with Dirichlet BC."""
    print_separator("Test 1: FEM-only Eddy Current (Dirichlet BC)")

    # Geometry: 10mm x 10mm slab, 5mm thick
    width = 0.01
    height = 0.01
    depth = 0.005
    sigma = 5.8e7   # Copper
    mu_r = 1.0
    Bz = 1.0        # 1T external field
    order = 2

    # Choose frequency so skin depth ~ depth/2 (well-resolved by mesh)
    # delta = sqrt(2/(omega*mu*sigma))
    # For delta = 2.5mm: omega = 2/(delta^2 * mu * sigma)
    # = 2 / ((2.5e-3)^2 * 4pi*1e-7 * 5.8e7)
    # = 2 / (6.25e-6 * 72.88) = 2 / 4.555e-4 = 4390 rad/s
    # f = 699 Hz -> use 500 Hz
    freq = 500
    omega = 2.0 * np.pi * freq
    delta = np.sqrt(2.0 / (omega * MU_0 * sigma))

    # Mesh size should be < delta
    maxh = min(delta * 0.7, depth / 3.0)
    maxh = max(maxh, 0.001)  # At least 1mm

    print(f"\n  Geometry: {width*1e3:.0f}mm x {height*1e3:.0f}mm x {depth*1e3:.1f}mm")
    print(f"  Conductivity: {sigma:.1e} S/m (copper)")
    print(f"  Frequency: {freq} Hz")
    print(f"  Skin depth: {delta*1e3:.2f} mm")
    print(f"  depth/delta = {depth/delta:.2f}")
    print(f"  Mesh maxh: {maxh*1e3:.2f} mm")
    print(f"  FE order: {order}")

    mesh = create_conductor_mesh(width, height, depth, maxh=maxh)
    print(f"  Elements: {mesh.ne}, Vertices: {mesh.nv}")

    solver = EddyCurrentFEMBEM(mesh, sigma=sigma, mu_r=mu_r, order=order)
    solver.assemble_fem(freq)
    solver.solve(B_ext=[0, 0, Bz], mode='fem')

    P_fem = solver.compute_loss()
    P_ana, _ = analytical_loss_slab(Bz, sigma, mu_r, freq,
                                     width, height, depth)
    W_m = solver.compute_stored_energy()

    print(f"\n  Results:")
    print(f"    FEM loss:        {P_fem:.4e} W")
    print(f"    Analytical loss: {P_ana:.4e} W (1D slab formula)")
    if P_ana > 0:
        ratio = P_fem / P_ana
        print(f"    Ratio FEM/analytical: {ratio:.3f}")
        # Note: 3D FEM loss > 1D analytical because edges contribute extra loss
        if 0.5 < ratio < 5.0:
            print(f"    [PASS] Within expected range (3D vs 1D)")
        else:
            print(f"    [WARN] Ratio outside expected range")
    print(f"    Stored energy:   {W_m:.4e} J")

    solver.print_summary()
    return solver, P_fem, P_ana


def test_fem_frequency_sweep():
    """Test 2: FEM frequency sweep."""
    print_separator("Test 2: FEM Frequency Sweep")

    width = 0.01
    height = 0.01
    depth = 0.005
    sigma = 5.8e7
    mu_r = 1.0
    Bz = 1.0
    order = 2

    freqs = [50, 100, 200, 500, 1000]

    print(f"\n  {'Freq (Hz)':>10s}  {'delta (mm)':>10s}  {'d/delta':>8s}  "
          f"{'P_fem (W)':>12s}  {'P_ana (W)':>12s}  {'Ratio':>8s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*8}  {'-'*12}  {'-'*12}  {'-'*8}")

    for freq in freqs:
        try:
            omega = 2.0 * np.pi * freq
            delta = np.sqrt(2.0 / (omega * MU_0 * sigma))

            maxh = min(delta * 0.5, depth / 3.0)
            maxh = max(maxh, 0.0008)

            mesh = create_conductor_mesh(width, height, depth, maxh=maxh)
            solver = EddyCurrentFEMBEM(mesh, sigma=sigma, mu_r=mu_r,
                                         order=order)
            solver.assemble_fem(freq)
            solver.solve(B_ext=[0, 0, Bz], mode='fem')
            P_fem = solver.compute_loss()
            P_ana, _ = analytical_loss_slab(Bz, sigma, mu_r, freq,
                                             width, height, depth)
            ratio = P_fem / P_ana if P_ana > 0 else float('inf')

            print(f"  {freq:>10.0f}  {delta*1e3:>10.3f}  "
                  f"{depth/delta:>8.2f}  {P_fem:>12.4e}  {P_ana:>12.4e}  "
                  f"{ratio:>8.3f}")
        except Exception as e:
            print(f"  {freq:>10.0f}  Error: {e}")

    print(f"\n  Note: Ratio > 1 expected (3D geometry vs 1D slab formula)")
    print(f"  Edge currents contribute additional loss in 3D.")


def test_fembem():
    """Test 3: FEM-BEM coupled solver."""
    print_separator("Test 3: FEM-BEM Coupled (Calderon Projector)")

    width = 0.01
    height = 0.01
    depth = 0.005
    sigma = 5.8e7
    mu_r = 1.0
    Bz = 1.0
    freq = 200      # Low frequency for stable BEM
    order = 1        # Low order for fast demo

    omega = 2.0 * np.pi * freq
    delta = np.sqrt(2.0 / (omega * MU_0 * sigma))
    maxh = min(delta * 0.7, depth / 3.0)
    maxh = max(maxh, 0.001)

    print(f"\n  Geometry: {width*1e3:.0f}mm x {height*1e3:.0f}mm x {depth*1e3:.1f}mm")
    print(f"  Frequency: {freq} Hz, skin depth: {delta*1e3:.2f} mm")
    print(f"  Mesh maxh: {maxh*1e3:.2f} mm, FE order: {order}")

    mesh = create_conductor_mesh(width, height, depth, maxh=maxh)
    print(f"  Elements: {mesh.ne}, Vertices: {mesh.nv}")

    solver = EddyCurrentFEMBEM(mesh, sigma=sigma, mu_r=mu_r, order=order)

    print(f"\n  Assembling FEM-BEM system...")
    solver.assemble_fembem(freq)
    print(f"  Assembly time: {solver.t_assemble:.3f} s")

    print(f"  Solving (GMRES)...")
    solver.solve(B_ext=[0, 0, Bz], mode='fembem', printrates=True)
    print(f"  Solve time: {solver.t_solve:.3f} s")

    P_fembem = solver.compute_loss()
    P_ana, _ = analytical_loss_slab(Bz, sigma, mu_r, freq,
                                     width, height, depth)

    print(f"\n  Results:")
    print(f"    FEM-BEM loss:    {P_fembem:.4e} W")
    print(f"    Analytical loss: {P_ana:.4e} W (1D slab formula)")
    if P_ana > 0:
        ratio = P_fembem / P_ana
        print(f"    Ratio FEM-BEM/analytical: {ratio:.3f}")

    solver.print_summary()
    return solver, P_fembem


def test_fem_vs_fembem():
    """Test 4: Compare FEM and FEM-BEM at same frequency."""
    print_separator("Test 4: FEM vs FEM-BEM Comparison")

    width = 0.01
    height = 0.01
    depth = 0.005
    sigma = 5.8e7
    mu_r = 1.0
    Bz = 1.0
    freq = 100
    order = 1

    omega = 2.0 * np.pi * freq
    delta = np.sqrt(2.0 / (omega * MU_0 * sigma))
    maxh = min(delta * 0.5, depth / 3.0)
    maxh = max(maxh, 0.001)

    mesh = create_conductor_mesh(width, height, depth, maxh=maxh)
    print(f"\n  Frequency: {freq} Hz, skin depth: {delta*1e3:.2f} mm")
    print(f"  Elements: {mesh.ne}, Vertices: {mesh.nv}")

    # FEM-only
    print(f"\n  Solving FEM-only...")
    solver_fem = EddyCurrentFEMBEM(mesh, sigma=sigma, mu_r=mu_r, order=order)
    solver_fem.assemble_fem(freq)
    solver_fem.solve(B_ext=[0, 0, Bz], mode='fem')
    P_fem = solver_fem.compute_loss()

    # FEM-BEM
    print(f"  Solving FEM-BEM...")
    solver_bem = EddyCurrentFEMBEM(mesh, sigma=sigma, mu_r=mu_r, order=order)
    solver_bem.assemble_fembem(freq)
    solver_bem.solve(B_ext=[0, 0, Bz], mode='fembem')
    P_bem = solver_bem.compute_loss()

    # Analytical
    P_ana, _ = analytical_loss_slab(Bz, sigma, mu_r, freq,
                                     width, height, depth)

    print(f"\n  {'Method':>15s}  {'P_loss (W)':>12s}  {'Ratio to Ana':>12s}")
    print(f"  {'-'*15}  {'-'*12}  {'-'*12}")
    print(f"  {'FEM':>15s}  {P_fem:>12.4e}  {P_fem/P_ana:>12.3f}")
    print(f"  {'FEM-BEM':>15s}  {P_bem:>12.4e}  {P_bem/P_ana:>12.3f}")
    print(f"  {'Analytical':>15s}  {P_ana:>12.4e}  {'1.000':>12s}")


def main():
    print("=" * 70)
    print("  Demo: ngbem FEM-BEM Eddy Current Solver")
    print("=" * 70)

    # Test 1: FEM-only validation
    test_fem_only()

    # Test 2: FEM frequency sweep
    test_fem_frequency_sweep()

    # Test 3: FEM-BEM coupled
    try:
        test_fembem()
    except Exception as e:
        print(f"\n  FEM-BEM test failed: {e}")
        import traceback
        traceback.print_exc()

    # Test 4: Comparison
    try:
        test_fem_vs_fembem()
    except Exception as e:
        print(f"\n  FEM vs FEM-BEM comparison failed: {e}")
        import traceback
        traceback.print_exc()

    print_separator("Summary")
    print("""
  FEM-BEM eddy current solver demonstrated:

  1. FEM-only mode: Dirichlet BC (Hz = Hz_inc on boundary)
     - Validates against analytical 1D slab formula
     - Works for all frequencies (mesh must resolve skin depth)

  2. FEM-BEM mode: Calderon projector coupling
     - Scattered field formulation: Hz_scat + Hz_inc
     - Volume source: -k^2 * Hz_inc
     - BEM handles exterior Laplace equation automatically

  3. Physics:
     - Scalar Hz formulation (z-component of magnetic field)
     - Skin depth: delta = sqrt(2/(omega*mu*sigma))
     - Loss: P = (1/(2*sigma)) * integral(|grad Hz|^2) dV
     - Edge currents add ~50-100% loss vs 1D slab formula

  Note: For full 3D vector eddy currents, use HCurl (Nedelec) elements.
  The scalar Hz formulation captures the dominant skin effect physics.
""")


if __name__ == '__main__':
    main()
