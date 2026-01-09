"""
RWG-EFIE 3D Induction Heating Demo

This example demonstrates the RWG-EFIE (Rao-Wilton-Glisson Electric Field
Integral Equation) solver for induction heating analysis with 3D surface
elements.

Unlike the segment-based PEEC approach, RWG-EFIE:
- Uses triangular surface elements
- Handles arbitrary 3D coil and workpiece geometries
- Accurately captures eddy current distribution on surfaces

Key Features Demonstrated:
1. Loop coil mesh generation
2. Spiral coil mesh generation
3. Plate workpiece mesh generation
4. Impedance calculation
5. Magnetic field computation

Author: Radia Development Team
Date: 2026-01-08
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src'))

import numpy as np

# Import RWG-EFIE solver
try:
    from radia.rwg_efie_solver import (
        RWGMesh,
        RWGEFIESolver,
        InductionHeatingSolver,
        create_induction_heating_model,
    )
    RWG_AVAILABLE = True
except ImportError as e:
    print(f"Warning: RWG-EFIE solver not available: {e}")
    RWG_AVAILABLE = False


def demo_loop_coil():
    """Demonstrate RWG mesh for a simple loop coil."""
    print("=" * 60)
    print("Demo 1: Loop Coil Mesh and Analysis")
    print("=" * 60)
    print()

    # Create loop coil mesh
    mesh = RWGMesh()
    mesh.create_loop_coil(
        center=[0, 0, 0],
        radius=0.05,           # 50mm radius
        normal=[0, 0, 1],      # Axis along z
        wire_radius=0.002,     # 2mm wire radius
        num_around=8,          # 8 panels around wire
        num_along=24           # 24 panels along loop
    )

    print(f"Loop Coil Mesh Statistics:")
    print(f"  Vertices: {len(mesh.vertices)}")
    print(f"  Triangles: {len(mesh.triangles)}")
    print(f"  Edges: {len(mesh.edges)}")
    print(f"  Interior edges (DOF): {mesh.num_interior_edges}")
    print(f"  Boundary edges: {len(mesh.get_boundary_edges())}")
    print(f"  Is closed surface: {mesh.is_closed()}")
    print()

    # Solve at different frequencies
    frequencies = [10e3, 50e3, 100e3, 200e3]
    print(f"{'Frequency [kHz]':>15} | {'R [mOhm]':>12} | {'L [uH]':>10} | {'Z [Ohm]':>15}")
    print("-" * 60)

    for freq in frequencies:
        solver = RWGEFIESolver(mesh)
        solver.set_frequency(freq)
        solver.set_conductivity(5.8e7)  # Copper

        # Find port edges (for closed loop, use arbitrary pair)
        interior_edges = [i for i, e in enumerate(mesh.edges) if not e.is_boundary]
        if interior_edges:
            solver.define_port([interior_edges[0]])

        solver.set_voltage_excitation(1.0)
        solver.solve()

        Z = solver.get_impedance()
        R = Z.real * 1e3  # mOhm
        L = Z.imag / (2 * np.pi * freq) * 1e6  # uH

        print(f"{freq/1e3:>15.0f} | {R:>12.4f} | {L:>10.3f} | {Z.real:>7.4f}+j{Z.imag:>6.4f}")

    print()


def demo_spiral_coil():
    """Demonstrate RWG mesh for a spiral coil (typical induction heating)."""
    print("=" * 60)
    print("Demo 2: Spiral Coil (Induction Heating Style)")
    print("=" * 60)
    print()

    # Create spiral coil mesh
    mesh = RWGMesh()
    mesh.create_spiral_coil(
        center=[0, 0, 0.02],    # 20mm above origin
        inner_radius=0.02,      # 20mm inner radius
        outer_radius=0.05,      # 50mm outer radius
        pitch=0.005,            # 5mm pitch
        num_turns=3,            # 3 turns
        axis=[0, 0, 1],         # Axis along z
        wire_radius=0.002,      # 2mm wire radius
        num_around=6            # 6 panels around wire
    )

    print(f"Spiral Coil Mesh Statistics:")
    print(f"  Vertices: {len(mesh.vertices)}")
    print(f"  Triangles: {len(mesh.triangles)}")
    print(f"  Edges: {len(mesh.edges)}")
    print(f"  Interior edges (DOF): {mesh.num_interior_edges}")
    print(f"  Boundary edges: {len(mesh.get_boundary_edges())}")
    print(f"  Is closed surface: {mesh.is_closed()}")
    print()

    # Solve at 50 kHz (typical induction heating frequency)
    freq = 50e3
    solver = RWGEFIESolver(mesh)
    solver.set_frequency(freq)
    solver.set_conductivity(5.8e7)

    # Use boundary edges as port
    boundary = mesh.get_boundary_edges()
    if boundary:
        solver.define_port(boundary[:2])
    else:
        interior_edges = [i for i, e in enumerate(mesh.edges) if not e.is_boundary]
        if interior_edges:
            solver.define_port([interior_edges[0]])

    solver.set_voltage_excitation(1.0)
    solver.solve()

    Z = solver.get_impedance()
    print(f"Analysis at {freq/1e3:.0f} kHz:")
    print(f"  Resistance: {Z.real*1e3:.4f} mOhm")
    print(f"  Inductance: {Z.imag/(2*np.pi*freq)*1e6:.3f} uH")
    print(f"  Impedance: {Z.real:.6f} + j{Z.imag:.6f} Ohm")
    print()


def demo_workpiece_mesh():
    """Demonstrate RWG mesh for workpiece (open surface)."""
    print("=" * 60)
    print("Demo 3: Workpiece Mesh (Open Surface)")
    print("=" * 60)
    print()

    # Rectangular plate workpiece
    plate_mesh = RWGMesh()
    plate_mesh.create_rectangular_plate(
        center=[0, 0, 0],
        Lx=0.1,                # 100mm x
        Ly=0.1,                # 100mm y
        normal=[0, 0, 1],
        nx=8,                  # 8 divisions in x
        ny=8                   # 8 divisions in y
    )

    print(f"Rectangular Plate Mesh:")
    print(f"  Vertices: {len(plate_mesh.vertices)}")
    print(f"  Triangles: {len(plate_mesh.triangles)}")
    print(f"  Interior edges: {plate_mesh.num_interior_edges}")
    print(f"  Boundary edges: {len(plate_mesh.get_boundary_edges())}")
    print(f"  Is closed: {plate_mesh.is_closed()}")
    print()

    # Circular disk workpiece
    disk_mesh = RWGMesh()
    disk_mesh.create_circular_disk(
        center=[0, 0, 0],
        radius=0.05,           # 50mm radius
        normal=[0, 0, 1],
        nr=5,                  # 5 radial divisions
        ntheta=16              # 16 angular divisions
    )

    print(f"Circular Disk Mesh:")
    print(f"  Vertices: {len(disk_mesh.vertices)}")
    print(f"  Triangles: {len(disk_mesh.triangles)}")
    print(f"  Interior edges: {disk_mesh.num_interior_edges}")
    print(f"  Boundary edges: {len(disk_mesh.get_boundary_edges())}")
    print(f"  Is closed: {disk_mesh.is_closed()}")
    print()


def demo_high_level_api():
    """Demonstrate high-level InductionHeatingSolver API."""
    print("=" * 60)
    print("Demo 4: High-Level Induction Heating API")
    print("=" * 60)
    print()

    # Create model using convenience function
    solver = create_induction_heating_model(
        coil_type='loop',
        workpiece_type=None,  # Coil only for this demo
        loop_radius=0.05,
        wire_radius=0.002,
        frequency=50e3,
        coil_conductivity=5.8e7,
        num_around=8
    )

    # Solve
    result = solver.solve(verbose=True)

    # Compute B field along axis
    print("B field on axis (above coil):")
    print(f"{'z [mm]':>10} | {'|B| [uT]':>12} | {'Bz_real [uT]':>12}")
    print("-" * 40)

    for z in [0.02, 0.03, 0.05, 0.07, 0.10]:
        B = solver.compute_B_field([0, 0, z])
        B_mag = np.abs(np.linalg.norm(B)) * 1e6  # uT
        Bz_real = B[2].real * 1e6
        print(f"{z*1000:>10.0f} | {B_mag:>12.3f} | {Bz_real:>12.3f}")

    print()


def demo_comparison_analytical():
    """Compare RWG-EFIE with analytical formula for loop coil."""
    print("=" * 60)
    print("Demo 5: Comparison with Analytical Formula")
    print("=" * 60)
    print()

    # Loop parameters
    R = 0.05      # 50mm radius
    a = 0.002     # 2mm wire radius

    # Analytical self-inductance (Neumann formula)
    # L = mu_0 * R * (ln(8R/a) - 2)
    mu_0 = 4 * np.pi * 1e-7
    L_analytical = mu_0 * R * (np.log(8 * R / a) - 2)

    print(f"Loop Coil: R = {R*1000:.0f}mm, wire radius = {a*1000:.1f}mm")
    print()
    print(f"Analytical formula: L = mu_0 * R * (ln(8R/a) - 2)")
    print(f"  L_analytical = {L_analytical*1e6:.4f} uH")
    print()

    # RWG-EFIE calculation at different mesh densities
    print("RWG-EFIE Results:")
    print(f"{'num_around':>12} | {'num_along':>10} | {'L [uH]':>12} | {'Error [%]':>10}")
    print("-" * 50)

    for num_around, num_along in [(6, 16), (8, 24), (10, 32)]:
        mesh = RWGMesh()
        mesh.create_loop_coil(
            center=[0, 0, 0],
            radius=R,
            normal=[0, 0, 1],
            wire_radius=a,
            num_around=num_around,
            num_along=num_along
        )

        solver = RWGEFIESolver(mesh)
        solver.set_frequency(1e6)  # 1 MHz (high enough for L >> R)
        solver.set_conductivity(5.8e7)

        interior_edges = [i for i, e in enumerate(mesh.edges) if not e.is_boundary]
        if interior_edges:
            solver.define_port([interior_edges[0]])

        solver.set_voltage_excitation(1.0)
        solver.solve()

        Z = solver.get_impedance()
        L_rwg = Z.imag / (2 * np.pi * 1e6)
        error = (L_rwg - L_analytical) / L_analytical * 100

        print(f"{num_around:>12} | {num_along:>10} | {L_rwg*1e6:>12.4f} | {error:>10.2f}")

    print()


def demo_induction_heating_setup():
    """Demonstrate typical induction heating coil-workpiece setup."""
    print("=" * 60)
    print("Demo 6: Induction Heating Coil-Workpiece Setup")
    print("=" * 60)
    print()

    print("System Configuration:")
    print("  Coil: 3-turn spiral, R_inner=20mm, R_outer=50mm")
    print("  Wire: 4mm x 2mm rectangular (approximated as 3mm radius)")
    print("  Workpiece: 100mm x 100mm steel plate")
    print("  Gap: 20mm (coil bottom to workpiece)")
    print("  Frequency: 50 kHz")
    print()

    # Create coil mesh
    coil_mesh = RWGMesh()
    coil_mesh.create_spiral_coil(
        center=[0, 0, 0.025],   # 25mm above workpiece (5mm gap to bottom of coil)
        inner_radius=0.02,
        outer_radius=0.05,
        pitch=0.005,
        num_turns=3,
        axis=[0, 0, 1],
        wire_radius=0.003,
        num_around=6
    )

    # Create workpiece mesh
    workpiece_mesh = RWGMesh()
    workpiece_mesh.create_rectangular_plate(
        center=[0, 0, 0],
        Lx=0.1,
        Ly=0.1,
        normal=[0, 0, 1],
        nx=10,
        ny=10
    )

    print(f"Mesh Statistics:")
    print(f"  Coil triangles: {len(coil_mesh.triangles)}")
    print(f"  Coil interior edges: {coil_mesh.num_interior_edges}")
    print(f"  Workpiece triangles: {len(workpiece_mesh.triangles)}")
    print(f"  Workpiece interior edges: {workpiece_mesh.num_interior_edges}")
    print()

    # Solve coil only (full coupling is future work)
    freq = 50e3
    solver = RWGEFIESolver(coil_mesh)
    solver.set_frequency(freq)
    solver.set_conductivity(5.8e7)  # Copper

    boundary = coil_mesh.get_boundary_edges()
    if boundary:
        solver.define_port(boundary[:2])

    solver.set_voltage_excitation(1.0)
    solver.solve()

    Z_coil = solver.get_impedance()
    L_coil = Z_coil.imag / (2 * np.pi * freq)
    R_coil = Z_coil.real

    # Compute B field at workpiece surface
    print(f"Coil Analysis at {freq/1e3:.0f} kHz:")
    print(f"  Coil inductance: {L_coil*1e6:.3f} uH")
    print(f"  Coil resistance: {R_coil*1e3:.4f} mOhm")
    print()

    # B field at workpiece center
    B_workpiece = solver.compute_B_field([0, 0, 0.005])  # 5mm above plate
    B_mag = np.abs(np.linalg.norm(B_workpiece))
    print(f"B field at workpiece (z=5mm above plate):")
    print(f"  |B| = {B_mag*1e3:.4f} mT")
    print()

    # Skin depth in steel at 50 kHz
    sigma_steel = 5e6  # S/m
    mu_r_steel = 100   # Relative permeability
    delta_steel = np.sqrt(2 / (2*np.pi*freq * 4*np.pi*1e-7 * mu_r_steel * sigma_steel)) * 1e3

    print(f"Steel Workpiece Properties:")
    print(f"  Conductivity: {sigma_steel/1e6:.1f} MS/m")
    print(f"  Relative permeability: {mu_r_steel}")
    print(f"  Skin depth at {freq/1e3:.0f} kHz: {delta_steel:.3f} mm")
    print()

    print("Note: Full coil-workpiece coupling is now available using the")
    print("      CoupledEFIESolver class. See demo_coupled_solver() below.")
    print()


def demo_coupled_solver():
    """Demonstrate coupled coil-workpiece EFIE solver."""
    print("=" * 60)
    print("Demo 7: Coupled Coil-Workpiece EFIE Solver")
    print("=" * 60)
    print()

    print("This demo solves the fully coupled coil-workpiece system:")
    print("  [Z_cc  Z_cw] [I_c]   [V_c]")
    print("  [Z_wc  Z_ww] [I_w] = [0  ]")
    print()

    # Use InductionHeatingSolver with coupled mode
    solver = create_induction_heating_model(
        coil_type='loop',
        workpiece_type='plate',
        coil_center=[0, 0, 0.03],       # 30mm above origin
        loop_radius=0.04,                # 40mm radius loop
        wire_radius=0.002,               # 2mm wire radius
        workpiece_center=[0, 0, 0],      # Workpiece at origin
        workpiece_Lx=0.08,               # 80mm x 80mm plate
        workpiece_Ly=0.08,
        frequency=50e3,                  # 50 kHz
        coil_conductivity=5.8e7,         # Copper
        workpiece_conductivity=5e6,      # Steel
        num_around=6,                    # Reduced mesh for faster demo
        nx=6,
        ny=6,
    )

    # Solve with coupling
    print("Solving coupled system (this may take a moment)...")
    result = solver.solve(verbose=True, coupled=True)

    # Compare with uncoupled (coil only)
    print("\nComparison: Coupled vs Uncoupled")
    print("-" * 40)

    result_uncoupled = solver.solve(verbose=False, coupled=False)

    print(f"Coil only:")
    print(f"  Z = {result_uncoupled['resistance']*1e3:.4f} + j{result_uncoupled['impedance_imag']*1e3:.4f} mOhm")
    print(f"  L = {result_uncoupled['inductance']*1e6:.4f} uH")
    print()

    print(f"With workpiece coupling:")
    print(f"  Z = {result['resistance']*1e3:.4f} + j{result['impedance_imag']*1e3:.4f} mOhm")
    print(f"  L = {result['inductance']*1e6:.4f} uH")
    print(f"  P_workpiece = {result['power_workpiece']:.4f} W (per 1V excitation)")
    print()

    # Explain the physics
    print("Physical interpretation:")
    print("  - Resistance increases due to eddy current losses in workpiece")
    print("  - Inductance decreases due to opposing field from workpiece currents")
    print("  - Power absorbed in workpiece = heating power")
    print()


def demo_hybrid_solver():
    """Demonstrate hybrid PEEC + RWG solver."""
    print("=" * 60)
    print("Demo 8: Hybrid PEEC + RWG Solver")
    print("=" * 60)
    print()

    print("The hybrid solver combines:")
    print("  - PEEC (wire segments) for the coil")
    print("  - RWG-EFIE (surface triangles) for the workpiece")
    print()
    print("This is efficient when the coil has simple geometry (wire loops/spirals)")
    print("but the workpiece has complex 3D geometry.")
    print()

    # Import hybrid solver
    try:
        from radia.rwg_efie_solver import create_hybrid_solver, HybridPEECRWGSolver
    except ImportError:
        print("Hybrid solver not available.")
        return

    # Create a hybrid model
    print("Creating hybrid model...")
    solver = create_hybrid_solver(
        coil_type='loop',
        workpiece_type='plate',
        coil_center=[0, 0, 0.02],    # 20mm above workpiece
        coil_radius=0.04,             # 40mm radius loop
        wire_radius=0.002,            # 2mm wire
        coil_n_segments=24,           # 24 PEEC segments
        workpiece_center=[0, 0, 0],
        workpiece_Lx=0.08,            # 80mm x 80mm plate
        workpiece_Ly=0.08,
        nx=6,                         # 6x6 mesh (72 triangles)
        ny=6,
        frequency=50e3,               # 50 kHz
    )

    print(f"  Coil: {len(solver.coil_segments)} PEEC segments")
    print(f"  Workpiece mesh: {len(solver.workpiece_mesh.triangles)} triangles")
    print()

    # Solve
    import time
    t_start = time.time()
    result = solver.solve(verbose=True)
    t_hybrid = time.time() - t_start

    print(f"Solution time: {t_hybrid:.2f} seconds")
    print()

    # Compare with full RWG solver
    print("Comparison with full RWG-EFIE solver...")

    solver_rwg = create_induction_heating_model(
        coil_type='loop',
        workpiece_type='plate',
        coil_center=[0, 0, 0.02],
        loop_radius=0.04,
        wire_radius=0.002,
        workpiece_center=[0, 0, 0],
        workpiece_Lx=0.08,
        workpiece_Ly=0.08,
        frequency=50e3,
        num_around=6,
        num_along=24,
        nx=6,
        ny=6,
    )

    t_start = time.time()
    result_rwg = solver_rwg.solve(verbose=False, coupled=True)
    t_rwg = time.time() - t_start

    print(f"\n  Full RWG-EFIE:")
    print(f"    Z = {result_rwg['resistance']*1e3:.4f} + j{result_rwg['impedance_imag']*1e3:.4f} mOhm")
    print(f"    Time: {t_rwg:.2f} s")
    print()
    print(f"  Hybrid PEEC+RWG:")
    print(f"    Z = {result['resistance']*1e3:.4f} + j{result['reactance']*1e3:.4f} mOhm")
    print(f"    Time: {t_hybrid:.2f} s")
    print()

    print("Note: Hybrid solver is faster when coil mesh would be dense.")
    print("      PEEC segments are computationally cheap compared to RWG triangles.")
    print()


def demo_pfft_solver():
    """Demonstrate pFFT-accelerated solver for large problems."""
    print("=" * 60)
    print("Demo 9: pFFT-Accelerated EFIE Solver")
    print("=" * 60)
    print()

    print("The pFFT (pre-corrected FFT) method provides O(N log N) matrix-vector")
    print("products instead of O(N^2), enabling solution of large problems.")
    print()

    # Import pFFT solver
    try:
        from radia.rwg_efie_solver import pFFTEFIESolver
    except ImportError:
        print("pFFT solver not available.")
        return

    # Create a moderate-sized mesh for demonstration
    print("Creating spiral coil mesh...")
    mesh = RWGMesh()
    mesh.create_spiral_coil(
        center=[0, 0, 0],
        inner_radius=0.02,
        outer_radius=0.05,
        pitch=0.005,
        num_turns=3,
        axis=[0, 0, 1],
        wire_radius=0.002,
        num_around=8
    )

    print(f"  Triangles: {len(mesh.triangles)}")
    print(f"  Interior edges (DOF): {mesh.num_interior_edges}")
    print()

    # Solve with pFFT
    solver = pFFTEFIESolver(mesh)
    solver.set_frequency(50e3)  # 50 kHz
    solver.set_conductivity(5.8e7)  # Copper

    # Find port
    interior_edges = [i for i, e in enumerate(mesh.edges) if not e.is_boundary]
    if interior_edges:
        solver.define_port([interior_edges[0]])

    solver.set_voltage_excitation(1.0)

    import time
    t_start = time.time()
    result = solver.solve(verbose=True)
    t_end = time.time()

    print(f"Solution time: {t_end - t_start:.2f} seconds")
    print()

    # Compare with direct solver
    print("Comparison with direct solver...")
    direct_solver = RWGEFIESolver(mesh)
    direct_solver.set_frequency(50e3)
    direct_solver.set_conductivity(5.8e7)
    if interior_edges:
        direct_solver.define_port([interior_edges[0]])
    direct_solver.set_voltage_excitation(1.0)

    t_start = time.time()
    direct_solver.solve()
    t_direct = time.time() - t_start

    Z_direct = direct_solver.get_impedance()
    Z_pfft = result['impedance']

    print(f"  Direct solver:")
    print(f"    Z = {Z_direct.real*1e3:.4f} + j{Z_direct.imag*1e3:.4f} mOhm")
    print(f"    Time: {t_direct:.2f} s")
    print()
    print(f"  pFFT solver:")
    print(f"    Z = {Z_pfft.real*1e3:.4f} + j{Z_pfft.imag*1e3:.4f} mOhm")
    print(f"    Time: {t_end - t_start:.2f} s")
    print()

    # Relative error
    error = abs(Z_pfft - Z_direct) / abs(Z_direct) * 100
    print(f"  Relative difference: {error:.2f}%")
    print()

    print("Note: pFFT becomes more advantageous for larger problems (N > 1000 DOF)")
    print("where the O(N log N) scaling provides significant speedup.")
    print()


if __name__ == "__main__":
    print()
    print("RWG-EFIE 3D Induction Heating Demo")
    print("=" * 60)
    print()

    if not RWG_AVAILABLE:
        print("RWG-EFIE solver not available. Please check dependencies.")
        sys.exit(1)

    # Run demos
    demo_loop_coil()
    demo_spiral_coil()
    demo_workpiece_mesh()
    demo_high_level_api()
    demo_comparison_analytical()
    demo_induction_heating_setup()
    demo_coupled_solver()
    demo_hybrid_solver()  # Hybrid PEEC + RWG solver
    demo_pfft_solver()    # pFFT-accelerated solver

    print("All demos completed!")
