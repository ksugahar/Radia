"""
ESIM + RWG-EFIE Nonlinear Coupled Solver Demo

This script demonstrates the integration of ESIM (Effective Surface Impedance Method)
with RWG-EFIE (Rao-Wilton-Glisson Electric Field Integral Equation) for analyzing
induction heating of ferromagnetic workpieces with nonlinear B-H curves.

The demo:
1. Creates an ESI table from a steel B-H curve
2. Sets up a coil mesh (source) and workpiece mesh (nonlinear target)
3. Runs the nonlinear coupled solver with iterative convergence
4. Computes power losses and impedance

Key Features:
- Nonlinear surface impedance Z(H) from cell problem solution
- Coupled coil-workpiece EFIE system
- Iterative solver with underrelaxation for nonlinear convergence
- Field-dependent material properties

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

try:
    from rwg_efie_solver import (
        RWGMesh,
        CoupledEFIESolver,
        NonlinearCoupledEFIESolver,
        create_nonlinear_solver,
        create_induction_heating_model_with_esim,
    )
    from esim_cell_problem import (
        ESIMCellProblemSolver,
        BHCurveInterpolator,
        ESITable,
        generate_esi_table_from_bh_curve,
    )
    RWG_AVAILABLE = True
except ImportError as e:
    print(f"RWG-EFIE solver not available: {e}")
    RWG_AVAILABLE = False


def print_header(title):
    """Print a formatted header."""
    print()
    print("=" * 70)
    print(title)
    print("=" * 70)


def print_subheader(title):
    """Print a formatted subheader."""
    print()
    print("-" * 50)
    print(title)
    print("-" * 50)


# Standard steel B-H curve (typical electrical steel)
STEEL_BH_CURVE = [
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


def demo_esi_table_generation():
    """Demonstrate ESI table generation from B-H curve."""
    print_header("1. ESI Table Generation from B-H Curve")

    sigma = 2e6  # S/m (hot steel)
    freq = 50000  # 50 kHz

    print(f"Material: Steel (nonlinear)")
    print(f"Conductivity: {sigma/1e6:.1f} MS/m")
    print(f"Frequency: {freq/1000:.0f} kHz")
    print()

    # Generate ESI table
    print("Generating ESI table by solving cell problems...")
    H_values = [100, 500, 1000, 2500, 5000, 10000]

    esi_table = ESITable(
        bh_curve=STEEL_BH_CURVE,
        sigma=sigma,
        frequency=freq,
        H_values=H_values
    )

    # Display table
    print()
    print(f"{'H0 [A/m]':>12} {'Re(Z) [mOhm]':>14} {'Im(Z) [mOhm]':>14} {'|Z| [mOhm]':>14}")
    print("-" * 60)

    for H in H_values:
        Z = esi_table.get_impedance(H)
        Z_mOhm = Z * 1000  # Convert to mOhm
        print(f"{H:>12.0f} {Z_mOhm.real:>14.4f} {Z_mOhm.imag:>14.4f} {abs(Z_mOhm):>14.4f}")

    print()
    print("Note: Re(Z) represents losses, Im(Z) represents field penetration")
    print("      Z increases with H due to reduced permeability at saturation")

    return esi_table


def demo_mesh_creation():
    """Demonstrate mesh creation for coil and workpiece."""
    print_header("2. Mesh Creation for Coupled Problem")

    # Create coil mesh (pancake coil above workpiece)
    print_subheader("2.1 Coil Mesh (Pancake Coil)")

    coil_mesh = RWGMesh()

    # Pancake coil parameters
    R_inner = 0.02   # 20 mm inner radius
    R_outer = 0.05   # 50 mm outer radius
    z_coil = 0.005   # 5 mm above workpiece
    n_turns = 5      # Number of turns (represented by rings)

    # Create multiple concentric rings
    vertices = []
    triangles = []
    n_theta = 16  # Angular divisions per ring

    for i_turn in range(n_turns):
        R = R_inner + (R_outer - R_inner) * i_turn / (n_turns - 1) if n_turns > 1 else (R_inner + R_outer) / 2

        # Add ring vertices
        ring_start = len(vertices)
        for j in range(n_theta):
            theta = 2 * np.pi * j / n_theta
            x = R * np.cos(theta)
            y = R * np.sin(theta)
            vertices.append([x, y, z_coil])

        # Add ring center
        center_idx = len(vertices)
        vertices.append([0, 0, z_coil])

        # Create triangles from center to ring
        for j in range(n_theta):
            j_next = (j + 1) % n_theta
            triangles.append([center_idx, ring_start + j, ring_start + j_next])

    coil_mesh.create_from_triangles(vertices, triangles)

    print(f"Coil mesh created:")
    print(f"  Inner radius: {R_inner*1000:.0f} mm")
    print(f"  Outer radius: {R_outer*1000:.0f} mm")
    print(f"  Height above workpiece: {z_coil*1000:.0f} mm")
    print(f"  Number of vertices: {coil_mesh.num_vertices}")
    print(f"  Number of triangles: {coil_mesh.num_triangles}")
    print(f"  Number of edges (basis functions): {coil_mesh.num_edges}")

    # Create workpiece mesh (rectangular plate)
    print_subheader("2.2 Workpiece Mesh (Rectangular Plate)")

    workpiece_mesh = RWGMesh()

    # Workpiece parameters
    Lx = 0.1   # 100 mm x-dimension
    Ly = 0.1   # 100 mm y-dimension
    z_work = 0.0  # At z=0
    nx = 8     # x divisions
    ny = 8     # y divisions

    workpiece_mesh.create_rectangular_plate(
        center=[0, 0, z_work],
        Lx=Lx, Ly=Ly,
        normal=[0, 0, 1],
        nx=nx, ny=ny
    )

    print(f"Workpiece mesh created:")
    print(f"  Dimensions: {Lx*1000:.0f} x {Ly*1000:.0f} mm")
    print(f"  Position: z = {z_work*1000:.0f} mm")
    print(f"  Divisions: {nx} x {ny}")
    print(f"  Number of vertices: {workpiece_mesh.num_vertices}")
    print(f"  Number of triangles: {workpiece_mesh.num_triangles}")
    print(f"  Number of edges (basis functions): {workpiece_mesh.num_edges}")

    return coil_mesh, workpiece_mesh


def demo_nonlinear_coupled_solver(coil_mesh, workpiece_mesh, esi_table):
    """Demonstrate nonlinear coupled EFIE solver with ESIM."""
    print_header("3. Nonlinear Coupled EFIE Solver")

    freq = 50000  # 50 kHz
    sigma = 2e6   # S/m

    print(f"Frequency: {freq/1000:.0f} kHz")
    print(f"Conductivity: {sigma/1e6:.1f} MS/m")
    print()

    # Create nonlinear coupled solver
    print("Creating NonlinearCoupledEFIESolver...")
    solver = NonlinearCoupledEFIESolver(
        coil_mesh=coil_mesh,
        workpiece_mesh=workpiece_mesh,
        esi_table=esi_table
    )

    # Configure solver
    solver.set_frequency(freq)
    solver.set_coil_conductivity(5.8e7)  # Copper coil

    # Nonlinear iteration parameters
    solver.set_max_nonlinear_iterations(20)
    solver.set_nonlinear_tolerance(1e-3)
    solver.set_relaxation_factor(0.5)  # Underrelaxation for stability

    print(f"Solver configuration:")
    print(f"  Max nonlinear iterations: 20")
    print(f"  Nonlinear tolerance: 1e-3")
    print(f"  Relaxation factor: 0.5")
    print()

    # Apply excitation (1V at coil port)
    V_exc = 1.0 + 0j  # 1V
    solver.set_voltage_excitation(V_exc)

    print(f"Excitation: V = {abs(V_exc):.1f} V")
    print()

    # Solve
    print("Solving nonlinear coupled system...")
    print("-" * 50)

    result = solver.solve()

    print("-" * 50)
    print()

    # Display results
    print_subheader("3.1 Convergence History")

    print(f"{'Iteration':>10} {'Residual':>14} {'Z_avg [mOhm]':>14}")
    print("-" * 45)

    for i, (res, Z_avg) in enumerate(zip(result['residual_history'], result['Z_history'])):
        Z_mOhm = Z_avg * 1000 if Z_avg else 0
        print(f"{i+1:>10} {res:>14.6f} {Z_mOhm:>14.4f}")

    print()
    print(f"Converged in {result['iterations']} iterations")
    print(f"Final residual: {result['final_residual']:.6e}")

    print_subheader("3.2 Impedance Results")

    Z_total = result['impedance']
    Z_coil = result['coil_impedance']
    Z_load = result['load_impedance']

    print(f"Total impedance:     Z = {Z_total.real*1000:.4f} + j{Z_total.imag*1000:.4f} mOhm")
    print(f"                     |Z| = {abs(Z_total)*1000:.4f} mOhm")
    print()
    print(f"Coil self-impedance: Z_coil = {Z_coil.real*1000:.4f} + j{Z_coil.imag*1000:.4f} mOhm")
    print(f"Load impedance:      Z_load = {Z_load.real*1000:.4f} + j{Z_load.imag*1000:.4f} mOhm")

    print_subheader("3.3 Power Analysis")

    I_total = result['current']
    P_total = result['power_total']
    P_coil = result['power_coil']
    P_workpiece = result['power_workpiece']

    print(f"Total current:    I = {abs(I_total):.4f} A (phase: {np.angle(I_total)*180/np.pi:.1f} deg)")
    print()
    print(f"Total power:      P_total = {P_total:.4f} W")
    print(f"Coil losses:      P_coil = {P_coil:.4f} W ({P_coil/P_total*100:.1f}%)")
    print(f"Workpiece power:  P_work = {P_workpiece:.4f} W ({P_workpiece/P_total*100:.1f}%)")
    print()
    print(f"Heating efficiency: {P_workpiece/P_total*100:.1f}%")

    return result


def demo_field_computation(solver, result):
    """Demonstrate field computation at observation points."""
    print_header("4. Field Computation")

    # Observation points (on workpiece surface)
    print_subheader("4.1 H-field on Workpiece Surface")

    obs_points = [
        [0, 0, 0],        # Center
        [0.02, 0, 0],     # 20 mm from center
        [0.04, 0, 0],     # 40 mm from center
        [0, 0.02, 0],     # 20 mm in y
        [0.03, 0.03, 0],  # Diagonal
    ]

    print(f"{'Point [mm]':>20} {'|H| [A/m]':>12} {'Phase [deg]':>12}")
    print("-" * 50)

    for pt in obs_points:
        H = solver.compute_H_field(pt, result['solution'])
        H_mag = np.sqrt(abs(H[0])**2 + abs(H[1])**2 + abs(H[2])**2)
        H_phase = np.angle(H[2]) * 180 / np.pi if abs(H[2]) > 0 else 0
        pt_mm = [p * 1000 for p in pt]
        print(f"({pt_mm[0]:5.0f}, {pt_mm[1]:5.0f}, {pt_mm[2]:5.0f}) {H_mag:>12.1f} {H_phase:>12.1f}")

    print_subheader("4.2 Surface Impedance Distribution")

    Z_distribution = result.get('Z_distribution', None)
    if Z_distribution is not None:
        Z_min = min(abs(z) for z in Z_distribution) * 1000
        Z_max = max(abs(z) for z in Z_distribution) * 1000
        Z_avg = np.mean([abs(z) for z in Z_distribution]) * 1000

        print(f"Surface impedance statistics:")
        print(f"  Min |Z|: {Z_min:.4f} mOhm")
        print(f"  Max |Z|: {Z_max:.4f} mOhm")
        print(f"  Avg |Z|: {Z_avg:.4f} mOhm")
        print(f"  Variation: {(Z_max - Z_min) / Z_avg * 100:.1f}%")
    else:
        print("Z distribution not available in result")


def demo_helper_functions():
    """Demonstrate helper functions for easy solver creation."""
    print_header("5. Helper Functions for Easy Usage")

    print_subheader("5.1 Using create_nonlinear_solver()")

    print("Example code:")
    print("""
    # Create meshes
    coil_mesh = RWGMesh()
    coil_mesh.create_circular_disk(...)

    workpiece_mesh = RWGMesh()
    workpiece_mesh.create_rectangular_plate(...)

    # Create solver with one function call
    solver = create_nonlinear_solver(
        coil_mesh=coil_mesh,
        workpiece_mesh=workpiece_mesh,
        bh_curve=STEEL_BH_CURVE,
        sigma=2e6,     # S/m
        frequency=50000  # Hz
    )

    solver.set_voltage_excitation(1.0)
    result = solver.solve()
    """)

    print_subheader("5.2 Using create_induction_heating_model_with_esim()")

    print("Example code:")
    print("""
    # All-in-one model creation
    model = create_induction_heating_model_with_esim(
        # Coil parameters
        coil_type='pancake',
        coil_radius=0.05,
        coil_height=0.005,
        coil_turns=5,

        # Workpiece parameters
        workpiece_type='plate',
        workpiece_size=[0.1, 0.1],
        workpiece_thickness=0.01,

        # Material parameters
        bh_curve=STEEL_BH_CURVE,
        sigma=2e6,
        frequency=50000,

        # Mesh parameters
        coil_mesh_density=16,
        workpiece_mesh_density=8,
    )

    solver = model['solver']
    solver.set_voltage_excitation(10.0)  # 10V
    result = solver.solve()

    print(f"Heating power: {result['power_workpiece']:.2f} W")
    """)


def main():
    """Main demo function."""
    print()
    print("*" * 70)
    print("*  ESIM + RWG-EFIE Nonlinear Coupled Solver Demo")
    print("*  Induction Heating Analysis with Nonlinear Ferromagnetic Materials")
    print("*" * 70)

    if not RWG_AVAILABLE:
        print()
        print("ERROR: RWG-EFIE solver is not available.")
        print("Please ensure scipy is installed: pip install scipy")
        return

    # Step 1: Generate ESI table
    esi_table = demo_esi_table_generation()

    # Step 2: Create meshes
    coil_mesh, workpiece_mesh = demo_mesh_creation()

    # Step 3: Run nonlinear solver
    result = demo_nonlinear_coupled_solver(coil_mesh, workpiece_mesh, esi_table)

    # Step 4: Show helper functions
    demo_helper_functions()

    print_header("Demo Complete")
    print()
    print("Key takeaways:")
    print("  1. ESIM provides field-dependent surface impedance Z(H)")
    print("  2. Nonlinear iteration converges in ~5-10 steps with underrelaxation")
    print("  3. Surface impedance varies spatially due to nonuniform H-field")
    print("  4. Helper functions simplify model creation")
    print()
    print("For more examples, see:")
    print("  - examples/induction_heating/demo_rwg_efie_3d.py")
    print("  - examples/induction_heating/esim_demo.py")


if __name__ == '__main__':
    main()
