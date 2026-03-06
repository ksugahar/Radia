"""
Loop-Star Decomposition Demo for Low-Frequency EFIE Stability

This script demonstrates the Loop-Star decomposition method for eliminating
low-frequency breakdown in the Electric Field Integral Equation (EFIE).

The Problem:
    Standard EFIE has the form: Z = j*omega*L + 1/(j*omega)*P
    At low frequencies (omega -> 0), the P/(j*omega) term becomes singular,
    causing numerical instability ("low-frequency breakdown").

The Solution:
    Loop-Star decomposition separates the current into:
    - Loop (solenoidal): div(J) = 0, represents circulating currents
    - Star (irrotational): curl(J) = 0, represents charge accumulation

    For MQS (Magneto-Quasi-Static) problems where displacement current is
    negligible, only the Loop component is needed, avoiding the singular P term.

Key Benefits:
    1. Stable at arbitrarily low frequencies
    2. Correct physics for MQS regime
    3. Reduced problem size (Loop DOF < RWG DOF)

Reference:
    G. Vecchi, "Loop-Star Decomposition of Basis Functions in the
    Discretization of the EFIE," IEEE TAP, 1999.

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
        RWGEFIESolver,
        LoopStarDecomposition,
        LoopStarEFIESolver,
        create_induction_heating_model,
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


def demo_loop_star_theory():
    """Explain the Loop-Star decomposition theory."""
    print_header("1. Loop-Star Decomposition Theory")

    print("""
The Low-Frequency Breakdown Problem:
------------------------------------

Standard EFIE impedance matrix:

    Z_mn = j*omega*L_mn + R_mn + P_mn/(j*omega)

where:
    L_mn = (mu_0/4*pi) * integral{ f_m . f_n / R } dA dA'  (vector potential)
    R_mn = skin effect resistance
    P_mn = (1/4*pi*eps_0) * integral{ div(f_m) div(f_n) / R } dA dA' (scalar potential)

At low frequencies (omega -> 0):
    - j*omega*L_mn -> 0  (inductive term vanishes)
    - P_mn/(j*omega) -> infinity  (capacitive term blows up!)

This causes numerical instability and poor conditioning.


The Loop-Star Solution:
-----------------------

Decompose RWG basis functions into:

    J = T_loop * J_loop + T_star * J_star

where:
    Loop basis: Solenoidal (div = 0), based on triangle circulations
    Star basis: Irrotational (curl = 0), based on vertex charges

Key insight: For Loop basis functions,
    div(Loop) = 0  =>  P_LL = T_loop^T * P * T_loop ~ 0

So the Loop-Loop block has NO singular P term!

For MQS problems:
    - div(J) = 0 (no charge accumulation)
    - Only Loop component is needed
    - System reduces from N_rwg to N_loop DOF
    - No low-frequency breakdown!
""")


def demo_decomposition_structure():
    """Demonstrate the structure of Loop-Star decomposition."""
    print_header("2. Loop-Star Decomposition Structure")

    # Create a simple mesh
    mesh = RWGMesh()
    mesh.create_rectangular_plate(
        center=[0, 0, 0],
        Lx=0.1, Ly=0.1,
        normal=[0, 0, 1],
        nx=4, ny=4
    )

    print(f"Test Mesh (4x4 rectangular plate):")
    print(f"  Vertices: {len(mesh.vertices)}")
    print(f"  Triangles: {len(mesh.triangles)}")
    print(f"  Edges: {len(mesh.edges)}")
    print(f"  Interior edges (RWG DOF): {mesh.num_interior_edges}")
    print()

    # Create Loop-Star decomposition
    decomp = LoopStarDecomposition(mesh)

    print(f"Loop-Star Decomposition:")
    print(f"  Loop DOF (= num_triangles): {decomp.num_loops}")
    print(f"  Star DOF (= num_vertices - 1): {decomp.num_stars}")
    print(f"  Total DOF: {decomp.num_loops + decomp.num_stars}")
    print()

    print(f"Transformation matrices:")
    print(f"  T_loop shape: {decomp.T_loop.shape}")
    print(f"  T_star shape: {decomp.T_star.shape}")
    print(f"  T_ls shape: {decomp.T_ls.shape}")
    print()

    # Verify orthogonality properties
    print_subheader("2.1 Orthogonality Properties")

    # T_loop^T * T_star should be zero (orthogonal bases)
    cross = decomp.T_loop.T @ decomp.T_star
    cross_norm = np.linalg.norm(cross.toarray())

    print(f"Loop-Star orthogonality:")
    print(f"  ||T_loop^T * T_star|| = {cross_norm:.6e}")
    print(f"  (Should be small, ideally zero)")
    print()

    # Verify that T_ls spans the RWG space
    # T_ls * T_ls^T should be close to identity (full rank)
    T_ls_dense = decomp.T_ls.toarray()
    rank = np.linalg.matrix_rank(T_ls_dense)

    print(f"Basis completeness:")
    print(f"  Rank of T_ls: {rank}")
    print(f"  Full rank would be: {min(decomp.num_rwg, decomp.num_loops + decomp.num_stars)}")

    return mesh, decomp


def demo_low_frequency_stability(mesh):
    """Demonstrate low-frequency stability of Loop-Star solver."""
    print_header("3. Low-Frequency Stability Comparison")

    print("Comparing standard RWG-EFIE vs Loop-Star at various frequencies...")
    print()

    frequencies = [1e6, 1e5, 1e4, 1e3, 1e2, 10]  # 1 MHz down to 10 Hz

    print(f"{'Frequency':>12} {'RWG Z (Ohm)':>20} {'LS Z (Ohm)':>20} {'Cond(Z_rwg)':>15}")
    print("-" * 75)

    for freq in frequencies:
        # Standard RWG solver
        rwg_solver = RWGEFIESolver(mesh)
        rwg_solver.set_frequency(freq)
        rwg_solver.set_conductivity(5.8e7)

        # Find port edges
        interior_edges = [i for i, e in enumerate(mesh.edges) if not e.is_boundary]
        if interior_edges:
            rwg_solver.define_port([interior_edges[0]])
        rwg_solver.set_voltage_excitation(1.0)

        # Assemble matrices
        rwg_solver.assemble_matrices()

        # Get condition number
        try:
            cond_rwg = np.linalg.cond(rwg_solver.Z_matrix)
        except:
            cond_rwg = float('inf')

        # Solve
        try:
            rwg_solver.solve()
            Z_rwg = rwg_solver.impedance
            Z_rwg_str = f"{Z_rwg.real:.4e} + j{Z_rwg.imag:.4e}"
        except:
            Z_rwg_str = "FAILED"

        # Loop-Star solver
        ls_solver = LoopStarEFIESolver(mesh)
        ls_solver.set_frequency(freq)
        ls_solver.set_conductivity(5.8e7)
        if interior_edges:
            ls_solver.define_port([interior_edges[0]])
        ls_solver.set_voltage_excitation(1.0)

        try:
            result_ls = ls_solver.solve()
            Z_ls = result_ls['impedance']
            Z_ls_str = f"{Z_ls.real:.4e} + j{Z_ls.imag:.4e}"
        except:
            Z_ls_str = "FAILED"

        # Format condition number
        if cond_rwg < 1e15:
            cond_str = f"{cond_rwg:.2e}"
        else:
            cond_str = "SINGULAR"

        print(f"{freq:>12.0e} {Z_rwg_str:>20} {Z_ls_str:>20} {cond_str:>15}")

    print()
    print("Observations:")
    print("  - Standard RWG becomes ill-conditioned at low frequencies")
    print("  - Loop-Star (MQS mode) remains stable even at 10 Hz")
    print("  - For MQS problems, Loop-Star is the correct physical model")


def demo_mqs_vs_fullwave():
    """Compare MQS (Loop-only) vs Full-wave (Loop+Star) solutions."""
    print_header("4. MQS vs Full-Wave Comparison")

    # Create mesh
    mesh = RWGMesh()
    mesh.create_rectangular_plate(
        center=[0, 0, 0],
        Lx=0.1, Ly=0.1,
        normal=[0, 0, 1],
        nx=4, ny=4
    )

    freq = 50e3  # 50 kHz

    print(f"Frequency: {freq/1e3:.0f} kHz")
    print()

    # Loop-Star solver in MQS mode
    print_subheader("4.1 MQS Mode (Loop-only)")

    ls_mqs = LoopStarEFIESolver(mesh)
    ls_mqs.set_frequency(freq)
    ls_mqs.set_conductivity(5.8e7)
    ls_mqs.use_mqs = True  # Use MQS approximation

    interior_edges = [i for i, e in enumerate(mesh.edges) if not e.is_boundary]
    if interior_edges:
        ls_mqs.define_port([interior_edges[0]])
    ls_mqs.set_voltage_excitation(1.0)

    result_mqs = ls_mqs.solve(verbose=True)

    print()
    print(f"MQS Results:")
    print(f"  DOF (Loop only): {result_mqs['dof_loop']}")
    print(f"  Impedance: {result_mqs['resistance']:.6f} + j{result_mqs['reactance']:.6f} Ohm")
    print(f"  Inductance: {result_mqs['inductance']*1e6:.3f} uH")

    # Loop-Star solver in Full-wave mode
    print_subheader("4.2 Full-Wave Mode (Loop + Star)")

    ls_full = LoopStarEFIESolver(mesh)
    ls_full.set_frequency(freq)
    ls_full.set_conductivity(5.8e7)
    ls_full.use_mqs = False  # Use full Loop+Star

    if interior_edges:
        ls_full.define_port([interior_edges[0]])
    ls_full.set_voltage_excitation(1.0)

    result_full = ls_full.solve(verbose=True)

    print()
    print(f"Full-Wave Results:")
    print(f"  DOF (Loop + Star): {result_full['dof_loop']} + {result_full['dof_star']} = {result_full['dof_loop'] + result_full['dof_star']}")
    print(f"  Impedance: {result_full['resistance']:.6f} + j{result_full['reactance']:.6f} Ohm")
    print(f"  Inductance: {result_full['inductance']*1e6:.3f} uH")

    # Comparison
    print_subheader("4.3 Comparison")

    Z_mqs = result_mqs['impedance']
    Z_full = result_full['impedance']
    diff = abs(Z_full - Z_mqs) / abs(Z_mqs) * 100 if abs(Z_mqs) > 0 else 0

    print(f"Impedance difference: {diff:.2f}%")
    print()
    print("At 50 kHz, displacement current is negligible, so MQS and full-wave")
    print("solutions should be nearly identical. MQS is faster due to fewer DOF.")


def demo_induction_heating_application():
    """Demonstrate Loop-Star for induction heating application."""
    print_header("5. Induction Heating Application")

    print("Typical induction heating operates at 1 kHz - 1 MHz.")
    print("At these frequencies, MQS (Loop-only) is the correct model.")
    print()

    # Create a spiral coil mesh
    mesh = RWGMesh()
    mesh.create_spiral_coil(
        center=[0, 0, 0],
        inner_radius=0.02,
        outer_radius=0.05,
        pitch=0.005,
        num_turns=3,
        axis=[0, 0, 1],
        wire_radius=0.002,
        num_around=6
    )

    print(f"Spiral Coil Mesh:")
    print(f"  Inner radius: 20 mm")
    print(f"  Outer radius: 50 mm")
    print(f"  Turns: 3")
    print(f"  Triangles: {len(mesh.triangles)}")
    print(f"  Interior edges: {mesh.num_interior_edges}")
    print()

    # Create Loop-Star decomposition
    decomp = LoopStarDecomposition(mesh)

    print(f"Loop-Star Decomposition:")
    print(f"  Loop DOF: {decomp.num_loops}")
    print(f"  Star DOF: {decomp.num_stars}")
    print(f"  DOF reduction: {mesh.num_interior_edges} -> {decomp.num_loops} ({100*(1-decomp.num_loops/mesh.num_interior_edges):.1f}% reduction)")
    print()

    # Solve at 50 kHz
    freq = 50e3

    ls_solver = LoopStarEFIESolver(mesh)
    ls_solver.set_frequency(freq)
    ls_solver.set_conductivity(5.8e7)
    ls_solver.use_mqs = True

    interior_edges = [i for i, e in enumerate(mesh.edges) if not e.is_boundary]
    if interior_edges:
        ls_solver.define_port([interior_edges[0]])
    ls_solver.set_voltage_excitation(1.0)

    print(f"Solving at {freq/1e3:.0f} kHz...")
    result = ls_solver.solve(verbose=True)

    print()
    print(f"Coil Self-Impedance:")
    print(f"  Z = {result['resistance']*1e3:.4f} + j{result['reactance']*1e3:.4f} mOhm")
    print(f"  L = {result['inductance']*1e6:.3f} uH")


def demo_usage_guide():
    """Show how to use Loop-Star decomposition."""
    print_header("6. Usage Guide")

    print("""
Basic Usage:
------------

from rwg_efie_solver import (
    RWGMesh,
    LoopStarDecomposition,
    LoopStarEFIESolver,
)

# Create mesh
mesh = RWGMesh()
mesh.create_rectangular_plate([0,0,0], 0.1, 0.1, [0,0,1], nx=8, ny=8)

# Option 1: Use LoopStarEFIESolver directly
solver = LoopStarEFIESolver(mesh)
solver.set_frequency(50e3)
solver.set_conductivity(5.8e7)
solver.use_mqs = True  # MQS mode (recommended for low frequencies)

solver.define_port([edge_indices])
solver.set_voltage_excitation(1.0)

result = solver.solve()
print(f"Z = {result['impedance']}")


# Option 2: Use decomposition directly
decomp = LoopStarDecomposition(mesh)

# Transform existing RWG matrix to Loop-Star basis
Z_ls = decomp.transform_to_loop_star(Z_rwg)

# Solve in Loop space only (MQS)
I_rwg = decomp.solve_mqs_loop_only(Z_rwg, V_rwg)


When to Use Loop-Star:
----------------------

1. Low frequencies (< 1 MHz): Always recommended
2. MQS problems: Induction heating, eddy currents, transformers
3. Large systems: Loop DOF < RWG DOF, so faster solve
4. Ill-conditioned problems: Better numerical stability

When NOT to Use Loop-Star:
--------------------------

1. High frequencies where displacement current matters
2. Antenna problems with radiation
3. Problems where charge accumulation is important
""")


def main():
    """Main demo function."""
    print()
    print("*" * 70)
    print("*  Loop-Star Decomposition Demo")
    print("*  Low-Frequency Stable EFIE for Induction Heating")
    print("*" * 70)

    if not RWG_AVAILABLE:
        print()
        print("ERROR: RWG-EFIE solver is not available.")
        print("Please ensure scipy is installed: pip install scipy")
        return

    # Step 1: Theory
    demo_loop_star_theory()

    # Step 2: Decomposition structure
    mesh, decomp = demo_decomposition_structure()

    # Step 3: Low-frequency stability
    demo_low_frequency_stability(mesh)

    # Step 4: MQS vs Full-wave
    demo_mqs_vs_fullwave()

    # Step 5: Induction heating application
    demo_induction_heating_application()

    # Step 6: Usage guide
    demo_usage_guide()

    print_header("Demo Complete")
    print()
    print("Key Takeaways:")
    print("  1. Standard EFIE breaks down at low frequencies (omega -> 0)")
    print("  2. Loop-Star decomposition separates inductive and capacitive effects")
    print("  3. MQS mode (Loop-only) is stable at any frequency")
    print("  4. For induction heating, MQS is the correct physical model")
    print("  5. Loop DOF < RWG DOF, so computation is faster")


if __name__ == '__main__':
    main()
