"""
Hybrid PEEC + RWG Solver Demo

This script demonstrates the hybrid PEEC (Partial Element Equivalent Circuit) + RWG
(Rao-Wilton-Glisson) solver for analyzing induction heating systems with:
- PEEC for coil (wire segments with thin-wire approximation)
- RWG for workpiece (surface triangles with full surface current)

The hybrid approach is optimal for induction heating where:
- Coil geometry is simple (circular loops, spirals)
- Workpiece requires accurate surface current modeling

Key Features:
- PEEC coil model with skin effect
- RWG workpiece model with ESIM surface impedance
- Mutual inductance coupling between coil and workpiece
- Efficient for typical induction heating configurations

Reference:
    A.E. Ruehli, "Equivalent Circuit Models for Three-Dimensional
    Multiconductor Systems," IEEE Trans. MTT, 1974

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
        HybridPEECRWGSolver,
        create_hybrid_solver,
    )
    from esim_cell_problem import ESITable
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


# Standard steel B-H curve
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


def create_pancake_coil_path(R_inner, R_outer, z, n_turns, n_points_per_turn=36):
    """Create path for pancake coil."""
    path = []

    for i_turn in range(n_turns):
        # Linear radius interpolation for spiral
        R_start = R_inner + (R_outer - R_inner) * i_turn / n_turns
        R_end = R_inner + (R_outer - R_inner) * (i_turn + 1) / n_turns

        for j in range(n_points_per_turn):
            t = j / n_points_per_turn
            theta = 2 * np.pi * t

            # Interpolate radius
            R = R_start + (R_end - R_start) * t

            x = R * np.cos(theta)
            y = R * np.sin(theta)

            path.append([x, y, z])

    # Close the path
    path.append(path[0])

    return path


def demo_peec_coil_model():
    """Demonstrate PEEC coil modeling."""
    print_header("1. PEEC Coil Model")

    # Pancake coil parameters
    R_inner = 0.02   # 20 mm
    R_outer = 0.05   # 50 mm
    z_coil = 0.005   # 5 mm above workpiece
    n_turns = 5

    # Wire parameters
    wire_radius = 0.002  # 2 mm
    sigma_copper = 5.8e7  # S/m

    print(f"Pancake Coil Parameters:")
    print(f"  Inner radius: {R_inner*1000:.0f} mm")
    print(f"  Outer radius: {R_outer*1000:.0f} mm")
    print(f"  Height: {z_coil*1000:.0f} mm")
    print(f"  Number of turns: {n_turns}")
    print(f"  Wire radius: {wire_radius*1000:.1f} mm")
    print(f"  Wire conductivity: {sigma_copper/1e7:.1f} x 10^7 S/m")
    print()

    # Create coil path
    coil_path = create_pancake_coil_path(R_inner, R_outer, z_coil, n_turns)

    print(f"Coil path points: {len(coil_path)}")
    print(f"Total wire length: {sum(np.linalg.norm(np.array(coil_path[i+1]) - np.array(coil_path[i])) for i in range(len(coil_path)-1))*1000:.1f} mm")
    print()

    # PEEC model info
    print_subheader("1.1 PEEC Model Description")

    print("""
PEEC represents the coil as a network of:
  - Partial self-inductances (L_ii) for each segment
  - Partial mutual inductances (L_ij) between segments
  - Segment resistances (R_i) including skin effect

The PEEC impedance matrix is:
  Z_ij = R_i*delta_ij + j*omega*L_ij

where:
  L_ij = (mu_0 / 4*pi) * integral{ dl_i . dl_j / |r_i - r_j| }

For thin wires, the Neumann formula is used:
  L_ij = (mu_0 / 4*pi) * integral_i integral_j { (dl_i . dl_j) / |r_i - r_j| }
""")

    return coil_path, wire_radius, sigma_copper


def demo_rwg_workpiece_model():
    """Demonstrate RWG workpiece modeling."""
    print_header("2. RWG Workpiece Model")

    # Workpiece parameters
    Lx = 0.1   # 100 mm
    Ly = 0.1   # 100 mm
    z_work = 0.0
    nx = 8
    ny = 8

    print(f"Workpiece Parameters:")
    print(f"  Dimensions: {Lx*1000:.0f} x {Ly*1000:.0f} mm")
    print(f"  Position: z = {z_work*1000:.0f} mm")
    print(f"  Mesh: {nx} x {ny} divisions")
    print()

    # Create workpiece mesh
    workpiece_mesh = RWGMesh()
    workpiece_mesh.create_rectangular_plate(
        center=[0, 0, z_work],
        Lx=Lx, Ly=Ly,
        normal=[0, 0, 1],
        nx=nx, ny=ny
    )

    print(f"RWG Mesh Statistics:")
    print(f"  Vertices: {workpiece_mesh.num_vertices}")
    print(f"  Triangles: {workpiece_mesh.num_triangles}")
    print(f"  Edges (RWG basis functions): {workpiece_mesh.num_edges}")
    print()

    print_subheader("2.1 RWG Model Description")

    print("""
RWG basis functions are edge-based vector functions:
  f_n(r) = (l_n / 2*A_+) * rho_n^+(r)  in triangle T+
         = (l_n / 2*A_-) * rho_n^-(r)  in triangle T-
         = 0                           elsewhere

The RWG-EFIE for workpiece with surface impedance Z_s:
  Z_s * J_s(r) + j*omega*A(r) + grad(Phi(r)) = E_inc(r)

where:
  A(r) = (mu_0/4*pi) * integral{ J_s(r') / |r-r'| dA' }
  Phi(r) = (1/4*pi*eps_0) * integral{ rho_s(r') / |r-r'| dA' }
  rho_s = -div(J_s) / j*omega  (surface charge from continuity)
""")

    return workpiece_mesh


def demo_hybrid_coupling():
    """Demonstrate hybrid PEEC-RWG coupling."""
    print_header("3. Hybrid PEEC-RWG Coupling")

    print_subheader("3.1 Coupling Mechanism")

    print("""
The hybrid system couples PEEC coil and RWG workpiece through:

1. Coil -> Workpiece (B-field excitation):
   E_inc(r_work) = -j*omega * A_coil(r_work)

   where A_coil is the vector potential from coil currents:
   A_coil(r) = (mu_0/4*pi) * integral{ I * dl / |r - r'| }

2. Workpiece -> Coil (reaction field):
   V_react_coil = -j*omega * integral{ I_coil . A_workpiece }

   where A_workpiece is from workpiece surface currents.

The coupled system in matrix form:
  [ Z_coil      M_cw   ] [ I_coil ]   [ V_exc ]
  [ M_wc        Z_work ] [ J_work ] = [ 0     ]

where:
  Z_coil = PEEC impedance matrix
  Z_work = RWG-EFIE impedance matrix with ESIM
  M_cw, M_wc = Mutual coupling matrices
""")

    print_subheader("3.2 Advantages of Hybrid Approach")

    print("""
Why Hybrid PEEC + RWG?

For Coil (PEEC):
  + Simple wire geometry -> accurate thin-wire model
  + Fewer DOF than full surface mesh
  + Skin effect easily included via surface impedance
  + Exact mutual inductance formulas available

For Workpiece (RWG):
  + Accurate surface current distribution
  + Handles complex geometries
  + Natural for ESIM integration
  + Correct edge singularity behavior

Combined Benefits:
  + Optimal DOF allocation (few for coil, many for workpiece)
  + Accurate mutual coupling
  + Compatible with nonlinear ESIM
  + Efficient for typical induction heating configurations
""")


def demo_hybrid_solver():
    """Demonstrate the hybrid PEEC+RWG solver."""
    print_header("4. Hybrid Solver Example")

    # Create coil
    coil_path, wire_radius, sigma_copper = create_pancake_coil_path(
        R_inner=0.02, R_outer=0.05, z=0.005, n_turns=5
    ), 0.002, 5.8e7

    coil_path = create_pancake_coil_path(0.02, 0.05, 0.005, 5)

    # Create workpiece mesh
    workpiece_mesh = RWGMesh()
    workpiece_mesh.create_rectangular_plate(
        center=[0, 0, 0],
        Lx=0.1, Ly=0.1,
        normal=[0, 0, 1],
        nx=8, ny=8
    )

    # Create ESI table for nonlinear workpiece
    sigma_steel = 2e6  # S/m
    freq = 50000       # 50 kHz

    print(f"Creating ESI table for nonlinear material...")
    esi_table = ESITable(
        bh_curve=STEEL_BH_CURVE,
        sigma=sigma_steel,
        frequency=freq
    )

    # Create hybrid solver
    print(f"Creating hybrid PEEC+RWG solver...")

    solver = HybridPEECRWGSolver(
        coil_path=coil_path,
        coil_wire_radius=wire_radius,
        coil_conductivity=sigma_copper,
        workpiece_mesh=workpiece_mesh,
        workpiece_esi_table=esi_table
    )

    solver.set_frequency(freq)

    print()
    print(f"Hybrid System Size:")
    print(f"  PEEC coil segments: {len(coil_path) - 1}")
    print(f"  RWG workpiece edges: {workpiece_mesh.num_edges}")
    print(f"  Total DOF: {len(coil_path) - 1 + workpiece_mesh.num_edges}")
    print()

    # Set excitation
    V_exc = 10.0  # 10V
    solver.set_voltage_excitation(V_exc)

    print(f"Excitation: V = {V_exc:.1f} V at {freq/1000:.0f} kHz")
    print()

    # Solve
    print("Solving hybrid system...")
    print("-" * 50)

    result = solver.solve()

    print("-" * 50)
    print()

    # Display results
    print_subheader("4.1 Impedance Results")

    Z_total = result['impedance']
    Z_coil = result['coil_self_impedance']
    Z_mutual = result['mutual_impedance']
    Z_work = result['workpiece_impedance']

    print(f"Total system impedance:")
    print(f"  Z_total = {Z_total.real*1000:.4f} + j{Z_total.imag*1000:.4f} mOhm")
    print(f"  |Z_total| = {abs(Z_total)*1000:.4f} mOhm")
    print()
    print(f"Component breakdown:")
    print(f"  Z_coil_self = {Z_coil.real*1000:.4f} + j{Z_coil.imag*1000:.4f} mOhm")
    print(f"  Z_mutual = {Z_mutual.real*1000:.4f} + j{Z_mutual.imag*1000:.4f} mOhm")
    print(f"  Z_workpiece = {Z_work.real*1000:.4f} + j{Z_work.imag*1000:.4f} mOhm")

    print_subheader("4.2 Power Analysis")

    I_coil = result['coil_current']
    P_total = result['power_total']
    P_coil = result['power_coil']
    P_work = result['power_workpiece']

    print(f"Coil current: I = {abs(I_coil):.4f} A")
    print()
    print(f"Power distribution:")
    print(f"  P_total = {P_total:.4f} W")
    print(f"  P_coil = {P_coil:.4f} W (coil losses)")
    print(f"  P_workpiece = {P_work:.4f} W (heating power)")
    print()
    print(f"Heating efficiency: {P_work/P_total*100:.1f}%")

    print_subheader("4.3 Coupling Analysis")

    k = result.get('coupling_coefficient', 0)
    M = result.get('mutual_inductance', 0)

    print(f"Mutual inductance: M = {M*1e6:.4f} uH")
    print(f"Coupling coefficient: k = {k:.4f}")

    return result


def demo_helper_function():
    """Demonstrate helper function for easy solver creation."""
    print_header("5. Helper Function Usage")

    print("""
The create_hybrid_solver() function provides a simple interface:

Example code:
-------------

from rwg_efie_solver import create_hybrid_solver

# Create solver with simple parameters
solver = create_hybrid_solver(
    # Coil parameters
    coil_type='pancake',
    coil_inner_radius=0.02,  # 20 mm
    coil_outer_radius=0.05,  # 50 mm
    coil_height=0.005,       # 5 mm
    coil_turns=5,
    wire_radius=0.002,       # 2 mm

    # Workpiece parameters
    workpiece_size=[0.1, 0.1],  # 100 x 100 mm
    workpiece_mesh_density=10,

    # Material parameters
    coil_conductivity=5.8e7,    # Copper
    workpiece_bh_curve=STEEL_BH_CURVE,
    workpiece_conductivity=2e6,  # Hot steel
    frequency=50000,             # 50 kHz
)

# Solve
solver.set_voltage_excitation(100.0)  # 100V
result = solver.solve()

print(f"Heating power: {result['power_workpiece']:.2f} W")
print(f"Efficiency: {result['power_workpiece']/result['power_total']*100:.1f}%")
""")


def main():
    """Main demo function."""
    print()
    print("*" * 70)
    print("*  Hybrid PEEC + RWG Solver Demo")
    print("*  Combining Wire and Surface Models for Induction Heating")
    print("*" * 70)

    if not RWG_AVAILABLE:
        print()
        print("ERROR: RWG-EFIE solver is not available.")
        print("Please ensure scipy is installed: pip install scipy")
        return

    # Step 1: PEEC coil model
    demo_peec_coil_model()

    # Step 2: RWG workpiece model
    demo_rwg_workpiece_model()

    # Step 3: Coupling mechanism
    demo_hybrid_coupling()

    # Step 4: Hybrid solver
    demo_hybrid_solver()

    # Step 5: Helper function
    demo_helper_function()

    print_header("Demo Complete")
    print()
    print("Key takeaways:")
    print("  1. PEEC is optimal for wire/coil modeling")
    print("  2. RWG is optimal for surface current modeling")
    print("  3. Hybrid combines strengths of both methods")
    print("  4. Compatible with nonlinear ESIM for ferromagnetic workpieces")


if __name__ == '__main__':
    main()
