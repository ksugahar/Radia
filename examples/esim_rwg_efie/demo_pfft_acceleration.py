"""
pFFT (Pre-corrected FFT) Acceleration Demo

This script demonstrates the pFFT acceleration for large-scale RWG-EFIE problems.
pFFT reduces the matrix-vector product complexity from O(N^2) to O(N log N),
enabling analysis of much larger problems.

The demo:
1. Compares direct vs pFFT matrix-vector products
2. Shows memory usage reduction
3. Demonstrates iterative solver with pFFT acceleration

Key Features:
- Pre-corrected FFT for near-field accuracy
- Hierarchical grid for optimal performance
- Compatible with nonlinear ESIM iteration

Reference:
    J.R. Phillips, J.K. White, "A precorrected-FFT method for electrostatic
    analysis of complicated 3-D structures," IEEE TCAD, 1997

Author: Radia Development Team
Date: 2026-01-08
"""

import sys
import os
import numpy as np
import time

# Add radia module path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

try:
    from rwg_efie_solver import (
        RWGMesh,
        RWGEFIESolver,
        pFFTAccelerator,
        pFFTEFIESolver,
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


def create_test_mesh(mesh_density=8):
    """Create a test mesh (rectangular plate)."""
    mesh = RWGMesh()

    # Plate parameters
    Lx = 0.1   # 100 mm
    Ly = 0.1   # 100 mm
    nx = mesh_density
    ny = mesh_density

    mesh.create_rectangular_plate(
        center=[0, 0, 0],
        Lx=Lx, Ly=Ly,
        normal=[0, 0, 1],
        nx=nx, ny=ny
    )

    return mesh


def demo_pfft_basics():
    """Demonstrate basic pFFT accelerator functionality."""
    print_header("1. pFFT Accelerator Basics")

    # Create test mesh
    mesh = create_test_mesh(mesh_density=10)

    print(f"Test mesh:")
    print(f"  Vertices: {mesh.num_vertices}")
    print(f"  Triangles: {mesh.num_triangles}")
    print(f"  Edges (DOF): {mesh.num_edges}")
    print()

    # Create pFFT accelerator
    print("Creating pFFT accelerator...")

    grid_size = 32  # FFT grid size
    near_threshold = 2.0  # Near-field threshold

    pfft = pFFTAccelerator(
        mesh=mesh,
        grid_size=grid_size,
        near_threshold=near_threshold
    )

    print(f"pFFT configuration:")
    print(f"  Grid size: {grid_size} x {grid_size} x {grid_size}")
    print(f"  Near-field threshold: {near_threshold}")
    print(f"  Memory usage: {pfft.memory_usage_mb:.2f} MB")
    print()

    # Compare with direct method
    print_subheader("1.1 Memory Comparison")

    N = mesh.num_edges
    dense_memory = N * N * 16 / (1024**2)  # Complex double = 16 bytes
    pfft_memory = pfft.memory_usage_mb

    print(f"Number of unknowns (N): {N}")
    print(f"Dense matrix memory: {dense_memory:.2f} MB (N^2 complex)")
    print(f"pFFT memory: {pfft_memory:.2f} MB")
    print(f"Memory reduction: {dense_memory / pfft_memory:.1f}x")

    return mesh, pfft


def demo_matvec_performance(mesh, pfft):
    """Compare matrix-vector product performance."""
    print_header("2. Matrix-Vector Product Performance")

    N = mesh.num_edges
    freq = 50000  # 50 kHz

    # Create random test vector
    np.random.seed(42)
    x = np.random.randn(N) + 1j * np.random.randn(N)

    # Direct method (for comparison)
    print_subheader("2.1 Direct Method (O(N^2))")

    print("Building full matrix...")
    t_start = time.time()

    solver_direct = RWGEFIESolver(mesh)
    solver_direct.set_frequency(freq)
    solver_direct.set_conductivity(5.8e7)
    solver_direct.assemble_matrices()

    t_build = time.time() - t_start
    print(f"Matrix assembly time: {t_build:.3f} s")

    # Time direct matvec
    t_start = time.time()
    for _ in range(10):
        y_direct = solver_direct.matvec(x)
    t_direct = (time.time() - t_start) / 10

    print(f"Direct matvec time: {t_direct*1000:.3f} ms")

    # pFFT method
    print_subheader("2.2 pFFT Method (O(N log N))")

    print("Setting up pFFT...")
    t_start = time.time()

    pfft.set_frequency(freq)
    pfft.setup()

    t_setup = time.time() - t_start
    print(f"pFFT setup time: {t_setup:.3f} s")

    # Time pFFT matvec
    t_start = time.time()
    for _ in range(10):
        y_pfft = pfft.matvec(x)
    t_pfft = (time.time() - t_start) / 10

    print(f"pFFT matvec time: {t_pfft*1000:.3f} ms")

    # Comparison
    print_subheader("2.3 Comparison")

    speedup = t_direct / t_pfft
    rel_error = np.linalg.norm(y_pfft - y_direct) / np.linalg.norm(y_direct)

    print(f"Speedup: {speedup:.2f}x")
    print(f"Relative error: {rel_error:.2e}")
    print()

    if rel_error < 1e-3:
        print("Result: pFFT provides accurate acceleration!")
    else:
        print("Warning: Accuracy may need tuning (adjust grid_size or near_threshold)")


def demo_pfft_solver():
    """Demonstrate pFFT-accelerated iterative solver."""
    print_header("3. pFFT-Accelerated Iterative Solver")

    # Create larger mesh to show scalability
    mesh = create_test_mesh(mesh_density=12)

    print(f"Larger test mesh:")
    print(f"  Vertices: {mesh.num_vertices}")
    print(f"  Triangles: {mesh.num_triangles}")
    print(f"  Edges (DOF): {mesh.num_edges}")
    print()

    # Create pFFT solver
    solver = pFFTEFIESolver(
        mesh=mesh,
        grid_size=64,
        near_threshold=1.5
    )

    freq = 50000  # 50 kHz
    sigma = 5.8e7  # Copper

    solver.set_frequency(freq)
    solver.set_conductivity(sigma)

    print(f"Frequency: {freq/1000:.0f} kHz")
    print(f"Conductivity: {sigma/1e7:.1f} x 10^7 S/m")
    print()

    # Set excitation
    V_exc = 1.0 + 0j
    solver.set_voltage_excitation(V_exc)

    # Configure iterative solver
    solver.set_tolerance(1e-6)
    solver.set_max_iterations(100)

    print("Solving with pFFT-accelerated GMRES...")
    print("-" * 50)

    t_start = time.time()
    result = solver.solve()
    t_solve = time.time() - t_start

    print("-" * 50)
    print()

    # Display results
    print_subheader("3.1 Solver Results")

    print(f"Solution time: {t_solve:.3f} s")
    print(f"Iterations: {result['iterations']}")
    print(f"Final residual: {result['residual']:.2e}")
    print()

    Z = result['impedance']
    print(f"Impedance: Z = {Z.real*1000:.4f} + j{Z.imag*1000:.4f} mOhm")
    print(f"           |Z| = {abs(Z)*1000:.4f} mOhm")


def demo_scaling():
    """Demonstrate pFFT scaling with problem size."""
    print_header("4. Scalability Analysis")

    print("Testing pFFT performance vs problem size...")
    print()

    densities = [6, 8, 10, 12]
    results = []

    print(f"{'Mesh':>8} {'N (DOF)':>10} {'Dense [MB]':>12} {'pFFT [MB]':>12} {'Ratio':>8}")
    print("-" * 60)

    for density in densities:
        mesh = create_test_mesh(mesh_density=density)
        N = mesh.num_edges

        # Dense memory estimate
        dense_mb = N * N * 16 / (1024**2)

        # pFFT memory
        grid_size = max(32, int(np.ceil(N / 4) * 2))  # Adaptive grid
        pfft = pFFTAccelerator(mesh, grid_size=min(grid_size, 128))
        pfft_mb = pfft.memory_usage_mb

        ratio = dense_mb / pfft_mb

        print(f"{density}x{density:>5} {N:>10} {dense_mb:>12.1f} {pfft_mb:>12.1f} {ratio:>8.1f}x")

        results.append({
            'density': density,
            'N': N,
            'dense_mb': dense_mb,
            'pfft_mb': pfft_mb,
            'ratio': ratio
        })

    print()
    print("As N increases, memory ratio improves from O(N^2) to O(N log N)")


def demo_integration_with_esim():
    """Show how pFFT integrates with ESIM for large nonlinear problems."""
    print_header("5. Integration with ESIM")

    print("""
The pFFT acceleration can be combined with ESIM for large-scale
nonlinear induction heating analysis:

Example code:
-------------

from rwg_efie_solver import (
    RWGMesh,
    pFFTAccelerator,
    NonlinearCoupledEFIESolver,
)
from esim_cell_problem import ESITable

# Create large workpiece mesh
workpiece_mesh = RWGMesh()
workpiece_mesh.create_rectangular_plate(
    center=[0, 0, 0],
    Lx=0.5, Ly=0.5,  # 500mm x 500mm
    nx=50, ny=50      # 5000+ triangles
)

# Create ESI table
esi_table = ESITable(
    bh_curve=STEEL_BH_CURVE,
    sigma=2e6,
    frequency=50000
)

# Create nonlinear solver with pFFT acceleration
solver = NonlinearCoupledEFIESolver(
    coil_mesh=coil_mesh,
    workpiece_mesh=workpiece_mesh,
    esi_table=esi_table,
    use_pfft=True,           # Enable pFFT
    pfft_grid_size=128,      # FFT grid
    pfft_near_threshold=1.5  # Near-field
)

solver.set_frequency(50000)
solver.set_voltage_excitation(100.0)  # 100V
result = solver.solve()

# pFFT enables solving 10,000+ DOF problems efficiently!
print(f"Solved {workpiece_mesh.num_edges} DOF in {result['solve_time']:.1f} s")
""")


def main():
    """Main demo function."""
    print()
    print("*" * 70)
    print("*  pFFT (Pre-corrected FFT) Acceleration Demo")
    print("*  Fast Matrix-Vector Products for Large-Scale RWG-EFIE")
    print("*" * 70)

    if not RWG_AVAILABLE:
        print()
        print("ERROR: RWG-EFIE solver is not available.")
        print("Please ensure scipy is installed: pip install scipy")
        return

    # Step 1: Basic pFFT demo
    mesh, pfft = demo_pfft_basics()

    # Step 2: Performance comparison
    demo_matvec_performance(mesh, pfft)

    # Step 3: pFFT solver
    demo_pfft_solver()

    # Step 4: Scaling analysis
    demo_scaling()

    # Step 5: ESIM integration
    demo_integration_with_esim()

    print_header("Demo Complete")
    print()
    print("Key takeaways:")
    print("  1. pFFT reduces memory from O(N^2) to O(N log N)")
    print("  2. Matrix-vector products are much faster for large N")
    print("  3. Accuracy is maintained through near-field correction")
    print("  4. Compatible with nonlinear ESIM iteration")


if __name__ == '__main__':
    main()
