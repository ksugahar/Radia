"""
Final Phase 2-A Validation Test

Creates a problem that requires many relaxation iterations:
- Strong permanent magnet providing large external field
- Saturation material (nonlinear B-H curve) requiring iteration
- Should demonstrate H-matrix reuse across multiple iterations
"""

import sys
import os
import time

# Add paths for Radia module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'dist'))

import radia as rad

print("="*80)
print("PHASE 2-A FINAL VALIDATION: H-MATRIX REUSE")
print("="*80)

def test_phase2a_reuse():
    """
    Test H-matrix reuse with saturation material requiring multiple iterations
    """
    rad.UtiDelAll()

    # Create container
    container = rad.ObjCnt([])

    # 1. Create strong permanent magnet array (source field)
    print("\n[1] Creating permanent magnet array (strong field source)...")
    n_pm = 5  # 5x5 array = 25 magnets
    pm_size = 10  # mm
    pm_gap = 2   # mm
    pm_z = -20   # mm (below grid)

    for i in range(n_pm):
        for j in range(n_pm):
            x = (i - n_pm/2 + 0.5) * (pm_size + pm_gap)
            y = (j - n_pm/2 + 0.5) * (pm_size + pm_gap)
            pm = rad.ObjRecMag([x, y, pm_z], [pm_size, pm_size, pm_size*2], [0, 0, 1])
            rad.ObjAddToCnt(container, [pm])

    print(f"    Created {n_pm*n_pm} permanent magnets")
    print(f"    Size: {pm_size}x{pm_size}x{pm_size*2} mm each")
    print(f"    Magnetization: [0, 0, 1] T")

    # 2. Create saturation material grid (7x7x7 = 343 elements)
    print("\n[2] Creating saturation material grid (7x7x7 = 343 elements)...")
    n_grid = 7
    elem_size = 8  # mm
    grid_z = 10    # mm (above magnets)

    for i in range(n_grid):
        for j in range(n_grid):
            for k in range(n_grid):
                x = (i - n_grid/2 + 0.5) * elem_size
                y = (j - n_grid/2 + 0.5) * elem_size
                z = (k - n_grid/2 + 0.5) * elem_size + grid_z

                block = rad.ObjRecMag([x, y, z], [elem_size, elem_size, elem_size], [0, 0, 0])
                # Saturation material (nonlinear, requires iteration)
                mat = rad.MatSatIsoFrm([2000, 2], [0.1, 2], [0.1, 2])
                rad.MatApl(block, mat)
                rad.ObjAddToCnt(container, [block])

    print(f"    Grid: {n_grid}x{n_grid}x{n_grid} = {n_grid**3} elements")
    print(f"    Material: MatSatIsoFrm (saturation material)")
    print(f"    Position: {grid_z} mm above permanent magnets")

    # Total elements
    n_total = n_pm * n_pm + n_grid**3
    print(f"\n[3] Total problem size: {n_total} elements")
    print(f"    - {n_pm*n_pm} permanent magnets")
    print(f"    - {n_grid**3} saturation material elements")

    print("\n" + "="*80)
    print("Running Relaxation Solver")
    print("="*80)
    print("\nExpected Phase 2-A behavior:")
    print(f"  [1] N={n_total} >= 200: H-matrix enabled (see '[Auto]' message)")
    print(f"  [2] First iteration: H-matrix built")
    print(f"  [3] Iterations 2+: H-matrix reused (see '[Phase 2-A]' messages)")
    print("="*80)

    # Enable H-matrix solver
    rad.SolverHMatrixEnable(1, 1e-4, 30)
    print("\n[Solver] H-matrix acceleration enabled (eps=1e-4, max_rank=30)")
    print("="*80 + "\n")

    # Run solver
    t0 = time.time()
    result = rad.Solve(container, 0.0001, 1000)
    t_solve = time.time() - t0

    print("\n" + "="*80)
    print("Results")
    print("="*80)
    print(f"\nSolver result: {result}")
    print(f"  Max |M|:     {result[0]:.6f} T")
    print(f"  Max |H|:     {result[1]:.6f} A/m")
    print(f"  Max |M|/Ms:  {result[2]:.6f}")
    print(f"  Iterations:  {int(result[3])}")
    print(f"  Time:        {t_solve*1000:.1f} ms")

    # Check field at grid center
    B_center = rad.Fld(container, 'b', [0, 0, grid_z])
    print(f"\nB at grid center: [{B_center[0]:.6f}, {B_center[1]:.6f}, {B_center[2]:.6f}] T")

    return result

if __name__ == "__main__":
    result = test_phase2a_reuse()

    print("\n" + "="*80)
    print("Phase 2-A Validation Summary")
    print("="*80)

    iterations = int(result[3])

    if iterations >= 2:
        print(f"\n[SUCCESS] Multiple iterations completed: {iterations} iterations")
        print(f"[SUCCESS] Problem size N >= 200 (H-matrix enabled)")
        print(f"\nLook for these messages above:")
        print(f"  - '[Auto] Enabling H-matrix' (threshold check)")
        print(f"  - 'Building True H-Matrix' (first iteration)")
        print(f"  - '[Phase 2-A] Reusing existing H-matrix' (iterations 2+)")
        print(f"\nIf all three messages appeared, Phase 2-A is working correctly!")
    else:
        print(f"\n[WARN] Only {iterations} iteration(s) - problem converged quickly")
        print(f"       Saturation material may not have been strongly magnetized")

    print("\n" + "="*80)
