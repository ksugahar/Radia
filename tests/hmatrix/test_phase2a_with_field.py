"""
Test Phase 2-A with External Field

This creates a realistic problem requiring multiple relaxation iterations:
- Permanent magnet providing external field
- Magnetizable material (high susceptibility) responding to field
- Should require 10-50 relaxation iterations
"""

import sys
import os
import time

# Add paths for Radia module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'dist'))

import radia as rad

print("="*80)
print("PHASE 2-A TEST: H-MATRIX REUSE WITH MAGNETIZABLE MATERIAL")
print("="*80)

def test_magnetizable_problem():
    """
    Create problem with:
    - Permanent magnet (source field)
    - Magnetizable material grid (high susceptibility)
    - Should require multiple iterations
    """
    print("\n" + "="*80)
    print("Creating Test Problem")
    print("="*80)

    rad.UtiDelAll()

    # Create container
    container = rad.ObjCnt([])

    # 1. Create permanent magnet (field source)
    print("\n[1] Creating permanent magnet (field source)...")
    pm = rad.ObjRecMag([0, 0, -30], [20, 20, 10], [0, 0, 1])  # Mz = 1 T
    rad.ObjAddToCnt(container, [pm])
    print(f"    Permanent magnet: 20x20x10 mm at z=-30 mm")
    print(f"    Magnetization: [0, 0, 1] T")

    # 2. Create magnetizable material grid (7x7x7 = 343 elements)
    print("\n[2] Creating magnetizable material grid (7x7x7 = 343 elements)...")
    n = 7
    elem_size = 8  # mm

    for i in range(n):
        for j in range(n):
            for k in range(n):
                x = (i - n/2 + 0.5) * elem_size
                y = (j - n/2 + 0.5) * elem_size
                z = (k - n/2 + 0.5) * elem_size + 20  # Above permanent magnet

                block = rad.ObjRecMag([x, y, z], [elem_size, elem_size, elem_size], [0, 0, 0])
                # High susceptibility material (soft iron)
                mat = rad.MatLin(999)  # chi = 999
                rad.MatApl(block, mat)
                rad.ObjAddToCnt(container, [block])

    print(f"    Grid: 7x7x7 = {n**3} elements")
    print(f"    Material: MatLin(999) - high susceptibility")
    print(f"    Position: Above permanent magnet")

    # 3. Calculate field from permanent magnet at grid center
    H_pm = rad.Fld(pm, 'h', [0, 0, 20])
    print(f"\n[3] External field from permanent magnet:")
    print(f"    H at grid center: [{H_pm[0]:.2f}, {H_pm[1]:.2f}, {H_pm[2]:.2f}] A/m")

    print("\n" + "="*80)
    print("Running Relaxation Solver (N=343 + 1 = 344 elements)")
    print("="*80)
    print("\nExpected behavior:")
    print("  - N >= 200: H-matrix threshold check passes")
    print("  - Multiple iterations required (magnetizable material)")
    print("  - Look for '[Auto]' and '[Phase 2-A]' messages below")
    print("="*80)

    # Enable H-matrix for solver (Phase 2-A will check threshold)
    rad.SolverHMatrixEnable(1, 1e-4, 30)  # enable, eps, max_rank
    print("\n[Enabled] H-matrix solver acceleration (eps=1e-4, max_rank=30)")
    print("="*80)

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
    print(f"  Iterations:  {result[3]}")
    print(f"  Time:        {t_solve*1000:.1f} ms")

    # Check magnetization in grid center
    B_center = rad.Fld(container, 'b', [0, 0, 20])
    print(f"\nB at grid center: [{B_center[0]:.6f}, {B_center[1]:.6f}, {B_center[2]:.6f}] T")

    return result

if __name__ == "__main__":
    result = test_magnetizable_problem()

    print("\n" + "="*80)
    print("Phase 2-A Validation")
    print("="*80)

    if result[3] > 1:
        print(f"\n[OK] Multiple iterations completed: {result[3]} iterations")
        print(f"[OK] Problem size: N >= 200 (H-matrix enabled)")
        print(f"\n  If '[Auto]' and '[Phase 2-A]' messages appeared above,")
        print(f"  then Phase 2-A H-matrix reuse is working correctly!")
    else:
        print(f"\n[WARN] Only {result[3]} iteration - problem converged too quickly")
        print(f"  (May need stronger permanent magnet or higher susceptibility)")

    print("\n" + "="*80)
