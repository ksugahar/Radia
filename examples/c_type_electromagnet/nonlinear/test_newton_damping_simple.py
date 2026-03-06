"""
Minimal test of Newton line search damping using cube model
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad
rad.FldUnits('m')
import time

mm = 1e-3  # 1 mm in meters

# Simple B-H curve
bh_data = [
    [0, 0], [100, 0.11], [500, 0.97], [1000, 1.79], [2000, 2.03],
    [5000, 2.24], [10000, 2.39], [20000, 2.53], [50000, 2.69], [100000, 2.81],
    [200000, 2.92], [500000, 3.06], [1000000, 3.17], [2000000, 3.28], [5000000, 3.42],
    [10000000, 3.53]
]

def create_cube(n=6):
    """Create simple cube with n^3 hexahedra"""
    mat = rad.MatSatIsoTab(bh_data)

    cube_size = 1000.0 * mm  # 1m
    dx = cube_size / n

    all_hexes = []
    for ix in range(n):
        for iy in range(n):
            for iz in range(n):
                x = ix * dx
                y = iy * dx
                z = iz * dx

                verts = [
                    [x, y, z], [x+dx, y, z], [x+dx, y+dx, z], [x, y+dx, z],
                    [x, y, z+dx], [x+dx, y, z+dx], [x+dx, y+dx, z+dx], [x, y+dx, z+dx]
                ]

                obj = rad.ObjHexahedron(verts, [0, 0, 0])
                rad.MatApl(obj, mat)
                all_hexes.append(obj)

    return rad.ObjCnt(all_hexes)

def run_test(name, use_newton, use_damping):
    print(f"\n{'='*50}")
    print(f"Test: {name}")
    print(f"{'='*50}")

    model = create_cube(n=6)

    # Add uniform external field
    ext_field = rad.ObjBckg(lambda p: [0, 0, 200000])  # 200 kA/m in Z direction
    model_with_field = rad.ObjCnt([model, ext_field])

    rad.SetSolver(2)  # HACApK
    rad.SetBiCGSTABTolerance(1e-4)
    rad.SetNewtonMethod(use_newton)

    if use_damping:
        rad.SetNewtonDamping(True, max_iter=5, min_omega=0.01)
        print("Damping: ENABLED")
    else:
        rad.SetNewtonDamping(False)
        print("Damping: DISABLED")

    # Build
    print("Building...")
    t0 = time.time()
    rad.BuildMatrix(model_with_field)
    print(f"Build: {time.time()-t0:.1f}s")

    # Solve
    print("Solving...")
    t0 = time.time()
    n_iter = rad.Solve(model_with_field, 1e-3, 200)
    t_solve = time.time() - t0

    # Get stats
    stats = rad.GetSolveStats()

    print(f"NL iterations: {n_iter}")
    print(f"Linear iterations: {stats['total_linear_iter']}")
    print(f"Time: {t_solve:.2f}s")

    return {'name': name, 'nl': n_iter, 'linear': stats['total_linear_iter'], 'time': t_solve}

def main():
    print("Newton Damping Test - Simple Cube 6x6x6")

    results = []

    # Test 1: Block Jacobi (Picard)
    r1 = run_test("Block Jacobi (Picard)", use_newton=False, use_damping=False)
    results.append(r1)

    # Test 2: Newton without damping
    r2 = run_test("Newton (no damping)", use_newton=True, use_damping=False)
    results.append(r2)

    # Test 3: Newton with damping
    r3 = run_test("Newton + Damping", use_newton=True, use_damping=True)
    results.append(r3)

    # Summary
    print(f"\n{'='*50}")
    print("SUMMARY")
    print(f"{'='*50}")
    for r in results:
        print(f"{r['name']:<25} NL={r['nl']:<4} Linear={r['linear']:<6} Time={r['time']:.1f}s")

    print(f"\nDamping effect: {r2['nl']} → {r3['nl']} NL iterations")
    print(f"  ({(r3['nl']/r2['nl']-1)*100:+.0f}% change)")

if __name__ == '__main__':
    main()
