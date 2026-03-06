"""
Quick test of Newton line search damping using existing C-type model

Tests 6x6x6 mesh with three configurations:
1. HACApK + Block Jacobi
2. HACApK + Newton (no damping)
3. HACApK + Newton + Damping
"""

import sys
sys.path.insert(0, r'S:\Radia\01_GitHub\src\radia')

import radia as rad
import time

# B-H curve (simplified, 20 points)
bh_data = [
    [0, 0], [100, 0.11], [500, 0.97], [1000, 1.79], [2000, 2.03],
    [5000, 2.24], [10000, 2.39], [20000, 2.53], [50000, 2.69], [100000, 2.81],
    [200000, 2.92], [500000, 3.06], [1000000, 3.17], [2000000, 3.28], [5000000, 3.42],
    [10000000, 3.53]
]

def create_simple_ctype(mesh_divs=6):
    """Create simplified C-type model using hexahedra"""

    mat = rad.MatSatIsoTab(bh_data)

    # Simple C-shape using hexahedra
    # Base dimensions (mm)
    width = 300.0
    depth = 300.0
    height = 300.0
    gap = 50.0

    all_objects = []

    # Create lower yoke (simple rectangular block)
    dx = width / mesh_divs
    dy = depth / mesh_divs
    dz = (height - gap) / (2 * mesh_divs)

    # Lower yoke hexahedra
    for ix in range(mesh_divs):
        for iy in range(mesh_divs):
            for iz in range(mesh_divs):
                x = ix * dx
                y = iy * dy
                z = iz * dz

                verts = [
                    [x, y, z],
                    [x+dx, y, z],
                    [x+dx, y+dy, z],
                    [x, y+dy, z],
                    [x, y, z+dz],
                    [x+dx, y, z+dz],
                    [x+dx, y+dy, z+dz],
                    [x, y+dy, z+dz]
                ]

                obj = rad.ObjHexahedron(verts, [0, 0, 0])
                rad.MatApl(obj, mat)
                all_objects.append(obj)

    yoke = rad.ObjCnt(all_objects)

    # Add simple coil
    coil = rad.ObjRaceTrk(
        [width/2, depth/2, height/2],
        [width/4, 0, 0],
        [0, depth/2, 0],
        10, 200, 'man', 'y'
    )

    model = rad.ObjCnt([yoke, coil])
    return model

def run_single_test(config_name, use_newton=False, use_damping=False):
    """Run single configuration"""

    print(f"\n{'='*60}")
    print(f"Configuration: {config_name}")
    print(f"{'='*60}")

    # Create model
    model = create_simple_ctype(mesh_divs=6)

    # Configure solver
    rad.SetSolver(2)  # HACApK
    rad.SetBiCGSTABTolerance(1e-4)
    rad.SetNewtonMethod(use_newton)

    if use_damping:
        rad.SetNewtonDamping(True, max_iter=5, min_omega=0.01)
        print("Newton damping: ENABLED")
    else:
        rad.SetNewtonDamping(False)
        print("Newton damping: DISABLED")

    # Build matrix
    print("Building matrix...")
    t0 = time.time()
    rad.BuildMatrix(model)
    t_build = time.time() - t0

    stats = rad.GetHACApKStats()
    print(f"Build time: {t_build:.2f}s")
    print(f"Matrix: {stats['matrix_size_mb']:.1f} MB, compression: {stats['compression_ratio']*100:.1f}%")

    # Solve
    print("Solving...")
    t0 = time.time()
    n_iter = rad.Solve(model, 1e-3, 200)
    t_solve = time.time() - t0

    # Results
    Bz = rad.Fld(model, 'bz', [0, 0, 0])
    solve_stats = rad.GetSolveStats()

    print(f"\nResults:")
    print(f"  NL iterations: {n_iter}")
    print(f"  Linear iterations: {solve_stats['total_linear_iter']}")
    print(f"  Avg linear/NL: {solve_stats['total_linear_iter']/n_iter:.1f}")
    print(f"  Bz: {Bz:.2f} mT")
    print(f"  Solve time: {t_solve:.2f}s")

    return {
        'config': config_name,
        'nl_iter': n_iter,
        'linear_iter': solve_stats['total_linear_iter'],
        'Bz': Bz,
        'time': t_solve
    }

def main():
    print("="*60)
    print("Quick Newton Damping Test (Simplified C-type 6x6x6)")
    print("="*60)

    results = []

    # Test 1: Baseline (Block Jacobi)
    r1 = run_single_test("HACApK+BlockJacobi", use_newton=False, use_damping=False)
    results.append(r1)

    # Test 2: Newton without damping
    r2 = run_single_test("HACApK+Newton", use_newton=True, use_damping=False)
    results.append(r2)

    # Test 3: Newton with damping
    r3 = run_single_test("HACApK+Newton+Damping", use_newton=True, use_damping=True)
    results.append(r3)

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"{'Config':<25} {'NL':<8} {'Linear':<10} {'Time':<8} {'Bz (mT)'}")
    print("-"*60)
    for r in results:
        print(f"{r['config']:<25} {r['nl_iter']:<8} {r['linear_iter']:<10} {r['time']:<8.1f} {r['Bz']:.2f}")

    # Analysis
    print(f"\n{'='*60}")
    print("ANALYSIS")
    print(f"{'='*60}")

    print(f"\nNewton (no damp) vs. Baseline:")
    print(f"  NL iter: {r1['nl_iter']} → {r2['nl_iter']} ({(r2['nl_iter']/r1['nl_iter']-1)*100:+.0f}%)")
    print(f"  Time: {r1['time']:.1f}s → {r2['time']:.1f}s ({(r2['time']/r1['time']-1)*100:+.0f}%)")

    print(f"\nNewton+Damping vs. Newton:")
    print(f"  NL iter: {r2['nl_iter']} → {r3['nl_iter']} ({(r3['nl_iter']/r2['nl_iter']-1)*100:+.0f}%)")
    print(f"  Time: {r2['time']:.1f}s → {r3['time']:.1f}s ({(r3['time']/r1['time']-1)*100:+.0f}%)")

    print(f"\n{'='*60}")
    print("Expected: Damping should reduce NL iterations from Newton")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
