"""
Sequential benchmark for C-type electromagnet: 6x6x6, 10x10x10, 20x20x20
Compares LU, BiCGSTAB, HACApK, and HACApK+Newton for paper Table 2

Run with: python -u benchmark_ctype_sequential.py
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad
rad.FldUnits('m')
import time
import json

mm = 1e-3  # 1 mm in meters

# B-H curve data (100 points, same as paper)
bh_data = [
    [0, 0], [10, 0.01], [20, 0.021], [50, 0.053], [100, 0.11], [150, 0.17],
    [200, 0.24], [250, 0.32], [300, 0.42], [350, 0.53], [400, 0.66],
    [450, 0.81], [500, 0.97], [550, 1.14], [600, 1.29], [650, 1.42],
    [700, 1.52], [750, 1.59], [800, 1.65], [850, 1.69], [900, 1.73],
    [950, 1.76], [1000, 1.79], [1100, 1.83], [1200, 1.87], [1300, 1.90],
    [1400, 1.92], [1500, 1.94], [1600, 1.96], [1700, 1.98], [1800, 2.00],
    [1900, 2.01], [2000, 2.03], [2200, 2.05], [2400, 2.07], [2600, 2.09],
    [2800, 2.11], [3000, 2.13], [3500, 2.16], [4000, 2.19], [4500, 2.22],
    [5000, 2.24], [6000, 2.28], [7000, 2.31], [8000, 2.34], [9000, 2.37],
    [10000, 2.39], [12000, 2.43], [14000, 2.46], [16000, 2.49], [18000, 2.51],
    [20000, 2.53], [25000, 2.57], [30000, 2.60], [35000, 2.63], [40000, 2.65],
    [45000, 2.67], [50000, 2.69], [60000, 2.72], [70000, 2.75], [80000, 2.77],
    [90000, 2.79], [100000, 2.81], [120000, 2.84], [140000, 2.86], [160000, 2.88],
    [180000, 2.90], [200000, 2.92], [250000, 2.95], [300000, 2.98], [350000, 3.00],
    [400000, 3.02], [450000, 3.04], [500000, 3.06], [600000, 3.09], [700000, 3.11],
    [800000, 3.13], [900000, 3.15], [1000000, 3.17], [1200000, 3.20], [1400000, 3.22],
    [1600000, 3.24], [1800000, 3.26], [2000000, 3.28], [2500000, 3.31], [3000000, 3.34],
    [3500000, 3.36], [4000000, 3.38], [4500000, 3.40], [5000000, 3.42], [6000000, 3.45],
    [7000000, 3.47], [8000000, 3.49], [9000000, 3.51], [10000000, 3.53]
]

def create_ctype_electromagnet(nx=6, ny=6, nz=6):
    """Create C-type electromagnet model with specified mesh density"""

    # Geometry parameters (mm -> meters)
    yoke_width = 304.8 * mm  # 12 inch
    yoke_height = 304.8 * mm
    yoke_depth = 288.0 * mm  # Reduced from 304.8 to avoid wedge
    gap_height = 50.0 * mm
    pole_width = 100.0 * mm

    # Coil parameters
    coil_current = 20000  # AT (Ampere-turns)
    coil_turns = 100
    current_per_turn = coil_current / coil_turns  # 200 A

    # Create material
    mat = rad.MatSatIsoTab(bh_data)

    # Create yoke parts (1/4 model with IMA symmetry)
    # Using subdivision for specified mesh density

    # Lower yoke (horizontal)
    lower_yoke = rad.ObjRecMag(
        [yoke_width/4, yoke_depth/4, (yoke_height - gap_height)/4],
        [yoke_width/2, yoke_depth/2, (yoke_height - gap_height)/2],
        [0, 0, 0]
    )
    rad.MatApl(lower_yoke, mat)
    rad.ObjDivMag(lower_yoke, [nx, ny, nz], 'Frame->Lab')

    # Upper pole
    upper_pole = rad.ObjRecMag(
        [pole_width/4, yoke_depth/4, yoke_height/2 + gap_height/4],
        [pole_width/2, yoke_depth/2, (yoke_height - gap_height)/2],
        [0, 0, 0]
    )
    rad.MatApl(upper_pole, mat)
    rad.ObjDivMag(upper_pole, [nx//2, ny, nz], 'Frame->Lab')

    # Vertical yoke
    vertical_yoke = rad.ObjRecMag(
        [yoke_width/2 - (yoke_width - pole_width)/4, yoke_depth/4, yoke_height/4 + gap_height/4],
        [(yoke_width - pole_width)/2, yoke_depth/2, (yoke_height - gap_height)/2],
        [0, 0, 0]
    )
    rad.MatApl(vertical_yoke, mat)
    rad.ObjDivMag(vertical_yoke, [nx//3, ny, nz*2], 'Frame->Lab')

    # Group yoke parts
    yoke = rad.ObjCnt([lower_yoke, upper_pole, vertical_yoke])

    # Create racetrack coil (approximation with two rectangular current loops)
    coil_height = gap_height + 20 * mm
    coil_width = pole_width + 20 * mm

    # Horizontal segments (main field contribution)
    coil_h1 = rad.ObjRaceTrk(
        [coil_width/4, yoke_depth/4, coil_height/2],
        [coil_width/2, 0, 0],
        [0, yoke_depth/2, 0],
        10, current_per_turn, 'man', 'y'
    )

    # Combine with yoke
    container = rad.ObjCnt([yoke, coil_h1])

    return container

def run_single_benchmark(mesh_size, solver_type, use_newton=False):
    """
    Run single benchmark configuration

    mesh_size: 6, 10, or 20
    solver_type: 0=LU, 1=BiCGSTAB, 2=HACApK
    use_newton: use Newton-Raphson (only for HACApK)
    """

    nx = ny = nz = mesh_size

    solver_names = {0: "LU", 1: "BiCGSTAB", 2: "HACApK"}
    solver_name = solver_names[solver_type]
    if use_newton and solver_type == 2:
        solver_name += "+Newton"

    print(f"\n{'='*70}")
    print(f"Mesh: {nx}x{ny}x{nz}, Solver: {solver_name}")
    print(f"{'='*70}")

    # Create model
    container = create_ctype_electromagnet(nx, ny, nz)

    # Configure solver
    rad.SetSolver(solver_type)
    if solver_type >= 1:
        rad.SetBiCGSTABTolerance(1e-4)

    # Configure Newton method
    if use_newton and solver_type == 2:
        rad.SetNewtonMethod(True)
        rad.SetNewtonDamping(False)  # No damping for benchmark consistency
        print("Newton-Raphson: ENABLED")
    else:
        rad.SetNewtonMethod(False)

    # Build system matrix with IMA symmetry (YZ and XY planes)
    print("Building interaction matrix with IMA symmetry...")
    t_build_start = time.time()
    rad.BuildMatrix(container, image="+x+z")  # 1/4 model
    t_build = time.time() - t_build_start

    # Get statistics
    if solver_type == 2:
        stats = rad.GetHACApKStats()
        print(f"Matrix build time: {t_build:.2f} s")
        print(f"Matrix size: {stats['matrix_size_mb']:.1f} MB")
        print(f"Compression ratio: {stats['compression_ratio']*100:.1f}%")
        print(f"Peak memory: {stats['peak_memory_mb']:.1f} MB")
        mem_mb = stats['peak_memory_mb']
        compression = stats['compression_ratio']
    else:
        print(f"Matrix build time: {t_build:.2f} s")
        mem_mb = None
        compression = None

    # Solve nonlinear system
    print(f"Solving nonlinear system (tol=1e-3, max_iter=200)...")
    t_solve_start = time.time()
    n_iter = rad.Solve(container, 1e-3, 200)
    t_solve = time.time() - t_solve_start

    # Get field at gap center
    Bz = rad.Fld(container, 'bz', [0, 0, 0])

    # Get solver statistics
    solve_stats = rad.GetSolveStats()

    print(f"\nResults:")
    print(f"  Nonlinear iterations: {n_iter}")
    print(f"  Total linear iterations: {solve_stats['total_linear_iter']}")
    print(f"  Avg linear iter/NL: {solve_stats['total_linear_iter']/n_iter:.1f}")
    print(f"  Bz at gap center: {Bz*1e3:.2f} mT ({Bz:.4f} T)")
    print(f"  Solve time: {t_solve:.2f} s")
    print(f"  Total time: {t_build + t_solve:.2f} s")

    result = {
        'mesh': f"{nx}x{ny}x{nz}",
        'solver': solver_name,
        'nl_iter': n_iter,
        'linear_iter': solve_stats['total_linear_iter'],
        'Bz_T': Bz,
        'build_time_s': t_build,
        'solve_time_s': t_solve,
        'total_time_s': t_build + t_solve,
        'converged': n_iter < 200
    }

    if mem_mb is not None:
        result['memory_mb'] = mem_mb
        result['compression_ratio'] = compression

    return result

def main():
    print("="*70)
    print("C-type Electromagnet Sequential Benchmark (6x6x6, 10x10x10, 20x20x20)")
    print("="*70)

    all_results = []

    mesh_sizes = [6, 10, 20]

    for mesh_size in mesh_sizes:
        print(f"\n{'#'*70}")
        print(f"# Mesh: {mesh_size}x{mesh_size}x{mesh_size}")
        print(f"{'#'*70}")

        # Test 1: LU (skip for 20x20x20 to save time)
        if mesh_size <= 10:
            print("\n--- LU Solver ---")
            r1 = run_single_benchmark(mesh_size, solver_type=0, use_newton=False)
            all_results.append(r1)
        else:
            print("\n--- LU Solver: SKIP (expected to fail due to memory) ---")

        # Test 2: BiCGSTAB (skip for 20x20x20 to save time)
        if mesh_size <= 10:
            print("\n--- BiCGSTAB Solver ---")
            r2 = run_single_benchmark(mesh_size, solver_type=1, use_newton=False)
            all_results.append(r2)
        else:
            print("\n--- BiCGSTAB Solver: SKIP (expected to fail due to memory) ---")

        # Test 3: HACApK
        print("\n--- HACApK Solver ---")
        r3 = run_single_benchmark(mesh_size, solver_type=2, use_newton=False)
        all_results.append(r3)

        # Test 4: HACApK + Newton
        print("\n--- HACApK+Newton Solver ---")
        r4 = run_single_benchmark(mesh_size, solver_type=2, use_newton=True)
        all_results.append(r4)

    # Summary table
    print("\n" + "="*70)
    print("SUMMARY TABLE")
    print("="*70)
    print(f"{'Mesh':<12} {'Solver':<18} {'NL':<6} {'Lin':<8} {'Time(s)':<10} {'Mem(MB)':<10} {'Bz(T)':<10}")
    print("-"*70)
    for r in all_results:
        mem_str = f"{r['memory_mb']:.1f}" if 'memory_mb' in r else "---"
        print(f"{r['mesh']:<12} {r['solver']:<18} {r['nl_iter']:<6} {r['linear_iter']:<8} "
              f"{r['solve_time_s']:<10.1f} {mem_str:<10} {r['Bz_T']:<10.4f}")

    # Save results
    output_file = 'ctype_benchmark_results.json'
    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to: {output_file}")

if __name__ == '__main__':
    main()
