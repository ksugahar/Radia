"""
Test Newton-Raphson line search damping on C-type electromagnet (6x6x6 mesh)

Compares three configurations:
1. HACApK + Block Jacobi (Picard)
2. HACApK + Newton (no damping)
3. HACApK + Newton + Line Search Damping

Expected: Damping reduces NL iterations while maintaining accuracy
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

def run_test(nx=6, ny=6, nz=6, use_newton=False, use_damping=False,
             max_iter=200, prec_tol=1e-4):
    """Run single test configuration"""

    print(f"\n{'='*70}")
    config_name = "HACApK"
    if use_newton:
        config_name += "+Newton"
        if use_damping:
            config_name += "+Damping"
    else:
        config_name += "+BlockJacobi"
    print(f"Configuration: {config_name}")
    print(f"Mesh: {nx}x{ny}x{nz}")
    print(f"{'='*70}")

    # Create model
    container = create_ctype_electromagnet(nx, ny, nz)

    # Configure solver
    rad.SetSolver(2)  # HACApK
    rad.SetBiCGSTABTolerance(prec_tol)

    # Configure Newton method
    rad.SetNewtonMethod(use_newton)

    # Configure Newton damping
    if use_damping:
        rad.SetNewtonDamping(True, max_iter=5, min_omega=0.01)
        print("Newton damping: ENABLED (max_iter=5, min_omega=0.01)")
    else:
        rad.SetNewtonDamping(False)
        print("Newton damping: DISABLED")

    # Build system matrix with IMA symmetry (YZ and XY planes)
    print("\nBuilding interaction matrix with IMA symmetry...")
    t_build_start = time.time()
    rad.BuildMatrix(container, image="+x+z")  # 1/4 model
    t_build = time.time() - t_build_start

    # Get HACApK statistics
    stats = rad.GetHACApKStats()
    print(f"Matrix build time: {t_build:.2f} s")
    print(f"Matrix size: {stats['matrix_size_mb']:.1f} MB")
    print(f"Compression ratio: {stats['compression_ratio']*100:.1f}%")
    print(f"Peak memory: {stats['peak_memory_mb']:.1f} MB")

    # Solve nonlinear system
    print(f"\nSolving nonlinear system (tol=1e-3, max_iter={max_iter})...")
    t_solve_start = time.time()
    n_iter = rad.Solve(container, 1e-3, max_iter)
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

    # Get Newton damping statistics if enabled
    if use_damping:
        damp_stats = rad.GetNewtonDampingStats()
        print(f"\nNewton Damping Stats:")
        print(f"  Enabled: {damp_stats['enabled']}")
        print(f"  Max iterations: {damp_stats['max_iter']}")
        print(f"  Min omega: {damp_stats['min_omega']}")

    return {
        'config': config_name,
        'mesh': f"{nx}x{ny}x{nz}",
        'dof': stats['total_dof'],
        'elements': stats['num_elements'],
        'nl_iter': n_iter,
        'linear_iter': solve_stats['total_linear_iter'],
        'avg_linear_per_nl': solve_stats['total_linear_iter']/n_iter,
        'Bz_T': Bz,
        'build_time_s': t_build,
        'solve_time_s': t_solve,
        'total_time_s': t_build + t_solve,
        'matrix_size_mb': stats['matrix_size_mb'],
        'compression_ratio': stats['compression_ratio'],
        'peak_memory_mb': stats['peak_memory_mb'],
        'converged': n_iter < max_iter
    }

def main():
    """Run all test configurations"""

    print("="*70)
    print("Testing Newton Line Search Damping on C-type Electromagnet (6x6x6)")
    print("="*70)

    results = []

    # Test 1: HACApK + Block Jacobi (baseline)
    print("\n" + "="*70)
    print("TEST 1: Baseline (Block Jacobi, no Newton)")
    print("="*70)
    r1 = run_test(nx=6, ny=6, nz=6, use_newton=False, use_damping=False)
    results.append(r1)

    # Test 2: HACApK + Newton (no damping)
    print("\n" + "="*70)
    print("TEST 2: Newton without damping")
    print("="*70)
    r2 = run_test(nx=6, ny=6, nz=6, use_newton=True, use_damping=False)
    results.append(r2)

    # Test 3: HACApK + Newton + Damping
    print("\n" + "="*70)
    print("TEST 3: Newton WITH line search damping")
    print("="*70)
    r3 = run_test(nx=6, ny=6, nz=6, use_newton=True, use_damping=True)
    results.append(r3)

    # Summary table
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"{'Configuration':<30} {'NL Iter':<10} {'Lin Iter':<10} {'Time (s)':<10} {'Bz (T)':<10}")
    print("-"*70)
    for r in results:
        print(f"{r['config']:<30} {r['nl_iter']:<10} {r['linear_iter']:<10} "
              f"{r['solve_time_s']:<10.1f} {r['Bz_T']:<10.4f}")

    # Performance analysis
    print("\n" + "="*70)
    print("PERFORMANCE ANALYSIS")
    print("="*70)

    baseline = results[0]
    newton_no_damp = results[1]
    newton_damp = results[2]

    print(f"\nNewton (no damping) vs. Baseline:")
    print(f"  NL iterations: {baseline['nl_iter']} → {newton_no_damp['nl_iter']} "
          f"({(newton_no_damp['nl_iter']/baseline['nl_iter']-1)*100:+.1f}%)")
    print(f"  Linear iterations: {baseline['linear_iter']} → {newton_no_damp['linear_iter']} "
          f"({(newton_no_damp['linear_iter']/baseline['linear_iter']-1)*100:+.1f}%)")
    print(f"  Time: {baseline['solve_time_s']:.1f}s → {newton_no_damp['solve_time_s']:.1f}s "
          f"({(newton_no_damp['solve_time_s']/baseline['solve_time_s']-1)*100:+.1f}%)")

    print(f"\nNewton (WITH damping) vs. Newton (no damping):")
    print(f"  NL iterations: {newton_no_damp['nl_iter']} → {newton_damp['nl_iter']} "
          f"({(newton_damp['nl_iter']/newton_no_damp['nl_iter']-1)*100:+.1f}%)")
    print(f"  Linear iterations: {newton_no_damp['linear_iter']} → {newton_damp['linear_iter']} "
          f"({(newton_damp['linear_iter']/newton_no_damp['linear_iter']-1)*100:+.1f}%)")
    print(f"  Time: {newton_no_damp['solve_time_s']:.1f}s → {newton_damp['solve_time_s']:.1f}s "
          f"({(newton_damp['solve_time_s']/newton_no_damp['solve_time_s']-1)*100:+.1f}%)")

    print(f"\nNewton (WITH damping) vs. Baseline:")
    print(f"  NL iterations: {baseline['nl_iter']} → {newton_damp['nl_iter']} "
          f"({(newton_damp['nl_iter']/baseline['nl_iter']-1)*100:+.1f}%)")
    print(f"  Linear iterations: {baseline['linear_iter']} → {newton_damp['linear_iter']} "
          f"({(newton_damp['linear_iter']/baseline['linear_iter']-1)*100:+.1f}%)")
    print(f"  Time: {baseline['solve_time_s']:.1f}s → {newton_damp['solve_time_s']:.1f}s "
          f"({(newton_damp['solve_time_s']/baseline['solve_time_s']-1)*100:+.1f}%)")

    # Accuracy check
    print(f"\nAccuracy (all Bz values should be within ~0.0001 T):")
    print(f"  Baseline: {baseline['Bz_T']:.4f} T")
    print(f"  Newton (no damp): {newton_no_damp['Bz_T']:.4f} T (Δ={abs(newton_no_damp['Bz_T']-baseline['Bz_T']):.4f})")
    print(f"  Newton (damping): {newton_damp['Bz_T']:.4f} T (Δ={abs(newton_damp['Bz_T']-baseline['Bz_T']):.4f})")

    # Save results
    output_file = 'newton_damping_6x6x6_results.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {output_file}")

    print("\n" + "="*70)
    print("Expected behavior:")
    print("  - Newton (no damp) may INCREASE NL iterations vs. baseline")
    print("  - Newton (WITH damp) should REDUCE NL iterations vs. Newton (no damp)")
    print("  - Best case: Newton+Damping achieves lower NL iter than baseline")
    print("="*70)

if __name__ == '__main__':
    main()
