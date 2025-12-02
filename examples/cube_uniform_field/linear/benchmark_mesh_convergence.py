#!/usr/bin/env python3
"""
Mesh Convergence Benchmark: Tetrahedra vs Hexahedra

Compare magnetization accuracy for different mesh sizes:
- Hexahedral mesh (using ObjDivMag subdivision)
- Tetrahedral mesh (using NGSolve/Netgen mesh import)

Test case: Linear magnetic cube in uniform external field
Analytical solution: M = chi * H0 / (1 + chi * N) where N = 1/3 for cube
"""

import sys
import os

# Ensure unbuffered output
sys.stdout.reconfigure(line_buffering=True)

# Add paths - build/Release first (radia.pyd), then src/python (Python utilities)
# Note: insert(0, ...) prepends, so add in reverse order
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/python'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build/Release'))

import numpy as np
import radia as rad
import time

# Physical constants
mu_0 = 4 * np.pi * 1e-7  # T/(A/m)

# Test parameters
CUBE_SIZE = 0.1  # 10cm cube (in meters)
MU_R = 100.0     # Relative permeability
CHI = MU_R - 1.0  # Susceptibility
B0 = 1.0         # Applied field (Tesla)
H0 = B0 / mu_0   # Applied H-field (A/m)

# Analytical solution for cube
N_CUBE = 1.0 / 3.0  # Demagnetization factor for cube
M_ANALYTICAL = CHI * H0 / (1 + CHI * N_CUBE)

print("=" * 80)
print("MESH CONVERGENCE BENCHMARK: Tetrahedra vs Hexahedra")
print("=" * 80)
print(f"Cube size: {CUBE_SIZE*1000:.0f} mm")
print(f"mu_r: {MU_R:.0f}")
print(f"chi: {CHI:.0f}")
print(f"B0: {B0:.1f} T")
print(f"H0: {H0:.0f} A/m")
print(f"N (demagnetization factor): {N_CUBE:.4f}")
print(f"M_analytical: {M_ANALYTICAL/1e6:.4f} MA/m")
print("=" * 80)


def uniform_B(pos):
    """Background field callback."""
    return [0, 0, B0]


def get_average_Mz(obj, n_elements):
    """Get average M_z from ObjM result, handling both single and subdivided objects.

    ObjM returns:
    - Single element: [[pos], [M]] -> M_z = M[1][2]
    - Subdivided/container: [[[pos],[M]], [[pos],[M]], ...] -> M_z = M[i][1][2]
    """
    M = rad.ObjM(obj)
    if n_elements == 1:
        return M[1][2]
    else:
        # Subdivided/container: M is list of [[pos], [M]] for each child
        return np.mean([m[1][2] for m in M])


def benchmark_hexahedra(subdivisions_list):
    """Benchmark hexahedral mesh with different subdivision levels."""
    results = []

    for subdivs in subdivisions_list:
        rad.UtiDelAll()
        rad.FldUnits('m')

        # Create cube with zero initial magnetization
        cube = rad.ObjRecMag([0, 0, 0], [CUBE_SIZE]*3, [0, 0, 0])

        # Subdivide
        if subdivs > 1:
            rad.ObjDivMag(cube, [subdivs, subdivs, subdivs])

        n_elements = subdivs ** 3

        # Apply material
        mat = rad.MatLin(CHI)
        rad.MatApl(cube, mat)

        # Background field
        bg = rad.ObjBckgCF(uniform_B)
        system = rad.ObjCnt([cube, bg])

        # Solve
        t0 = time.time()
        result = rad.Solve(system, 0.0001, 1000)
        solve_time = time.time() - t0

        # Get magnetization (average over all sub-elements)
        M_z = get_average_Mz(cube, n_elements)

        # Compute error
        error_pct = abs(M_z - M_ANALYTICAL) / M_ANALYTICAL * 100

        results.append({
            'subdivs': subdivs,
            'n_elements': n_elements,
            'M_z': M_z,
            'error_pct': error_pct,
            'solve_time': solve_time,
            'iterations': result[0] if isinstance(result, (list, tuple)) else result
        })

        print(f"  {subdivs}x{subdivs}x{subdivs} ({n_elements:4d} elem): "
              f"M_z = {M_z/1e6:.4f} MA/m, error = {error_pct:.2f}%, "
              f"time = {solve_time:.3f}s")

    return results


def benchmark_tetrahedra_ngsolve(maxh_list):
    """Benchmark tetrahedral mesh using NGSolve/Netgen mesh generation.

    Uses Method 9 (LU direct solver) because relaxation methods don't work
    well for tetrahedral elements. This matches ELF_MAGIC's approach.
    """
    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh
        from netgen_mesh_import import netgen_mesh_to_radia, TETRA_FACES
    except ImportError as e:
        print(f"  NGSolve/Netgen not available: {e}")
        return []

    results = []
    hs = CUBE_SIZE / 2

    # Use original Radia method for tetrahedra
    rad.SolverTetraMethod(0)

    for maxh in maxh_list:
        rad.UtiDelAll()
        rad.FldUnits('m')

        # Create cube geometry using OCC
        cube_solid = Box(Pnt(-hs, -hs, -hs), Pnt(hs, hs, hs))
        cube_solid.mat('magnetic')
        cube_geo = OCCGeometry(cube_solid)

        # Generate mesh
        try:
            ngmesh = cube_geo.GenerateMesh(maxh=maxh)
            mesh = Mesh(ngmesh)
            n_elements = mesh.ne
        except Exception as e:
            print(f"  maxh={maxh}: Mesh generation failed: {e}")
            continue

        # Import to Radia
        try:
            mag_obj = netgen_mesh_to_radia(mesh,
                                           material={'magnetization': [0, 0, 0]},
                                           units='m',
                                           material_filter='magnetic')
        except Exception as e:
            print(f"  maxh={maxh}: Radia import failed: {e}")
            continue

        # Apply material
        mat = rad.MatLin(CHI)
        rad.MatApl(mag_obj, mat)

        # Background field
        bg = rad.ObjBckgCF(uniform_B)
        system = rad.ObjCnt([mag_obj, bg])

        # Solve using Method 9 (LU direct solver)
        # Relaxation methods don't converge well for tetrahedra
        t0 = time.time()
        try:
            # Build interaction matrix
            intrc = rad.RlxPre(system, system)
            # Use Method 9 (LU direct solver)
            result = rad.RlxMan(intrc, 9, 1, 1.0)
            solve_time = time.time() - t0
        except Exception as e:
            print(f"  maxh={maxh}: Solve failed: {e}")
            continue

        # Get magnetization (average over all elements)
        M_z = get_average_Mz(mag_obj, n_elements)

        # Check for NaN
        if np.isnan(M_z):
            print(f"  maxh={maxh:.3f} ({n_elements:4d} elem): NaN result")
            continue

        # Compute error
        error_pct = abs(M_z - M_ANALYTICAL) / M_ANALYTICAL * 100

        results.append({
            'maxh': maxh,
            'n_elements': n_elements,
            'M_z': M_z,
            'error_pct': error_pct,
            'solve_time': solve_time,
            'iterations': 1  # LU solver is direct - 1 iteration
        })

        print(f"  maxh={maxh:.3f} ({n_elements:4d} elem): "
              f"M_z = {M_z/1e6:.4f} MA/m, error = {error_pct:.2f}%, "
              f"time = {solve_time:.3f}s")

    return results


def benchmark_tetrahedra_manual(n_subdivs_list):
    """Benchmark tetrahedral mesh by subdividing a cube manually into 5 tetrahedra.

    Uses Method 9 (LU direct solver) because relaxation methods don't work
    well for tetrahedral elements. This matches ELF_MAGIC's approach.
    """
    results = []
    hs = CUBE_SIZE / 2

    # Use original Radia method for tetrahedra
    rad.SolverTetraMethod(0)

    # 5-tetrahedra decomposition of unit cube (scaled to actual size)
    # This is a standard decomposition that preserves symmetry
    TETRA_FACES = [[1, 3, 2], [1, 2, 4], [2, 3, 4], [3, 1, 4]]

    def create_cube_tetrahedra(cx, cy, cz, size):
        """Create 5 tetrahedra from a cube centered at (cx, cy, cz)."""
        h = size / 2
        # Cube vertices
        v = [
            [cx - h, cy - h, cz - h],  # 0: ---
            [cx + h, cy - h, cz - h],  # 1: +--
            [cx - h, cy + h, cz - h],  # 2: -+-
            [cx + h, cy + h, cz - h],  # 3: ++-
            [cx - h, cy - h, cz + h],  # 4: --+
            [cx + h, cy - h, cz + h],  # 5: +-+
            [cx - h, cy + h, cz + h],  # 6: -++
            [cx + h, cy + h, cz + h],  # 7: +++
        ]

        # 5 tetrahedra decomposition (Freudenthal triangulation)
        tetra_indices = [
            [0, 1, 3, 5],  # Tetra 1
            [0, 3, 2, 6],  # Tetra 2
            [0, 5, 4, 6],  # Tetra 3
            [3, 5, 6, 7],  # Tetra 4
            [0, 3, 5, 6],  # Tetra 5 (center)
        ]

        tetrahedra = []
        for indices in tetra_indices:
            verts = [v[i] for i in indices]
            tetrahedra.append(verts)

        return tetrahedra

    for n_subdivs in n_subdivs_list:
        rad.UtiDelAll()
        rad.FldUnits('m')

        # Create subdivided tetrahedra
        cell_size = CUBE_SIZE / n_subdivs
        all_tetra_objs = []

        for ix in range(n_subdivs):
            for iy in range(n_subdivs):
                for iz in range(n_subdivs):
                    cx = -hs + cell_size/2 + ix * cell_size
                    cy = -hs + cell_size/2 + iy * cell_size
                    cz = -hs + cell_size/2 + iz * cell_size

                    tetras = create_cube_tetrahedra(cx, cy, cz, cell_size)
                    for verts in tetras:
                        obj = rad.ObjPolyhdr(verts, TETRA_FACES, [0, 0, 0])
                        all_tetra_objs.append(obj)

        n_elements = len(all_tetra_objs)

        # Group all tetrahedra
        mag_obj = rad.ObjCnt(all_tetra_objs)

        # Apply material
        mat = rad.MatLin(CHI)
        rad.MatApl(mag_obj, mat)

        # Background field
        bg = rad.ObjBckgCF(uniform_B)
        system = rad.ObjCnt([mag_obj, bg])

        # Solve using Method 9 (LU direct solver)
        # Relaxation methods don't converge well for tetrahedra
        t0 = time.time()
        try:
            # Build interaction matrix
            intrc = rad.RlxPre(system, system)
            # Use Method 9 (LU direct solver)
            result = rad.RlxMan(intrc, 9, 1, 1.0)
            solve_time = time.time() - t0
        except Exception as e:
            print(f"  {n_subdivs}^3 x 5 ({n_elements:4d} elem): Solve failed: {e}")
            continue

        # Get magnetization (average over all elements)
        M_z = get_average_Mz(mag_obj, n_elements)

        # Check for NaN
        if np.isnan(M_z):
            print(f"  {n_subdivs}^3 x 5 ({n_elements:4d} elem): NaN result")
            continue

        # Compute error
        error_pct = abs(M_z - M_ANALYTICAL) / M_ANALYTICAL * 100

        results.append({
            'n_subdivs': n_subdivs,
            'n_elements': n_elements,
            'M_z': M_z,
            'error_pct': error_pct,
            'solve_time': solve_time,
            'iterations': 1  # LU solver is direct - 1 iteration
        })

        print(f"  {n_subdivs}^3 x 5 ({n_elements:4d} elem): "
              f"M_z = {M_z/1e6:.4f} MA/m, error = {error_pct:.2f}%, "
              f"time = {solve_time:.3f}s")

    return results


def print_results_table(hex_results, tetra_results, tetra_ngsolve_results):
    """Print formatted results table."""
    print("\n" + "=" * 80)
    print("RESULTS SUMMARY")
    print("=" * 80)
    print(f"Analytical M_z = {M_ANALYTICAL/1e6:.4f} MA/m")
    print()

    # Hexahedral results
    print("HEXAHEDRAL MESH (ObjRecMag + ObjDivMag)")
    print("-" * 70)
    print(f"{'Subdivisions':<15} {'Elements':<10} {'M_z (MA/m)':<12} {'Error (%)':<12} {'Time (s)':<10}")
    print("-" * 70)
    for r in hex_results:
        print(f"{r['subdivs']}x{r['subdivs']}x{r['subdivs']:<10} {r['n_elements']:<10} "
              f"{r['M_z']/1e6:<12.4f} {r['error_pct']:<12.2f} {r['solve_time']:<10.3f}")
    print()

    # Tetrahedral results (manual)
    if tetra_results:
        print("TETRAHEDRAL MESH (Manual 5-tetra decomposition)")
        print("-" * 70)
        print(f"{'Subdivisions':<15} {'Elements':<10} {'M_z (MA/m)':<12} {'Error (%)':<12} {'Time (s)':<10}")
        print("-" * 70)
        for r in tetra_results:
            print(f"{r['n_subdivs']}^3 x 5{'':<8} {r['n_elements']:<10} "
                  f"{r['M_z']/1e6:<12.4f} {r['error_pct']:<12.2f} {r['solve_time']:<10.3f}")
        print()

    # Tetrahedral results (NGSolve)
    if tetra_ngsolve_results:
        print("TETRAHEDRAL MESH (NGSolve/Netgen)")
        print("-" * 70)
        print(f"{'maxh (m)':<15} {'Elements':<10} {'M_z (MA/m)':<12} {'Error (%)':<12} {'Time (s)':<10}")
        print("-" * 70)
        for r in tetra_ngsolve_results:
            print(f"{r['maxh']:<15.3f} {r['n_elements']:<10} "
                  f"{r['M_z']/1e6:<12.4f} {r['error_pct']:<12.2f} {r['solve_time']:<10.3f}")
        print()

    print("=" * 80)


def main():
    print("\n" + "-" * 80)
    print("1. HEXAHEDRAL MESH BENCHMARK")
    print("-" * 80)
    hex_subdivs = [1, 2, 3, 4, 5]
    hex_results = benchmark_hexahedra(hex_subdivs)

    print("\n" + "-" * 80)
    print("2. TETRAHEDRAL MESH BENCHMARK (Manual 5-tetra decomposition)")
    print("-" * 80)
    tetra_subdivs = [1, 2, 3]
    tetra_results = benchmark_tetrahedra_manual(tetra_subdivs)

    print("\n" + "-" * 80)
    print("3. TETRAHEDRAL MESH BENCHMARK (NGSolve/Netgen)")
    print("-" * 80)
    maxh_list = [0.1, 0.08, 0.06, 0.05, 0.04]  # in meters (100mm to 40mm)
    tetra_ngsolve_results = benchmark_tetrahedra_ngsolve(maxh_list)

    # Print summary table
    print_results_table(hex_results, tetra_results, tetra_ngsolve_results)


if __name__ == '__main__':
    main()
