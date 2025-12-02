#!/usr/bin/env python3
"""
External Field Benchmark: Tetrahedra vs Hexahedra

Compare EXTERNAL field accuracy for different mesh types and sizes.
This is the proper metric for MMM (Magnetic Moment Method) validation,
as MMM computes field accurately outside magnetic materials.

Test case: Linear magnetic cube in uniform external field
Reference: Fine hexahedral mesh (8x8x8 = 512 elements)
"""

import sys
import os

# Ensure unbuffered output
sys.stdout.reconfigure(line_buffering=True)

# Add paths - build/Release first (radia.pyd), then src/python (Python utilities)
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
print("EXTERNAL FIELD BENCHMARK: Tetrahedra vs Hexahedra")
print("=" * 80)
print(f"Cube size: {CUBE_SIZE*1000:.0f} mm")
print(f"mu_r: {MU_R:.0f}")
print(f"chi: {CHI:.0f}")
print(f"B0: {B0:.1f} T")
print(f"H0: {H0:.0f} A/m")
print(f"M_analytical (for reference): {M_ANALYTICAL/1e6:.4f} MA/m")
print("=" * 80)
print()
print("NOTE: External field (outside cube) is the proper metric for MMM validation.")
print("      Internal magnetization distribution may differ between methods,")
print("      but external field should converge to same values.")
print()


def uniform_B(pos):
    """Background field callback."""
    return [0, 0, B0]


# Observation points OUTSIDE the cube (cube has vertices at +/- 0.05m)
OBS_POINTS = [
    ([0.10, 0, 0], "+x, r=0.10m"),
    ([0.12, 0, 0], "+x, r=0.12m"),
    ([0.15, 0, 0], "+x, r=0.15m"),
    ([0, 0, 0.10], "+z, r=0.10m"),
    ([0, 0, 0.12], "+z, r=0.12m"),
    ([0, 0, 0.15], "+z, r=0.15m"),
    ([0.10, 0.10, 0], "xy diagonal"),
]


def get_reference_fields(n_ref=8):
    """Compute reference field values using fine hexahedral mesh."""
    rad.UtiDelAll()
    rad.FldUnits('m')

    cube_ref = rad.ObjRecMag([0, 0, 0], [CUBE_SIZE]*3, [0, 0, 0])
    rad.ObjDivMag(cube_ref, [n_ref, n_ref, n_ref])
    mat = rad.MatLin(CHI)
    rad.MatApl(cube_ref, mat)
    bg = rad.ObjBckgCF(uniform_B)
    system_ref = rad.ObjCnt([cube_ref, bg])
    rad.Solve(system_ref, 0.00001, 2000)

    H_ref = {}
    for pt, desc in OBS_POINTS:
        H_ref[tuple(pt)] = np.array(rad.Fld(system_ref, 'h', pt))

    return H_ref, n_ref**3


def benchmark_hexahedra(subdivisions_list, H_ref):
    """Benchmark hexahedral mesh with different subdivision levels."""
    results = []

    for subdivs in subdivisions_list:
        rad.UtiDelAll()
        rad.FldUnits('m')

        cube = rad.ObjRecMag([0, 0, 0], [CUBE_SIZE]*3, [0, 0, 0])
        if subdivs > 1:
            rad.ObjDivMag(cube, [subdivs, subdivs, subdivs])

        n_elements = subdivs ** 3
        mat = rad.MatLin(CHI)
        rad.MatApl(cube, mat)

        bg = rad.ObjBckgCF(uniform_B)
        system = rad.ObjCnt([cube, bg])

        t0 = time.time()
        result = rad.Solve(system, 0.00001, 2000)
        solve_time = time.time() - t0

        # Compute external field errors
        errors = []
        for pt, desc in OBS_POINTS:
            H = np.array(rad.Fld(system, 'h', pt))
            H_r = H_ref[tuple(pt)]
            err = np.linalg.norm(H - H_r) / np.linalg.norm(H_r) * 100
            errors.append(err)

        avg_error = np.mean(errors)
        max_error = max(errors)

        results.append({
            'subdivs': subdivs,
            'n_elements': n_elements,
            'avg_error': avg_error,
            'max_error': max_error,
            'solve_time': solve_time,
        })

        print(f"  {subdivs}x{subdivs}x{subdivs} ({n_elements:4d} elem): "
              f"Avg error = {avg_error:.2f}%, Max = {max_error:.2f}%, "
              f"time = {solve_time:.3f}s")

    return results


def benchmark_tetrahedra_manual(n_subdivs_list, H_ref):
    """Benchmark tetrahedral mesh by subdividing a cube manually into 5 tetrahedra."""
    results = []
    hs = CUBE_SIZE / 2

    rad.SolverTetraMethod(0)
    TETRA_FACES = [[1, 3, 2], [1, 2, 4], [2, 3, 4], [3, 1, 4]]

    def create_cube_tetrahedra(cx, cy, cz, size):
        h = size / 2
        v = [
            [cx - h, cy - h, cz - h],
            [cx + h, cy - h, cz - h],
            [cx - h, cy + h, cz - h],
            [cx + h, cy + h, cz - h],
            [cx - h, cy - h, cz + h],
            [cx + h, cy - h, cz + h],
            [cx - h, cy + h, cz + h],
            [cx + h, cy + h, cz + h],
        ]
        tetra_indices = [[0,1,3,5], [0,3,2,6], [0,5,4,6], [3,5,6,7], [0,3,5,6]]
        return [[v[i] for i in indices] for indices in tetra_indices]

    for n_subdivs in n_subdivs_list:
        rad.UtiDelAll()
        rad.FldUnits('m')

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
        mag_obj = rad.ObjCnt(all_tetra_objs)
        mat = rad.MatLin(CHI)
        rad.MatApl(mag_obj, mat)

        bg = rad.ObjBckgCF(uniform_B)
        system = rad.ObjCnt([mag_obj, bg])

        t0 = time.time()
        try:
            result = rad.Solve(system, 0.00001, 2000)
            solve_time = time.time() - t0
        except Exception as e:
            print(f"  {n_subdivs}^3 x 5 ({n_elements:4d} elem): Solve failed: {e}")
            continue

        # Compute external field errors
        errors = []
        for pt, desc in OBS_POINTS:
            H = np.array(rad.Fld(system, 'h', pt))
            H_r = H_ref[tuple(pt)]
            err = np.linalg.norm(H - H_r) / np.linalg.norm(H_r) * 100
            errors.append(err)

        avg_error = np.mean(errors)
        max_error = max(errors)

        results.append({
            'n_subdivs': n_subdivs,
            'n_elements': n_elements,
            'avg_error': avg_error,
            'max_error': max_error,
            'solve_time': solve_time,
        })

        print(f"  {n_subdivs}^3 x 5 ({n_elements:4d} elem): "
              f"Avg error = {avg_error:.2f}%, Max = {max_error:.2f}%, "
              f"time = {solve_time:.3f}s")

    return results


def benchmark_tetrahedra_ngsolve(maxh_list, H_ref):
    """Benchmark tetrahedral mesh using NGSolve/Netgen mesh generation."""
    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh
        from netgen_mesh_import import netgen_mesh_to_radia
    except ImportError as e:
        print(f"  NGSolve/Netgen not available: {e}")
        return []

    results = []
    hs = CUBE_SIZE / 2

    rad.SolverTetraMethod(0)

    for maxh in maxh_list:
        rad.UtiDelAll()
        rad.FldUnits('m')

        cube_solid = Box(Pnt(-hs, -hs, -hs), Pnt(hs, hs, hs))
        cube_solid.mat('magnetic')
        cube_geo = OCCGeometry(cube_solid)

        try:
            ngmesh = cube_geo.GenerateMesh(maxh=maxh)
            mesh = Mesh(ngmesh)
            n_elements = mesh.ne
        except Exception as e:
            print(f"  maxh={maxh}: Mesh generation failed: {e}")
            continue

        try:
            mag_obj = netgen_mesh_to_radia(mesh,
                                           material={'magnetization': [0, 0, 0]},
                                           units='m',
                                           material_filter='magnetic')
        except Exception as e:
            print(f"  maxh={maxh}: Radia import failed: {e}")
            continue

        mat = rad.MatLin(CHI)
        rad.MatApl(mag_obj, mat)

        bg = rad.ObjBckgCF(uniform_B)
        system = rad.ObjCnt([mag_obj, bg])

        t0 = time.time()
        try:
            result = rad.Solve(system, 0.00001, 2000)
            solve_time = time.time() - t0
        except Exception as e:
            print(f"  maxh={maxh}: Solve failed: {e}")
            continue

        # Compute external field errors
        errors = []
        for pt, desc in OBS_POINTS:
            H = np.array(rad.Fld(system, 'h', pt))
            H_r = H_ref[tuple(pt)]
            err = np.linalg.norm(H - H_r) / np.linalg.norm(H_r) * 100
            errors.append(err)

        avg_error = np.mean(errors)
        max_error = max(errors)

        results.append({
            'maxh': maxh,
            'n_elements': n_elements,
            'avg_error': avg_error,
            'max_error': max_error,
            'solve_time': solve_time,
        })

        print(f"  maxh={maxh:.3f} ({n_elements:4d} elem): "
              f"Avg error = {avg_error:.2f}%, Max = {max_error:.2f}%, "
              f"time = {solve_time:.3f}s")

    return results


def print_results_table(hex_results, tetra_results, tetra_ngsolve_results):
    """Print formatted results table."""
    print("\n" + "=" * 80)
    print("RESULTS SUMMARY")
    print("=" * 80)
    print("External field error vs 8x8x8 hexahedral reference")
    print()

    print("HEXAHEDRAL MESH (ObjRecMag + ObjDivMag)")
    print("-" * 60)
    print(f"{'Subdivisions':<15} {'Elements':<10} {'Avg Err (%)':<12} {'Max Err (%)':<12}")
    print("-" * 60)
    for r in hex_results:
        print(f"{r['subdivs']}x{r['subdivs']}x{r['subdivs']:<10} {r['n_elements']:<10} "
              f"{r['avg_error']:<12.2f} {r['max_error']:<12.2f}")
    print()

    if tetra_results:
        print("TETRAHEDRAL MESH (Manual 5-tetra decomposition)")
        print("-" * 60)
        print(f"{'Subdivisions':<15} {'Elements':<10} {'Avg Err (%)':<12} {'Max Err (%)':<12}")
        print("-" * 60)
        for r in tetra_results:
            print(f"{r['n_subdivs']}^3 x 5{'':<8} {r['n_elements']:<10} "
                  f"{r['avg_error']:<12.2f} {r['max_error']:<12.2f}")
        print()

    if tetra_ngsolve_results:
        print("TETRAHEDRAL MESH (NGSolve/Netgen)")
        print("-" * 60)
        print(f"{'maxh (m)':<15} {'Elements':<10} {'Avg Err (%)':<12} {'Max Err (%)':<12}")
        print("-" * 60)
        for r in tetra_ngsolve_results:
            print(f"{r['maxh']:<15.3f} {r['n_elements']:<10} "
                  f"{r['avg_error']:<12.2f} {r['max_error']:<12.2f}")
        print()

    print("=" * 80)


def main():
    # Build reference
    print("-" * 80)
    print("Building reference: 8x8x8 hexahedral mesh...")
    print("-" * 80)
    H_ref, n_ref = get_reference_fields(8)
    print(f"  Reference mesh: {n_ref} elements")
    print()

    print("-" * 80)
    print("1. HEXAHEDRAL MESH BENCHMARK")
    print("-" * 80)
    hex_subdivs = [1, 2, 3, 4, 5, 6]
    hex_results = benchmark_hexahedra(hex_subdivs, H_ref)

    print("\n" + "-" * 80)
    print("2. TETRAHEDRAL MESH BENCHMARK (Manual 5-tetra decomposition)")
    print("-" * 80)
    tetra_subdivs = [1, 2, 3]
    tetra_results = benchmark_tetrahedra_manual(tetra_subdivs, H_ref)

    print("\n" + "-" * 80)
    print("3. TETRAHEDRAL MESH BENCHMARK (NGSolve/Netgen)")
    print("-" * 80)
    maxh_list = [0.1, 0.08, 0.06, 0.05, 0.04]
    tetra_ngsolve_results = benchmark_tetrahedra_ngsolve(maxh_list, H_ref)

    # Print summary table
    print_results_table(hex_results, tetra_results, tetra_ngsolve_results)


if __name__ == '__main__':
    main()
