#!/usr/bin/env python
"""
Linear Material Benchmark: Tetrahedral (Radia) vs NGSolve FEM

Compares Radia tetrahedral element results with NGSolve FEM solution
for a linear magnetic material in uniform external field.

Evaluation points are OUTSIDE the cube (air region).

Problem Setup:
- Cube: 0.1m x 0.1m x 0.1m, centered at origin
- Material: Linear isotropic, mu_r = 100
- External field: B0 = [0, 0, 1e-6] T (uniform along z-axis)
- Evaluation: Points OUTSIDE the cube

Author: Radia Development Team
Date: 2025-12-03
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/python'))

import numpy as np
import time

# Physical constants
MU_0 = 4 * np.pi * 1e-7  # T/(A/m)

# Problem parameters
CUBE_SIZE = 0.1      # 0.1 m cube (100 mm)
CUBE_HALF = 0.05     # half size
MU_R = 100           # Relative permeability
CHI = MU_R - 1       # Susceptibility = 99
B0_mag = 1e-6        # 1 uT external field
B0 = np.array([0.0, 0.0, B0_mag])
H0 = B0 / MU_0       # External H field


def create_evaluation_points():
    """Create evaluation points outside the cube."""
    points = []
    descriptions = []

    # Points on z-axis (above cube)
    for z in [0.06, 0.08, 0.10, 0.15, 0.20]:
        points.append([0.0, 0.0, z])
        descriptions.append(f'z={z:.2f}m (z-axis)')

    # Points on x-axis (beside cube)
    for x in [0.06, 0.08, 0.10, 0.15, 0.20]:
        points.append([x, 0.0, 0.0])
        descriptions.append(f'x={x:.2f}m (x-axis)')

    # Points on diagonal
    for r in [0.10, 0.15, 0.20]:
        d = r / np.sqrt(3)
        points.append([d, d, d])
        descriptions.append(f'r={r:.2f}m (diagonal)')

    return points, descriptions


def run_radia_hexahedral(n_div):
    """Run Radia with hexahedral mesh."""
    import radia as rad
    rad.FldUnits('m')
    rad.UtiDelAll()

    # Create cube
    cube = rad.ObjRecMag([0, 0, 0], [CUBE_SIZE, CUBE_SIZE, CUBE_SIZE], [0, 0, 0])
    rad.ObjDivMag(cube, [n_div, n_div, n_div])

    # Apply material
    mat = rad.MatLin(CHI)
    rad.MatApl(cube, mat)

    # Background field
    background = rad.ObjBckg([B0[0], B0[1], B0[2]])
    container = rad.ObjCnt([cube, background])

    # Solve
    start = time.time()
    res = rad.Solve(container, 0.0001, 10000, 9)  # Method 9 for stability
    solve_time = time.time() - start

    return container, cube, res, solve_time


def run_radia_tetrahedral(max_h):
    """Run Radia with tetrahedral mesh from Netgen."""
    import radia as rad
    rad.FldUnits('m')
    rad.UtiDelAll()
    rad.SolverTetraMethod(0)  # Method 0: polygon-based

    try:
        from netgen.occ import Box, OCCGeometry
        from ngsolve import Mesh
        from netgen_mesh_import import netgen_mesh_to_radia
    except ImportError as e:
        print(f"[SKIP] Tetrahedral: {e}")
        return None, None, None, 0

    # Create mesh
    cube_geom = Box((-CUBE_HALF, -CUBE_HALF, -CUBE_HALF),
                    (CUBE_HALF, CUBE_HALF, CUBE_HALF))
    geo = OCCGeometry(cube_geom)
    ngmesh = Mesh(geo.GenerateMesh(maxh=max_h))
    n_elem = ngmesh.ne

    # Import to Radia
    cube = netgen_mesh_to_radia(ngmesh,
                                 material={'magnetization': [0, 0, 0]},
                                 units='m',
                                 verbose=False)

    # Apply material
    mat = rad.MatLin(CHI)
    rad.MatApl(cube, mat)

    # Background field
    background = rad.ObjBckg([B0[0], B0[1], B0[2]])
    container = rad.ObjCnt([cube, background])

    # Solve
    start = time.time()
    res = rad.Solve(container, 0.001, 10000, 9)  # Method 9
    solve_time = time.time() - start

    return container, cube, res, solve_time, n_elem


def run_ngsolve_fem(max_h):
    """Run NGSolve FEM for comparison."""
    try:
        import ngsolve
        from ngsolve import Mesh, HCurl, HDiv, BilinearForm, GridFunction
        from ngsolve import curl, dx, CoefficientFunction, x, y, z, BND
        from netgen.occ import Box, OCCGeometry, Glue
    except ImportError as e:
        print(f"[SKIP] NGSolve: {e}")
        return None

    # Create geometry: cube inside larger air domain
    cube = Box((-CUBE_HALF, -CUBE_HALF, -CUBE_HALF),
               (CUBE_HALF, CUBE_HALF, CUBE_HALF))
    air_size = 0.5  # Air domain extends to +/- 0.5m
    air = Box((-air_size, -air_size, -air_size),
              (air_size, air_size, air_size))

    cube.mat("iron")
    air.mat("air")

    # Combine geometries
    geo = OCCGeometry(Glue([cube, air - cube]))
    mesh = Mesh(geo.GenerateMesh(maxh=max_h))

    # Material coefficients
    mu_iron = MU_0 * MU_R
    mu_air = MU_0

    mu_cf = mesh.MaterialCF({"iron": mu_iron, "air": mu_air}, default=mu_air)
    nu_cf = 1.0 / mu_cf

    # Magnetostatic formulation: curl(nu * curl A) = 0
    # with boundary condition for applied field
    fes = HCurl(mesh, order=2, dirichlet=".*")

    # External field: A = 0.5 * B0 x r = 0.5 * [B0_z*y - B0_y*z, B0_x*z - B0_z*x, ...]
    # For B0 = [0, 0, B0_z]: A = [-0.5*B0_z*y, 0.5*B0_z*x, 0]
    A_ext = CoefficientFunction((-0.5*B0_mag*y, 0.5*B0_mag*x, 0))

    u, v = fes.TnT()

    a = BilinearForm(fes)
    a += nu_cf * curl(u) * curl(v) * dx
    a.Assemble()

    # Apply boundary conditions
    gfu = GridFunction(fes)
    gfu.Set(A_ext, BND)

    # Solve
    r = gfu.vec.CreateVector()
    r.data = -a.mat * gfu.vec

    # Use direct solver
    inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    gfu.vec.data += inv * r

    # Compute B = curl(A)
    B_gf = GridFunction(HDiv(mesh, order=1))
    B_gf.Set(curl(gfu))

    return mesh, gfu, B_gf


def evaluate_field_radia(container, points):
    """Evaluate B field at points using Radia."""
    import radia as rad

    B_values = []
    for pt in points:
        B = rad.Fld(container, 'b', pt)
        B_values.append(np.array(B))

    return B_values


def evaluate_field_ngsolve(mesh, B_gf, points):
    """Evaluate B field at points using NGSolve."""
    from ngsolve import Mesh

    B_values = []
    for pt in points:
        mp = mesh(pt[0], pt[1], pt[2])
        B = B_gf(mp)
        B_values.append(np.array([B[0], B[1], B[2]]))

    return B_values


def main():
    import radia as rad

    print("=" * 80)
    print("Linear Material Benchmark: Tetrahedral (Radia) vs NGSolve FEM")
    print("=" * 80)

    print(f"\nProblem Setup:")
    print(f"  Cube size:   {CUBE_SIZE*1000:.0f} mm (centered at origin)")
    print(f"  mu_r:        {MU_R}")
    print(f"  B0:          {B0_mag*1e6:.1f} uT (along z-axis)")
    print(f"  H0:          {H0[2]:.4f} A/m")

    # Evaluation points
    points, descriptions = create_evaluation_points()

    results = {}

    # 1. Hexahedral reference (8x8x8)
    print("\n" + "-" * 80)
    print("Reference: Hexahedral 8x8x8 (Method 9)")
    print("-" * 80)

    container_hex, cube_hex, res_hex, time_hex = run_radia_hexahedral(8)
    print(f"  Elements: 512")
    print(f"  Iterations: {res_hex[3]:.0f}")
    print(f"  Time: {time_hex:.3f} s")

    B_hex = evaluate_field_radia(container_hex, points)
    results['hex_8x8x8'] = {'B': B_hex, 'n_elem': 512}

    # 2. Tetrahedral meshes
    print("\n" + "-" * 80)
    print("Tetrahedral Meshes (Netgen import)")
    print("-" * 80)

    for max_h in [0.06, 0.04, 0.03]:
        rad.UtiDelAll()
        ret = run_radia_tetrahedral(max_h)
        if ret[0] is None:
            continue
        container_tet, cube_tet, res_tet, time_tet, n_elem = ret

        print(f"\n  maxh={max_h}:")
        print(f"    Elements: {n_elem}")
        print(f"    Iterations: {res_tet[3]:.0f}")
        print(f"    Time: {time_tet:.3f} s")

        B_tet = evaluate_field_radia(container_tet, points)
        results[f'tet_maxh{max_h}'] = {'B': B_tet, 'n_elem': n_elem}

    # 3. NGSolve FEM reference
    print("\n" + "-" * 80)
    print("NGSolve FEM Reference")
    print("-" * 80)

    try:
        ngsolve_result = run_ngsolve_fem(max_h=0.05)
        if ngsolve_result is not None:
            ng_mesh, ng_gfu, ng_B = ngsolve_result
            print(f"  Mesh elements: {ng_mesh.ne}")

            B_ng = evaluate_field_ngsolve(ng_mesh, ng_B, points)
            results['ngsolve'] = {'B': B_ng, 'n_elem': ng_mesh.ne}
    except Exception as e:
        print(f"  [SKIP] NGSolve failed: {e}")

    # 4. Comparison table
    print("\n" + "=" * 80)
    print("RESULTS: B-field at External Points")
    print("=" * 80)

    # Use hex_8x8x8 as reference
    B_ref = results['hex_8x8x8']['B']

    print(f"\n{'Point':<25} {'Hex 8^3':>12} ", end='')
    for key in results:
        if key != 'hex_8x8x8':
            print(f"{key:>12} ", end='')
    print()
    print("-" * 100)

    errors = {key: [] for key in results if key != 'hex_8x8x8'}

    for i, (pt, desc) in enumerate(zip(points, descriptions)):
        Bz_ref = B_ref[i][2]
        pt_str = f"({pt[0]:.2f},{pt[1]:.2f},{pt[2]:.2f})"
        print(f"{pt_str:<25} {Bz_ref*1e9:>12.4f} ", end='')  # nT

        for key in results:
            if key == 'hex_8x8x8':
                continue
            Bz = results[key]['B'][i][2]
            print(f"{Bz*1e9:>12.4f} ", end='')  # nT

            # Calculate error
            if abs(Bz_ref) > 1e-15:
                err = abs(Bz - Bz_ref) / abs(Bz_ref) * 100
            else:
                err = 0
            errors[key].append(err)
        print()

    print("-" * 100)
    print("(Bz in nT)")

    # Error summary
    print("\n" + "=" * 80)
    print("ERROR SUMMARY (vs Hexahedral 8x8x8 reference)")
    print("=" * 80)

    print(f"\n{'Mesh Type':<25} {'Elements':>10} {'Avg Error':>12} {'Max Error':>12}")
    print("-" * 65)

    for key in errors:
        n_elem = results[key]['n_elem']
        avg_err = np.mean(errors[key])
        max_err = np.max(errors[key])
        print(f"{key:<25} {n_elem:>10} {avg_err:>11.2f}% {max_err:>11.2f}%")

    # Return summary for README
    return results, errors


if __name__ == "__main__":
    results, errors = main()
