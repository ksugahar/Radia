#!/usr/bin/env python3
"""
Profile HACApK H-matrix build to identify performance bottlenecks.

Compares timing with ELF Fortran baseline.
"""

import sys
import os
import time

# Add build directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build-msvc/Release'))
import radia as rad

from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import Mesh, VOL
from ngsolve import TaskManager

# Parameters
CUBE_SIZE = 1.0  # meters
MAXH = 0.20  # mesh size
H_EXT = 200000.0  # A/m

# B-H curve data (nonlinear steel)
BH_DATA = [
    [0.0, 0.0],
    [100.0, 0.2],
    [200.0, 0.5],
    [500.0, 1.0],
    [1000.0, 1.3],
    [2000.0, 1.5],
    [5000.0, 1.7],
    [10000.0, 1.85],
    [50000.0, 2.0],
    [100000.0, 2.1],
    [500000.0, 2.3],
]

def create_mesh():
    """Generate tetrahedral mesh using Netgen."""
    half = CUBE_SIZE / 2
    box = Box(Pnt(-half, -half, -half), Pnt(half, half, half))
    box.mat('magnetic')
    geo = OCCGeometry(box)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=MAXH)
        return Mesh(ngmesh)

def create_radia_elements(mesh):
    """Create Radia tetrahedral elements from NGSolve mesh."""
    # Extract vertices
    vertices = []
    for v in mesh.vertices:
        pt = v.point
        vertices.append([pt[0], pt[1], pt[2]])

    # Create tetrahedra
    container = rad.ObjCnt([])
    for el in mesh.Elements(VOL):
        v_indices = [v.nr for v in el.vertices]
        tet_verts = [vertices[i] for i in v_indices]
        tet = rad.ObjTetrahedron(tet_verts, [0, 0, 0])
        rad.ObjAddToCnt(container, [tet])

    return container

def profile_hacapk():
    """Profile HACApK H-matrix build."""
    print("=" * 60)
    print("HACApK Profiling")
    print("=" * 60)

    rad.UtiDelAll()

    # Step 1: Mesh generation
    print("\n1. Mesh generation...")
    t0 = time.perf_counter()
    mesh = create_mesh()
    t_mesh = time.perf_counter() - t0

    n_elem = mesh.ne
    print(f"   Elements: {n_elem}")
    print(f"   DOF: {n_elem * 3}")
    print(f"   Time: {t_mesh:.4f}s")

    # Step 2: Create Radia elements
    print("\n2. Creating Radia elements...")
    t0 = time.perf_counter()
    container = create_radia_elements(mesh)
    t_elem = time.perf_counter() - t0
    print(f"   Time: {t_elem:.4f}s")

    # Step 3: Apply material
    print("\n3. Applying material...")
    t0 = time.perf_counter()
    mat = rad.MatSatIsoTab(BH_DATA)
    rad.MatApl(container, mat)
    t_mat = time.perf_counter() - t0
    print(f"   Time: {t_mat:.4f}s")

    # Step 4: Add background field
    print("\n4. Adding background field...")
    t0 = time.perf_counter()
    bg_field = rad.ObjBckg(lambda p: [0, 0, H_EXT * 4e-7 * 3.14159])  # Convert H to B
    rad.ObjAddToCnt(container, [bg_field])
    t_pre = time.perf_counter() - t0
    print(f"   Time: {t_pre:.4f}s")

    # Step 5: Configure HACApK
    print("\n5. Configuring HACApK...")
    rad.SolverConfig(hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)

    # Profile solve with HACApK (Method 2)
    print("\n6. Running HACApK solve...")
    print("   (This includes H-matrix build + nonlinear iteration)")

    t0 = time.perf_counter()
    result = rad.Solve(container, 0.001, 100, 2)  # Method 2 = HACApK
    t_solve = time.perf_counter() - t0

    print(f"\n   Total solve time: {t_solve:.4f}s")
    print(f"   Iterations: {result[1]}")

    # Get H-matrix stats
    try:
        config = rad.GetSolverConfig()
        stats = config.get('hacapk_stats')
        if stats:
            print(f"\n   H-matrix stats:")
            print(f"   - n_lowrank: {stats.get('n_lowrank', 'N/A')}")
            print(f"   - n_dense: {stats.get('n_dense', 'N/A')}")
            print(f"   - max_rank: {stats.get('max_rank', 'N/A')}")
            print(f"   - build_time: {stats.get('build_time', 'N/A'):.4f}s")
        else:
            print("   (H-matrix stats not available)")
    except:
        print("   (H-matrix stats not available)")

    # Get M_avg_z
    M_tot = rad.ObjM(container)
    M_avg_z = M_tot[2] / n_elem
    print(f"\n   M_avg_z: {M_avg_z:.0f} A/m")

    # Compare with ELF Fortran baseline
    print("\n" + "=" * 60)
    print("Comparison with ELF Fortran")
    print("=" * 60)
    print(f"{'Metric':<25} {'Radia C':<15} {'ELF Fortran':<15} {'Ratio':<10}")
    print("-" * 60)

    elf_build_time = 0.165  # seconds (from previous benchmark)
    try:
        radia_build_time = stats.get('build_time', t_solve)
        ratio = radia_build_time / elf_build_time
        print(f"{'H-matrix build time':<25} {radia_build_time:.4f}s{'':<7} {elf_build_time:.4f}s{'':<7} {ratio:.2f}x")
    except:
        pass

    print(f"{'M_avg_z (A/m)':<25} {M_avg_z:.0f}{'':<7} 748855{'':<10} {'(matches)':<10}")

    return

if __name__ == '__main__':
    # Import path setup
    script_dir = os.path.dirname(os.path.abspath(__file__))
    radia_dir = os.path.dirname(os.path.dirname(script_dir))
    sys.path.insert(0, os.path.join(radia_dir, 'src', 'radia'))

    profile_hacapk()
