#!/usr/bin/env python
"""
test_mesh_import.py - Test mesh import functionality for tetrahedral and hexahedral elements

This test compares:
1. Tetrahedral elements from Netgen mesh (ObjPolyhdr with TETRA_FACES)
2. Hexahedral elements via ObjPolyhdr (general polyhedron MSC)
3. Hexahedral elements via ObjRecMag (analytical formula)

The tests verify that different mesh import methods produce consistent results.

Author: Radia Development Team
Created: 2025-12-05
"""

import sys
import os
import unittest
import numpy as np

# Add Radia build path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/radia'))

import radia as rad

# Import face topologies and mesh utilities
from netgen_mesh_import import TETRA_FACES, HEX_FACES, extract_elements

# Try to import Netgen and NGSolve
try:
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen.meshing import MeshingParameters
    from ngsolve import Mesh as NGSolveMesh
    NETGEN_AVAILABLE = True
except ImportError:
    NETGEN_AVAILABLE = False
    print("Warning: Netgen/NGSolve not available, some tests will be skipped")

MU_0 = 4 * np.pi * 1e-7


def create_cube_hex_mesh(center, size, n_div=2):
    """
    Create a cube mesh using hexahedral elements.

    Parameters
    ----------
    center : list
        [x, y, z] center of cube
    size : list
        [dx, dy, dz] dimensions
    n_div : int
        Subdivisions per axis (total elements = n_div^3)

    Returns
    -------
    tuple
        (nodes, hexahedra) - nodes array and element connectivity
    """
    cx, cy, cz = center
    dx, dy, dz = size

    # Create grid of nodes
    nodes = []
    node_map = {}

    for k in range(n_div + 1):
        for j in range(n_div + 1):
            for i in range(n_div + 1):
                x = cx - dx/2 + i * dx / n_div
                y = cy - dy/2 + j * dy / n_div
                z = cz - dz/2 + k * dz / n_div
                node_map[(i, j, k)] = len(nodes)
                nodes.append([x, y, z])

    # Create hexahedra
    hexahedra = []

    for k in range(n_div):
        for j in range(n_div):
            for i in range(n_div):
                # 8 vertices of hexahedron (standard ordering)
                v = [
                    node_map[(i, j, k)],         # 0
                    node_map[(i+1, j, k)],       # 1
                    node_map[(i+1, j+1, k)],     # 2
                    node_map[(i, j+1, k)],       # 3
                    node_map[(i, j, k+1)],       # 4
                    node_map[(i+1, j, k+1)],     # 5
                    node_map[(i+1, j+1, k+1)],   # 6
                    node_map[(i, j+1, k+1)],     # 7
                ]
                hexahedra.append(v)

    return nodes, hexahedra


def create_netgen_tet_mesh(center, size, maxh=0.3):
    """
    Create a tetrahedral mesh using Netgen.

    Parameters
    ----------
    center : list
        [x, y, z] center of cube
    size : list or float
        [dx, dy, dz] dimensions or single value for cube
    maxh : float
        Maximum element size

    Returns
    -------
    tuple
        (nodes, tetrahedra) - nodes array and element connectivity
    """
    if not NETGEN_AVAILABLE:
        raise ImportError("Netgen is not available")

    cx, cy, cz = center
    if isinstance(size, (int, float)):
        dx = dy = dz = float(size)
    else:
        dx, dy, dz = size

    # Create box geometry
    p1 = Pnt(cx - dx/2, cy - dy/2, cz - dz/2)
    p2 = Pnt(cx + dx/2, cy + dy/2, cz + dz/2)
    box = Box(p1, p2)

    # Generate mesh - need to wrap in OCCGeometry first
    geo = OCCGeometry(box)
    mp = MeshingParameters(maxh=maxh)
    netgen_mesh = geo.GenerateMesh(mp)

    # Wrap in NGSolve Mesh and use extract_elements for correct indexing
    ngsolve_mesh = NGSolveMesh(netgen_mesh)
    elements, _ = extract_elements(ngsolve_mesh)

    # Build nodes list from element vertices (deduplicated)
    nodes_dict = {}
    tetrahedra = []
    for el_data in elements:
        tet_indices = []
        for v in el_data['vertices']:
            v_tuple = tuple(v)
            if v_tuple not in nodes_dict:
                nodes_dict[v_tuple] = len(nodes_dict)
            tet_indices.append(nodes_dict[v_tuple])
        tetrahedra.append(tet_indices)

    # Convert nodes dict to list
    nodes = [None] * len(nodes_dict)
    for v, idx in nodes_dict.items():
        nodes[idx] = list(v)

    return nodes, tetrahedra


class TestNetgenMeshImport(unittest.TestCase):
    """Test mesh import using Netgen tetrahedral mesh."""

    def setUp(self):
        """Set up test fixtures."""
        rad.FldUnits('m')
        rad.UtiDelAll()

    @unittest.skipIf(not NETGEN_AVAILABLE, "Netgen not available")
    def test_netgen_tet_field(self):
        """Test field from Netgen tetrahedral mesh."""
        rad.UtiDelAll()

        center = [0.0, 0.0, 0.0]
        size = 1.0
        magnetization = [0.0, 0.0, 1.0e6]  # 1 MA/m in z

        # Create mesh
        nodes, tets = create_netgen_tet_mesh(center, size, maxh=0.5)
        print(f"\nNetgen mesh: {len(nodes)} nodes, {len(tets)} tetrahedra")

        # Create Radia objects
        polyhedra = []
        for tet_indices in tets:
            tet_verts = [nodes[i] for i in tet_indices]
            obj = rad.ObjPolyhdr(tet_verts, TETRA_FACES, magnetization)
            polyhedra.append(obj)

        container = rad.ObjCnt(polyhedra)

        # Test field at point outside
        test_point = [0.0, 0.0, 1.0]
        B = rad.Fld(container, 'b', test_point)

        print(f"B at {test_point}: {B}")
        print(f"Bz: {B[2]:.6e} T")

        # Field should be non-zero and primarily in z direction
        self.assertGreater(abs(B[2]), 0.0, "Bz should be non-zero")

    @unittest.skipIf(not NETGEN_AVAILABLE, "Netgen not available")
    def test_netgen_tet_demagnetization(self):
        """Test demagnetization factor for cube with Netgen tetrahedral mesh."""
        rad.UtiDelAll()

        center = [0.0, 0.0, 0.0]
        size = 1.0
        magnetization = [0.0, 0.0, 1.0]  # Unit magnetization

        # Create mesh with different densities
        for maxh in [0.4, 0.3, 0.25]:
            rad.UtiDelAll()

            nodes, tets = create_netgen_tet_mesh(center, size, maxh=maxh)

            polyhedra = []
            for tet_indices in tets:
                tet_verts = [nodes[i] for i in tet_indices]
                obj = rad.ObjPolyhdr(tet_verts, TETRA_FACES, magnetization)
                polyhedra.append(obj)

            container = rad.ObjCnt(polyhedra)

            # Get H at center
            H = rad.Fld(container, 'h', center)
            N_zz = -H[2] / magnetization[2]

            print(f"\nNetgen maxh={maxh}: {len(tets)} tetrahedra")
            print(f"  H at center: {H}")
            print(f"  N_zz computed: {N_zz:.6f}")
            print(f"  N_zz theoretical: {1/3:.6f}")
            print(f"  Error: {abs(N_zz - 1/3) * 100:.2f}%")

        # Final test with finest mesh should be accurate
        self.assertLess(abs(N_zz - 1/3), 0.05,
                        f"Demagnetization factor error exceeds 5%")

    @unittest.skipIf(not NETGEN_AVAILABLE, "Netgen not available")
    def test_netgen_vs_recmag_comparison(self):
        """Compare Netgen tet mesh vs ObjRecMag for same geometry."""
        center = [0.0, 0.0, 0.0]
        size = 1.0
        magnetization = [0.0, 0.0, 1.0e6]

        # Reference: Fine ObjRecMag mesh
        rad.UtiDelAll()
        cube_ref = rad.ObjRecMag(center, [size, size, size], magnetization)
        rad.ObjDivMag(cube_ref, [10, 10, 10])  # 1000 elements
        H_ref = rad.Fld(cube_ref, 'h', center)

        # Netgen mesh
        rad.UtiDelAll()
        nodes, tets = create_netgen_tet_mesh(center, size, maxh=0.25)

        polyhedra = []
        for tet_indices in tets:
            tet_verts = [nodes[i] for i in tet_indices]
            obj = rad.ObjPolyhdr(tet_verts, TETRA_FACES, magnetization)
            polyhedra.append(obj)

        cube_netgen = rad.ObjCnt(polyhedra)
        H_netgen = rad.Fld(cube_netgen, 'h', center)

        print(f"\nNetgen vs ObjRecMag comparison:")
        print(f"  ObjRecMag (10x10x10, 1000 elements): Hz = {H_ref[2]:.6e}")
        print(f"  Netgen ({len(tets)} tets): Hz = {H_netgen[2]:.6e}")

        if abs(H_ref[2]) > 1e-10:
            diff_pct = abs(H_netgen[2] - H_ref[2]) / abs(H_ref[2]) * 100
            print(f"  Difference: {diff_pct:.2f}%")
            self.assertLess(diff_pct, 10.0, "Netgen should match analytical within 10%")


class TestHexMeshImport(unittest.TestCase):
    """Test hexahedral mesh import."""

    def setUp(self):
        rad.FldUnits('m')
        rad.UtiDelAll()

    def test_single_hexahedron_polyhdr_vs_recmag(self):
        """Compare ObjPolyhdr (MSC) vs ObjRecMag (analytical) for single hexahedron."""

        center = [0.5, 0.5, 0.5]
        size = [1.0, 1.0, 1.0]
        magnetization = [0.0, 0.0, 1.0e6]  # 1 MA/m in z

        # Method 1: ObjPolyhdr with HEX_FACES (MSC)
        rad.UtiDelAll()
        vertices = [
            [0.0, 0.0, 0.0],  # 0
            [1.0, 0.0, 0.0],  # 1
            [1.0, 1.0, 0.0],  # 2
            [0.0, 1.0, 0.0],  # 3
            [0.0, 0.0, 1.0],  # 4
            [1.0, 0.0, 1.0],  # 5
            [1.0, 1.0, 1.0],  # 6
            [0.0, 1.0, 1.0],  # 7
        ]
        hex_polyhdr = rad.ObjPolyhdr(vertices, HEX_FACES, magnetization)

        # Method 2: ObjRecMag (analytical)
        rad.UtiDelAll()
        hex_recmag = rad.ObjRecMag(center, size, magnetization)

        # Compare fields at test points
        test_points = [
            [0.5, 0.5, 2.0],   # +z
            [0.5, 0.5, -1.0],  # -z
            [2.0, 0.5, 0.5],   # +x
            [0.5, 2.0, 0.5],   # +y
        ]

        print("\nComparing ObjPolyhdr (MSC) vs ObjRecMag (analytical):")
        print(f"{'Point':^25} | {'Polyhdr Bz':^15} | {'RecMag Bz':^15} | {'Diff %':^10}")
        print("-" * 75)

        max_diff_pct = 0.0

        for pt in test_points:
            # Get fields
            rad.UtiDelAll()
            hex_p = rad.ObjPolyhdr(vertices, HEX_FACES, magnetization)
            B_polyhdr = rad.Fld(hex_p, 'b', pt)

            rad.UtiDelAll()
            hex_r = rad.ObjRecMag(center, size, magnetization)
            B_recmag = rad.Fld(hex_r, 'b', pt)

            # Compare
            if abs(B_recmag[2]) > 1e-10:
                diff_pct = abs(B_polyhdr[2] - B_recmag[2]) / abs(B_recmag[2]) * 100
            else:
                diff_pct = 0.0

            max_diff_pct = max(max_diff_pct, diff_pct)

            print(f"({pt[0]:.1f}, {pt[1]:.1f}, {pt[2]:.1f})".ljust(25) +
                  f" | {B_polyhdr[2]:+.6e} | {B_recmag[2]:+.6e} | {diff_pct:.4f}%")

        # Allow up to 1% difference (numerical precision)
        self.assertLess(max_diff_pct, 1.0,
                        f"Max difference {max_diff_pct:.4f}% exceeds 1% tolerance")

    def test_cube_hex_mesh_demagnetization(self):
        """Test demagnetization factor for cube discretized with hexahedra."""
        rad.UtiDelAll()

        center = [0.0, 0.0, 0.0]
        size = [1.0, 1.0, 1.0]
        n_div = 2  # 2^3 = 8 hexahedra

        nodes, hexahedra = create_cube_hex_mesh(center, size, n_div)

        # Create hexahedra with unit Mz
        magnetization = [0.0, 0.0, 1.0]

        polyhedra = []
        for hex_indices in hexahedra:
            hex_verts = [nodes[i] for i in hex_indices]
            obj = rad.ObjPolyhdr(hex_verts, HEX_FACES, magnetization)
            polyhedra.append(obj)

        container = rad.ObjCnt(polyhedra)

        # Get H at center
        H = rad.Fld(container, 'h', center)

        N_zz = -H[2] / magnetization[2]

        print(f"\nHexahedral mesh (ObjPolyhdr) demagnetization test:")
        print(f"  N divisions: {n_div} (total {len(hexahedra)} hexahedra)")
        print(f"  H at center: {H}")
        print(f"  N_zz computed: {N_zz:.6f}")
        print(f"  N_zz theoretical: {1/3:.6f}")
        print(f"  Error: {abs(N_zz - 1/3) * 100:.2f}%")

        self.assertLess(abs(N_zz - 1/3), 0.05)

    def test_recmag_vs_polyhdr_mesh(self):
        """Compare ObjRecMag+ObjDivMag vs ObjPolyhdr mesh for same geometry."""

        center = [0.0, 0.0, 0.0]
        size = [1.0, 1.0, 1.0]
        n_div = 3
        magnetization = [0.0, 0.0, 1.0e6]

        # Method 1: ObjRecMag + ObjDivMag (analytical formula)
        rad.UtiDelAll()
        cube_analytical = rad.ObjRecMag(center, size, magnetization)
        rad.ObjDivMag(cube_analytical, [n_div, n_div, n_div])

        H_analytical = rad.Fld(cube_analytical, 'h', center)

        # Method 2: ObjPolyhdr mesh (MSC)
        rad.UtiDelAll()
        nodes, hexahedra = create_cube_hex_mesh(center, size, n_div)

        polyhedra = []
        for hex_indices in hexahedra:
            hex_verts = [nodes[i] for i in hex_indices]
            obj = rad.ObjPolyhdr(hex_verts, HEX_FACES, magnetization)
            polyhedra.append(obj)

        cube_mesh = rad.ObjCnt(polyhedra)
        H_mesh = rad.Fld(cube_mesh, 'h', center)

        print(f"\nObjRecMag+ObjDivMag vs ObjPolyhdr mesh comparison:")
        print(f"  N divisions: {n_div} ({n_div**3} elements)")
        print(f"  H_analytical at center: {H_analytical}")
        print(f"  H_mesh at center: {H_mesh}")

        # Compare Hz
        if abs(H_analytical[2]) > 1e-10:
            diff_pct = abs(H_mesh[2] - H_analytical[2]) / abs(H_analytical[2]) * 100
        else:
            diff_pct = 0.0

        print(f"  Difference: {diff_pct:.4f}%")

        # Should be very close (< 1%)
        self.assertLess(diff_pct, 1.0)


class TestMeshImportSolver(unittest.TestCase):
    """Test mesh import with solver for soft magnetic materials."""

    def setUp(self):
        rad.FldUnits('m')
        rad.UtiDelAll()

    @unittest.skipIf(not NETGEN_AVAILABLE, "Netgen not available")
    def test_netgen_tet_mesh_with_material(self):
        """Test Netgen tetrahedral mesh with linear material and solver."""
        rad.UtiDelAll()

        center = [0.0, 0.0, 0.0]
        size = 1.0

        # Create tetrahedral mesh
        nodes, tets = create_netgen_tet_mesh(center, size, maxh=0.4)
        print(f"\nNetgen mesh with material: {len(tets)} tetrahedra")

        polyhedra = []
        for tet_indices in tets:
            tet_verts = [nodes[i] for i in tet_indices]
            obj = rad.ObjPolyhdr(tet_verts, TETRA_FACES, [0, 0, 0])  # Zero initial M
            polyhedra.append(obj)

        container = rad.ObjCnt(polyhedra)

        # Apply linear material (mu_r = 1000)
        mat = rad.MatLin(999.0)  # chi = mu_r - 1
        rad.MatApl(container, mat)

        # External field
        H_ext = 1000.0  # A/m
        B_ext = MU_0 * H_ext
        ext_field = rad.ObjBckg([0, 0, B_ext])
        grp = rad.ObjCnt([container, ext_field])

        # Solve
        result = rad.Solve(grp, 0.001, 100, 1)  # BiCGSTAB

        # Get magnetization
        all_M = rad.ObjM(container)
        M_list = [m[1] for m in all_M]
        M_avg_z = np.mean([m[2] for m in M_list])

        print(f"  mu_r: 1000")
        print(f"  H_ext: {H_ext} A/m")
        print(f"  M_avg_z: {M_avg_z:.0f} A/m")
        print(f"  Iterations: {result[3]:.0f}")

        # M should be positive (induced by external field)
        self.assertGreater(M_avg_z, 0, "Magnetization should be positive")

    def test_hex_mesh_polyhdr_vs_recmag_solver(self):
        """Compare ObjPolyhdr vs ObjRecMag hex mesh with solver."""

        center = [0.0, 0.0, 0.0]
        size = [1.0, 1.0, 1.0]
        n_div = 3
        mu_r = 1000.0
        H_ext = 1000.0  # A/m
        B_ext = MU_0 * H_ext

        # Method 1: ObjPolyhdr mesh
        rad.UtiDelAll()
        nodes, hexahedra = create_cube_hex_mesh(center, size, n_div)

        polyhedra = []
        for hex_indices in hexahedra:
            hex_verts = [nodes[i] for i in hex_indices]
            obj = rad.ObjPolyhdr(hex_verts, HEX_FACES, [0, 0, 0])
            polyhedra.append(obj)

        cube_polyhdr = rad.ObjCnt(polyhedra)
        mat = rad.MatLin(mu_r - 1)
        rad.MatApl(cube_polyhdr, mat)
        ext = rad.ObjBckg([0, 0, B_ext])
        grp_polyhdr = rad.ObjCnt([cube_polyhdr, ext])

        result_polyhdr = rad.Solve(grp_polyhdr, 0.001, 100, 1)
        M_polyhdr = rad.ObjM(cube_polyhdr)
        M_avg_polyhdr = np.mean([m[1][2] for m in M_polyhdr])

        # Method 2: ObjRecMag + ObjDivMag
        rad.UtiDelAll()
        cube_recmag = rad.ObjRecMag(center, size, [0, 0, 0])
        rad.ObjDivMag(cube_recmag, [n_div, n_div, n_div])
        mat2 = rad.MatLin(mu_r - 1)
        rad.MatApl(cube_recmag, mat2)
        ext2 = rad.ObjBckg([0, 0, B_ext])
        grp_recmag = rad.ObjCnt([cube_recmag, ext2])

        result_recmag = rad.Solve(grp_recmag, 0.001, 100, 1)
        M_recmag = rad.ObjM(cube_recmag)
        M_avg_recmag = np.mean([m[1][2] for m in M_recmag])

        print(f"\nObjPolyhdr vs ObjRecMag solver comparison:")
        print(f"  N: {n_div} ({n_div**3} elements)")
        print(f"  mu_r: {mu_r}")
        print(f"  H_ext: {H_ext} A/m")
        print(f"  M_avg_z (Polyhdr): {M_avg_polyhdr:.0f} A/m")
        print(f"  M_avg_z (RecMag): {M_avg_recmag:.0f} A/m")

        diff_pct = abs(M_avg_polyhdr - M_avg_recmag) / abs(M_avg_recmag) * 100
        print(f"  Difference: {diff_pct:.2f}%")

        # Should be close (< 5%)
        self.assertLess(diff_pct, 5.0)


class TestMethodComparison(unittest.TestCase):
    """Test comparison between different field computation methods."""

    def setUp(self):
        rad.FldUnits('m')
        rad.UtiDelAll()

    @unittest.skipIf(not NETGEN_AVAILABLE, "Netgen not available")
    def test_tet_vs_hex_accuracy(self):
        """Compare accuracy of tetrahedral (Netgen) vs hexahedral mesh."""

        center = [0.0, 0.0, 0.0]
        size = 1.0
        magnetization = [0.0, 0.0, 1.0]

        # Reference: Fine ObjRecMag mesh (analytical, high accuracy)
        rad.UtiDelAll()
        cube_ref = rad.ObjRecMag(center, [size, size, size], magnetization)
        rad.ObjDivMag(cube_ref, [10, 10, 10])  # 1000 elements
        H_ref = rad.Fld(cube_ref, 'h', center)
        N_ref = -H_ref[2] / magnetization[2]

        print(f"\nTet vs Hex mesh accuracy comparison:")
        print(f"  Reference (ObjRecMag 10x10x10): N_zz = {N_ref:.6f}")
        print()

        results = []

        # Hexahedral mesh (ObjPolyhdr)
        for n_div in [2, 3, 4]:
            rad.UtiDelAll()
            nodes, hexs = create_cube_hex_mesh(center, [size, size, size], n_div)
            hex_objs = []
            for hex_idx in hexs:
                hex_verts = [nodes[i] for i in hex_idx]
                obj = rad.ObjPolyhdr(hex_verts, HEX_FACES, magnetization)
                hex_objs.append(obj)
            cube_hex = rad.ObjCnt(hex_objs)
            H_hex = rad.Fld(cube_hex, 'h', center)
            N_hex = -H_hex[2] / magnetization[2]
            err_hex = abs(N_hex - N_ref) / N_ref * 100

            # ObjRecMag - analytical
            rad.UtiDelAll()
            cube_ana = rad.ObjRecMag(center, [size, size, size], magnetization)
            rad.ObjDivMag(cube_ana, [n_div, n_div, n_div])
            H_ana = rad.Fld(cube_ana, 'h', center)
            N_ana = -H_ana[2] / magnetization[2]
            err_ana = abs(N_ana - N_ref) / N_ref * 100

            results.append({
                'n_div': n_div,
                'n_hex': len(hexs),
                'N_hex': N_hex,
                'N_ana': N_ana,
                'err_hex': err_hex,
                'err_ana': err_ana,
            })

        # Netgen tetrahedral mesh
        tet_results = []
        for maxh in [0.5, 0.4, 0.3]:
            rad.UtiDelAll()
            nodes, tets = create_netgen_tet_mesh(center, size, maxh=maxh)
            tet_objs = []
            for tet_idx in tets:
                tet_verts = [nodes[i] for i in tet_idx]
                obj = rad.ObjPolyhdr(tet_verts, TETRA_FACES, magnetization)
                tet_objs.append(obj)
            cube_tet = rad.ObjCnt(tet_objs)
            H_tet = rad.Fld(cube_tet, 'h', center)
            N_tet = -H_tet[2] / magnetization[2]
            err_tet = abs(N_tet - N_ref) / N_ref * 100

            tet_results.append({
                'maxh': maxh,
                'n_tet': len(tets),
                'N_tet': N_tet,
                'err_tet': err_tet,
            })

        # Print comparison table - Hexahedral
        print("Hexahedral mesh (ObjPolyhdr vs ObjRecMag):")
        print(f"{'N':^5} | {'#Hex':^6} | {'N_hex':^10} | {'N_ana':^10} | {'Err_hex%':^10} | {'Err_ana%':^10}")
        print("-" * 65)
        for r in results:
            print(f"{r['n_div']:^5} | {r['n_hex']:^6} | {r['N_hex']:^10.6f} | {r['N_ana']:^10.6f} | {r['err_hex']:^10.4f} | {r['err_ana']:^10.4f}")

        print()
        print("Tetrahedral mesh (Netgen):")
        print(f"{'maxh':^6} | {'#Tet':^6} | {'N_tet':^10} | {'Err_tet%':^10}")
        print("-" * 40)
        for r in tet_results:
            print(f"{r['maxh']:^6.2f} | {r['n_tet']:^6} | {r['N_tet']:^10.6f} | {r['err_tet']:^10.4f}")

        # Verify errors are reasonable
        self.assertLess(results[0]['err_hex'], 10.0)
        self.assertLess(tet_results[-1]['err_tet'], 10.0)


if __name__ == '__main__':
    print("=" * 70)
    print("Radia Mesh Import Tests")
    print("=" * 70)
    print()
    print("Testing tetrahedral and hexahedral mesh import functionality")
    print("Comparing ObjPolyhdr (MSC) vs ObjRecMag (analytical)")
    print()

    unittest.main(verbosity=2)
