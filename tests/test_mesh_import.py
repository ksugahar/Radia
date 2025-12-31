#!/usr/bin/env python
"""
test_mesh_import.py - Test mesh import functionality for tetrahedral and hexahedral elements

This test verifies:
1. Tetrahedral elements from Netgen mesh (ObjTetrahedron)
2. Hexahedral elements (ObjHexahedron)

The tests verify that mesh import methods produce correct results.

Author: Radia Development Team
Created: 2025-12-05
Updated: 2025-12-30 - Use ObjTetrahedron and ObjHexahedron APIs
"""

import sys
import os
import unittest
import numpy as np

# Add Radia build path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/radia'))

import radia as rad

# Import mesh utilities
from netgen_mesh_import import extract_elements

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
            obj = rad.ObjTetrahedron(tet_verts, magnetization)
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
                obj = rad.ObjTetrahedron(tet_verts, magnetization)
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


class TestHexMeshImport(unittest.TestCase):
    """Test hexahedral mesh import via ObjHexahedron."""

    def setUp(self):
        rad.FldUnits('m')
        rad.UtiDelAll()

    def test_single_hexahedron_field(self):
        """Test field from single ObjHexahedron."""

        magnetization = [0.0, 0.0, 1.0e6]  # 1 MA/m in z

        # ObjHexahedron (MSC)
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
        hex_polyhdr = rad.ObjHexahedron(vertices, magnetization)

        # Test field at points
        test_points = [
            [0.5, 0.5, 2.0],   # +z
            [0.5, 0.5, -1.0],  # -z
            [2.0, 0.5, 0.5],   # +x
            [0.5, 2.0, 0.5],   # +y
        ]

        print("\nObjHexahedron field test:")
        print(f"{'Point':^25} | {'Bz':^15}")
        print("-" * 45)

        for pt in test_points:
            B = rad.Fld(hex_polyhdr, 'b', pt)
            print(f"({pt[0]:.1f}, {pt[1]:.1f}, {pt[2]:.1f})".ljust(25) +
                  f" | {B[2]:+.6e}")

            # Field should be non-zero
            self.assertNotEqual(B[2], 0.0, f"Bz should be non-zero at {pt}")

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
            obj = rad.ObjHexahedron(hex_verts, magnetization)
            polyhedra.append(obj)

        container = rad.ObjCnt(polyhedra)

        # Get H at center
        H = rad.Fld(container, 'h', center)

        N_zz = -H[2] / magnetization[2]

        print(f"\nHexahedral mesh (ObjHexahedron) demagnetization test:")
        print(f"  N divisions: {n_div} (total {len(hexahedra)} hexahedra)")
        print(f"  H at center: {H}")
        print(f"  N_zz computed: {N_zz:.6f}")
        print(f"  N_zz theoretical: {1/3:.6f}")
        print(f"  Error: {abs(N_zz - 1/3) * 100:.2f}%")

        self.assertLess(abs(N_zz - 1/3), 0.05)


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
            obj = rad.ObjTetrahedron(tet_verts, [0, 0, 0])  # Zero initial M
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

    def test_hex_mesh_solver(self):
        """Test ObjHexahedron mesh with solver."""

        center = [0.0, 0.0, 0.0]
        size = [1.0, 1.0, 1.0]
        n_div = 3
        mu_r = 1000.0
        H_ext = 1000.0  # A/m
        B_ext = MU_0 * H_ext

        # ObjHexahedron mesh
        rad.UtiDelAll()
        nodes, hexahedra = create_cube_hex_mesh(center, size, n_div)

        polyhedra = []
        for hex_indices in hexahedra:
            hex_verts = [nodes[i] for i in hex_indices]
            obj = rad.ObjHexahedron(hex_verts, [0, 0, 0])
            polyhedra.append(obj)

        cube_hex = rad.ObjCnt(polyhedra)
        mat = rad.MatLin(mu_r - 1)
        rad.MatApl(cube_hex, mat)
        ext = rad.ObjBckg([0, 0, B_ext])
        grp_hex = rad.ObjCnt([cube_hex, ext])

        result_hex = rad.Solve(grp_hex, 0.001, 100, 1)
        M_hex = rad.ObjM(cube_hex)
        M_avg_hex = np.mean([m[1][2] for m in M_hex])

        print(f"\nObjHexahedron mesh solver test:")
        print(f"  N: {n_div} ({n_div**3} elements)")
        print(f"  mu_r: {mu_r}")
        print(f"  H_ext: {H_ext} A/m")
        print(f"  M_avg_z: {M_avg_hex:.0f} A/m")

        # M should be positive (induced by external field)
        self.assertGreater(M_avg_hex, 0, "Magnetization should be positive")


class TestMethodComparison(unittest.TestCase):
    """Test comparison between tetrahedral and hexahedral mesh accuracy."""

    def setUp(self):
        rad.FldUnits('m')
        rad.UtiDelAll()

    @unittest.skipIf(not NETGEN_AVAILABLE, "Netgen not available")
    def test_tet_vs_hex_accuracy(self):
        """Compare accuracy of tetrahedral (Netgen) vs hexahedral mesh."""

        center = [0.0, 0.0, 0.0]
        size = 1.0
        magnetization = [0.0, 0.0, 1.0]

        # Theoretical demagnetization factor for cube: N = 1/3
        N_theoretical = 1.0 / 3.0

        print(f"\nTet vs Hex mesh accuracy comparison:")
        print(f"  Theoretical N_zz: {N_theoretical:.6f}")
        print()

        # Hexahedral mesh (ObjHexahedron)
        hex_results = []
        for n_div in [2, 3, 4]:
            rad.UtiDelAll()
            nodes, hexs = create_cube_hex_mesh(center, [size, size, size], n_div)
            hex_objs = []
            for hex_idx in hexs:
                hex_verts = [nodes[i] for i in hex_idx]
                obj = rad.ObjHexahedron(hex_verts, magnetization)
                hex_objs.append(obj)
            cube_hex = rad.ObjCnt(hex_objs)
            H_hex = rad.Fld(cube_hex, 'h', center)
            N_hex = -H_hex[2] / magnetization[2]
            err_hex = abs(N_hex - N_theoretical) / N_theoretical * 100

            hex_results.append({
                'n_div': n_div,
                'n_hex': len(hexs),
                'N_hex': N_hex,
                'err_hex': err_hex,
            })

        # Netgen tetrahedral mesh
        tet_results = []
        for maxh in [0.5, 0.4, 0.3]:
            rad.UtiDelAll()
            nodes, tets = create_netgen_tet_mesh(center, size, maxh=maxh)
            tet_objs = []
            for tet_idx in tets:
                tet_verts = [nodes[i] for i in tet_idx]
                obj = rad.ObjTetrahedron(tet_verts, magnetization)
                tet_objs.append(obj)
            cube_tet = rad.ObjCnt(tet_objs)
            H_tet = rad.Fld(cube_tet, 'h', center)
            N_tet = -H_tet[2] / magnetization[2]
            err_tet = abs(N_tet - N_theoretical) / N_theoretical * 100

            tet_results.append({
                'maxh': maxh,
                'n_tet': len(tets),
                'N_tet': N_tet,
                'err_tet': err_tet,
            })

        # Print comparison table - Hexahedral
        print("Hexahedral mesh (ObjHexahedron):")
        print(f"{'N':^5} | {'#Hex':^6} | {'N_hex':^10} | {'Err%':^10}")
        print("-" * 40)
        for r in hex_results:
            print(f"{r['n_div']:^5} | {r['n_hex']:^6} | {r['N_hex']:^10.6f} | {r['err_hex']:^10.4f}")

        print()
        print("Tetrahedral mesh (Netgen):")
        print(f"{'maxh':^6} | {'#Tet':^6} | {'N_tet':^10} | {'Err%':^10}")
        print("-" * 40)
        for r in tet_results:
            print(f"{r['maxh']:^6.2f} | {r['n_tet']:^6} | {r['N_tet']:^10.6f} | {r['err_tet']:^10.4f}")

        # Verify errors are reasonable
        self.assertLess(hex_results[0]['err_hex'], 10.0)
        self.assertLess(tet_results[-1]['err_tet'], 10.0)


if __name__ == '__main__':
    print("=" * 70)
    print("Radia Mesh Import Tests")
    print("=" * 70)
    print()
    print("Testing tetrahedral and hexahedral mesh import functionality")
    print("Using ObjTetrahedron and ObjHexahedron APIs")
    print()

    unittest.main(verbosity=2)
