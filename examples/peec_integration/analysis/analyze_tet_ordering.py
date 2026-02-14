#!/usr/bin/env python
"""
Analyze Netgen tetrahedral vertex ordering to determine correct mapping to Gmsh.
"""

import sys
import numpy as np

# Use custom NGSolve with the fix
NGSOLVE_PATH = r'S:\NGSolve\01_GitHub\install_ksugahar\lib\site-packages'
sys.path.insert(0, NGSOLVE_PATH)

from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import Mesh, VOL

# Create simple tetrahedral mesh
box = Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))
geo = OCCGeometry(box)
mesh = Mesh(geo.GenerateMesh(maxh=2.0))  # Very coarse to get few elements

print("=" * 60)
print("NETGEN TETRAHEDRAL VERTEX ORDERING ANALYSIS")
print("=" * 60)

# Analyze first tetrahedral element
for el in mesh.Elements(VOL):
    if len(el.vertices) == 4:  # Tetrahedron
        print("\nFirst tetrahedral element found:")

        # Get vertex coordinates
        coords = []
        for i, v in enumerate(el.vertices):
            vertex = mesh.vertices[v.nr]
            pt = vertex.point
            coords.append([pt[0], pt[1], pt[2]])
            print(f"  Netgen[{i}]: ({pt[0]:7.4f}, {pt[1]:7.4f}, {pt[2]:7.4f})")

        # Compute volume using scalar triple product
        v0 = np.array(coords[0])
        v1 = np.array(coords[1])
        v2 = np.array(coords[2])
        v3 = np.array(coords[3])

        edge1 = v1 - v0
        edge2 = v2 - v0
        edge3 = v3 - v0

        volume = np.dot(edge1, np.cross(edge2, edge3))

        print(f"\n  Scalar triple product: {volume:.6f}")
        if volume > 0:
            print("  -> POSITIVE (correct orientation)")
        else:
            print("  -> NEGATIVE (wrong orientation)")

        # Gmsh standard tetrahedral vertex ordering:
        # https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
        # Gmsh tet: vertices ordered so that (v1-v0) x (v2-v0) points towards v3
        #
        # Standard Gmsh tet ordering (looking from above):
        #        3
        #       /|\
        #      / | \
        #     /  |  \
        #    0---1---2
        #
        # Netgen might use different ordering. Let's try different permutations.

        print("\n" + "=" * 60)
        print("TESTING DIFFERENT PERMUTATIONS:")
        print("=" * 60)

        # Current tetGmsh: {0,1,2,3,4,5,8,6,7,10,9}
        # For 4-node tet, only first 4 matter: {0,1,2,3} = identity

        # Try different permutations
        permutations = {
            "Identity {0,1,2,3}": [0, 1, 2, 3],
            "Swap 1<->2 {0,2,1,3}": [0, 2, 1, 3],
            "Swap 2<->3 {0,1,3,2}": [0, 1, 3, 2],
            "Rotate {1,2,3,0}": [1, 2, 3, 0],
            "Reverse {3,2,1,0}": [3, 2, 1, 0],
        }

        for name, perm in permutations.items():
            p_coords = [coords[perm[i]] for i in range(4)]
            p_v0 = np.array(p_coords[0])
            p_v1 = np.array(p_coords[1])
            p_v2 = np.array(p_coords[2])
            p_v3 = np.array(p_coords[3])

            p_edge1 = p_v1 - p_v0
            p_edge2 = p_v2 - p_v0
            p_edge3 = p_v3 - p_v0

            p_volume = np.dot(p_edge1, np.cross(p_edge2, p_edge3))

            sign = "POSITIVE" if p_volume > 0 else "NEGATIVE"
            print(f"  {name:25s}: {p_volume:10.6f} ({sign})")

        break  # Only analyze first tet

print("\n" + "=" * 60)
print("GMSH REFERENCE:")
print("=" * 60)
print("Gmsh tet vertex ordering (from gmsh.info):")
print("  Vertices 0,1,2,3 ordered so that:")
print("  (v1-v0) x (v2-v0) points towards v3")
print("  => Scalar triple product > 0")
