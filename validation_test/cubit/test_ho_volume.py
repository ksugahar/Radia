"""
Verify order=2 mesh export quality by comparing volumes.

Creates a sphere (r=0.05 m, V_exact = 5.236e-4 m^3), exports as
order=1 and order=2 GMSH .msh, reads both back, computes volume
from tet elements, and reports the error vs analytical.

If order=2 projection is working, TET10 volume should be closer
to the analytical value than TET4.

Usage:
  python validation_test/cubit/test_ho_volume.py
"""

import sys
import os
import math
import numpy as np

# --- Setup paths ---
_test_dir = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.dirname(os.path.dirname(_test_dir))
sys.path.insert(0, os.path.join(_repo_root, 'src', 'radia'))
from install_panels import find_cubit_bin

_cubit_path = find_cubit_bin()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)

_plugin_dir = os.path.join(os.path.dirname(_cubit_path), 'bin', 'plugins') \
    if _cubit_path else None
if _plugin_dir and os.path.isdir(_plugin_dir):
    os.environ['CUBIT_PLUGIN_DIR'] = _plugin_dir

import cubit
cubit.init(['cubit', '-nojournal', '-batch',
            '-commandplugindir', _plugin_dir or ''])

# --- Output directory ---
OUT_DIR = os.path.join(_test_dir, 'ho_volume_test')
os.makedirs(OUT_DIR, exist_ok=True)

RADIUS = 0.05  # meters
V_EXACT = (4.0 / 3.0) * math.pi * RADIUS**3


def build_sphere_mesh(mesh_size=0.02):
    """Create a meshed sphere in Cubit."""
    cubit.cmd("reset")
    cubit.cmd(f"create sphere radius {RADIUS}")
    cubit.cmd("volume 1 scheme tetmesh")
    cubit.cmd(f"volume 1 size {mesh_size}")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add tet all")
    cubit.cmd('block 1 name "sphere"')
    n_tet = len(cubit.get_block_tets(1))
    n_nodes = cubit.get_node_count()
    print(f"Sphere mesh: {n_tet} tets, {n_nodes} nodes")
    return n_tet, n_nodes


def parse_gmsh_v41(filename):
    """Parse GMSH v4.1 file and return nodes dict and element connectivity.

    v4.1 format (block-structured):
      $Nodes
      <numEntityBlocks numNodes minTag maxTag>
      for each block:
        <entityDim entityTag parametric numNodesInBlock>
        nodeTag1
        nodeTag2
        ...
        x1 y1 z1
        x2 y2 z2
        ...
      $EndNodes

      $Elements (analogous, one block per entity/type)
    """
    nodes = {}
    elements = []

    with open(filename, 'r') as f:
        lines = [ln.strip() for ln in f.readlines()]

    i = 0
    n = len(lines)
    while i < n:
        ln = lines[i]
        if ln == '$Nodes':
            i += 1
            hdr = lines[i].split()
            n_blocks = int(hdr[0])
            i += 1
            for _ in range(n_blocks):
                bh = lines[i].split()
                n_in_block = int(bh[3])
                i += 1
                tag_ids = [int(lines[i + k]) for k in range(n_in_block)]
                i += n_in_block
                for k in range(n_in_block):
                    parts = lines[i + k].split()
                    nodes[tag_ids[k]] = (
                        float(parts[0]), float(parts[1]), float(parts[2]))
                i += n_in_block
            while i < n and lines[i] != '$EndNodes':
                i += 1
        elif ln == '$Elements':
            i += 1
            hdr = lines[i].split()
            n_blocks = int(hdr[0])
            i += 1
            for _ in range(n_blocks):
                bh = lines[i].split()
                etype = int(bh[2])
                n_in_block = int(bh[3])
                i += 1
                for k in range(n_in_block):
                    parts = lines[i + k].split()
                    conn = [int(p) for p in parts[1:]]
                    elements.append((etype, conn))
                i += n_in_block
            while i < n and lines[i] != '$EndElements':
                i += 1
        i += 1

    print(f"  Parsed {filename}: {len(nodes)} nodes, {len(elements)} elements")
    return nodes, elements


# Back-compat alias: older callers (and any user scripts that imported
# the v2.2 parser) keep working against v4.1-only output.
parse_gmsh_v22 = parse_gmsh_v41


def tet4_volume(p0, p1, p2, p3):
    """Volume of a linear tetrahedron (4 vertices)."""
    v1 = np.array(p1) - np.array(p0)
    v2 = np.array(p2) - np.array(p0)
    v3 = np.array(p3) - np.array(p0)
    return abs(np.dot(v1, np.cross(v2, v3))) / 6.0


def tet10_volume_gauss(nodes_dict, conn):
    """
    Volume of a TET10 element using 4-point Gauss quadrature.

    TET10 shape functions (in barycentric coords L1,L2,L3,L4):
      N0 = L1*(2*L1-1)  (vertex 0)
      N1 = L2*(2*L2-1)  (vertex 1)
      N2 = L3*(2*L3-1)  (vertex 2)
      N3 = L4*(2*L4-1)  (vertex 3)
      N4 = 4*L1*L2      (edge 0-1)
      N5 = 4*L2*L3      (edge 1-2)
      N6 = 4*L1*L3      (edge 0-2)
      N7 = 4*L1*L4      (edge 0-3)
      N8 = 4*L2*L4      (edge 1-3)
      N9 = 4*L3*L4      (edge 2-3)

    where L4 = 1 - L1 - L2 - L3.
    """
    # 4-point Gauss quadrature for tetrahedra (degree 2)
    a = 0.1381966011250105
    b = 0.5854101966249685
    gauss_pts = [
        (a, a, a, b),
        (a, a, b, a),
        (a, b, a, a),
        (b, a, a, a),
    ]
    gauss_w = 1.0 / 24.0  # each weight = 1/24 (reference tet volume = 1/6)

    # Get node coordinates
    coords = []
    for nid in conn:
        coords.append(np.array(nodes_dict[nid]))

    total_vol = 0.0

    for L1, L2, L3, L4 in gauss_pts:
        # Shape functions
        N = np.array([
            L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1), L4*(2*L4-1),
            4*L1*L2, 4*L2*L3, 4*L1*L3, 4*L1*L4, 4*L2*L4, 4*L3*L4
        ])

        # Shape function derivatives w.r.t. (L1, L2, L3)
        # dN/dL1:
        dNdL1 = np.array([
            4*L1-1, 0, 0, -(4*L4-1),
            4*L2, 0, 4*L3, 4*(L4-L1), -4*L2, -4*L3
        ])
        # dN/dL2:
        dNdL2 = np.array([
            0, 4*L2-1, 0, -(4*L4-1),
            4*L1, 4*L3, 0, -4*L1, 4*(L4-L2), -4*L3
        ])
        # dN/dL3:
        dNdL3 = np.array([
            0, 0, 4*L3-1, -(4*L4-1),
            0, 4*L2, 4*L1, -4*L1, -4*L2, 4*(L4-L3)
        ])

        # Jacobian: dx/dLi
        dxdL1 = sum(dNdL1[k] * coords[k] for k in range(10))
        dxdL2 = sum(dNdL2[k] * coords[k] for k in range(10))
        dxdL3 = sum(dNdL3[k] * coords[k] for k in range(10))

        J = np.array([dxdL1, dxdL2, dxdL3])
        detJ = abs(np.linalg.det(J))

        total_vol += gauss_w * detJ

    return total_vol


def compute_volume(nodes_dict, elements):
    """Compute total volume from tet elements (type 4=TET4, type 11=TET10)."""
    vol = 0.0
    n_tet4 = 0
    n_tet10 = 0
    n_other = 0

    for etype, conn in elements:
        if etype == 4:  # TET4
            p = [nodes_dict[c] for c in conn[:4]]
            vol += tet4_volume(p[0], p[1], p[2], p[3])
            n_tet4 += 1
        elif etype == 11:  # TET10
            vol += tet10_volume_gauss(nodes_dict, conn)
            n_tet10 += 1
        else:
            n_other += 1  # skip non-tet (e.g., tri boundary)

    print(f"  Volume elements: TET4={n_tet4}, TET10={n_tet10}, other={n_other}")
    return vol


def check_midnode_displacement(nodes_dict, elements, radius):
    """Check if TET10 mid-nodes on surface are projected onto the sphere."""
    max_disp = 0.0
    n_surface_mid = 0
    n_interior_mid = 0

    for etype, conn in elements:
        if etype != 11:
            continue
        # Mid-nodes are conn[4:10]
        vertices = conn[:4]
        midnodes = conn[4:10]
        # Edge table for TET10: same as in HighOrderMesh.cpp
        edges = [(0,1),(1,2),(0,2),(0,3),(1,3),(2,3)]

        for idx, (i, j) in enumerate(edges):
            mid_id = midnodes[idx]
            v0 = np.array(nodes_dict[vertices[i]])
            v1 = np.array(nodes_dict[vertices[j]])
            mid = np.array(nodes_dict[mid_id])

            # Linear midpoint
            linear_mid = 0.5 * (v0 + v1)
            disp = np.linalg.norm(mid - linear_mid)

            # Check if both endpoints are on the sphere surface
            r0 = np.linalg.norm(v0)
            r1 = np.linalg.norm(v1)
            r_mid = np.linalg.norm(mid)

            if abs(r0 - radius) < 0.002 and abs(r1 - radius) < 0.002:
                # Both endpoints on surface
                n_surface_mid += 1
                if disp > max_disp:
                    max_disp = disp
            else:
                n_interior_mid += 1

    print(f"  Surface mid-nodes: {n_surface_mid}, "
          f"Interior mid-nodes: {n_interior_mid}")
    print(f"  Max displacement of surface mid-nodes from linear: {max_disp:.6e} m")
    return max_disp, n_surface_mid


def main():
    print("=" * 70)
    print("Order=2 GMSH Export: Volume Accuracy Test")
    print(f"Sphere radius = {RADIUS} m, V_exact = {V_EXACT:.6e} m^3")
    print("=" * 70)

    for mesh_size in [0.02, 0.01]:
        print(f"\n--- Mesh size = {mesh_size} ---")
        n_tet, n_nodes = build_sphere_mesh(mesh_size)

        # Export order=1
        f1 = os.path.join(OUT_DIR, f"sphere_order1_sz{mesh_size}.msh")
        cubit.cmd(f'export gmsh "{f1}" order 1 overwrite')

        # Export order=2
        f2 = os.path.join(OUT_DIR, f"sphere_order2_sz{mesh_size}.msh")
        cubit.cmd(f'export gmsh "{f2}" order 2 overwrite')

        # Parse and compute volumes
        print("\nOrder 1 (TET4):")
        nodes1, elems1 = parse_gmsh_v41(f1)
        vol1 = compute_volume(nodes1, elems1)
        err1 = (vol1 - V_EXACT) / V_EXACT * 100.0
        print(f"  V = {vol1:.6e} m^3, error = {err1:+.2f}%")

        print("\nOrder 2 (TET10):")
        nodes2, elems2 = parse_gmsh_v41(f2)
        vol2 = compute_volume(nodes2, elems2)
        err2 = (vol2 - V_EXACT) / V_EXACT * 100.0
        print(f"  V = {vol2:.6e} m^3, error = {err2:+.2f}%")

        # Check if mid-nodes are actually displaced
        print("\nMid-node displacement check:")
        max_disp, n_surf = check_midnode_displacement(nodes2, elems2, RADIUS)

        # Compare node counts
        print(f"\nNode count: order1={len(nodes1)}, order2={len(nodes2)}")
        if len(nodes2) > len(nodes1):
            print(f"  -> {len(nodes2) - len(nodes1)} mid-nodes added (OK)")
        else:
            print("  -> WARNING: No additional nodes! order=2 may not be working")

        # Verdict
        print(f"\nVolume improvement: order1 err={err1:+.2f}%, order2 err={err2:+.2f}%")
        if abs(err2) < abs(err1) * 0.5:
            print("  -> PASS: order=2 significantly improves volume accuracy")
        elif abs(err2) < abs(err1):
            print("  -> MARGINAL: order=2 improves but not dramatically")
        elif max_disp < 1e-15:
            print("  -> FAIL: mid-nodes NOT projected (displacement = 0)")
        else:
            print("  -> FAIL: order=2 does not improve volume accuracy")

    print("\n" + "=" * 70)


if __name__ == "__main__":
    main()
