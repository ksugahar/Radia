"""
Verify order=2 mesh export quality across ALL formats by comparing volumes.

Sphere (r=0.05 m, V_exact = 5.236e-4 m^3):
  - GMSH v2.2:  order 1/2, TET4/TET10
  - Nastran:    order 1/2, CTETRA/CTETRA(10)
  - VTK:        order 1/2, cell type 10/24

If order=2 projection is working, TET10 volume should be closer
to the analytical value than TET4 in ALL formats.

Usage:
  python tests/cubit/test_ho_volume_all_formats.py
"""

import sys
import os
import math
import re
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

RADIUS = 0.05
V_EXACT = (4.0 / 3.0) * math.pi * RADIUS**3


# ================================================================
# Cubit model builder
# ================================================================
def build_sphere_mesh(mesh_size=0.015):
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
    return n_tet, n_nodes


# ================================================================
# TET10 Gauss quadrature volume
# ================================================================
def tet4_volume(p0, p1, p2, p3):
    v1 = np.array(p1) - np.array(p0)
    v2 = np.array(p2) - np.array(p0)
    v3 = np.array(p3) - np.array(p0)
    return abs(np.dot(v1, np.cross(v2, v3))) / 6.0


def tet10_volume_gauss(coords):
    """Volume of TET10 using 4-point Gauss quadrature.
    coords: list of 10 arrays, shape (10,3)."""
    a = 0.1381966011250105
    b = 0.5854101966249685
    gauss_pts = [(a,a,a,b), (a,a,b,a), (a,b,a,a), (b,a,a,a)]

    vol = 0.0
    for L1, L2, L3, L4 in gauss_pts:
        N = np.array([
            L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1), L4*(2*L4-1),
            4*L1*L2, 4*L2*L3, 4*L1*L3, 4*L1*L4, 4*L2*L4, 4*L3*L4
        ])
        dNdL1 = np.array([
            4*L1-1, 0, 0, -(4*L4-1),
            4*L2, 0, 4*L3, 4*(L4-L1), -4*L2, -4*L3
        ])
        dNdL2 = np.array([
            0, 4*L2-1, 0, -(4*L4-1),
            4*L1, 4*L3, 0, -4*L1, 4*(L4-L2), -4*L3
        ])
        dNdL3 = np.array([
            0, 0, 4*L3-1, -(4*L4-1),
            0, 4*L2, 4*L1, -4*L1, -4*L2, 4*(L4-L3)
        ])

        dxdL1 = sum(dNdL1[k] * coords[k] for k in range(10))
        dxdL2 = sum(dNdL2[k] * coords[k] for k in range(10))
        dxdL3 = sum(dNdL3[k] * coords[k] for k in range(10))

        J = np.array([dxdL1, dxdL2, dxdL3])
        vol += abs(np.linalg.det(J)) / 24.0

    return vol


# ================================================================
# GMSH v2.2 parser
# ================================================================
def parse_gmsh_volume(filename):
    """Parse GMSH v2.2 or v4.1 and compute volume from TET4/TET10."""
    with open(filename, 'r') as f:
        content = f.read()

    # Detect version
    if '$MeshFormat' in content:
        fmt_line = content.split('$MeshFormat\n')[1].split('\n')[0]
        version = fmt_line.split()[0]
    else:
        version = '2.2'

    if version.startswith('4'):
        return _parse_gmsh_v41_volume(content)
    else:
        return _parse_gmsh_v22_volume(content)


def _parse_gmsh_v22_volume(content):
    """Parse GMSH v2.2 format."""
    nodes = {}
    elements = []
    lines = content.split('\n')

    section = None
    n_nodes = 0
    n_elems = 0

    for line in lines:
        line = line.strip()
        if line.startswith('$'):
            if line == '$Nodes':
                section = 'nodes'
            elif line == '$Elements':
                section = 'elements'
            elif line.startswith('$End'):
                section = None
            continue

        if section == 'nodes':
            if n_nodes == 0:
                n_nodes = int(line)
            else:
                parts = line.split()
                if len(parts) >= 4:
                    nid = int(parts[0])
                    nodes[nid] = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
        elif section == 'elements':
            if n_elems == 0:
                n_elems = int(line)
            else:
                parts = line.split()
                if len(parts) >= 4:
                    etype = int(parts[1])
                    ntags = int(parts[2])
                    conn = [int(p) for p in parts[3 + ntags:]]
                    if etype in (4, 11):
                        elements.append((etype, conn))

    return _compute_tet_volume(nodes, elements)


def _parse_gmsh_v41_volume(content):
    """Parse GMSH v4.1 format."""
    nodes = {}
    elements = []
    lines = content.split('\n')

    # Parse $Nodes
    i = 0
    while i < len(lines):
        if lines[i].strip() == '$Nodes':
            i += 1
            header = lines[i].strip().split()
            num_blocks = int(header[0])
            total_nodes = int(header[1])
            i += 1

            for _ in range(num_blocks):
                block_header = lines[i].strip().split()
                # entityDim entityTag parametric numNodes
                num_block_nodes = int(block_header[3])
                i += 1

                # Read node tags
                node_ids = []
                for _ in range(num_block_nodes):
                    node_ids.append(int(lines[i].strip()))
                    i += 1

                # Read node coordinates
                for j, nid in enumerate(node_ids):
                    parts = lines[i].strip().split()
                    nodes[nid] = np.array([float(parts[0]), float(parts[1]), float(parts[2])])
                    i += 1

            break
        i += 1

    # Parse $Elements
    i = 0
    while i < len(lines):
        if lines[i].strip() == '$Elements':
            i += 1
            header = lines[i].strip().split()
            num_blocks = int(header[0])
            i += 1

            for _ in range(num_blocks):
                block_header = lines[i].strip().split()
                # entityDim entityTag elementType numElements
                etype = int(block_header[2])
                num_block_elems = int(block_header[3])
                i += 1

                for _ in range(num_block_elems):
                    parts = lines[i].strip().split()
                    # elementTag node1 node2 ...
                    conn = [int(p) for p in parts[1:]]
                    if etype in (4, 11):
                        elements.append((etype, conn))
                    i += 1

            break
        i += 1

    return _compute_tet_volume(nodes, elements)


def _compute_tet_volume(nodes, elements):
    """Compute total volume from tet elements."""
    vol = 0.0
    n_tet4 = n_tet10 = 0
    for etype, conn in elements:
        if etype == 4:
            p = [nodes[c] for c in conn[:4]]
            vol += tet4_volume(p[0], p[1], p[2], p[3])
            n_tet4 += 1
        elif etype == 11:
            coords = [nodes[c] for c in conn[:10]]
            vol += tet10_volume_gauss(coords)
            n_tet10 += 1
    return vol, len(nodes), n_tet4, n_tet10


# ================================================================
# Nastran BDF parser
# ================================================================
def parse_nastran_volume(filename):
    """Parse Nastran BDF and compute volume from CTETRA / CTETRA(10)."""
    nodes = {}
    elements = []

    with open(filename, 'r') as f:
        content = f.read()

    # Parse GRID* (large-field: 2 lines per node)
    grid_pattern = re.compile(
        r'GRID\*\s+(\d+)\s+\d+\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)\s*\n'
        r'\*\s+([\d.eE+\-]+)',
        re.MULTILINE
    )
    for m in grid_pattern.finditer(content):
        nid = int(m.group(1))
        x, y, z = float(m.group(2)), float(m.group(3)), float(m.group(4))
        nodes[nid] = np.array([x, y, z])

    # Also try small-field GRID (8-char fields)
    if len(nodes) == 0:
        for line in content.split('\n'):
            if line.startswith('GRID') and not line.startswith('GRID*'):
                parts = line.split()
                if len(parts) >= 6:
                    nid = int(parts[1])
                    x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                    nodes[nid] = np.array([x, y, z])

    # Parse CTETRA cards
    # Format: CTETRA eid pid n1 n2 n3 n4 [n5 n6+]
    #         +      n7 n8 n9 n10
    # The '+' suffix on last field means continuation follows
    def parse_nastran_int(s):
        """Parse integer, stripping any continuation marker (+)."""
        return int(s.rstrip('+'))

    lines = content.split('\n')
    i = 0
    n_tet4 = n_tet10 = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('CTETRA'):
            parts = line.split()
            conn = [parse_nastran_int(p) for p in parts[3:] if p.strip()]

            # Check for continuation line(s)
            while i + 1 < len(lines) and lines[i + 1].startswith('+'):
                i += 1
                cparts = lines[i].split()
                conn.extend([parse_nastran_int(p) for p in cparts[1:] if p.strip()])

            if len(conn) >= 10:
                coords = [nodes[c] for c in conn[:10]]
                vol = tet10_volume_gauss(coords)
                elements.append(('TET10', vol))
                n_tet10 += 1
            elif len(conn) >= 4:
                p = [nodes[c] for c in conn[:4]]
                vol = tet4_volume(p[0], p[1], p[2], p[3])
                elements.append(('TET4', vol))
                n_tet4 += 1
        i += 1

    total_vol = sum(v for _, v in elements)
    return total_vol, len(nodes), n_tet4, n_tet10


# ================================================================
# VTK Legacy parser
# ================================================================
def parse_vtk_volume(filename):
    """Parse VTK Legacy and compute volume from TET4 (type 10) / TET10 (type 24)."""
    nodes = {}
    cells = []
    cell_types = []

    with open(filename, 'r') as f:
        lines = f.readlines()

    section = None
    n_nodes = 0
    n_cells = 0
    node_idx = 0
    cell_idx = 0

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith('POINTS'):
            parts = line.split()
            n_nodes = int(parts[1])
            section = 'points'
            node_idx = 0
            i += 1
            continue
        elif line.startswith('CELLS'):
            parts = line.split()
            n_cells = int(parts[1])
            section = 'cells'
            cell_idx = 0
            i += 1
            continue
        elif line.startswith('CELL_TYPES'):
            section = 'cell_types'
            i += 1
            continue

        if section == 'points' and node_idx < n_nodes:
            parts = line.split()
            # VTK can have multiple points per line
            for j in range(0, len(parts), 3):
                if node_idx < n_nodes and j + 2 < len(parts):
                    nodes[node_idx] = np.array([
                        float(parts[j]), float(parts[j+1]), float(parts[j+2])
                    ])
                    node_idx += 1
        elif section == 'cells' and cell_idx < n_cells:
            parts = line.split()
            nn = int(parts[0])
            conn = [int(p) for p in parts[1:nn+1]]
            cells.append(conn)
            cell_idx += 1
        elif section == 'cell_types':
            parts = line.split()
            for p in parts:
                try:
                    cell_types.append(int(p))
                except ValueError:
                    pass

        i += 1

    vol = 0.0
    n_tet4 = n_tet10 = 0
    for ci in range(min(len(cells), len(cell_types))):
        ct = cell_types[ci]
        conn = cells[ci]
        if ct == 10 and len(conn) >= 4:  # VTK_TETRA
            p = [nodes[c] for c in conn[:4]]
            vol += tet4_volume(p[0], p[1], p[2], p[3])
            n_tet4 += 1
        elif ct == 24 and len(conn) >= 10:  # VTK_QUADRATIC_TETRA
            coords = [nodes[c] for c in conn[:10]]
            vol += tet10_volume_gauss(coords)
            n_tet10 += 1

    return vol, len(nodes), n_tet4, n_tet10


# ================================================================
# Main test
# ================================================================
def main():
    print("=" * 70)
    print("Order=2 Volume Accuracy: ALL Formats (GMSH, Nastran, VTK)")
    print(f"Sphere r={RADIUS} m, V_exact={V_EXACT:.6e} m^3")
    print("=" * 70)

    mesh_size = 0.015
    n_tet, n_nodes = build_sphere_mesh(mesh_size)
    print(f"Mesh: {n_tet} tets, {n_nodes} nodes, size={mesh_size}")

    results = []

    # Define test cases: (format_name, ext, export_cmd_template, parser)
    test_cases = [
        ("GMSH v2.2", ".msh",
         'export gmsh "{f}" order {o} version 2 overwrite',
         parse_gmsh_volume),
        ("GMSH v4.1", ".msh",
         'export gmsh "{f}" order {o} version 4 overwrite',
         parse_gmsh_volume),
        ("Nastran", ".bdf",
         'export nastran "{f}" order {o} dimension 3 overwrite',
         parse_nastran_volume),
        ("VTK", ".vtk",
         'export vtk "{f}" order {o} dimension 3 overwrite',
         parse_vtk_volume),
    ]

    all_pass = True

    for fmt_name, ext, cmd_template, parser in test_cases:
        print(f"\n{'─' * 50}")
        print(f"  {fmt_name}")
        print(f"{'─' * 50}")

        for order in [1, 2]:
            fname = os.path.join(OUT_DIR, f"vol_test_{fmt_name.replace(' ','_')}_o{order}{ext}")
            cmd = cmd_template.format(f=fname, o=order)
            cubit.cmd(cmd)

            vol, nn, n4, n10 = parser(fname)
            err = (vol - V_EXACT) / V_EXACT * 100.0

            label = f"order={order}"
            tet_info = f"TET4={n4}" if order == 1 else f"TET10={n10}"
            print(f"  {label}: V={vol:.6e}, err={err:+.4f}%, "
                  f"nodes={nn}, {tet_info}")

            results.append({
                'format': fmt_name, 'order': order,
                'volume': vol, 'error_pct': err,
                'nodes': nn, 'n_tet4': n4, 'n_tet10': n10
            })

        # Compare order 1 vs 2
        r1 = results[-2]
        r2 = results[-1]
        ratio = abs(r1['error_pct']) / abs(r2['error_pct']) if abs(r2['error_pct']) > 1e-10 else float('inf')

        if abs(r2['error_pct']) < abs(r1['error_pct']) * 0.5:
            verdict = "PASS"
        elif abs(r2['error_pct']) < abs(r1['error_pct']):
            verdict = "MARGINAL"
        else:
            verdict = "FAIL"
            all_pass = False

        print(f"  Improvement: {ratio:.0f}x -> [{verdict}]")

    # Cross-format consistency check
    print(f"\n{'=' * 70}")
    print("Cross-Format Consistency Check")
    print(f"{'=' * 70}")

    for order in [1, 2]:
        vols = [(r['format'], r['volume']) for r in results if r['order'] == order]
        if len(vols) < 2:
            continue
        ref_fmt, ref_vol = vols[0]
        print(f"\n  Order={order} (reference: {ref_fmt}, V={ref_vol:.6e}):")
        for fmt, vol in vols[1:]:
            diff_pct = (vol - ref_vol) / ref_vol * 100.0
            status = "OK" if abs(diff_pct) < 0.01 else "MISMATCH"
            if status == "MISMATCH":
                all_pass = False
            print(f"    {fmt:<15} V={vol:.6e}  diff={diff_pct:+.6f}%  [{status}]")

    # Summary table
    print(f"\n{'=' * 70}")
    print(f"{'Format':<15} {'Order 1 err':>12} {'Order 2 err':>12} {'Ratio':>8} {'Verdict':>8}")
    print(f"{'─' * 70}")
    for i in range(0, len(results), 2):
        r1 = results[i]
        r2 = results[i + 1]
        ratio = abs(r1['error_pct']) / abs(r2['error_pct']) if abs(r2['error_pct']) > 1e-10 else float('inf')
        verdict = "PASS" if abs(r2['error_pct']) < abs(r1['error_pct']) * 0.5 else "FAIL"
        print(f"{r1['format']:<15} {r1['error_pct']:>+11.4f}% {r2['error_pct']:>+11.4f}% {ratio:>7.0f}x {verdict:>8}")
    print(f"{'=' * 70}")

    if all_pass:
        print("\nALL PASS: order=2 projection and cross-format consistency verified.")
    else:
        print("\nFAILURE: see details above.")

    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
