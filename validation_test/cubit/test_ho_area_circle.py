"""
Verify order=2 mesh export for 2D circle: area accuracy.

Circle (r=0.05 m, A_exact = pi*r^2 = 7.854e-3 m^2).
Exports as GMSH/Nastran/VTK order=1 and order=2, computes area
from TRI3/TRI6 elements, compares with analytical.

Usage:
  python validation_test/cubit/test_ho_area_circle.py
"""

import sys
import os
import math
import numpy as np

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

OUT_DIR = os.path.join(_test_dir, 'ho_area_test')
os.makedirs(OUT_DIR, exist_ok=True)

RADIUS = 0.05
A_EXACT = math.pi * RADIUS**2


def build_circle_mesh(mesh_size=0.015):
    cubit.cmd("reset")
    cubit.cmd(f"create surface circle radius {RADIUS} zplane")
    cubit.cmd("surface 1 scheme trimesh")
    cubit.cmd(f"surface 1 size {mesh_size}")
    cubit.cmd("mesh surface 1")
    cubit.cmd("block 1 add tri all")
    cubit.cmd('block 1 name "circle"')
    n_tri = len(cubit.get_block_tris(1))
    n_nodes = cubit.get_node_count()
    print(f"Circle mesh: {n_tri} tris, {n_nodes} nodes")
    return n_tri, n_nodes


def tri3_area(p0, p1, p2):
    v1 = np.array(p1) - np.array(p0)
    v2 = np.array(p2) - np.array(p0)
    return 0.5 * np.linalg.norm(np.cross(v1, v2))


def tri6_area_gauss(coords):
    """Area of TRI6 using 3-point Gauss quadrature."""
    # Reference triangle: (L1, L2), L3 = 1-L1-L2
    # 3-point rule (degree 2)
    gauss_pts = [(1/6, 1/6), (2/3, 1/6), (1/6, 2/3)]
    gauss_w = 1/6  # each weight (reference area = 1/2)

    area = 0.0
    for L1, L2 in gauss_pts:
        L3 = 1 - L1 - L2
        # TRI6 shape functions
        N = np.array([
            L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1),
            4*L1*L2, 4*L2*L3, 4*L3*L1
        ])
        # dN/dL1
        dNdL1 = np.array([4*L1-1, 0, -(4*L3-1), 4*L2, -4*L2, 4*(L3-L1)])
        # dN/dL2
        dNdL2 = np.array([0, 4*L2-1, -(4*L3-1), 4*L1, 4*(L3-L2), -4*L1])

        dxdL1 = sum(dNdL1[k] * coords[k] for k in range(6))
        dxdL2 = sum(dNdL2[k] * coords[k] for k in range(6))

        # Cross product magnitude (for 3D surface area)
        cross = np.cross(dxdL1, dxdL2)
        area += gauss_w * np.linalg.norm(cross)

    return area


def parse_gmsh_area(filename):
    """Parse GMSH v2.2 or v4.1, compute area from TRI3 (type 2) / TRI6 (type 9)."""
    with open(filename, 'r') as f:
        content = f.read()
    lines = content.split('\n')

    # Detect version
    version = '2.2'
    for line in lines:
        s = line.strip()
        if s.startswith('4.1') or s.startswith('4.0'):
            version = '4.1'
            break
        elif s.startswith('2.'):
            version = '2.2'
            break

    nodes = {}
    elements = []

    if version.startswith('4'):
        # v4.1 parser
        i = 0
        while i < len(lines):
            if lines[i].strip() == '$Nodes':
                i += 1
                header = lines[i].strip().split()
                num_blocks = int(header[0])
                i += 1
                for _ in range(num_blocks):
                    bh = lines[i].strip().split()
                    nb = int(bh[3])
                    i += 1
                    nids = [int(lines[i+j].strip()) for j in range(nb)]
                    i += nb
                    for j, nid in enumerate(nids):
                        p = lines[i].strip().split()
                        nodes[nid] = np.array([float(p[0]), float(p[1]), float(p[2])])
                        i += 1
                break
            i += 1
        i = 0
        while i < len(lines):
            if lines[i].strip() == '$Elements':
                i += 1
                header = lines[i].strip().split()
                num_blocks = int(header[0])
                i += 1
                for _ in range(num_blocks):
                    bh = lines[i].strip().split()
                    etype = int(bh[2])
                    ne = int(bh[3])
                    i += 1
                    for _ in range(ne):
                        p = lines[i].strip().split()
                        conn = [int(x) for x in p[1:]]
                        if etype in (2, 9):
                            elements.append((etype, conn))
                        i += 1
                break
            i += 1
    else:
        # v2.2 parser
        section = None
        n_count = [0]
        e_count = [0]
        for line in lines:
            s = line.strip()
            if s.startswith('$'):
                section = s if not s.startswith('$End') else None
                continue
            if section == '$Nodes':
                if n_count[0] == 0:
                    n_count[0] = int(s)
                else:
                    p = s.split()
                    if len(p) >= 4:
                        nodes[int(p[0])] = np.array([float(p[1]), float(p[2]), float(p[3])])
            elif section == '$Elements':
                if e_count[0] == 0:
                    e_count[0] = int(s)
                else:
                    p = s.split()
                    if len(p) >= 4:
                        etype = int(p[1])
                        ntags = int(p[2])
                        conn = [int(x) for x in p[3+ntags:]]
                        if etype in (2, 9):
                            elements.append((etype, conn))

    return _compute_tri_area(nodes, elements)


def parse_nastran_area(filename):
    """Parse Nastran BDF, compute area from CTRIA3/CTRIA6."""
    import re
    with open(filename, 'r') as f:
        content = f.read()

    nodes = {}
    elements = []

    # GRID*
    grid_pat = re.compile(
        r'GRID\*\s+(\d+)\s+\d+\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)\s*\n'
        r'\*\s+([\d.eE+\-]+)', re.MULTILINE)
    for m in grid_pat.finditer(content):
        nid = int(m.group(1))
        nodes[nid] = np.array([float(m.group(2)), float(m.group(3)), float(m.group(4))])

    lines = content.split('\n')
    for line in lines:
        if line.startswith('CTRIA3'):
            p = line.split()
            conn = [int(x) for x in p[3:6]]
            elements.append((3, conn))
        elif line.startswith('CTRIA6'):
            p = line.split()
            conn = [int(x) for x in p[3:9]]
            elements.append((6, conn))

    return _compute_tri_area(nodes, elements, nastran=True)


def parse_vtk_area(filename):
    """Parse VTK legacy, compute area from TRI (type 5) / TRI6 (type 22)."""
    with open(filename, 'r') as f:
        lines = f.readlines()

    nodes = {}
    cells = []
    cell_types = []
    section = None
    n_nodes = 0
    n_cells = 0
    ni = 0
    ci = 0

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('POINTS'):
            n_nodes = int(line.split()[1])
            section = 'points'
            ni = 0
            i += 1
            continue
        elif line.startswith('CELLS'):
            n_cells = int(line.split()[1])
            section = 'cells'
            ci = 0
            i += 1
            continue
        elif line.startswith('CELL_TYPES'):
            section = 'cell_types'
            i += 1
            continue

        if section == 'points' and ni < n_nodes:
            parts = line.split()
            for j in range(0, len(parts), 3):
                if ni < n_nodes and j+2 < len(parts):
                    nodes[ni] = np.array([float(parts[j]), float(parts[j+1]), float(parts[j+2])])
                    ni += 1
        elif section == 'cells' and ci < n_cells:
            parts = line.split()
            nn = int(parts[0])
            cells.append([int(p) for p in parts[1:nn+1]])
            ci += 1
        elif section == 'cell_types':
            for p in line.split():
                try:
                    cell_types.append(int(p))
                except ValueError:
                    pass
        i += 1

    elements = []
    for ci in range(min(len(cells), len(cell_types))):
        ct = cell_types[ci]
        if ct == 5:  # VTK_TRIANGLE
            elements.append((3, cells[ci]))
        elif ct == 22:  # VTK_QUADRATIC_TRIANGLE
            elements.append((6, cells[ci]))

    return _compute_tri_area(nodes, elements, vtk=True)


def _compute_tri_area(nodes, elements, nastran=False, vtk=False):
    area = 0.0
    n_tri3 = n_tri6 = 0
    for etype, conn in elements:
        if (etype == 2 or etype == 3) and not vtk:
            # GMSH type 2 or Nastran TRI3
            coords = [nodes[c] for c in conn[:3]]
            area += tri3_area(coords[0], coords[1], coords[2])
            n_tri3 += 1
        elif etype == 3 and vtk:
            coords = [nodes[c] for c in conn[:3]]
            area += tri3_area(coords[0], coords[1], coords[2])
            n_tri3 += 1
        elif (etype == 9 or etype == 6) and not vtk:
            # GMSH type 9 or Nastran TRI6
            coords = [nodes[c] for c in conn[:6]]
            area += tri6_area_gauss(coords)
            n_tri6 += 1
        elif etype == 6 and vtk:
            coords = [nodes[c] for c in conn[:6]]
            area += tri6_area_gauss(coords)
            n_tri6 += 1
    return area, len(nodes), n_tri3, n_tri6


def main():
    print("=" * 70)
    print("Order=2 Area Accuracy: Circle (2D surface mesh)")
    print(f"Circle r={RADIUS} m, A_exact={A_EXACT:.6e} m^2")
    print("=" * 70)

    mesh_size = 0.012
    n_tri, n_nodes = build_circle_mesh(mesh_size)

    test_cases = [
        ("GMSH v4.1", ".msh",
         'export gmsh "{f}" order {o} dimension 2 overwrite',
         parse_gmsh_area),
        ("Nastran", ".bdf",
         'export jmag_nastran "{f}" order {o} dimension 2 overwrite',
         parse_nastran_area),
        ("VTK", ".vtk",
         'export vtk "{f}" order {o} dimension 2 overwrite',
         parse_vtk_area),
    ]

    results = []
    all_pass = True

    for fmt_name, ext, cmd_template, parser in test_cases:
        print(f"\n{'=' * 50}")
        print(f"  {fmt_name}")
        print(f"{'=' * 50}")

        for order in [1, 2]:
            fname = os.path.join(OUT_DIR, f"circle_{fmt_name.replace(' ','_')}_o{order}{ext}")
            cmd = cmd_template.format(f=fname, o=order)
            cubit.cmd(cmd)

            area, nn, n3, n6 = parser(fname)
            err = (area - A_EXACT) / A_EXACT * 100.0
            label = f"order={order}"
            info = f"TRI3={n3}" if order == 1 else f"TRI6={n6}"
            print(f"  {label}: A={area:.6e}, err={err:+.4f}%, nodes={nn}, {info}")
            results.append({
                'format': fmt_name, 'order': order,
                'area': area, 'error_pct': err,
            })

        r1 = results[-2]
        r2 = results[-1]
        ratio = abs(r1['error_pct']) / abs(r2['error_pct']) if abs(r2['error_pct']) > 1e-10 else float('inf')
        verdict = "PASS" if abs(r2['error_pct']) < abs(r1['error_pct']) * 0.5 else "FAIL"
        if verdict == "FAIL":
            all_pass = False
        print(f"  Improvement: {ratio:.0f}x -> [{verdict}]")

    # Summary
    print(f"\n{'=' * 70}")
    print(f"{'Format':<15} {'Order 1 err':>12} {'Order 2 err':>12} {'Ratio':>8} {'Verdict':>8}")
    print(f"{'-' * 70}")
    for i in range(0, len(results), 2):
        r1 = results[i]
        r2 = results[i+1]
        ratio = abs(r1['error_pct']) / abs(r2['error_pct']) if abs(r2['error_pct']) > 1e-10 else float('inf')
        verdict = "PASS" if abs(r2['error_pct']) < abs(r1['error_pct']) * 0.5 else "FAIL"
        print(f"{r1['format']:<15} {r1['error_pct']:>+11.4f}% {r2['error_pct']:>+11.4f}% {ratio:>7.0f}x {verdict:>8}")
    print(f"{'=' * 70}")

    if all_pass:
        print("\nALL PASS: order=2 circle area accuracy verified.")
    else:
        print("\nFAILURE detected.")

    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
