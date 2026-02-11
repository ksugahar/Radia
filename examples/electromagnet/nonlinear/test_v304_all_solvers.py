#!/usr/bin/env python
"""
Test: V304 EIEM2 mesh with all solvers (LU, BiCGSTAB, HACApK).

Verifies:
1. Hex-only quarter (74 hex) + full coil + IMA (+x-z) with all 3 solvers
2. Hex-only full model (4x74=296 hex, mirrored) + no IMA (LU)
3. Hex+wedge quarter (74 hex + 2 wedge) + full coil + IMA (+x-z) (LU)
4. IMA correctness: quarter+IMA vs full model
5. Comparison with both ELF references (hex-only and hex+wedge)

V304 quarter mesh: X=[0,65], Y=[-20,202.5], Z=[0,105] mm
ELF reference (hex-only):   Bz = -963.23 mT at origin
ELF reference (hex+wedge):  Bz = -961.47 mT at origin
"""

import sys
import os
import time
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(work_dir)))
sys.path.insert(0, os.path.join(repo_root, 'src'))
sys.path.insert(0, work_dir)

import radia as rad
from coil_model import create_racetrack_coil

# ELF mesh paths
ELF_V304_BASE = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\nonlinear_20000AT\ELF_MMB8T_EIEM2_V304"
ELF_HEX_ONLY = os.path.join(ELF_V304_BASE, "MMB8T_only")
ELF_HEX_WEDGE = ELF_V304_BASE  # includes hex + wedge

scale = 0.001  # mm to m
SOLVER_NAMES = {0: "LU", 1: "BiCGSTAB", 2: "HACApK"}


def load_elf_geometry(path, load_wedge=False):
    """Load ELF geometry from .meg file.

    Returns (nodes, hex_elements, wedge_elements).
    hex_elements: list of (elem_id, [8 node_ids])
    wedge_elements: list of (elem_id, [6 node_ids]) if load_wedge=True
    """
    nodes = {}
    hex_elements = []
    wedge_elements = []
    with open(os.path.join(path, "ELF_magic.meg"), 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('MGR1'):
                parts = line.split()
                node_id = int(parts[1])
                x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                nodes[node_id] = np.array([x, y, z])
            elif line.startswith('MMB8T'):
                parts = line.split()
                elem_id = int(parts[1])
                node_ids = [int(parts[i]) for i in range(4, 12)]
                hex_elements.append((elem_id, node_ids))
            elif load_wedge and line.startswith('MMB6T'):
                parts = line.split()
                elem_id = int(parts[1])
                node_ids = [int(parts[i]) for i in range(4, 10)]
                wedge_elements.append((elem_id, node_ids))
    return nodes, hex_elements, wedge_elements


def load_bh_curve(filepath):
    """Load B-H curve from text file."""
    data = np.loadtxt(filepath, comments='#')
    return data.tolist()


def load_elf_field_at_origin(path):
    """Load ELF field at (0,0,0) from .mag file."""
    with open(os.path.join(path, "ELF_magic.mag"), 'r') as f:
        for line in f:
            if line.startswith('M3GB') and not line.startswith('UNIT'):
                parts = line.split()
                if len(parts) >= 8:
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    if abs(x) < 1e-3 and abs(y) < 1e-3 and abs(z) < 1e-3:
                        Bx = float(parts[5])
                        By = float(parts[6])
                        Bz = float(parts[7])
                        return np.array([Bx, By, Bz])
    return None


def mirror_hex(verts, mx=False, mz=False):
    """Mirror hex vertices with correct vertex reordering."""
    mirrored = [[-v[0] if mx else v[0], v[1], -v[2] if mz else v[2]] for v in verts]
    n_ref = int(mx) + int(mz)
    if n_ref == 0 or n_ref == 2:
        return mirrored
    elif mx:
        return [mirrored[1], mirrored[0], mirrored[3], mirrored[2],
                mirrored[5], mirrored[4], mirrored[7], mirrored[6]]
    else:  # mz only
        return [mirrored[5], mirrored[4], mirrored[7], mirrored[6],
                mirrored[1], mirrored[0], mirrored[3], mirrored[2]]


def run_solver_quarter_ima(nodes, hex_elements, bh_data, solver_method,
                           wedge_elements=None):
    """Run quarter model + IMA and return results."""
    rad.UtiDelAll()
    rad.FldUnits('m')

    mat = rad.MatSatIsoTab(bh_data)
    all_objects = []

    # Hex elements
    for elem_id, node_ids in hex_elements:
        verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
                 for nid in node_ids]
        obj = rad.ObjHexahedron(verts, [0, 0, 0])
        rad.MatApl(obj, mat)
        all_objects.append(obj)

    # Wedge elements
    if wedge_elements:
        for elem_id, node_ids in wedge_elements:
            verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
                     for nid in node_ids]
            obj = rad.ObjWedge(verts, [0, 0, 0])
            rad.MatApl(obj, mat)
            all_objects.append(obj)

    yoke = rad.ObjCnt(all_objects)
    coil = create_racetrack_coil(20000.0)
    model = rad.ObjCnt([yoke, coil])

    if solver_method == 2:
        rad.SetHACApKParams(1e-4, 10, 2.0)

    t_start = time.time()
    result = rad.Solve(model, 0.001, 100, solver_method, image='+x-z')
    t_solve = time.time() - t_start

    radia_B = np.array(rad.Fld(model, 'b', [0, 0, 0]))

    return {
        'B': radia_B,
        't_solve': t_solve,
        'n_elements': len(all_objects),
        'iterations': int(result[2]) if len(result) > 2 else 0,
    }


def run_full_model(nodes, hex_elements, bh_data):
    """Run full model (mirrored quarter hex, no IMA) with LU solver."""
    rad.UtiDelAll()
    rad.FldUnits('m')

    hex_objects = []
    for elem_id, node_ids in hex_elements:
        verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
                 for nid in node_ids]
        for (mx, mz) in [(False, False), (True, False), (False, True), (True, True)]:
            mv = mirror_hex(verts, mx, mz)
            mat = rad.MatSatIsoTab(bh_data)
            obj = rad.ObjHexahedron(mv, [0, 0, 0])
            rad.MatApl(obj, mat)
            hex_objects.append(obj)

    yoke = rad.ObjCnt(hex_objects)
    coil = create_racetrack_coil(20000.0)
    model = rad.ObjCnt([yoke, coil])

    t_start = time.time()
    result = rad.Solve(model, 0.001, 100, 0)  # LU, no IMA
    t_solve = time.time() - t_start

    radia_B = np.array(rad.Fld(model, 'b', [0, 0, 0]))

    return {
        'B': radia_B,
        't_solve': t_solve,
        'n_elements': len(hex_objects),
        'iterations': int(result[2]) if len(result) > 2 else 0,
    }


# ============================================================
# Main
# ============================================================
print("=" * 70)
print("Test: V304 EIEM2 with All Solvers + IMA (+x-z)")
print("=" * 70)

# Load hex-only geometry
nodes_ho, hex_elems_ho, _ = load_elf_geometry(ELF_HEX_ONLY)
print(f"Hex-only quarter: {len(hex_elems_ho)} hex, DOF: {len(hex_elems_ho) * 6}")

# Load hex+wedge geometry
nodes_hw, hex_elems_hw, wedge_elems_hw = load_elf_geometry(ELF_HEX_WEDGE, load_wedge=True)
print(f"Hex+wedge quarter: {len(hex_elems_hw)} hex + {len(wedge_elems_hw)} wedge"
      f" = {len(hex_elems_hw) + len(wedge_elems_hw)} elem")

# Show coordinate ranges
all_coords = np.array(list(nodes_ho.values()))
print(f"  X: [{np.min(all_coords[:,0]):.1f}, {np.max(all_coords[:,0]):.1f}] mm")
print(f"  Y: [{np.min(all_coords[:,1]):.1f}, {np.max(all_coords[:,1]):.1f}] mm")
print(f"  Z: [{np.min(all_coords[:,2]):.1f}, {np.max(all_coords[:,2]):.1f}] mm")

# Load B-H curve
bh_file = os.path.join(work_dir, "BH.txt")
bh_data = load_bh_curve(bh_file)
print(f"B-H curve: {len(bh_data)} points")

# Load ELF references
elf_B_ho = load_elf_field_at_origin(ELF_HEX_ONLY)
elf_B_hw = load_elf_field_at_origin(ELF_HEX_WEDGE)
if elf_B_ho is not None:
    print(f"ELF reference (hex-only):  Bz = {elf_B_ho[2]*1000:.2f} mT")
if elf_B_hw is not None:
    print(f"ELF reference (hex+wedge): Bz = {elf_B_hw[2]*1000:.2f} mT")

# ============================================================
# Part 1: Hex-only quarter + IMA with all solvers
# ============================================================
print("\n" + "=" * 70)
print("Part 1: Hex-only quarter (74 elem) + full coil + IMA +x-z")
print("=" * 70)

results_ho = []
for solver_method in [0, 1, 2]:
    solver_name = SOLVER_NAMES[solver_method]
    print(f"\nTesting {solver_name} (Method {solver_method})...")

    try:
        result = run_solver_quarter_ima(nodes_ho, hex_elems_ho, bh_data, solver_method)

        rel_diff_ho = 0
        rel_diff_hw = 0
        if elf_B_ho is not None:
            rel_diff_ho = abs(result['B'][2] - elf_B_ho[2]) / abs(elf_B_ho[2]) * 100
        if elf_B_hw is not None:
            rel_diff_hw = abs(result['B'][2] - elf_B_hw[2]) / abs(elf_B_hw[2]) * 100
        status = "PASS" if rel_diff_ho < 5.0 else "WARN" if rel_diff_ho < 10.0 else "FAIL"

        print(f"  Bz = {result['B'][2]*1000:.2f} mT")
        print(f"  Time = {result['t_solve']:.2f} s")
        print(f"  vs ELF hex-only:  {rel_diff_ho:.2f}%")
        print(f"  vs ELF hex+wedge: {rel_diff_hw:.2f}%")

        results_ho.append({
            'solver': solver_name,
            'method': solver_method,
            'Bz_mT': result['B'][2] * 1000,
            'B': result['B'],
            't_solve': result['t_solve'],
            'diff_ho': rel_diff_ho,
            'diff_hw': rel_diff_hw,
            'status': status,
        })

    except Exception as e:
        print(f"  ERROR: {e}")
        results_ho.append({'solver': solver_name, 'method': solver_method, 'status': 'ERROR'})

# ============================================================
# Part 2: Hex+wedge quarter + IMA (LU only)
# ============================================================
print("\n" + "=" * 70)
print("Part 2: Hex+wedge quarter (74 hex + 2 wedge) + full coil + IMA +x-z")
print("=" * 70)

result_hw = None
try:
    result_hw = run_solver_quarter_ima(nodes_hw, hex_elems_hw, bh_data, 0,
                                       wedge_elements=wedge_elems_hw)
    diff_hw_elf = 0
    if elf_B_hw is not None:
        diff_hw_elf = abs(result_hw['B'][2] - elf_B_hw[2]) / abs(elf_B_hw[2]) * 100

    print(f"  Elements: {result_hw['n_elements']} (74 hex + 2 wedge)")
    print(f"  Bz = {result_hw['B'][2]*1000:.2f} mT")
    print(f"  Time = {result_hw['t_solve']:.2f} s")
    if elf_B_hw is not None:
        print(f"  vs ELF hex+wedge: {diff_hw_elf:.2f}%")
except Exception as e:
    print(f"  ERROR: {e}")

# ============================================================
# Part 3: Hex-only full model (mirrored, no IMA)
# ============================================================
print("\n" + "=" * 70)
print("Part 3: Hex-only full model (4x74=296 elem) + full coil, no IMA")
print("=" * 70)

full_result = run_full_model(nodes_ho, hex_elems_ho, bh_data)

full_diff_ho = 0
full_diff_hw = 0
if elf_B_ho is not None:
    full_diff_ho = abs(full_result['B'][2] - elf_B_ho[2]) / abs(elf_B_ho[2]) * 100
if elf_B_hw is not None:
    full_diff_hw = abs(full_result['B'][2] - elf_B_hw[2]) / abs(elf_B_hw[2]) * 100

print(f"  Elements: {full_result['n_elements']}")
print(f"  Bz = {full_result['B'][2]*1000:.2f} mT")
print(f"  Time = {full_result['t_solve']:.2f} s")
print(f"  vs ELF hex-only:  {full_diff_ho:.2f}%")
print(f"  vs ELF hex+wedge: {full_diff_hw:.2f}%")

# ============================================================
# Part 4: IMA vs Full model comparison
# ============================================================
print("\n" + "=" * 70)
print("Part 4: IMA vs Full Model Comparison")
print("=" * 70)

lu_result = next((r for r in results_ho if r['method'] == 0 and r.get('B') is not None), None)
if lu_result is not None:
    ima_vs_full = np.linalg.norm(lu_result['B'] - full_result['B'])
    rel_ima_full = ima_vs_full / np.linalg.norm(full_result['B']) * 100
    print(f"  Hex-only full (LU):       Bz = {full_result['B'][2]*1000:.2f} mT")
    print(f"  Hex-only quarter+IMA (LU): Bz = {lu_result['Bz_mT']:.2f} mT")
    print(f"  Vector difference:         {rel_ima_full:.4f}%")
    if rel_ima_full < 0.01:
        print(f"  --> IMA is CORRECT (matches full model)")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"\n{'Model':<40} {'Solver':<10} {'Bz(mT)':<10} {'vs ELF-HO':<10} {'vs ELF-HW':<10}")
print("-" * 80)

for r in results_ho:
    if r.get('Bz_mT') is not None:
        print(f"{'Hex-only quarter+IMA (74)':<40} {r['solver']:<10} {r['Bz_mT']:<10.2f} "
              f"{r['diff_ho']:<10.2f} {r['diff_hw']:<10.2f}")

print(f"{'Hex-only full model (296)':<40} {'LU':<10} {full_result['B'][2]*1000:<10.2f} "
      f"{full_diff_ho:<10.2f} {full_diff_hw:<10.2f}")

if result_hw is not None and result_hw.get('B') is not None:
    print(f"{'Hex+wedge quarter+IMA (76)':<40} {'LU':<10} {result_hw['B'][2]*1000:<10.2f} "
          f"{'--':<10} {diff_hw_elf:<10.2f}")

print(f"\nELF References:")
if elf_B_ho is not None:
    print(f"  Hex-only (MMB8T_only): Bz = {elf_B_ho[2]*1000:.2f} mT")
if elf_B_hw is not None:
    print(f"  Hex+wedge (V304):      Bz = {elf_B_hw[2]*1000:.2f} mT")
