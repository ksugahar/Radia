#!/usr/bin/env python
"""
Test: 1x1x1 EIEM2 mesh with all solvers (LU, BiCGSTAB, HACApK).
Verifies quarter model + IMA (+x-z) for nonlinear materials.
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

# ELF mesh path
ELF_PATH = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\nonlinear_20000AT\ELF_MMB8T_EIEM2_1x1x1\quater"

scale = 0.001  # mm to m
SOLVER_NAMES = {0: "LU", 1: "BiCGSTAB", 2: "HACApK"}


def load_elf_geometry(path):
    """Load ELF geometry from .meg file."""
    nodes = {}
    elements = []
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
                elements.append((elem_id, node_ids))
    return nodes, elements


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
                    if abs(x) < 1e-6 and abs(y) < 1e-6 and abs(z) < 1e-6:
                        Bx = float(parts[5])
                        By = float(parts[6])
                        Bz = float(parts[7])
                        return np.array([Bx, By, Bz])
    return None


def run_solver(nodes, elements, bh_data, solver_method):
    """Run solver and return results."""
    rad.UtiDelAll()
    rad.FldUnits('m')

    # Create nonlinear material
    mat = rad.MatSatIsoTab(bh_data)

    # Create yoke geometry
    hex_objects = []
    for elem_id, node_ids in elements:
        verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
                 for nid in node_ids]
        hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
        rad.MatApl(hex_obj, mat)
        hex_objects.append(hex_obj)

    yoke = rad.ObjCnt(hex_objects)

    # Create coil
    coil = create_racetrack_coil(20000.0)

    # Combine model
    model = rad.ObjCnt([yoke, coil])

    # Set HACApK params if needed
    if solver_method == 2:
        rad.SetHACApKParams(1e-4, 10, 2.0)

    # Solve
    t_start = time.time()
    result = rad.Solve(model, 0.0001, 100, solver_method, image='+x-z')
    t_solve = time.time() - t_start

    # Get field at origin
    radia_B = np.array(rad.Fld(model, 'b', [0, 0, 0]))

    return {
        'B': radia_B,
        't_solve': t_solve,
        'max_M': result[0],
        'iterations': int(result[2]) if len(result) > 2 else 0,
    }


print("=" * 70)
print("Test: EIEM2_1x1x1 with All Solvers + IMA (+x-z)")
print("=" * 70)

# Load geometry
nodes, elements = load_elf_geometry(ELF_PATH)
print(f"Elements: {len(elements)}, DOF: {len(elements) * 6}")

# Load B-H curve
bh_file = os.path.join(work_dir, "BH.txt")
bh_data = load_bh_curve(bh_file)
print(f"B-H curve: {len(bh_data)} points")

# Load ELF reference
elf_B = load_elf_field_at_origin(ELF_PATH)
if elf_B is not None:
    print(f"ELF reference Bz: {elf_B[2]*1000:.2f} mT")

print("\n" + "-" * 70)
results = []

# Test each solver
for solver_method in [0, 1, 2]:
    solver_name = SOLVER_NAMES[solver_method]
    print(f"\nTesting {solver_name} (Method {solver_method})...")

    try:
        result = run_solver(nodes, elements, bh_data, solver_method)

        if elf_B is not None:
            diff_Bz = abs(result['B'][2] - elf_B[2])
            rel_diff = abs(diff_Bz / elf_B[2]) * 100
            status = "PASS" if rel_diff < 1.0 else "WARN" if rel_diff < 5.0 else "FAIL"
        else:
            rel_diff = 0
            status = "N/A"

        print(f"  Bz = {result['B'][2]*1000:.2f} mT")
        print(f"  Time = {result['t_solve']:.2f} s")
        print(f"  Iterations = {result['iterations']}")
        print(f"  vs ELF: {rel_diff:.4f}% [{status}]")

        results.append({
            'solver': solver_name,
            'Bz_mT': result['B'][2] * 1000,
            't_solve': result['t_solve'],
            'iterations': result['iterations'],
            'rel_diff_pct': rel_diff,
            'status': status,
        })

    except Exception as e:
        print(f"  ERROR: {e}")
        results.append({
            'solver': solver_name,
            'status': 'ERROR',
        })

# Summary
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"\n{'Solver':<10} {'Bz(mT)':<12} {'Time(s)':<10} {'Iter':<8} {'Diff%':<10} {'Status'}")
print("-" * 60)

for r in results:
    if r['status'] != 'ERROR':
        print(f"{r['solver']:<10} {r['Bz_mT']:<12.2f} {r['t_solve']:<10.2f} "
              f"{r['iterations']:<8} {r['rel_diff_pct']:<10.4f} {r['status']}")
    else:
        print(f"{r['solver']:<10} ERROR")

pass_count = sum(1 for r in results if r['status'] == 'PASS')
print(f"\n{pass_count}/3 tests PASSED (<1% difference)")

if all(r['status'] == 'PASS' for r in results):
    print("\n*** ALL SOLVERS VERIFIED ***")
