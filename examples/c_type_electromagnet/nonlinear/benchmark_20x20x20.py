#!/usr/bin/env python
"""
Benchmark: 20x20x20 C-type electromagnet with LU, BiCGSTAB, and HACApK solvers.

Compares solvers on EIEM2 20x20x20 mesh (27600 elements, 165600 DOF).
Dense matrix requires ~209 GB, so LU and BiCGSTAB will fail on 128 GB RAM.
Only HACApK with H-matrix compression can handle this scale.
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

# ELF 20x20x20 mesh
ELF_20x20x20 = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\nonlinear_20000AT\ELF_MMB8T_EIEM2_20x20x20"

scale = 0.001  # mm to m


def load_elf_geometry(path):
    """Load ELF geometry from .meg file."""
    nodes = {}
    hex_elements = []
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
    return nodes, hex_elements


def load_bh_curve(filepath):
    """Load B-H curve from text file."""
    data = np.loadtxt(filepath, comments='#')
    return data.tolist()


def run_solver(nodes, hex_elements, bh_data, solver_method, solver_name,
               use_ima=True):
    """Run solver and return results with timing."""
    rad.UtiDelAll()
    rad.FldUnits('m')

    mat = rad.MatSatIsoTab(bh_data)
    all_objects = []

    for elem_id, node_ids in hex_elements:
        verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale,
                  nodes[nid][2] * scale] for nid in node_ids]
        obj = rad.ObjHexahedron(verts, [0, 0, 0])
        rad.MatApl(obj, mat)
        all_objects.append(obj)

    yoke = rad.ObjCnt(all_objects)
    coil = create_racetrack_coil(20000.0)
    model = rad.ObjCnt([yoke, coil])

    if solver_method == 2:
        rad.SetHACApKParams(1e-4, 10, 2.0)

    ima = '+x-z' if use_ima else ''

    t_start = time.time()
    if ima:
        result = rad.Solve(model, 0.001, 100, solver_method, image=ima)
    else:
        result = rad.Solve(model, 0.001, 100, solver_method)
    t_solve = time.time() - t_start

    # Get H-matrix stats if HACApK was used
    hacapk_stats = None
    if solver_method == 2:
        try:
            hacapk_stats = rad.GetHACApKStats()
        except Exception:
            pass

    B = rad.Fld(model, 'b', [0, 0, 0])

    # Dense matrix memory estimate: N^2 * 8 bytes (double)
    dof = len(all_objects) * 6
    dense_mem_MB = dof * dof * 8 / (1024 * 1024)

    return {
        'name': solver_name,
        'method': solver_method,
        'n_elements': len(all_objects),
        'dof': dof,
        't_solve': t_solve,
        'dense_mem_MB': dense_mem_MB,
        'Bz': B[2],
        'Bz_mT': B[2] * 1000,
        'max_M': result[0],
        'rel_err': result[1],
        'hacapk_stats': hacapk_stats,
    }


# ============================================================
# Main
# ============================================================
print("=" * 70)
print("Benchmark: 20x20x20 C-type Electromagnet (Nonlinear 20000 AT)")
print("=" * 70)

# Load geometry
nodes, hex_elements = load_elf_geometry(ELF_20x20x20)
print(f"Mesh: 20x20x20 EIEM2")
print(f"Elements: {len(hex_elements)}")
dof = len(hex_elements) * 6
print(f"DOF: {dof}")

# Load B-H curve
bh_file = os.path.join(work_dir, "BH.txt")
bh_data = load_bh_curve(bh_file)
print(f"B-H curve: {len(bh_data)} points")

# Dense matrix memory for reference
dense_mem = dof * dof * 8 / (1024 * 1024)
dense_mem_GB = dense_mem / 1024
print(f"Dense matrix memory: {dense_mem:.0f} MB ({dense_mem_GB:.1f} GB)")
print(f"  ({dof}x{dof} doubles)")
print(f"  NOTE: Exceeds 128 GB RAM - LU and BiCGSTAB will fail")

# Run HACApK first (the one that will work), then try LU and BiCGSTAB
# HACApK: IMA now works via kernel-level mirror contributions (no dense IMA matrix)
# LU/BiCGSTAB: will fail because dense matrix exceeds RAM
use_ima = True
print(f"Image symmetry: {'enabled' if use_ima else 'DISABLED'}")

results = []
for method, name in [(2, 'HACApK'), (0, 'LU'), (1, 'BiCGSTAB')]:
    print(f"\n{'='*50}")
    print(f"Solver: {name} (Method {method})")
    print(f"{'='*50}")

    if method in (0, 1):
        print(f"  Dense matrix requires {dense_mem_GB:.1f} GB > 128 GB RAM")
        print(f"  Attempting anyway to demonstrate failure...")

    try:
        r = run_solver(nodes, hex_elements, bh_data, method, name,
                       use_ima=use_ima)
        results.append(r)

        print(f"  Elements: {r['n_elements']}, DOF: {r['dof']}")
        print(f"  Bz = {r['Bz_mT']:.2f} mT")
        print(f"  Time: {r['t_solve']:.2f} s")
        print(f"  Dense matrix memory: {r['dense_mem_MB']:.0f} MB")

        if r['hacapk_stats']:
            s = r['hacapk_stats']
            print(f"\n  H-matrix statistics:")
            print(f"    DOF: {s.get('n_dof', 'N/A')}")
            print(f"    Leaves: {s.get('n_leaves', 'N/A')}")
            print(f"    Low-rank blocks: {s.get('n_lowrank', 'N/A')}")
            print(f"    Dense blocks: {s.get('n_dense', 'N/A')}")
            print(f"    Max rank: {s.get('max_rank', 'N/A')}")
            print(f"    Compression ratio: {s.get('compression', 'N/A')}")
            print(f"    H-matrix memory: {s.get('memory_mb', 'N/A')} MB")
            print(f"    Dense memory: {s.get('dense_memory_mb', 'N/A')} MB")
            print(f"    Build time: {s.get('build_time', 'N/A')} s")
            print(f"    Linear iterations: {s.get('linear_iterations', 'N/A')}")
    except Exception as e:
        print(f"  FAILED: {e}")
        results.append({
            'name': name, 'method': method,
            'Bz_mT': 0, 't_solve': 0, 'dense_mem_MB': dense_mem,
            'hacapk_stats': None, 'dof': dof,
            'failed': True, 'error': str(e),
        })

# Summary
print(f"\n{'='*70}")
print(f"SUMMARY: 20x20x20 C-type Electromagnet ({len(hex_elements)} elem, {dof} DOF)")
print(f"{'='*70}")

print(f"\n{'Solver':<12} {'Bz(mT)':<12} {'Time(s)':<12} {'Memory'}")
print("-" * 56)
for r in results:
    if r.get('failed'):
        print(f"{r['name']:<12} {'FAILED':<12} {'--':<12} "
              f"Dense={r['dense_mem_MB']/1024:.1f} GB (>128 GB)")
    else:
        mem_str = ""
        if r['hacapk_stats'] and r['hacapk_stats'].get('memory_mb', 0) > 0:
            hm = r['hacapk_stats']['memory_mb']
            mem_str = f"H-mat={hm:.0f} MB"
        else:
            mem_str = f"Dense={r['dense_mem_MB']/1024:.1f} GB"
        print(f"{r['name']:<12} {r['Bz_mT']:<12.2f} {r['t_solve']:<12.1f} {mem_str}")

# HACApK compression summary
hacapk_r = [r for r in results if r['name'] == 'HACApK' and not r.get('failed')]
if hacapk_r:
    r = hacapk_r[0]
    s = r.get('hacapk_stats', {})
    if s and s.get('memory_mb', 0) > 0:
        hm = s['memory_mb']
        dm = s.get('dense_memory_mb', dense_mem)
        print(f"\nH-matrix compression: {hm:.0f} MB / {dm:.0f} MB = {hm/dm*100:.1f}%")
        print(f"Memory reduction: {(1-hm/dm)*100:.0f}%")
