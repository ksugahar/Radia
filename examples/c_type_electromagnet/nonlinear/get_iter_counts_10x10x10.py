#!/usr/bin/env python
"""Get nonlinear iteration counts for 10x10x10 C-type electromagnet."""
import sys, os, time
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.join(work_dir, '..', '..', '..')
sys.path.insert(0, os.path.join(repo_root, 'src'))
sys.path.insert(0, work_dir)

import radia as rad
from coil_model import create_racetrack_coil

scale = 0.001

def load_elf_geometry(path):
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
    data = np.loadtxt(filepath, comments='#')
    return data.tolist()

def run_and_get_stats(nodes, hex_elements, bh_data, solver_method, solver_name):
    rad.UtiDelAll()
    rad.FldUnits('m')
    mat = rad.MatSatIsoTab(bh_data)
    all_objects = []
    for elem_id, node_ids in hex_elements:
        verts = [[nodes[nid][0]*scale, nodes[nid][1]*scale, nodes[nid][2]*scale] for nid in node_ids]
        obj = rad.ObjHexahedron(verts, [0, 0, 0])
        rad.MatApl(obj, mat)
        all_objects.append(obj)
    yoke = rad.ObjCnt(all_objects)
    coil = create_racetrack_coil(20000.0)
    model = rad.ObjCnt([yoke, coil])
    if solver_method == 2:
        rad.SetHACApKParams(1e-4, 10, 2.0)
    t0 = time.time()
    result = rad.Solve(model, 0.001, 100, solver_method, image='+x-z')
    t_solve = time.time() - t0
    B = rad.Fld(model, 'b', [0, 0, 0])
    stats = rad.GetSolveStats()
    return {
        'name': solver_name,
        'Bz_mT': B[2] * 1000,
        't_solve': t_solve,
        'nonl_iter': stats.get('nonl_iterations', 'N/A') if stats else 'N/A',
        'linear_iter': stats.get('linear_iterations', 'N/A') if stats else 'N/A',
    }

bh_file = os.path.join(work_dir, "BH.txt")
bh_data = load_bh_curve(bh_file)

# 10x10x10
ELF_10x10x10 = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\nonlinear_20000AT\ELF_MMB8T_EIEM2_10x10x10"
nodes, hex_elements = load_elf_geometry(ELF_10x10x10)
print(f"=== 10x10x10 (Elements: {len(hex_elements)}, DOF: {len(hex_elements)*6}) ===")

for method, name in [(0, 'LU'), (1, 'BiCGSTAB'), (2, 'HACApK')]:
    print(f"Running {name}...", flush=True)
    r = run_and_get_stats(nodes, hex_elements, bh_data, method, name)
    print(f"  {name}: Bz={r['Bz_mT']:.2f} mT, Time={r['t_solve']:.1f}s, "
          f"NL_iter={r['nonl_iter']}, Linear_iter={r['linear_iter']}")
