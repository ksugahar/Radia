#!/usr/bin/env python
"""Test Newton-Raphson + Block Jacobi on 6x6x6 C-type electromagnet."""
import sys, os, time
import numpy as np

work_dir = r"S:\Radia\01_GitHub\examples\electromagnet\nonlinear"
repo_root = r"S:\Radia\01_GitHub"
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

def run_test(nodes, hex_elements, bh_data, solver_method, solver_name, use_newton=False, max_iter=100):
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

    # Enable/disable Newton
    rad.SetNewtonMethod(use_newton)

    t0 = time.time()
    try:
        result = rad.Solve(model, 0.001, max_iter, solver_method, image='+x-z')
        converged = True
    except RuntimeError as e:
        converged = False
        print(f"  WARNING: {e}")
    t_solve = time.time() - t0

    B = rad.Fld(model, 'b', [0, 0, 0])
    stats = rad.GetSolveStats()

    return {
        'name': solver_name,
        'Bz_mT': B[2] * 1000,
        't_solve': t_solve,
        'nonl_iter': stats.get('nonl_iterations', 'N/A') if stats else 'N/A',
        'linear_iter': stats.get('linear_iterations', 'N/A') if stats else 'N/A',
        'converged': converged,
    }

bh_file = os.path.join(work_dir, "BH.txt")
bh_data = load_bh_curve(bh_file)

# 6x6x6
ELF_6x6x6 = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\nonlinear_20000AT\ELF_MMB8T_EIEM2_6x6x6"
nodes, hex_elements = load_elf_geometry(ELF_6x6x6)
print(f"=== 6x6x6 (Elements: {len(hex_elements)}, DOF: {len(hex_elements)*6}) ===")
print()

# Test configurations
configs = [
    (2, 'HACApK+Picard (baseline)', False, 100),
    (2, 'HACApK+Newton', True, 200),
]

for method, name, newton, max_iter in configs:
    print(f"Running {name}...", flush=True)
    r = run_test(nodes, hex_elements, bh_data, method, name, use_newton=newton, max_iter=max_iter)
    status = "CONVERGED" if r['converged'] else "NOT CONVERGED"
    print(f"  {name}: Bz={r['Bz_mT']:.2f} mT, Time={r['t_solve']:.1f}s, "
          f"NL_iter={r['nonl_iter']}, Linear_iter={r['linear_iter']}, {status}")
    print()

# Restore default
rad.SetNewtonMethod(False)
print("Done.")
