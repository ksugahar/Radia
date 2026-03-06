#!/usr/bin/env python
"""Test: HACApK with IMA on 20x20x20 C-type mesh with extended iterations."""
import sys, os, time
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(work_dir)))
sys.path.insert(0, os.path.join(repo_root, 'src'))
sys.path.insert(0, work_dir)

import radia as rad
from coil_model import create_racetrack_coil

ELF_20x20x20 = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\nonlinear_20000AT\ELF_MMB8T_EIEM2_20x20x20"
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
                node_ids = [int(parts[i]) for i in range(4, 12)]
                hex_elements.append((int(parts[1]), node_ids))
    return nodes, hex_elements

nodes, hex_elements = load_elf_geometry(ELF_20x20x20)
bh_data = np.loadtxt(os.path.join(work_dir, "BH.txt"), comments='#').tolist()

print(f"20x20x20: {len(hex_elements)} elements, {len(hex_elements)*6} DOF")

rad.UtiDelAll()
rad.FldUnits('m')
mat = rad.MatSatIsoTab(bh_data)
all_objects = []
for elem_id, node_ids in hex_elements:
    verts = [[nodes[nid][0]*scale, nodes[nid][1]*scale, nodes[nid][2]*scale]
             for nid in node_ids]
    obj = rad.ObjHexahedron(verts, [0, 0, 0])
    rad.MatApl(obj, mat)
    all_objects.append(obj)

yoke = rad.ObjCnt(all_objects)
coil = create_racetrack_coil(20000.0)
model = rad.ObjCnt([yoke, coil])

rad.SetHACApKParams(1e-4, 10, 2.0)

# Run 200 iterations for convergence
print(f"Solving with HACApK + IMA (+x-z), maxiter=200, prec=0.001...")
t0 = time.time()
try:
    result = rad.Solve(model, 0.001, 200, 2, image='+x-z')
    t1 = time.time()
    B = rad.Fld(model, 'b', [0, 0, 0])
    stats = rad.GetHACApKStats()

    print(f"\nSUCCESS:")
    print(f"  Bz = {B[2]*1000:.2f} mT")
    print(f"  Time = {t1-t0:.1f} s")
    print(f"  Result: {result}")
    if stats:
        print(f"  DOF: {stats.get('n_dof', 'N/A')}")
        print(f"  H-matrix memory: {stats.get('memory_mb', 'N/A'):.1f} MB")
        print(f"  Compression: {stats.get('compression', 'N/A')}")
        print(f"  Build time: {stats.get('build_time', 'N/A'):.1f} s")
        print(f"  Linear iterations: {stats.get('linear_iterations', 'N/A')}")
except Exception as e:
    t1 = time.time()
    print(f"\nFAILED after {t1-t0:.1f} s: {e}")
    # Try to get field anyway
    try:
        B = rad.Fld(model, 'b', [0, 0, 0])
        print(f"  Bz (partial) = {B[2]*1000:.2f} mT")
    except:
        pass
    try:
        stats = rad.GetHACApKStats()
        if stats:
            print(f"  H-matrix memory: {stats.get('memory_mb', 'N/A'):.1f} MB")
    except:
        pass
