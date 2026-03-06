#!/usr/bin/env python
"""
Test quarter model with NONLINEAR material but WITHOUT IMA.
This isolates whether the issue is with nonlinear material or with IMA.
"""

import sys
import os
import time
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
nonlinear_dir = os.path.dirname(work_dir)
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(nonlinear_dir)))
sys.path.insert(0, os.path.join(repo_root, 'src'))
sys.path.insert(0, nonlinear_dir)

import radia as rad

ELF_QUARTER = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\nonlinear_20000AT\ELF_MMB8T_EIEM2_1x1x1\quater"
scale = 0.001

# Import the standard coil model
from coil_model import create_racetrack_coil


def load_elf_geometry(path):
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
    data = np.loadtxt(filepath, comments='#')
    return data.tolist()


def test_nonlinear(nodes, elements, bh_data, use_ima=False):
    """Test quarter model with nonlinear material."""
    rad.UtiDelAll()

    # NONLINEAR material
    mat = rad.MatSatIsoTab(bh_data)

    hex_objects = []
    for elem_id, node_ids in elements:
        verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
                 for nid in node_ids]
        hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
        rad.MatApl(hex_obj, mat)
        hex_objects.append(hex_obj)

    yoke = rad.ObjCnt(hex_objects)
    coil = create_racetrack_coil(20000.0)
    model = rad.ObjCnt([yoke, coil])

    t_start = time.time()
    if use_ima:
        result = rad.Solve(model, 0.001, 100, 0, image='+x-z')
    else:
        result = rad.Solve(model, 0.001, 100, 0)
    t_solve = time.time() - t_start

    B = np.array(rad.Fld(model, 'b', [0, 0, 0]))

    # Correct result interpretation:
    # result[0] = MisfitM (convergence measure)
    # result[1] = MaxModM (Max |M|)
    # result[2] = MaxModH
    # result[3] = ActualIterNum
    return {
        'misfit_M': result[0],
        'max_M': result[1],
        'max_H': result[2],
        'iterations': int(result[3]),
        'time': t_solve,
        'Bz': B[2],
    }


print("=" * 70)
print("Quarter Model with NONLINEAR Material")
print("=" * 70)

nodes, elements = load_elf_geometry(ELF_QUARTER)
bh_file = os.path.join(nonlinear_dir, "BH.txt")
bh_data = load_bh_curve(bh_file)

print(f"Geometry: {len(elements)} elements, {len(elements) * 6} DOF")
print(f"Material: {len(bh_data)} B-H points")
print()

# Test WITHOUT IMA
print("-" * 70)
print("Test 1: Quarter model WITHOUT IMA (quarter geometry, no mirrors)")
print("-" * 70)
result_no_ima = test_nonlinear(nodes, elements, bh_data, use_ima=False)
print(f"  Misfit M:   {result_no_ima['misfit_M']:.6f}")
print(f"  Max |M|:    {result_no_ima['max_M']:.2f} A/m")
print(f"  Max |H|:    {result_no_ima['max_H']:.2f} A/m")
print(f"  Iterations: {result_no_ima['iterations']}")
print(f"  Time:       {result_no_ima['time']:.3f} s")
print(f"  Bz:         {result_no_ima['Bz']*1000:.2f} mT")
print()

# Test WITH IMA
print("-" * 70)
print("Test 2: Quarter model WITH IMA (+x-z)")
print("-" * 70)
result_ima = test_nonlinear(nodes, elements, bh_data, use_ima=True)
print(f"  Misfit M:   {result_ima['misfit_M']:.6f}")
print(f"  Max |M|:    {result_ima['max_M']:.2f} A/m")
print(f"  Max |H|:    {result_ima['max_H']:.2f} A/m")
print(f"  Iterations: {result_ima['iterations']}")
print(f"  Time:       {result_ima['time']:.3f} s")
print(f"  Bz:         {result_ima['Bz']*1000:.2f} mT")
print()

# Compare
print("=" * 70)
print("COMPARISON (Nonlinear Material)")
print("=" * 70)
print(f"{'Test':<25} {'Max |M| (A/m)':<15} {'Iter':<8} {'Bz (mT)':<12}")
print("-" * 60)
print(f"{'No IMA (quarter only)':<25} {result_no_ima['max_M']:<15.2f} {result_no_ima['iterations']:<8} {result_no_ima['Bz']*1000:<12.2f}")
print(f"{'IMA +x-z':<25} {result_ima['max_M']:<15.2f} {result_ima['iterations']:<8} {result_ima['Bz']*1000:<12.2f}")
