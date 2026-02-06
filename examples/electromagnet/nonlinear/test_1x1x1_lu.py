#!/usr/bin/env python
"""
Simple test: 1x1x1 EIEM2 mesh with LU solver.
Verifies basic quarter model + IMA (+x-z) for nonlinear materials.
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


print("=" * 60)
print("Test: EIEM2_1x1x1 with LU Solver + IMA (+x-z)")
print("=" * 60)

# Load geometry
nodes, elements = load_elf_geometry(ELF_PATH)
print(f"Elements: {len(elements)}")
print(f"DOF: {len(elements) * 6}")

# Load B-H curve
bh_file = os.path.join(work_dir, "BH.txt")
bh_data = load_bh_curve(bh_file)
print(f"B-H curve: {len(bh_data)} points")

# Load ELF reference
elf_B = load_elf_field_at_origin(ELF_PATH)
if elf_B is not None:
    print(f"\nELF reference at (0,0,0):")
    print(f"  Bz = {elf_B[2]*1000:.2f} mT")

# Build Radia model
print("\nBuilding Radia model...")
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
print(f"Yoke: {len(hex_objects)} elements")

# Create coil
coil = create_racetrack_coil(20000.0)
print("Coil: Racetrack 20000 AT")

# Combine model
model = rad.ObjCnt([yoke, coil])

# Solve with LU
print("\nSolving with LU (Method 0)...")
print("  Image: +x-z")
print("  Precision: 0.0001")

t_start = time.time()
result = rad.Solve(model, 0.0001, 100, 0, image='+x-z')
t_solve = time.time() - t_start

print(f"\nSolve completed:")
print(f"  Max |M|: {result[0]:.0f} A/m")
print(f"  Iterations: {int(result[2])}")
print(f"  Time: {t_solve:.2f} s")

# Get field at origin
radia_B = np.array(rad.Fld(model, 'b', [0, 0, 0]))
print(f"\nRadia field at (0,0,0):")
print(f"  Bx = {radia_B[0]*1000:.4f} mT")
print(f"  By = {radia_B[1]*1000:.4f} mT")
print(f"  Bz = {radia_B[2]*1000:.4f} mT")

# Compare
if elf_B is not None:
    diff_Bz = abs(radia_B[2] - elf_B[2])
    rel_diff = abs(diff_Bz / elf_B[2]) * 100

    print("\n" + "=" * 60)
    print("COMPARISON")
    print("=" * 60)
    print(f"ELF Bz:   {elf_B[2]*1000:.2f} mT")
    print(f"Radia Bz: {radia_B[2]*1000:.2f} mT")
    print(f"Diff:     {diff_Bz*1000:.2f} mT ({rel_diff:.4f}%)")

    if rel_diff < 1.0:
        print("\n*** PASS (<1%) ***")
    elif rel_diff < 5.0:
        print(f"\n*** WARN (1-5%) ***")
    else:
        print(f"\n*** FAIL (>5%) ***")
