#!/usr/bin/env python
"""
Verify ELF vs Radia for nonlinear electromagnet using LU solver (Method 0).

Model: C-type electromagnet with saturable iron
- 13 hexahedral elements (quarter model with MIMA X + MIMA -Z)
- Nonlinear B-H curve material
- 20000 AT racetrack coil

ELF Reference: Bz = -994.16 mT at (0,0,0)

LU Solver (Method 0):
- Dense direct solver using LAPACK dgesv
- Best for small problems (N < 500 elements)
- Guaranteed convergence for well-posed problems
"""

import sys
import os
import time
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
quater_dir = os.path.dirname(work_dir)               # .../quater
nonlinear_dir = os.path.dirname(quater_dir)          # .../nonlinear
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(nonlinear_dir)))  # repo root
sys.path.insert(0, os.path.join(repo_root, 'src'))
sys.path.insert(0, nonlinear_dir)  # For coil_model_quarter.py

import radia as rad
from coil_model_quarter import create_racetrack_coil_quarter

ELF_NONLINEAR = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\nonlinear_20000AT\ELF_MMB8T_EIEM2_1x1x1\quater"

if not os.path.isdir(ELF_NONLINEAR):
    print(f"Skipping: ELF reference data not found at {ELF_NONLINEAR}")
    sys.exit(0)

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
            if line.startswith('M3GB'):
                parts = line.split()
                if len(parts) >= 8:
                    point_id = int(parts[1])
                    if point_id == 1:
                        Bx = float(parts[5])
                        By = float(parts[6])
                        Bz = float(parts[7])
                        return np.array([Bx, By, Bz])
    return None


print("=" * 70)
print("ELF vs Radia Nonlinear Electromagnet - LU Solver")
print("=" * 70)
print("Solver:   LU (Method 0) - Dense direct solver")
print("Image:    Quarter model (+x-z)")
print("Material: Nonlinear B-H curve")
print("Current:  20000 AT")

# Load geometry
nodes, elements = load_elf_geometry(ELF_NONLINEAR)
n_elem = len(elements)
n_dof = n_elem * 6
print(f"\nGeometry: {n_elem} hexahedral elements, {n_dof} DOF")

# Load B-H curve
bh_file = os.path.join(nonlinear_dir, "BH.txt")
bh_data = load_bh_curve(bh_file)
print(f"Material: {len(bh_data)} B-H data points")

# Load ELF reference field
elf_B = load_elf_field_at_origin(ELF_NONLINEAR)
if elf_B is not None:
    print(f"\nELF reference at (0,0,0):")
    print(f"  Bz = {elf_B[2]*1000:.2f} mT")

# Create Radia model
print("\n" + "-" * 70)
print("Building Radia Model")
print("-" * 70)

rad.UtiDelAll()

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

# Create racetrack coil (20000 AT) - QUARTER model for IMA symmetry
print("Coil: Racetrack QUARTER (x >= 0), 20000 AT")
coil = create_racetrack_coil_quarter(20000.0)

# Combine model
model = rad.ObjCnt([yoke, coil])

# Solve with LU solver
print("\n" + "-" * 70)
print("Solving with LU Solver (Method 0)")
print("-" * 70)
print("  Precision: 0.0001")
print("  Max iterations: 100")
print("  Image: '+x-z' (quarter model)")

t_start = time.time()
result = rad.Solve(model, 0.0001, 100, 0, image='+x-z')
t_solve = time.time() - t_start

print(f"\nSolve Result:")
print(f"  Max |M|: {result[0]:.2f} A/m")
print(f"  Iterations: {int(result[2])}")
print(f"  Time: {t_solve:.3f} s")

# Compute field at origin
radia_B = np.array(rad.Fld(model, 'b', [0, 0, 0]))
print(f"\nRadia field at (0,0,0):")
print(f"  Bx = {radia_B[0]*1000:.4f} mT")
print(f"  By = {radia_B[1]*1000:.4f} mT")
print(f"  Bz = {radia_B[2]*1000:.4f} mT")
print(f"  |B| = {np.linalg.norm(radia_B)*1000:.4f} mT")

# Compare with ELF
if elf_B is not None:
    diff_Bz = abs(radia_B[2] - elf_B[2])
    rel_diff = abs(diff_Bz / elf_B[2]) * 100

    print("\n" + "=" * 70)
    print("COMPARISON")
    print("=" * 70)
    print(f"ELF Bz:    {elf_B[2]*1000:.2f} mT")
    print(f"Radia Bz:  {radia_B[2]*1000:.2f} mT")
    print(f"Difference: {diff_Bz*1000:.2f} mT ({rel_diff:.2f}%)")

    if rel_diff < 5.0:
        print("\n*** PASS: Field within 5% ***")
        status = "PASS"
    elif rel_diff < 10.0:
        print(f"\n*** ACCEPTABLE: Field within 10% ***")
        status = "ACCEPTABLE"
    else:
        print(f"\n*** NEEDS REVIEW: Field differs by {rel_diff:.2f}% ***")
        status = "NEEDS REVIEW"

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY - LU Solver (Method 0)")
    print("=" * 70)
    print(f"Elements:     {n_elem}")
    print(f"DOF:          {n_dof}")
    print(f"Iterations:   {int(result[2])}")
    print(f"Solve time:   {t_solve:.3f} s")
    print(f"ELF Bz:       {elf_B[2]*1000:.2f} mT")
    print(f"Radia Bz:     {radia_B[2]*1000:.2f} mT")
    print(f"Error:        {rel_diff:.2f}%")
    print(f"Status:       {status}")
