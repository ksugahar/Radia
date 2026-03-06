#!/usr/bin/env python
"""
Verify ELF vs Radia for nonlinear electromagnet - BiCGSTAB solver.

IMPORTANT: Uses MIMA emulation (solve quarter, then 4x iron field)
to match ELF's MIMA behavior.

Model: C-type electromagnet with saturable iron
- 13 hexahedral elements (quarter yoke at x>=0, z>=0)
- Nonlinear B-H curve material
- 20000 AT full racetrack coil excitation

ELF Reference (with MIMA X + MIMA -Z): Bz = -994.16 mT at (0,0,0)
"""

import sys
import os
import time
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(work_dir)
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(parent_dir)))
sys.path.insert(0, os.path.join(repo_root, 'src'))
sys.path.insert(0, parent_dir)

import radia as rad
from coil_model import create_racetrack_coil

ELF_NONLINEAR = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\nonlinear_20000AT\ELF_MMB8T_EIEM2_1x1x1"

if not os.path.isdir(ELF_NONLINEAR):
    print(f"Skipping: ELF reference data not found at {ELF_NONLINEAR}")
    sys.exit(0)

scale = 0.001  # mm to m


def load_elf_geometry(path):
    """Load ELF yoke geometry from .meg file (MMB8T elements only)."""
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
                        return np.array([float(parts[5]), float(parts[6]), float(parts[7])])
    return None


def create_yoke_model(nodes, elements, mat):
    """Create yoke geometry from ELF nodes and elements."""
    hex_objects = []
    for elem_id, node_ids in elements:
        verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
                 for nid in node_ids]
        hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
        rad.MatApl(hex_obj, mat)
        hex_objects.append(hex_obj)
    return rad.ObjCnt(hex_objects), hex_objects


print("=" * 70)
print("ELF vs Radia Nonlinear Electromagnet - BiCGSTAB Solver")
print("=" * 70)
print("Solver:   BiCGSTAB (Method 1)")
print("Method:   MIMA emulation (solve quarter, then 4x iron field)")
print()

# Load geometry
nodes, elements = load_elf_geometry(ELF_NONLINEAR)
n_elem = len(elements)
n_dof = n_elem * 6
print(f"Geometry: {n_elem} hexahedral elements, {n_dof} DOF")

# Load B-H curve
bh_file = os.path.join(parent_dir, "BH.txt")
bh_data = load_bh_curve(bh_file)
print(f"Material: {len(bh_data)} B-H data points")

# Load ELF reference field
elf_B = load_elf_field_at_origin(ELF_NONLINEAR)
print(f"\nELF Reference at (0,0,0): Bz = {elf_B[2]*1000:.2f} mT")

# Build model
print("\n" + "-" * 70)
print("Building Radia Model")
print("-" * 70)

rad.UtiDelAll()

mat = rad.MatSatIsoTab(bh_data)
yoke, hex_objects = create_yoke_model(nodes, elements, mat)
coil = create_racetrack_coil(20000.0)
model = rad.ObjCnt([yoke, coil])

B_coil = np.array(rad.Fld(coil, 'b', [0, 0, 0]))
print(f"Coil-only Bz: {B_coil[2]*1000:.4f} mT")

# Solve with BiCGSTAB
print("\n" + "-" * 70)
print("Solving with BiCGSTAB (Method 1) - NO IMA")
print("-" * 70)

t_start = time.time()
result = rad.Solve(model, 0.0001, 1000, 1)  # Method 1 = BiCGSTAB
t_solve = time.time() - t_start

print(f"Iterations: {int(result[2])}")
print(f"Time: {t_solve:.3f} s")

# Compute field with MIMA emulation
B_total_quarter = np.array(rad.Fld(model, 'b', [0, 0, 0]))
B_iron_quarter = B_total_quarter[2] - B_coil[2]
B_mima_Bz = B_coil[2] + 4 * B_iron_quarter

print(f"\nMIMA emulation:")
print(f"  Quarter iron Bz: {B_iron_quarter*1000:.4f} mT")
print(f"  Total Bz (4x):   {B_mima_Bz*1000:.4f} mT")

# Compare with ELF
rel_error = abs(B_mima_Bz - elf_B[2]) / abs(elf_B[2]) * 100

print("\n" + "=" * 70)
print("COMPARISON WITH ELF")
print("=" * 70)
print(f"ELF Bz:     {elf_B[2]*1000:.2f} mT")
print(f"Radia Bz:   {B_mima_Bz*1000:.2f} mT")
print(f"Error:      {rel_error:.2f}%")

if rel_error < 5.0:
    print("\n*** PASS: Field within 5% of ELF ***")
elif rel_error < 10.0:
    print("\n*** ACCEPTABLE: Field within 10% ***")
else:
    print(f"\n*** NEEDS REVIEW: Field differs by {rel_error:.2f}% ***")
