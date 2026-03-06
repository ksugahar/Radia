#!/usr/bin/env python
"""
Verify ELF vs Radia for nonlinear electromagnet - LU solver.

IMPORTANT: ELF's MIMA works differently than Radia's IMA!
- ELF MIMA: Solve quarter model, then 4x field contribution
- Radia IMA: Modify interaction matrix (gives different M values)

This script uses the MIMA-emulation approach for proper ELF comparison:
1. Solve quarter model without IMA
2. Apply 4x factor to iron field manually

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
nonlinear_dir = os.path.dirname(parent_dir)
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(nonlinear_dir)))
sys.path.insert(0, os.path.join(repo_root, 'src'))
sys.path.insert(0, nonlinear_dir)

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
                        Bx = float(parts[5])
                        By = float(parts[6])
                        Bz = float(parts[7])
                        return np.array([Bx, By, Bz])
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


def get_magnetization_stats(hex_objects):
    """Get magnetization statistics from hex objects."""
    M_values = []
    for hex_obj in hex_objects:
        M_dict = rad.ObjM(hex_obj)
        if 'magnetization' in M_dict:
            Mx, My, Mz = M_dict['magnetization']
            M_mag = np.sqrt(Mx**2 + My**2 + Mz**2)
            M_values.append((Mx, My, Mz, M_mag))
    return M_values


print("=" * 70)
print("ELF vs Radia Nonlinear Electromagnet - LU Solver")
print("=" * 70)
print("Method: MIMA emulation (solve quarter, then 4x iron field)")
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
print(f"  H range: {bh_data[0][0]:.1f} to {bh_data[-1][0]:.1f} A/m")
print(f"  B range: {bh_data[0][1]:.4f} to {bh_data[-1][1]:.4f} T")

# Load ELF reference field
elf_B = load_elf_field_at_origin(ELF_NONLINEAR)
print(f"\nELF Reference at (0,0,0) [with MIMA X + MIMA -Z]:")
print(f"  Bz = {elf_B[2]*1000:.2f} mT")

# ============================================================================
# Build and solve model
# ============================================================================
print("\n" + "-" * 70)
print("Building Radia Model")
print("-" * 70)

rad.UtiDelAll()
rad.FldUnits('m')

mat = rad.MatSatIsoTab(bh_data)
yoke, hex_objects = create_yoke_model(nodes, elements, mat)
print(f"Yoke: {len(hex_objects)} hexahedral elements")

# Full coil (matches ELF's MIMA-mirrored MCL8T elements)
coil = create_racetrack_coil(20000.0)
print("Coil: Full racetrack (20000 AT)")

# Check coil-only field
B_coil = np.array(rad.Fld(coil, 'b', [0, 0, 0]))
print(f"\nCoil-only field at (0,0,0):")
print(f"  Bz = {B_coil[2]*1000:.4f} mT")

model = rad.ObjCnt([yoke, coil])

# ============================================================================
# Solve WITHOUT IMA (MIMA emulation)
# ============================================================================
print("\n" + "-" * 70)
print("Solving with LU (Method 0) - NO IMA")
print("-" * 70)
print("  MIMA emulation: Solve quarter model, then apply 4x factor")

t_start = time.time()
result = rad.Solve(model, 0.0001, 100, 0)  # NO IMA
t_solve = time.time() - t_start

print(f"\nSolve result:")
print(f"  Max|M|: {result[0]:.2f} A/m (from solver)")
print(f"  Avg|M|: {result[1]:.2f} A/m (from solver)")
print(f"  Iterations: {int(result[2])}")
print(f"  Time: {t_solve:.3f} s")

# Get actual magnetization stats
M_stats = get_magnetization_stats(hex_objects)
if M_stats:
    M_mags = [m[3] for m in M_stats]
    print(f"\nElement magnetization (actual):")
    print(f"  Max|M|: {max(M_mags):.1f} A/m")
    print(f"  Min|M|: {min(M_mags):.1f} A/m")
    print(f"  Avg|M|: {np.mean(M_mags):.1f} A/m")

# ============================================================================
# Compute field with MIMA emulation
# ============================================================================
print("\n" + "-" * 70)
print("Field at Origin (MIMA emulation)")
print("-" * 70)

B_total_quarter = np.array(rad.Fld(model, 'b', [0, 0, 0]))
B_iron_quarter = B_total_quarter[2] - B_coil[2]

# Apply MIMA factor: Coil + 4x Iron (for +x-z symmetry)
B_mima_Bz = B_coil[2] + 4 * B_iron_quarter

print(f"Quarter model field:")
print(f"  Coil Bz:       {B_coil[2]*1000:.4f} mT")
print(f"  Iron Bz (1x):  {B_iron_quarter*1000:.4f} mT")
print()
print(f"MIMA emulation (4x iron):")
print(f"  Iron Bz (4x):  {4*B_iron_quarter*1000:.4f} mT")
print(f"  Total Bz:      {B_mima_Bz*1000:.4f} mT")

# ============================================================================
# Comparison with ELF
# ============================================================================
print("\n" + "=" * 70)
print("COMPARISON WITH ELF")
print("=" * 70)

diff_Bz = B_mima_Bz - elf_B[2]
rel_error = abs(diff_Bz / elf_B[2]) * 100

print(f"ELF Bz:     {elf_B[2]*1000:.2f} mT")
print(f"Radia Bz:   {B_mima_Bz*1000:.2f} mT")
print(f"Difference: {diff_Bz*1000:.2f} mT")
print(f"Error:      {rel_error:.2f}%")

if rel_error < 5.0:
    status = "PASS"
    print(f"\n*** PASS: Field within 5% of ELF ***")
elif rel_error < 10.0:
    status = "ACCEPTABLE"
    print(f"\n*** ACCEPTABLE: Field within 10% ***")
else:
    status = "NEEDS REVIEW"
    print(f"\n*** NEEDS REVIEW: Field differs by {rel_error:.2f}% ***")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 70)
print("SUMMARY - LU Solver with MIMA Emulation")
print("=" * 70)
print(f"Elements:     {n_elem}")
print(f"DOF:          {n_dof}")
print(f"Iterations:   {int(result[2])}")
print(f"Solve time:   {t_solve:.3f} s")
print(f"ELF Bz:       {elf_B[2]*1000:.2f} mT")
print(f"Radia Bz:     {B_mima_Bz*1000:.2f} mT")
print(f"Error:        {rel_error:.2f}%")
print(f"Status:       {status}")

print("\nNote: MIMA emulation = solve quarter model + 4x iron field")
print("      This matches ELF's MIMA behavior better than Radia's IMA.")
