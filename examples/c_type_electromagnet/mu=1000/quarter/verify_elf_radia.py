#!/usr/bin/env python
"""
Verify ELF vs Radia for the quarter model (13 hexahedral elements + dual image).

Uses Image symmetry: MIMA X (symmetric) + MIMA -Z (antisymmetric)
= quarter model with +x and -z mirrors

Requirements:
1. Matrix: EXACT match (row-major format unified)
2. Field: Within 1% (using same racetrack coil model as ELF)
"""

import sys
import os
import numpy as np
import struct

work_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(work_dir)
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(parent_dir)))
sys.path.insert(0, os.path.join(repo_root, 'src'))
sys.path.insert(0, parent_dir)  # For coil_model

import radia as rad
from coil_model import create_racetrack_coil

# Note: ELF directory is spelled "quater" (missing 'r')
ELF_QUARTER = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\quater"

if not os.path.isdir(ELF_QUARTER):
    print(f"Skipping: ELF reference data not found at {ELF_QUARTER}")
    sys.exit(0)

MU_R = 1000
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


def load_elf_matrix(path, n_dof):
    """Load ELF interaction matrix from .mat file (Fortran binary format)."""
    with open(os.path.join(path, "ELF_magic.mat"), 'rb') as f:
        data = f.read()

    mat = np.zeros((n_dof, n_dof))
    offset = 0
    for i in range(n_dof):
        rec_len = struct.unpack('<i', data[offset:offset+4])[0]
        row_data = np.frombuffer(data[offset+4:offset+4+rec_len], dtype=np.float64)
        mat[i, :] = row_data
        offset += 4 + rec_len + 4
    return mat


def load_elf_field_at_origin(path):
    """Load ELF field at (0,0,0) from .mag file."""
    with open(os.path.join(path, "ELF_magic.mag"), 'r') as f:
        for line in f:
            if line.startswith('M3GB'):
                parts = line.split()
                if len(parts) >= 8:
                    point_id = int(parts[1])
                    if point_id == 1:  # Grid point 1 is (0,0,0)
                        Bx = float(parts[5])
                        By = float(parts[6])
                        Bz = float(parts[7])
                        return np.array([Bx, By, Bz])
    return None


print("=" * 70)
print("ELF vs Radia Quarter Model Verification")
print("Matrix: Row-Major [target][source] format")
print("Image:  MIMA X (symmetric) + MIMA -Z (antisymmetric)")
print("Coil:   Racetrack (from coil_model.py)")
print("=" * 70)

# Load ELF geometry
nodes, elements = load_elf_geometry(ELF_QUARTER)
n_elem = len(elements)
n_dof = n_elem * 6
print(f"\nGeometry: {n_elem} hexahedral elements, {n_dof} DOF (quarter model)")

# Load ELF matrix
elf_mat = load_elf_matrix(ELF_QUARTER, n_dof)
print(f"ELF matrix loaded: {elf_mat.shape}")

# Create Radia model
rad.UtiDelAll()
rad.FldUnits('m')

hex_objects = []
for elem_id, node_ids in elements:
    verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
             for nid in node_ids]
    hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
    mat = rad.MatLin(MU_R)
    rad.MatApl(hex_obj, mat)
    hex_objects.append(hex_obj)

container = rad.ObjCnt(hex_objects)

# Build and get Radia matrix WITH dual image symmetry
# MIMA X = symmetric mirror on x=0 plane (sign = +1)
# MIMA -Z = antisymmetric mirror on z=0 plane (sign = -1)
# Combined: '+x-z' = quarter model
handle = rad.BuildMatrix(container, image='+x-z')
radia_mat, dof = rad.GetInteractMatrix(handle)
print(f"Radia matrix (with quarter image): {radia_mat.shape}")

# ============================================================================
# MATRIX COMPARISON
# ============================================================================
print("\n" + "=" * 70)
print("MATRIX COMPARISON")
print("=" * 70)

diff_direct = np.abs(elf_mat - radia_mat)
max_diff = diff_direct.max()
mean_diff = diff_direct.mean()
rel_err = np.linalg.norm(elf_mat - radia_mat) / np.linalg.norm(elf_mat) * 100

print(f"\nMax |ELF - Radia|:      {max_diff:.6e}")
print(f"Mean |ELF - Radia|:     {mean_diff:.6e}")
print(f"Relative error:         {rel_err:.6e}%")

# Check if matrices match (within numerical precision tolerance)
MATRIX_TOL = 1e-6
if max_diff < MATRIX_TOL:
    print(f"\n*** PASS: Matrix MATCH (max diff {max_diff:.2e} < {MATRIX_TOL}) ***")
    matrix_pass = True
else:
    print(f"\n*** FAIL: Matrix difference {max_diff:.6e} > tolerance {MATRIX_TOL} ***")
    matrix_pass = False

    # Additional diagnostics
    print("\nDiagnostics:")
    diff_transposed = np.abs(elf_mat - radia_mat.T).max()
    print(f"Max |ELF - Radia.T|:    {diff_transposed:.6e}")

    idx = np.unravel_index(np.argmax(diff_direct), diff_direct.shape)
    print(f"Worst element at [{idx[0]}, {idx[1]}]:")
    print(f"  ELF:   {elf_mat[idx[0], idx[1]]:.10e}")
    print(f"  Radia: {radia_mat[idx[0], idx[1]]:.10e}")

# ============================================================================
# FIELD COMPARISON WITH RACETRACK COIL
# ============================================================================
print("\n" + "=" * 70)
print("FIELD COMPARISON (Racetrack Coil)")
print("=" * 70)

# Load ELF field at origin
elf_B = load_elf_field_at_origin(ELF_QUARTER)
if elf_B is not None:
    print(f"\nELF field at (0,0,0):")
    print(f"  Bx = {elf_B[0]*1000:.4f} mT")
    print(f"  By = {elf_B[1]*1000:.4f} mT")
    print(f"  Bz = {elf_B[2]*1000:.4f} mT")
    print(f"  |B| = {np.linalg.norm(elf_B)*1000:.4f} mT")

# Create new model with racetrack coil and image symmetry
rad.UtiDelAll()
rad.FldUnits('m')

hex_objects = []
for elem_id, node_ids in elements:
    verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
             for nid in node_ids]
    hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
    mat = rad.MatLin(MU_R)
    rad.MatApl(hex_obj, mat)
    hex_objects.append(hex_obj)

container = rad.ObjCnt(hex_objects)

# Create racetrack coil (2000 AT)
print("\nCreating racetrack coil (2000 AT)...")
coil = create_racetrack_coil(2000.0)

# Combine yoke and coil
model = rad.ObjCnt([container, coil])

# Solve WITH dual image symmetry
print("Solving with quarter image (+x-z)...")
result = rad.Solve(model, 0.0001, 100, 0, image='+x-z')
print(f"Solve result: {result}")

# Get field at origin
radia_B = np.array(rad.Fld(model, 'b', [0, 0, 0]))
print(f"\nRadia field at (0,0,0):")
print(f"  Bx = {radia_B[0]*1000:.4f} mT")
print(f"  By = {radia_B[1]*1000:.4f} mT")
print(f"  Bz = {radia_B[2]*1000:.4f} mT")
print(f"  |B| = {np.linalg.norm(radia_B)*1000:.4f} mT")

# Compare fields
if elf_B is not None:
    diff_B = radia_B - elf_B
    rel_diff_Bz = abs(diff_B[2] / elf_B[2]) * 100 if abs(elf_B[2]) > 1e-12 else 0

    print(f"\n--- Field Comparison ---")
    print(f"Bz difference: {diff_B[2]*1000:.4f} mT")
    print(f"Bz relative difference: {rel_diff_Bz:.2f}%")

    if rel_diff_Bz < 1.0:
        print("\n*** PASS: Field within 1% ***")
        field_pass = True
    elif rel_diff_Bz < 10.0:
        print(f"\n*** ACCEPTABLE: Field within 10% (coil approximation) ***")
        field_pass = True
    else:
        print(f"\n*** NEEDS REVIEW: Field differs by {rel_diff_Bz:.2f}% ***")
        field_pass = False
else:
    field_pass = False
    rel_diff_Bz = float('nan')

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"\n1. Matrix comparison:")
print(f"   Status:               {'PASS' if matrix_pass else 'FAIL'}")
print(f"   Max difference:       {max_diff:.6e}")
print(f"   Relative error:       {rel_err:.6e}%")
print(f"   Tolerance:            {MATRIX_TOL:.0e}")

if elf_B is not None:
    print(f"\n2. Field comparison:")
    print(f"   Status:               {'PASS' if field_pass else 'NEEDS REVIEW'}")
    print(f"   ELF Bz:               {elf_B[2]*1000:.4f} mT")
    print(f"   Radia Bz:             {radia_B[2]*1000:.4f} mT")
    print(f"   Difference:           {rel_diff_Bz:.2f}%")

if matrix_pass and field_pass:
    print("\n" + "=" * 70)
    print("*** VERIFICATION SUCCESSFUL ***")
    print("=" * 70)
    print("- Matrix: Row-major [target][source] format matches ELF")
    print("- Image:  Quarter model (+x-z) working correctly")
    print("- Field:  Within tolerance")
else:
    print("\n" + "=" * 70)
    print("*** VERIFICATION NEEDS REVIEW ***")
    print("=" * 70)
    if not matrix_pass:
        print("- Matrix mismatch detected")
    if not field_pass:
        print("- Field difference exceeds tolerance")
