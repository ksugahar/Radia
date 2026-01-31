#!/usr/bin/env python
"""
Compare field results between Radia and reference values.
This is the ultimate test - if the fields match, the Image API is working correctly.
"""

import sys
import os
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(work_dir))))
sys.path.insert(0, os.path.join(repo_root, 'src'))

import radia as rad

# === ELF paths ===
ELF_QUARTER_PATH = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\quater"
ELF_FULL_PATH = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\full"

MU_R = 1000
scale = 0.001  # mm to m

print("=" * 70)
print("Field Comparison: Radia Quarter Model vs Full Model")
print("=" * 70)

# === Load quarter model geometry ===
nodes_quarter = {}
elements_quarter = []
with open(os.path.join(ELF_QUARTER_PATH, "ELF_magic.meg"), 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('MGR1'):
            parts = line.split()
            node_id = int(parts[1])
            x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
            nodes_quarter[node_id] = np.array([x, y, z])
        elif line.startswith('MMB8T'):
            parts = line.split()
            elem_id = int(parts[1])
            node_ids = [int(parts[i]) for i in range(4, 12)]
            elements_quarter.append((elem_id, node_ids))

# === Load full model geometry ===
nodes_full = {}
elements_full = []
with open(os.path.join(ELF_FULL_PATH, "ELF_magic.meg"), 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('MGR1'):
            parts = line.split()
            node_id = int(parts[1])
            x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
            nodes_full[node_id] = np.array([x, y, z])
        elif line.startswith('MMB8T'):
            parts = line.split()
            elem_id = int(parts[1])
            node_ids = [int(parts[i]) for i in range(4, 12)]
            elements_full.append((elem_id, node_ids))

print(f"Quarter model: {len(elements_quarter)} elements")
print(f"Full model: {len(elements_full)} elements")

# === Create Radia full model ===
def create_full_model():
    rad.UtiDelAll()
    rad.FldUnits('m')

    hex_objects = []
    for elem_id, node_ids in elements_full:
        verts = [[nodes_full[nid][0] * scale, nodes_full[nid][1] * scale, nodes_full[nid][2] * scale]
                 for nid in node_ids]
        hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
        mat = rad.MatLin(MU_R)
        rad.MatApl(hex_obj, mat)
        hex_objects.append(hex_obj)

    container = rad.ObjCnt(hex_objects)

    # Background field: B_z = 0.1 T
    bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])
    model = rad.ObjCnt([container, bkg])

    return model

# === Create Radia quarter model with Image ===
def create_quarter_model_with_image():
    rad.UtiDelAll()
    rad.FldUnits('m')

    hex_objects = []
    for elem_id, node_ids in elements_quarter:
        verts = [[nodes_quarter[nid][0] * scale, nodes_quarter[nid][1] * scale, nodes_quarter[nid][2] * scale]
                 for nid in node_ids]
        hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
        mat = rad.MatLin(MU_R)
        rad.MatApl(hex_obj, mat)
        hex_objects.append(hex_obj)

    container = rad.ObjCnt(hex_objects)

    # Background field
    bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])
    model = rad.ObjCnt([container, bkg])

    return model

# === Test 1: Full model (reference) ===
print("\n--- Test 1: Full Model (52 elements) ---")
model_full = create_full_model()
rad.Solve(model_full, 0.0001, 100, 0)

# Field at various points
test_points = [
    [0, 0, 0],
    [0, 0.05, 0],
    [0, 0.10, 0],
    [0.05, 0, 0],
    [0, 0, 0.05],
]

print("\nFields at test points (Full model):")
B_full = []
for pt in test_points:
    B = rad.Fld(model_full, 'b', pt)
    Bz_mT = B[2] * 1000
    B_full.append(Bz_mT)
    print(f"  Point {pt}: Bz = {Bz_mT:.4f} mT")

# === Test 2: Quarter model with Image ===
print("\n--- Test 2: Quarter Model with Image (13 elements) ---")
model_quarter = create_quarter_model_with_image()
# New API: Solve with image parameter for IMA symmetry
rad.Solve(model_quarter, 0.0001, 100, 0, image='+x-z')

print("\nFields at test points (Quarter model with Image):")
B_quarter = []
for pt in test_points:
    B = rad.Fld(model_quarter, 'b', pt)
    Bz_mT = B[2] * 1000
    B_quarter.append(Bz_mT)
    print(f"  Point {pt}: Bz = {Bz_mT:.4f} mT")

# === Test 3: Quarter model WITHOUT Image ===
print("\n--- Test 3: Quarter Model WITHOUT Image (13 elements) ---")
rad.UtiDelAll()
rad.FldUnits('m')

hex_objects = []
for elem_id, node_ids in elements_quarter:
    verts = [[nodes_quarter[nid][0] * scale, nodes_quarter[nid][1] * scale, nodes_quarter[nid][2] * scale]
             for nid in node_ids]
    hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
    mat = rad.MatLin(MU_R)
    rad.MatApl(hex_obj, mat)
    hex_objects.append(hex_obj)

container = rad.ObjCnt(hex_objects)
bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])
model_no_image = rad.ObjCnt([container, bkg])
rad.Solve(model_no_image, 0.0001, 100, 0)

print("\nFields at test points (Quarter model, NO Image):")
B_no_image = []
for pt in test_points:
    B = rad.Fld(model_no_image, 'b', pt)
    Bz_mT = B[2] * 1000
    B_no_image.append(Bz_mT)
    print(f"  Point {pt}: Bz = {Bz_mT:.4f} mT")

# === Summary ===
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"\n{'Point':<20} {'Full (ref)':<15} {'Quarter+Image':<15} {'No Image':<15} {'Diff (%)':<12}")
print("-" * 70)
for i, pt in enumerate(test_points):
    diff_pct = abs(B_quarter[i] - B_full[i]) / abs(B_full[i]) * 100 if B_full[i] != 0 else 0
    print(f"{str(pt):<20} {B_full[i]:<15.4f} {B_quarter[i]:<15.4f} {B_no_image[i]:<15.4f} {diff_pct:<12.4f}")

print("\n" + "=" * 70)

# Check if Image gives significant improvement
max_diff_with_image = max(abs(B_quarter[i] - B_full[i]) / abs(B_full[i]) * 100 for i in range(len(test_points)) if B_full[i] != 0)
max_diff_no_image = max(abs(B_no_image[i] - B_full[i]) / abs(B_full[i]) * 100 for i in range(len(test_points)) if B_full[i] != 0)

print(f"Max error with Image:    {max_diff_with_image:.4f}%")
print(f"Max error without Image: {max_diff_no_image:.4f}%")

if max_diff_with_image < 0.1:
    print("\nPASS: Image API gives accurate field results")
elif max_diff_with_image < max_diff_no_image:
    print(f"\nPARTIAL: Image improves accuracy but still {max_diff_with_image:.2f}% error")
else:
    print("\nFAIL: Image API does not improve accuracy")
