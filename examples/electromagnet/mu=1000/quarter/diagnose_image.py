#!/usr/bin/env python
"""
Diagnose Image API: Check if magnetization values change with Image.

Uses the new Solve(image=...) API.
"""

import sys
import os
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(work_dir))))
sys.path.insert(0, os.path.join(repo_root, 'src'))

import radia as rad

ELF_QUARTER_PATH = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\quater"
ELF_FULL_PATH = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\full"

MU_R = 1000
scale = 0.001

print("=" * 70)
print("Diagnosing Image API (New Solve API)")
print("=" * 70)

# Load quarter geometry
nodes = {}
elements = []
with open(os.path.join(ELF_QUARTER_PATH, "ELF_magic.meg"), 'r') as f:
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

print(f"Loaded {len(elements)} elements from ELF quarter model")

def create_model():
    """Create model without solving."""
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
    bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])
    model = rad.ObjCnt([container, bkg])

    return model, hex_objects

# === Test 1: No Image ===
print("\n--- Test 1: No Image ---")
model1, hexes1 = create_model()
result1 = rad.Solve(model1, 0.0001, 100, 0)  # No image
print(f"  Solve iterations: {result1}")

M1 = rad.ObjM(hexes1[0])
print(f"  Element 1 magnetization: {M1}")

B1_origin = rad.Fld(model1, 'b', [0, 0, 0])
print(f"  B at origin: {[b*1000 for b in B1_origin]} mT")

# === Test 2: With Image (+x-z) using new API ===
print("\n--- Test 2: With Image (+x-z) - NEW API ---")
model2, hexes2 = create_model()
# New API: Solve with image parameter
result2 = rad.Solve(model2, 0.0001, 100, 0, image='+x-z')
print(f"  Solve iterations: {result2}")

M2 = rad.ObjM(hexes2[0])
print(f"  Element 1 magnetization: {M2}")

B2_origin = rad.Fld(model2, 'b', [0, 0, 0])
print(f"  B at origin: {[b*1000 for b in B2_origin]} mT")

# === Test 3: With Image (+x only) ===
print("\n--- Test 3: With Image (+x only) ---")
model3, hexes3 = create_model()
result3 = rad.Solve(model3, 0.0001, 100, 0, image='+x')
print(f"  Solve iterations: {result3}")

M3 = rad.ObjM(hexes3[0])
print(f"  Element 1 magnetization: {M3}")

B3_origin = rad.Fld(model3, 'b', [0, 0, 0])
print(f"  B at origin: {[b*1000 for b in B3_origin]} mT")

# === Test 4: With Image (-z only) ===
print("\n--- Test 4: With Image (-z only) ---")
model4, hexes4 = create_model()
result4 = rad.Solve(model4, 0.0001, 100, 0, image='-z')
print(f"  Solve iterations: {result4}")

M4 = rad.ObjM(hexes4[0])
print(f"  Element 1 magnetization: {M4}")

B4_origin = rad.Fld(model4, 'b', [0, 0, 0])
print(f"  B at origin: {[b*1000 for b in B4_origin]} mT")

# === Test 5: Using BuildMatrix API ===
print("\n--- Test 5: Using BuildMatrix API ---")
model5, hexes5 = create_model()
# Build matrix first
handle = rad.BuildMatrix(model5, image='+x-z')
print(f"  BuildMatrix handle: {handle}")
# Get matrix for verification
matrix, dof = rad.GetInteractMatrix(handle)
print(f"  Matrix shape: {matrix.shape}, DOF: {dof}")
# Solve (uses cached matrix)
result5 = rad.Solve(model5, 0.0001, 100, 0)
print(f"  Solve iterations: {result5}")

M5 = rad.ObjM(hexes5[0])
print(f"  Element 1 magnetization: {M5}")

# === Summary ===
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print("\nMagnetization comparison (Element 1):")
print(f"  No Image:       M = {M1}")
print(f"  +x only:        M = {M3}")
print(f"  -z only:        M = {M4}")
print(f"  +x-z (Solve):   M = {M2}")
print(f"  +x-z (Build):   M = {M5}")

print("\nField at origin comparison:")
print(f"  No Image:       Bz = {B1_origin[2]*1000:.4f} mT")
print(f"  +x only:        Bz = {B3_origin[2]*1000:.4f} mT")
print(f"  -z only:        Bz = {B4_origin[2]*1000:.4f} mT")
print(f"  +x-z (Solve):   Bz = {B2_origin[2]*1000:.4f} mT")

# Extract magnetization values for comparison
def get_mag(M):
    """Extract magnetization tuple from ObjM result."""
    return np.array(M['magnetization'])

M1_val = get_mag(M1)
M2_val = get_mag(M2)
M5_val = get_mag(M5)

# Check if magnetization changes
print("\n--- Verification ---")
if np.allclose(M1_val, M2_val, rtol=1e-3):
    print("WARNING: Magnetization NOT changed by Image API!")
else:
    print("OK: Magnetization changed by Image API")
    print(f"  No Image Mz: {M1_val[2]:.1f} A/m")
    print(f"  +x-z Mz:     {M2_val[2]:.1f} A/m")
    print(f"  Difference:  {abs(M1_val[2] - M2_val[2]):.1f} A/m")

# Check consistency between Solve(image=...) and BuildMatrix + Solve
# Note: BuildMatrix creates IMA matrix but Solve() without image param
# creates its own matrix internally, so they won't match
print("\nNote: BuildMatrix + Solve (without image) uses different matrix")
print(f"  BuildMatrix+Solve Mz: {M5_val[2]:.1f} A/m")
print(f"  Solve(image) Mz:      {M2_val[2]:.1f} A/m")
