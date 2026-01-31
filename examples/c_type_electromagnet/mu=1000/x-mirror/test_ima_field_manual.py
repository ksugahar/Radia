#!/usr/bin/env python
"""
Test IMA field by manually computing contributions.

This script tests the IMA field computation by:
1. Solving with IMA to get sigma values
2. Computing field manually from sigma values (without IMA context)
3. Comparing with automatic IMA field computation
"""

import sys
import os
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(work_dir)
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(parent_dir)))
sys.path.insert(0, os.path.join(repo_root, 'src'))

import radia as rad

ELF_X_MIRROR = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\x-mirror"

MU_R = 1000
scale = 0.001

N_TURNS = 1000
CURRENT = 2.0
COIL_AT = N_TURNS * CURRENT
GAP_CENTER = [0.0, 0.055, 0.0275]


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


print("=" * 70)
print("Test IMA Field - Manual Computation")
print("=" * 70)

# Load geometry
nodes, elements = load_elf_geometry(ELF_X_MIRROR)
n_elem = len(elements)
print(f"X-mirror model: {n_elem} elements")

# Calculate background field
H_ext = COIL_AT / (2 * 0.055)
B_ext_y = 4e-7 * np.pi * H_ext
print(f"External field: B_y = {B_ext_y*1e3:.2f} mT")

# Test 1: Full model (explicit duplication, no IMA)
print("\n--- Test 1: Full model (explicit duplication) ---")

rad.UtiDelAll()
rad.FldUnits('m')

hex_objects_orig = []
for elem_id, node_ids in elements:
    verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
             for nid in node_ids]
    hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
    mat = rad.MatLin(MU_R)
    rad.MatApl(hex_obj, mat)
    hex_objects_orig.append(hex_obj)

hex_objects_mirror = []
for elem_id, node_ids in elements:
    verts_mirror = [[-nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
                     for nid in node_ids]
    hex_obj = rad.ObjHexahedron(verts_mirror, [0, 0, 0])
    mat = rad.MatLin(MU_R)
    rad.MatApl(hex_obj, mat)
    hex_objects_mirror.append(hex_obj)

iron_full = rad.ObjCnt(hex_objects_orig + hex_objects_mirror)
bkg_full = rad.ObjBckg(lambda p: [0, B_ext_y, 0])
container_full = rad.ObjCnt([iron_full, bkg_full])

print("Solving full model...")
result_full = rad.Solve(container_full, 0.0001, 1000, 0)
print(f"Solve result: {result_full}")

B_full = rad.Fld(container_full, 'b', GAP_CENTER)
print(f"B at gap: Bx={B_full[0]*1e3:.2f}, By={B_full[1]*1e3:.2f}, Bz={B_full[2]*1e3:.2f} mT")

# Compute field from original elements only (should be half of total)
B_orig_only = rad.Fld(rad.ObjCnt(hex_objects_orig), 'b', GAP_CENTER)
print(f"B from original elements only: By={B_orig_only[1]*1e3:.2f}, Bz={B_orig_only[2]*1e3:.2f} mT")

# Test 2: IMA model
print("\n--- Test 2: IMA model ---")

rad.UtiDelAll()
rad.FldUnits('m')

hex_objects_ima = []
for elem_id, node_ids in elements:
    verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
             for nid in node_ids]
    hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
    mat = rad.MatLin(MU_R)
    rad.MatApl(hex_obj, mat)
    hex_objects_ima.append(hex_obj)

iron_ima = rad.ObjCnt(hex_objects_ima)
bkg_ima = rad.ObjBckg(lambda p: [0, B_ext_y, 0])
container_ima = rad.ObjCnt([iron_ima, bkg_ima])

print("Solving IMA model...")
result_ima = rad.Solve(container_ima, 0.0001, 1000, 0, image='+x')
print(f"Solve result: {result_ima}")

# Compute field WITHOUT IMA context (just from the elements as-is)
B_ima_direct = rad.Fld(iron_ima, 'b', GAP_CENTER)
print(f"B from IMA elements (no IMA ctx): By={B_ima_direct[1]*1e3:.2f}, Bz={B_ima_direct[2]*1e3:.2f} mT")

# Compute field WITH IMA context (includes mirror contributions)
B_ima_full = rad.Fld(container_ima, 'b', GAP_CENTER)
print(f"B from IMA elements (with IMA ctx): By={B_ima_full[1]*1e3:.2f}, Bz={B_ima_full[2]*1e3:.2f} mT")

# Summary
print("\n" + "=" * 70)
print("Summary")
print("=" * 70)
print(f"Full model (52 elem):       By={B_full[1]*1e3:.2f}, Bz={B_full[2]*1e3:.2f} mT")
print(f"Full orig only (26 elem):   By={B_orig_only[1]*1e3:.2f}, Bz={B_orig_only[2]*1e3:.2f} mT")
print(f"IMA direct (26 elem):       By={B_ima_direct[1]*1e3:.2f}, Bz={B_ima_direct[2]*1e3:.2f} mT")
print(f"IMA with ctx (26 elem):     By={B_ima_full[1]*1e3:.2f}, Bz={B_ima_full[2]*1e3:.2f} mT")

# Expected: IMA direct should equal orig-only, IMA with ctx should equal full
print(f"\nIMA direct vs orig-only: dBy={abs(B_ima_direct[1]-B_orig_only[1])*1e3:.2f} mT")
print(f"IMA with ctx vs full:    dBy={abs(B_ima_full[1]-B_full[1])*1e3:.2f} mT")
