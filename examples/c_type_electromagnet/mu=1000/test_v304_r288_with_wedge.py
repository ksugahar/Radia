#!/usr/bin/env python
"""
Test V304 and R288 with wedge (MMB6T) element support - LINEAR material (mu=1000).

Compares Radia results against ELF reference:
  S:/ELF_MAGIC/2020_03_07_CEFC_2020/model_C-Type/mu=1000/ELF_MMB8T_EIEM2_V304
  S:/ELF_MAGIC/2020_03_07_CEFC_2020/model_C-Type/mu=1000/ELF_MMB8T_EIEM2_R288

Tests:
  [1] Full model (explicit 4x mirroring, no IMA)
  [2] Quarter model with IMA (+x-z)
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

ELF_BASE = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000"
MU_R = 1000
scale = 0.001  # mm to m


def load_elf_geometry(path):
    """Load ELF geometry from .meg file.

    Supports both MMB8T (8-node hexahedra) and MMB6T (6-node wedge/prism) elements.
    """
    nodes = {}
    hex_elements = []
    wedge_elements = []

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
                hex_elements.append(node_ids)
            elif line.startswith('MMB6T'):
                parts = line.split()
                node_ids = [int(parts[i]) for i in range(4, 10)]
                wedge_elements.append(node_ids)

    return nodes, hex_elements, wedge_elements


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


def mirror_vertices(verts, x_mirror=False, z_mirror=False):
    """Mirror vertices across x=0 and/or z=0 planes."""
    result = []
    for v in verts:
        x = -v[0] if x_mirror else v[0]
        z = -v[2] if z_mirror else v[2]
        result.append([x, v[1], z])
    return result


def run_test_full_model(nodes, hex_elements, wedge_elements):
    """Run test with FULL MODEL (explicit 4x mirroring, no IMA)."""
    rad.UtiDelAll()
    rad.FldUnits('m')

    mat = rad.MatLin(MU_R)
    all_objects = []

    hex_base = [[[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
                 for nid in nids] for nids in hex_elements]
    wedge_base = [[[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
                   for nid in nids] for nids in wedge_elements]

    for x_mir, z_mir in [(False, False), (True, False), (False, True), (True, True)]:
        for verts in hex_base:
            obj = rad.ObjHexahedron(mirror_vertices(verts, x_mir, z_mir), [0, 0, 0])
            rad.MatApl(obj, mat)
            all_objects.append(obj)
        for verts in wedge_base:
            obj = rad.ObjWedge(mirror_vertices(verts, x_mir, z_mir), [0, 0, 0])
            rad.MatApl(obj, mat)
            all_objects.append(obj)

    yoke = rad.ObjCnt(all_objects)
    coil = create_racetrack_coil(2000.0)
    model = rad.ObjCnt([yoke, coil])

    print(f"    Full model: {len(all_objects)} elements = 4 x ({len(hex_elements)} hex + {len(wedge_elements)} wedge)")

    t_start = time.time()
    result = rad.Solve(model, 0.0001, 100, 0)
    t_solve = time.time() - t_start

    B = np.array(rad.Fld(model, 'b', [0, 0, 0]))
    return B, t_solve


def run_test_ima(nodes, hex_elements, wedge_elements):
    """Run test with IMA (quarter model)."""
    rad.UtiDelAll()
    rad.FldUnits('m')

    mat = rad.MatLin(MU_R)
    all_objects = []

    for nids in hex_elements:
        verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
                 for nid in nids]
        obj = rad.ObjHexahedron(verts, [0, 0, 0])
        rad.MatApl(obj, mat)
        all_objects.append(obj)

    for nids in wedge_elements:
        verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
                 for nid in nids]
        obj = rad.ObjWedge(verts, [0, 0, 0])
        rad.MatApl(obj, mat)
        all_objects.append(obj)

    yoke = rad.ObjCnt(all_objects)
    coil = create_racetrack_coil(2000.0)
    model = rad.ObjCnt([yoke, coil])

    print(f"    IMA model: {len(hex_elements)} hex + {len(wedge_elements)} wedge = {len(all_objects)} elements")

    t_start = time.time()
    result = rad.Solve(model, 0.0001, 100, 0, image='+x-z')
    t_solve = time.time() - t_start

    B = np.array(rad.Fld(model, 'b', [0, 0, 0]))
    return B, t_solve


# Main
print("=" * 70)
print("V304 and R288 with Wedge - LINEAR material (mu_r=%d)" % MU_R)
print("ELF ref: %s" % ELF_BASE)
print("=" * 70)

for mesh_name in ["ELF_MMB8T_EIEM2_V304", "ELF_MMB8T_EIEM2_R288"]:
    print(f"\n{'='*70}")
    print(f"Testing: {mesh_name}")
    print("=" * 70)

    path = os.path.join(ELF_BASE, mesh_name)
    nodes, hex_elements, wedge_elements = load_elf_geometry(path)
    elf_B = load_elf_field_at_origin(path)

    print(f"  Hexahedra (MMB8T): {len(hex_elements)}")
    print(f"  Wedges (MMB6T): {len(wedge_elements)}")
    print(f"  Total elements: {len(hex_elements) + len(wedge_elements)}")
    print(f"  ELF Bz: {elf_B[2]*1000:.2f} mT")

    # Test full model
    print("\n  [1] Full Model (explicit 4x mirroring, no IMA):")
    B_full, t_full = run_test_full_model(nodes, hex_elements, wedge_elements)
    diff_full = abs(B_full[2] - elf_B[2]) / abs(elf_B[2]) * 100
    print(f"    Bz = {B_full[2]*1000:.2f} mT, time = {t_full:.2f}s")
    print(f"    vs ELF: {diff_full:.4f}%")

    # Test IMA
    print("\n  [2] Quarter + IMA (+x-z):")
    B_ima, t_ima = run_test_ima(nodes, hex_elements, wedge_elements)
    diff_ima = abs(B_ima[2] - elf_B[2]) / abs(elf_B[2]) * 100
    print(f"    Bz = {B_ima[2]*1000:.2f} mT, time = {t_ima:.2f}s")
    print(f"    vs ELF: {diff_ima:.4f}%")

    # Compare
    diff_full_ima = abs(B_full[2] - B_ima[2]) / abs(B_full[2]) * 100
    print(f"\n  Full vs IMA: {diff_full_ima:.4f}%")

    if diff_full < 5.0:
        print(f"  => Radia matches ELF within {diff_full:.2f}%")
    else:
        print(f"  => Discrepancy with ELF: {diff_full:.2f}%")

print(f"\n{'='*70}")
print("Test completed.")
print("=" * 70)
