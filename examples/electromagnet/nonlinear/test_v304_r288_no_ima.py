#!/usr/bin/env python
"""
Test V304 and R288 WITHOUT IMA - build full model by explicit mirroring.
This helps identify if the issue is IMA-specific or geometry-specific.
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

ELF_BASE = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\nonlinear_20000AT"
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


def mirror_vertices(verts, x_mirror=False, z_mirror=False):
    """Mirror vertices across x=0 and/or z=0 planes."""
    result = []
    for v in verts:
        x = -v[0] if x_mirror else v[0]
        z = -v[2] if z_mirror else v[2]
        result.append([x, v[1], z])
    return result


def run_test_full_model(mesh_name, nodes, elements, bh_data):
    """Run test with FULL MODEL (explicit 4x mirroring, no IMA)."""
    rad.UtiDelAll()
    rad.FldUnits('m')

    mat = rad.MatSatIsoTab(bh_data)
    hex_objects = []

    # Get base vertices for each element
    base_verts_list = []
    for elem_id, node_ids in elements:
        verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
                 for nid in node_ids]
        base_verts_list.append(verts)

    # Create 4 copies: original + x-mirror + z-mirror + x-z-mirror
    mirror_configs = [
        (False, False),  # Original (x>=0, z>=0)
        (True, False),   # x-mirror (x<=0, z>=0)
        (False, True),   # z-mirror (x>=0, z<=0)
        (True, True),    # x-z-mirror (x<=0, z<=0)
    ]

    for x_mir, z_mir in mirror_configs:
        for base_verts in base_verts_list:
            verts = mirror_vertices(base_verts, x_mir, z_mir)
            hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
            rad.MatApl(hex_obj, mat)
            hex_objects.append(hex_obj)

    yoke = rad.ObjCnt(hex_objects)

    # Create coil (full model - 4 quadrants)
    coil = create_racetrack_coil(20000.0)

    model = rad.ObjCnt([yoke, coil])

    print(f"    Full model: {len(hex_objects)} elements (4 x {len(elements)})")

    # Solve WITHOUT IMA
    t_start = time.time()
    result = rad.Solve(model, 0.0001, 100, 0)  # LU, no IMA
    t_solve = time.time() - t_start

    B = np.array(rad.Fld(model, 'b', [0, 0, 0]))
    return B, t_solve


def run_test_ima(mesh_name, nodes, elements, bh_data):
    """Run test with IMA (quarter model)."""
    rad.UtiDelAll()
    rad.FldUnits('m')

    mat = rad.MatSatIsoTab(bh_data)
    hex_objects = []

    for elem_id, node_ids in elements:
        verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
                 for nid in node_ids]
        hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
        rad.MatApl(hex_obj, mat)
        hex_objects.append(hex_obj)

    yoke = rad.ObjCnt(hex_objects)

    # Create coil
    coil = create_racetrack_coil(20000.0)

    model = rad.ObjCnt([yoke, coil])

    print(f"    IMA model: {len(hex_objects)} elements")

    # Solve WITH IMA
    t_start = time.time()
    result = rad.Solve(model, 0.0001, 100, 0, image='+x-z')
    t_solve = time.time() - t_start

    B = np.array(rad.Fld(model, 'b', [0, 0, 0]))
    return B, t_solve


# Main
print("=" * 70)
print("Comparison: Full Model (no IMA) vs IMA")
print("=" * 70)

bh_file = os.path.join(work_dir, "BH.txt")
bh_data = load_bh_curve(bh_file)

for mesh_name in ["ELF_MMB8T_EIEM2_V304", "ELF_MMB8T_EIEM2_R288"]:
    print(f"\n{'='*70}")
    print(f"Testing: {mesh_name}")
    print("=" * 70)

    path = os.path.join(ELF_BASE, mesh_name)
    nodes, elements = load_elf_geometry(path)
    elf_B = load_elf_field_at_origin(path)

    print(f"  Elements (quarter): {len(elements)}")
    print(f"  ELF Bz: {elf_B[2]*1000:.2f} mT")

    # Test full model (no IMA)
    print("\n  [1] Full Model (explicit 4x mirroring, no IMA):")
    try:
        B_full, t_full = run_test_full_model(mesh_name, nodes, elements, bh_data)
        diff_full = abs(B_full[2] - elf_B[2]) / abs(elf_B[2]) * 100
        print(f"    Bz = {B_full[2]*1000:.2f} mT, time = {t_full:.2f}s")
        print(f"    vs ELF: {diff_full:.4f}%")
    except Exception as e:
        print(f"    ERROR: {e}")
        B_full = None

    # Test IMA
    print("\n  [2] Quarter + IMA (+x-z):")
    try:
        B_ima, t_ima = run_test_ima(mesh_name, nodes, elements, bh_data)
        diff_ima = abs(B_ima[2] - elf_B[2]) / abs(elf_B[2]) * 100
        print(f"    Bz = {B_ima[2]*1000:.2f} mT, time = {t_ima:.2f}s")
        print(f"    vs ELF: {diff_ima:.4f}%")
    except Exception as e:
        print(f"    ERROR: {e}")
        B_ima = None

    # Compare
    if B_full is not None and B_ima is not None:
        diff_full_ima = abs(B_full[2] - B_ima[2]) / abs(B_full[2]) * 100
        print(f"\n  Full vs IMA: {diff_full_ima:.4f}%")

        if diff_full_ima < 0.1:
            print("  => IMA is correct (matches full model)")
            print("  => Discrepancy with ELF is geometry/model issue, not IMA")
        else:
            print("  => IMA produces different result from full model!")
            print("  => IMA implementation issue detected")
