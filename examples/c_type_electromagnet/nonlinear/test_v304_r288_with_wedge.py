#!/usr/bin/env python
"""
Test V304 and R288 with wedge (MMB6T) element support.
This test loads both MMB8T (hexahedra) and MMB6T (wedge) elements from ELF mesh.
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
    """Load ELF geometry from .meg file.

    Supports both MMB8T (8-node hexahedra) and MMB6T (6-node wedge/prism) elements.
    """
    nodes = {}
    hex_elements = []  # 8-node elements
    wedge_elements = []  # 6-node elements

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
                node_ids = [int(parts[i]) for i in range(4, 12)]  # 8 nodes
                hex_elements.append((elem_id, node_ids))
            elif line.startswith('MMB6T'):
                parts = line.split()
                elem_id = int(parts[1])
                node_ids = [int(parts[i]) for i in range(4, 10)]  # 6 nodes
                wedge_elements.append((elem_id, node_ids))

    return nodes, hex_elements, wedge_elements


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


def mirror_hex_vertices(verts, x_mirror=False, z_mirror=False):
    """Mirror hexahedron vertices across x=0 and/or z=0 planes."""
    result = []
    for v in verts:
        x = -v[0] if x_mirror else v[0]
        z = -v[2] if z_mirror else v[2]
        result.append([x, v[1], z])
    return result


def mirror_wedge_vertices(verts, x_mirror=False, z_mirror=False):
    """Mirror wedge vertices across x=0 and/or z=0 planes."""
    result = []
    for v in verts:
        x = -v[0] if x_mirror else v[0]
        z = -v[2] if z_mirror else v[2]
        result.append([x, v[1], z])
    return result


def run_test_full_model(mesh_name, nodes, hex_elements, wedge_elements, bh_data):
    """Run test with FULL MODEL (explicit 4x mirroring, no IMA)."""
    rad.UtiDelAll()
    rad.FldUnits('m')

    mat = rad.MatSatIsoTab(bh_data)
    all_objects = []

    # Get base vertices for hexahedra
    hex_base_verts = []
    for elem_id, node_ids in hex_elements:
        verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
                 for nid in node_ids]
        hex_base_verts.append(verts)

    # Get base vertices for wedges
    wedge_base_verts = []
    for elem_id, node_ids in wedge_elements:
        verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
                 for nid in node_ids]
        wedge_base_verts.append(verts)

    # Create 4 copies: original + x-mirror + z-mirror + x-z-mirror
    mirror_configs = [
        (False, False),  # Original (x>=0, z>=0)
        (True, False),   # x-mirror (x<=0, z>=0)
        (False, True),   # z-mirror (x>=0, z<=0)
        (True, True),    # x-z-mirror (x<=0, z<=0)
    ]

    for x_mir, z_mir in mirror_configs:
        # Add hexahedra
        for base_verts in hex_base_verts:
            verts = mirror_hex_vertices(base_verts, x_mir, z_mir)
            obj = rad.ObjHexahedron(verts, [0, 0, 0])
            rad.MatApl(obj, mat)
            all_objects.append(obj)

        # Add wedges
        for base_verts in wedge_base_verts:
            verts = mirror_wedge_vertices(base_verts, x_mir, z_mir)
            obj = rad.ObjWedge(verts, [0, 0, 0])
            rad.MatApl(obj, mat)
            all_objects.append(obj)

    yoke = rad.ObjCnt(all_objects)

    # Create coil (full model - 4 quadrants)
    coil = create_racetrack_coil(20000.0)

    model = rad.ObjCnt([yoke, coil])

    n_total = len(hex_elements) + len(wedge_elements)
    print(f"    Full model: {len(all_objects)} elements = 4 x ({len(hex_elements)} hex + {len(wedge_elements)} wedge)")

    # Solve WITHOUT IMA
    t_start = time.time()
    result = rad.Solve(model, 0.0001, 100, 0)  # LU, no IMA
    t_solve = time.time() - t_start

    B = np.array(rad.Fld(model, 'b', [0, 0, 0]))
    return B, t_solve


def run_test_ima(mesh_name, nodes, hex_elements, wedge_elements, bh_data):
    """Run test with IMA (quarter model)."""
    rad.UtiDelAll()
    rad.FldUnits('m')

    mat = rad.MatSatIsoTab(bh_data)
    all_objects = []

    # Add hexahedra
    for elem_id, node_ids in hex_elements:
        verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
                 for nid in node_ids]
        obj = rad.ObjHexahedron(verts, [0, 0, 0])
        rad.MatApl(obj, mat)
        all_objects.append(obj)

    # Add wedges
    for elem_id, node_ids in wedge_elements:
        verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
                 for nid in node_ids]
        obj = rad.ObjWedge(verts, [0, 0, 0])
        rad.MatApl(obj, mat)
        all_objects.append(obj)

    yoke = rad.ObjCnt(all_objects)

    # Create coil
    coil = create_racetrack_coil(20000.0)

    model = rad.ObjCnt([yoke, coil])

    print(f"    IMA model: {len(hex_elements)} hex + {len(wedge_elements)} wedge = {len(all_objects)} elements")

    # Solve WITH IMA
    t_start = time.time()
    result = rad.Solve(model, 0.0001, 100, 0, image='+x-z')
    t_solve = time.time() - t_start

    B = np.array(rad.Fld(model, 'b', [0, 0, 0]))
    return B, t_solve


# Main
print("=" * 70)
print("V304 and R288 with Wedge (MMB6T) Element Support")
print("=" * 70)

bh_file = os.path.join(work_dir, "BH.txt")
bh_data = load_bh_curve(bh_file)

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

    # Test full model (no IMA)
    print("\n  [1] Full Model (explicit 4x mirroring, no IMA):")
    try:
        B_full, t_full = run_test_full_model(mesh_name, nodes, hex_elements, wedge_elements, bh_data)
        diff_full = abs(B_full[2] - elf_B[2]) / abs(elf_B[2]) * 100
        print(f"    Bz = {B_full[2]*1000:.2f} mT, time = {t_full:.2f}s")
        print(f"    vs ELF: {diff_full:.4f}%")
    except Exception as e:
        import traceback
        print(f"    ERROR: {e}")
        traceback.print_exc()
        B_full = None

    # Test IMA
    print("\n  [2] Quarter + IMA (+x-z):")
    try:
        B_ima, t_ima = run_test_ima(mesh_name, nodes, hex_elements, wedge_elements, bh_data)
        diff_ima = abs(B_ima[2] - elf_B[2]) / abs(elf_B[2]) * 100
        print(f"    Bz = {B_ima[2]*1000:.2f} mT, time = {t_ima:.2f}s")
        print(f"    vs ELF: {diff_ima:.4f}%")
    except Exception as e:
        import traceback
        print(f"    ERROR: {e}")
        traceback.print_exc()
        B_ima = None

    # Compare
    if B_full is not None and B_ima is not None:
        diff_full_ima = abs(B_full[2] - B_ima[2]) / abs(B_full[2]) * 100
        print(f"\n  Full vs IMA: {diff_full_ima:.4f}%")

        if diff_full < 1.0:
            print(f"  => SUCCESS! Radia matches ELF within {diff_full:.4f}%")
        else:
            print(f"  => Discrepancy with ELF: {diff_full:.4f}%")

print("\n" + "=" * 70)
print("Test completed.")
print("=" * 70)
