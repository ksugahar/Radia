#!/usr/bin/env python
"""
Verify ELF vs Radia with FULL coil model (not quarter).

The original verify_elf_radia.py uses a quarter coil, but Radia's IMA
only mirrors magnetic materials, not current sources. This causes a
~13 mT difference in the coil field contribution.

Solution: Create a full coil model that accounts for all 4 quadrants
by explicitly adding the mirrored coil contributions.
"""

import sys
import os
import time
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
quater_dir = os.path.dirname(work_dir)
nonlinear_dir = os.path.dirname(quater_dir)
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(nonlinear_dir)))
sys.path.insert(0, os.path.join(repo_root, 'src'))
sys.path.insert(0, nonlinear_dir)

import radia as rad
from coil_model_quarter import create_racetrack_coil_quarter

ELF_NONLINEAR = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\nonlinear_20000AT\ELF_MMB8T_EIEM2_1x1x1\quater"

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


def create_coil_segment(center, size, j_vec):
    """Create a single rectangular current segment."""
    return rad.ObjRecCur(center, size, j_vec)


def create_full_coil(current_at=20000.0, n_arc_seg=4):
    """
    Create full coil by modeling all 4 quadrants explicitly.

    ELF uses MIMA X and MIMA -Z which mirrors the coil:
    - MIMA X: x-mirror, same current direction
    - MIMA -Z: z-mirror, reversed current direction (antisymmetric)

    For a coil segment at (x, y, z) with current (jx, jy, jz):
    - x-mirror: (-x, y, z) with current (-jx, jy, jz) [flip x, flip jx]
    - z-mirror: (x, y, -z) with current (-jx, -jy, jz) [flip z, reverse current = flip jx,jy]
    - xz-mirror: (-x, y, -z) with current (jx, -jy, jz) [both mirrors]
    """
    import numpy as np

    # Coil dimensions from ELF (in meters)
    x_inner = 0.030
    x_outer = 0.065
    width = x_outer - x_inner
    x_center = (x_inner + x_outer) / 2

    z_bottom = 0.0
    z_top = 0.0525
    height = z_top - z_bottom
    z_center = (z_bottom + z_top) / 2

    y_bottom = 0.100
    y_top = 0.1625
    straight_len = y_top - y_bottom
    y_mid = (y_bottom + y_top) / 2

    r_inner = 0.005
    r_outer = 0.040
    r_mean = (r_inner + r_outer) / 2
    arc_center_x = 0.025

    y_cap_top = 0.185
    y_cap_bottom = 0.0775

    cross_section = width * height
    j = current_at / cross_section

    coil_objects = []

    # Helper to add segment with all 4 symmetry copies
    # For racetrack coil:
    # - x-mirror: this is the return path, so jy should flip sign
    # - z-mirror: coil is symmetric in z, same current
    def add_segment_4(cx, cy, cz, sx, sy, sz, jx, jy, jz):
        # Original (x >= 0, z >= 0)
        coil_objects.append(create_coil_segment([cx, cy, cz], [sx, sy, sz], [jx, jy, jz]))
        # x-mirror (x <= 0, z >= 0): return path of the loop, jy flips
        coil_objects.append(create_coil_segment([-cx, cy, cz], [sx, sy, sz], [-jx, -jy, jz]))
        # z-mirror (x >= 0, z <= 0): coil symmetric in z, same current
        coil_objects.append(create_coil_segment([cx, cy, -cz], [sx, sy, sz], [jx, jy, jz]))
        # xz-mirror (x <= 0, z <= 0): both
        coil_objects.append(create_coil_segment([-cx, cy, -cz], [sx, sy, sz], [-jx, -jy, jz]))

    # Y-straight section
    add_segment_4(x_center, y_mid, z_center, width, straight_len, height, 0, j, 0)

    # Top arcs
    for i in range(n_arc_seg):
        phi1 = np.radians(0 + 90 * i / n_arc_seg)
        phi2 = np.radians(0 + 90 * (i + 1) / n_arc_seg)
        phi_mid = (phi1 + phi2) / 2
        seg_x = arc_center_x + r_mean * np.cos(phi_mid)
        seg_y = y_top + r_mean * np.sin(phi_mid)
        arc_len = r_mean * abs(phi2 - phi1)
        jx = -j * np.sin(phi_mid)
        jy = j * np.cos(phi_mid)
        add_segment_4(seg_x, seg_y, z_center, width, arc_len, height, jx, jy, 0)

    # Top cap (x direction, from center to +x)
    cap_x_center = arc_center_x / 2
    cap_x_len = arc_center_x
    add_segment_4(cap_x_center, y_cap_top, z_center, cap_x_len, width, height, -j, 0, 0)

    # Bottom arcs
    for i in range(n_arc_seg):
        phi1 = np.radians(-90 + 90 * i / n_arc_seg)
        phi2 = np.radians(-90 + 90 * (i + 1) / n_arc_seg)
        phi_mid = (phi1 + phi2) / 2
        seg_x = arc_center_x + r_mean * np.cos(phi_mid)
        seg_y = y_bottom + r_mean * np.sin(phi_mid)
        arc_len = r_mean * abs(phi2 - phi1)
        jx = -j * np.sin(phi_mid)
        jy = j * np.cos(phi_mid)
        add_segment_4(seg_x, seg_y, z_center, width, arc_len, height, jx, jy, 0)

    # Bottom cap
    add_segment_4(cap_x_center, y_cap_bottom, z_center, cap_x_len, width, height, j, 0, 0)

    return rad.ObjCnt(coil_objects)


print("=" * 70)
print("ELF vs Radia with FULL Coil Model")
print("=" * 70)
print("This version uses full coil (4 quadrants) to match ELF's MIMA behavior.")
print("Solver:   LU (Method 0)")
print("Image:    Quarter model (+x-z) for yoke only")
print("Coil:     FULL model (all 4 quadrants)")
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
print("Building Radia Model (FULL COIL)")
print("-" * 70)

rad.UtiDelAll()
rad.FldUnits('m')

# Create nonlinear material
mat = rad.MatSatIsoTab(bh_data)

# Create yoke geometry (quarter only - will be mirrored by IMA)
hex_objects = []
for elem_id, node_ids in elements:
    verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
             for nid in node_ids]
    hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
    rad.MatApl(hex_obj, mat)
    hex_objects.append(hex_obj)

yoke = rad.ObjCnt(hex_objects)
print(f"Yoke: {len(hex_objects)} elements (quarter, IMA will mirror)")

# Create FULL racetrack coil (all 4 quadrants)
print("Coil: Racetrack FULL (4 quadrants), 20000 AT")
coil = create_full_coil(20000.0)

# Check full coil field
B_coil_full = np.array(rad.Fld(coil, 'b', [0, 0, 0]))
print(f"\nFull coil field at (0,0,0):")
print(f"  Bx = {B_coil_full[0]*1000:.4f} mT")
print(f"  By = {B_coil_full[1]*1000:.4f} mT")
print(f"  Bz = {B_coil_full[2]*1000:.4f} mT")

# Combine model
model = rad.ObjCnt([yoke, coil])

# Solve with LU solver
print("\n" + "-" * 70)
print("Solving with LU Solver (Method 0)")
print("-" * 70)
print("  Precision: 0.0001")
print("  Max iterations: 100")
print("  Image: '+x-z' (yoke only)")

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
    print("COMPARISON (Full Coil)")
    print("=" * 70)
    print(f"ELF Bz:    {elf_B[2]*1000:.2f} mT")
    print(f"Radia Bz:  {radia_B[2]*1000:.2f} mT")
    print(f"Difference: {diff_Bz*1000:.2f} mT ({rel_diff:.2f}%)")

    if rel_diff < 1.0:
        print(f"\n*** EXCELLENT: Field within 1% ***")
        status = "EXCELLENT"
    elif rel_diff < 5.0:
        print("\n*** PASS: Field within 5% ***")
        status = "PASS"
    else:
        print(f"\n*** NEEDS REVIEW: Field differs by {rel_diff:.2f}% ***")
        status = "NEEDS REVIEW"

    print("\n" + "=" * 70)
    print("SUMMARY - Full Coil Model")
    print("=" * 70)
    print(f"Coil Bz (full): {B_coil_full[2]*1000:.2f} mT")
    print(f"ELF Bz:         {elf_B[2]*1000:.2f} mT")
    print(f"Radia Bz:       {radia_B[2]*1000:.2f} mT")
    print(f"Error:          {rel_diff:.2f}%")
    print(f"Status:         {status}")
