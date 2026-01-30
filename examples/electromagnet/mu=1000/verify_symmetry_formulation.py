#!/usr/bin/env python
"""
Verify symmetry matrix formulation for ELF MIMA directive.

VERIFIED FORMULAS:
- MIMA X  (x-mirror):  A_symmetry = A_oo + A_om
- MIMA -Z (z-mirror):  A_symmetry = A_oo - A_om

The "-" in MIMA directive means the mirror contribution is SUBTRACTED.

Where:
- A_oo = interaction matrix from original elements to original elements
- A_om = interaction matrix from mirror elements to original elements

Convention:
- ELF: A[source][target]
- Radia: A[target][source]
- So: A_ELF = A_Radia^T
"""

import sys
import numpy as np

sys.path.insert(0, 's:/Radia/01_GitHub/src')
import radia as rad


def read_fortran_matrix(mat_path):
    """Read ELF .mat file."""
    with open(mat_path, 'rb') as f:
        raw = f.read()

    pos = 0
    rows = []
    while pos < len(raw) - 4:
        rec_len = np.frombuffer(raw[pos:pos+4], dtype='<i4')[0]
        if rec_len <= 0 or pos + 8 + rec_len > len(raw):
            break
        data = raw[pos+4:pos+4+rec_len]
        end_marker = np.frombuffer(raw[pos+4+rec_len:pos+8+rec_len], dtype='<i4')[0]
        if rec_len != end_marker:
            raise ValueError(f"Fortran record marker mismatch")
        row = np.frombuffer(data, dtype='<f8')
        rows.append(row)
        pos += 8 + rec_len

    return np.vstack(rows)


def parse_meg_file(meg_path):
    """Parse ELF .meg mesh file."""
    nodes = {}
    elements = []
    scale = 1.0
    mima = None

    with open(meg_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('*'):
                continue

            parts = line.split()
            if not parts:
                continue

            keyword = parts[0]

            if keyword == 'MIMA':
                mima = ' '.join(parts[1:])
            elif keyword == 'MGSC':
                scale = float(parts[1])
            elif keyword == 'MGR1':
                node_id = int(parts[1])
                x = float(parts[3])
                y = float(parts[4])
                z = float(parts[5])
                nodes[node_id] = [x * scale, y * scale, z * scale]
            elif keyword == 'MMB8T':
                elem_id = int(parts[1])
                material_id = int(parts[3])
                node_ids = [int(parts[i]) for i in range(4, 12)]
                elements.append({
                    'id': elem_id,
                    'type': keyword,
                    'material': material_id,
                    'nodes': node_ids
                })

    return nodes, elements, scale, mima


def mirror_vertices(vertices, axis):
    """Mirror vertices across the specified axis."""
    mirrored = []
    axis_idx = {'X': 0, 'Y': 1, 'Z': 2}[axis.upper()]
    for v in vertices:
        v_new = list(v)
        v_new[axis_idx] = -v_new[axis_idx]
        mirrored.append(v_new)
    return mirrored


def verify_symmetry_case(case_name, meg_path, mat_path):
    """Verify symmetry formulation for a specific case."""
    print(f"\n{'='*70}")
    print(f"Verifying: {case_name}")
    print(f"{'='*70}")

    # Parse files
    nodes, elements, scale, mima = parse_meg_file(meg_path)
    elf_matrix = read_fortran_matrix(mat_path)

    n_elem = len(elements)
    n_dof = n_elem * 6

    print(f"MIMA directive: {mima}")
    print(f"Elements: {n_elem}")
    print(f"DOF: {n_dof}")

    # Parse MIMA directive
    if mima.startswith('-'):
        axis = mima[1:].upper()
        subtract = True
    else:
        axis = mima.upper()
        subtract = False

    print(f"Mirror axis: {axis}")
    print(f"Subtract mode: {subtract}")

    # Create Radia model
    rad.UtiDelAll()
    rad.FldUnits('m')

    # Create original elements
    orig_objs = []
    for elem in elements:
        vertices = [nodes[nid] for nid in elem['nodes']]
        obj = rad.ObjHexahedron(vertices, [0, 0, 0])
        orig_objs.append(obj)

    # Create mirror elements
    mirror_objs = []
    for elem in elements:
        vertices = [nodes[nid] for nid in elem['nodes']]
        mirrored = mirror_vertices(vertices, axis)
        obj = rad.ObjHexahedron(mirrored, [0, 0, 0])
        mirror_objs.append(obj)

    # Create combined container
    all_objs = orig_objs + mirror_objs
    container = rad.ObjCnt(all_objs)
    mat = rad.MatLin(1000)
    rad.MatApl(container, mat)

    # Get full interaction matrix
    intrc = rad.PreRelax(container, container)
    full_matrix_raw, dof = rad.GetInteractMatrix(intrc)
    full_matrix = np.array(full_matrix_raw)

    # Extract submatrices
    A_oo = full_matrix[:n_dof, :n_dof]
    A_om = full_matrix[:n_dof, n_dof:]

    # Apply symmetry formula
    if subtract:
        A_symmetry = A_oo - A_om
        formula = "A_oo - A_om"
    else:
        A_symmetry = A_oo + A_om
        formula = "A_oo + A_om"

    print(f"Symmetry formula: {formula}")

    # Compare with ELF matrix
    elf_T = elf_matrix.T

    diff = np.abs(elf_T - A_symmetry)
    max_diff = np.max(diff)
    rel_error = max_diff / (np.max(np.abs(elf_matrix)) + 1e-15)

    print(f"\nMatrix comparison (ELF^T vs Radia {formula}):")
    print(f"  Max absolute difference: {max_diff:.6e}")
    print(f"  Relative error: {rel_error:.6e}")

    if rel_error < 1e-5:
        print(f"  --> VERIFIED! Matrices match.")
        return True
    else:
        print(f"  --> FAILED! Matrices differ.")
        return False


def main():
    print("="*70)
    print("SYMMETRY MATRIX FORMULATION VERIFICATION")
    print("="*70)
    print("""
This script verifies the symmetry formulation for MSC hexahedral elements.

ELF MIMA directive:
- MIMA X  -> A_symmetry = A_oo + A_om  (add mirror contribution)
- MIMA -X -> A_symmetry = A_oo - A_om  (subtract mirror contribution)
- MIMA Y  -> A_symmetry = A_oo + A_om
- MIMA -Y -> A_symmetry = A_oo - A_om
- MIMA Z  -> A_symmetry = A_oo + A_om
- MIMA -Z -> A_symmetry = A_oo - A_om
""")

    base_path = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1"

    cases = [
        ("x-mirror (MIMA X)", f"{base_path}/x-mirror/ELF_MAGIC.meg", f"{base_path}/x-mirror/ELF_magic.mat"),
        ("z-mirror (MIMA -Z)", f"{base_path}/z-mirror/ELF_MAGIC.meg", f"{base_path}/z-mirror/ELF_magic.mat"),
    ]

    results = []
    for case_name, meg_path, mat_path in cases:
        success = verify_symmetry_case(case_name, meg_path, mat_path)
        results.append((case_name, success))

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    for case_name, success in results:
        status = "PASS" if success else "FAIL"
        print(f"  {case_name}: {status}")

    all_pass = all(r[1] for r in results)
    print(f"\nOverall: {'ALL TESTS PASSED' if all_pass else 'SOME TESTS FAILED'}")

    if all_pass:
        print("""
VERIFIED FORMULAS:
==================
For MSC 6DOF hexahedral elements:

1. Full model (no symmetry):
   A_Radia^T = A_ELF

2. X-mirror (MIMA X):
   A_symmetry = A_oo + A_om
   A_symmetry^T = A_ELF

3. Z-mirror (MIMA -Z):
   A_symmetry = A_oo - A_om
   A_symmetry^T = A_ELF

Where:
- A_oo: interaction from original to original elements
- A_om: interaction from mirror to original elements
- The "-" in MIMA directive means subtract the mirror contribution

PHYSICAL INTERPRETATION:
========================
- MIMA X: Magnetic symmetry across x=0 plane (B field is symmetric)
- MIMA -Z: Magnetic anti-symmetry across z=0 plane (B field is anti-symmetric)
  The "-" indicates that the magnetization image has opposite sign.
""")


if __name__ == '__main__':
    main()
