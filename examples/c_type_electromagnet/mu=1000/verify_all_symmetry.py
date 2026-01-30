#!/usr/bin/env python
"""
Complete verification of ELF MIMA symmetry formulations.

VERIFIED FORMULAS:
==================

1. Full model (52 elements, no symmetry):
   A_Radia^T = A_ELF

2. X-mirror (26 elements, MIMA X):
   A_symmetry = A_oo + A_o_mx
   A_symmetry^T = A_ELF

3. Z-mirror (26 elements, MIMA -Z):
   A_symmetry = A_oo - A_o_mz
   A_symmetry^T = A_ELF

4. Quarter model (13 elements, MIMA X + MIMA -Z):
   A_symmetry = A_oo + A_o_mx - A_o_mz - A_o_mxz
   A_symmetry^T = A_ELF

SIGN RULE:
==========
- MIMA X  (no sign): ADD mirror contribution
- MIMA -Z (minus sign): SUBTRACT mirror contribution
- Combined: Multiply signs
  * X and -Z -> + * - = - for mxz term

CONVENTION:
===========
- ELF: A[source][target]
- Radia: A[target][source]
- Therefore: A_ELF = A_Radia^T
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
    mima_list = []

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
                mima_list.append(' '.join(parts[1:]))
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

    return nodes, elements, scale, mima_list


def mirror_vertices(vertices, axes):
    """Mirror vertices across specified axes."""
    mirrored = []
    for v in vertices:
        v_new = list(v)
        if 'X' in axes.upper():
            v_new[0] = -v_new[0]
        if 'Z' in axes.upper():
            v_new[2] = -v_new[2]
        mirrored.append(v_new)
    return mirrored


def verify_full_model(meg_path, mat_path):
    """Verify full model (no symmetry)."""
    nodes, elements, scale, mima_list = parse_meg_file(meg_path)
    elf_matrix = read_fortran_matrix(mat_path)

    n_elem = len(elements)
    n_dof = n_elem * 6

    rad.UtiDelAll()
    rad.FldUnits('m')

    hex_objs = []
    for elem in elements:
        vertices = [nodes[nid] for nid in elem['nodes']]
        obj = rad.ObjHexahedron(vertices, [0, 0, 0])
        hex_objs.append(obj)

    container = rad.ObjCnt(hex_objs)
    mat = rad.MatLin(1000)
    rad.MatApl(container, mat)

    intrc = rad.PreRelax(container, container)
    radia_matrix_raw, dof = rad.GetInteractMatrix(intrc)
    radia_matrix = np.array(radia_matrix_raw)

    elf_T = elf_matrix.T
    error = np.max(np.abs(elf_T - radia_matrix))
    rel_error = error / (np.max(np.abs(elf_matrix)) + 1e-15)

    return n_elem, error, rel_error


def verify_single_symmetry(meg_path, mat_path):
    """Verify single symmetry case (x-mirror or z-mirror)."""
    nodes, elements, scale, mima_list = parse_meg_file(meg_path)
    elf_matrix = read_fortran_matrix(mat_path)

    n_elem = len(elements)
    n_dof = n_elem * 6

    mima = mima_list[0]
    if mima.startswith('-'):
        axis = mima[1:].upper()
        subtract = True
    else:
        axis = mima.upper()
        subtract = False

    rad.UtiDelAll()
    rad.FldUnits('m')

    orig_objs = []
    for elem in elements:
        vertices = [nodes[nid] for nid in elem['nodes']]
        obj = rad.ObjHexahedron(vertices, [0, 0, 0])
        orig_objs.append(obj)

    mirror_objs = []
    for elem in elements:
        vertices = [nodes[nid] for nid in elem['nodes']]
        mirrored = mirror_vertices(vertices, axis)
        obj = rad.ObjHexahedron(mirrored, [0, 0, 0])
        mirror_objs.append(obj)

    all_objs = orig_objs + mirror_objs
    container = rad.ObjCnt(all_objs)
    mat = rad.MatLin(1000)
    rad.MatApl(container, mat)

    intrc = rad.PreRelax(container, container)
    full_matrix_raw, dof = rad.GetInteractMatrix(intrc)
    full_matrix = np.array(full_matrix_raw)

    A_oo = full_matrix[:n_dof, :n_dof]
    A_om = full_matrix[:n_dof, n_dof:]

    if subtract:
        A_symmetry = A_oo - A_om
    else:
        A_symmetry = A_oo + A_om

    elf_T = elf_matrix.T
    error = np.max(np.abs(elf_T - A_symmetry))
    rel_error = error / (np.max(np.abs(elf_matrix)) + 1e-15)

    return n_elem, mima, error, rel_error


def verify_quarter_model(meg_path, mat_path):
    """Verify quarter model (double symmetry)."""
    nodes, elements, scale, mima_list = parse_meg_file(meg_path)
    elf_matrix = read_fortran_matrix(mat_path)

    n_elem = len(elements)
    n_dof = n_elem * 6

    rad.UtiDelAll()
    rad.FldUnits('m')

    orig_objs = []
    for elem in elements:
        vertices = [nodes[nid] for nid in elem['nodes']]
        obj = rad.ObjHexahedron(vertices, [0, 0, 0])
        orig_objs.append(obj)

    mx_objs = []
    for elem in elements:
        vertices = [nodes[nid] for nid in elem['nodes']]
        mirrored = mirror_vertices(vertices, 'X')
        obj = rad.ObjHexahedron(mirrored, [0, 0, 0])
        mx_objs.append(obj)

    mz_objs = []
    for elem in elements:
        vertices = [nodes[nid] for nid in elem['nodes']]
        mirrored = mirror_vertices(vertices, 'Z')
        obj = rad.ObjHexahedron(mirrored, [0, 0, 0])
        mz_objs.append(obj)

    mxz_objs = []
    for elem in elements:
        vertices = [nodes[nid] for nid in elem['nodes']]
        mirrored = mirror_vertices(vertices, 'XZ')
        obj = rad.ObjHexahedron(mirrored, [0, 0, 0])
        mxz_objs.append(obj)

    all_objs = orig_objs + mx_objs + mz_objs + mxz_objs
    container = rad.ObjCnt(all_objs)
    mat = rad.MatLin(1000)
    rad.MatApl(container, mat)

    intrc = rad.PreRelax(container, container)
    full_matrix_raw, dof = rad.GetInteractMatrix(intrc)
    full_matrix = np.array(full_matrix_raw)

    A_oo = full_matrix[:n_dof, :n_dof]
    A_o_mx = full_matrix[:n_dof, n_dof:2*n_dof]
    A_o_mz = full_matrix[:n_dof, 2*n_dof:3*n_dof]
    A_o_mxz = full_matrix[:n_dof, 3*n_dof:4*n_dof]

    # MIMA X -> +mx, MIMA -Z -> -mz, combined -> -mxz
    A_symmetry = A_oo + A_o_mx - A_o_mz - A_o_mxz

    elf_T = elf_matrix.T
    error = np.max(np.abs(elf_T - A_symmetry))
    rel_error = error / (np.max(np.abs(elf_matrix)) + 1e-15)

    return n_elem, mima_list, error, rel_error


def main():
    print("="*70)
    print("COMPLETE SYMMETRY VERIFICATION")
    print("ELF vs Radia MSC 6DOF Hexahedral Interaction Matrices")
    print("="*70)

    base_path = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1"

    results = []

    # 1. Full model
    print("\n1. Full Model (no symmetry)")
    print("-" * 50)
    n_elem, err, rel = verify_full_model(
        f"{base_path}/full/ELF_MAGIC.meg",
        f"{base_path}/full/ELF_magic.mat"
    )
    status = "PASS" if rel < 1e-5 else "FAIL"
    print(f"   Elements: {n_elem}, DOF: {n_elem*6}")
    print(f"   Formula: A_Radia^T = A_ELF")
    print(f"   Max error: {err:.6e}, Relative: {rel:.6e}")
    print(f"   Status: {status}")
    results.append(("Full (52 elem)", status, rel))

    # 2. X-mirror
    print("\n2. X-Mirror (MIMA X)")
    print("-" * 50)
    n_elem, mima, err, rel = verify_single_symmetry(
        f"{base_path}/x-mirror/ELF_MAGIC.meg",
        f"{base_path}/x-mirror/ELF_magic.mat"
    )
    status = "PASS" if rel < 1e-5 else "FAIL"
    print(f"   Elements: {n_elem}, DOF: {n_elem*6}")
    print(f"   MIMA: {mima}")
    print(f"   Formula: A_symmetry = A_oo + A_o_mx")
    print(f"   Max error: {err:.6e}, Relative: {rel:.6e}")
    print(f"   Status: {status}")
    results.append(("X-mirror (26 elem)", status, rel))

    # 3. Z-mirror
    print("\n3. Z-Mirror (MIMA -Z)")
    print("-" * 50)
    n_elem, mima, err, rel = verify_single_symmetry(
        f"{base_path}/z-mirror/ELF_MAGIC.meg",
        f"{base_path}/z-mirror/ELF_magic.mat"
    )
    status = "PASS" if rel < 1e-5 else "FAIL"
    print(f"   Elements: {n_elem}, DOF: {n_elem*6}")
    print(f"   MIMA: {mima}")
    print(f"   Formula: A_symmetry = A_oo - A_o_mz")
    print(f"   Max error: {err:.6e}, Relative: {rel:.6e}")
    print(f"   Status: {status}")
    results.append(("Z-mirror (26 elem)", status, rel))

    # 4. Quarter model
    print("\n4. Quarter Model (MIMA X + MIMA -Z)")
    print("-" * 50)
    n_elem, mima_list, err, rel = verify_quarter_model(
        f"{base_path}/quater/ELF_MAGIC.meg",
        f"{base_path}/quater/ELF_magic.mat"
    )
    status = "PASS" if rel < 1e-5 else "FAIL"
    print(f"   Elements: {n_elem}, DOF: {n_elem*6}")
    print(f"   MIMA: {mima_list}")
    print(f"   Formula: A_symmetry = A_oo + A_o_mx - A_o_mz - A_o_mxz")
    print(f"   Max error: {err:.6e}, Relative: {rel:.6e}")
    print(f"   Status: {status}")
    results.append(("Quarter (13 elem)", status, rel))

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"\n{'Case':<25} {'Status':<10} {'Relative Error':<15}")
    print("-" * 50)
    for case, status, rel in results:
        print(f"{case:<25} {status:<10} {rel:.2e}")

    all_pass = all(r[1] == "PASS" for r in results)
    print("\n" + "="*70)
    if all_pass:
        print("ALL TESTS PASSED!")
        print("="*70)
        print("""
VERIFIED SYMMETRY FORMULAS:
===========================

For MSC 6DOF hexahedral elements, the ELF MIMA directive specifies how to
combine interaction matrices from original and mirrored elements:

1. No symmetry:
   A_ELF = A_Radia^T

2. Single axis mirror (MIMA X or MIMA Y or MIMA Z):
   A_symmetry = A_oo + A_o_mirror
   A_ELF = A_symmetry^T

3. Single axis anti-mirror (MIMA -X or MIMA -Y or MIMA -Z):
   A_symmetry = A_oo - A_o_mirror
   A_ELF = A_symmetry^T

4. Double symmetry (MIMA X + MIMA -Z):
   A_symmetry = A_oo + A_o_mx - A_o_mz - A_o_mxz
   A_ELF = A_symmetry^T
   (Signs multiply: X gives +, -Z gives -, combined gives -)

PHYSICAL INTERPRETATION:
========================
- MIMA X: Magnetic field is symmetric across x=0 (Bx symmetric, By/Bz anti-symmetric)
- MIMA -Z: Magnetic field is anti-symmetric across z=0 (Bz symmetric, Bx/By anti-symmetric)
- The "-" sign indicates anti-symmetric boundary conditions

CONVENTION:
===========
- ELF matrix: A[source_element][target_element]
- Radia matrix: A[target_element][source_element]
- Therefore: A_ELF = A_Radia^T (full matrix transpose)
""")
    else:
        print("SOME TESTS FAILED!")
        print("="*70)


if __name__ == '__main__':
    main()
