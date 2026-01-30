#!/usr/bin/env python
"""
Verify quarter model symmetry matrix formulation.

Quarter model uses two MIMA directives:
- MIMA X  (add mirror across x=0)
- MIMA -Z (subtract mirror across z=0)

For double symmetry, we need 4 copies:
1. Original (o)
2. X-mirrored (mx)
3. Z-mirrored (mz)
4. X-and-Z-mirrored (mxz)

The interaction formula should be:
A_symmetry = A_o,o + A_o,mx - A_o,mz - A_o,mxz

Where:
- A_o,mx contributes positively (from MIMA X)
- A_o,mz contributes negatively (from MIMA -Z)
- A_o,mxz combines both: x-add then z-subtract = negative
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
    """
    Mirror vertices across specified axes.

    Args:
        vertices: List of [x, y, z] coordinates
        axes: String like 'X', 'Z', 'XZ' indicating which axes to mirror
    """
    mirrored = []
    for v in vertices:
        v_new = list(v)
        if 'X' in axes.upper():
            v_new[0] = -v_new[0]
        if 'Y' in axes.upper():
            v_new[1] = -v_new[1]
        if 'Z' in axes.upper():
            v_new[2] = -v_new[2]
        mirrored.append(v_new)
    return mirrored


def main():
    meg_path = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\quater\ELF_MAGIC.meg"
    mat_path = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\quater\ELF_magic.mat"

    print("="*70)
    print("Quarter Model Verification (MIMA X + MIMA -Z)")
    print("="*70)

    # Parse files
    nodes, elements, scale, mima_list = parse_meg_file(meg_path)
    elf_matrix = read_fortran_matrix(mat_path)

    n_elem = len(elements)
    n_dof = n_elem * 6

    print(f"MIMA directives: {mima_list}")
    print(f"Elements: {n_elem}")
    print(f"DOF: {n_dof}")
    print(f"ELF matrix shape: {elf_matrix.shape}")

    # Create Radia model with all 4 copies
    rad.UtiDelAll()
    rad.FldUnits('m')

    # Original elements (o)
    orig_objs = []
    for elem in elements:
        vertices = [nodes[nid] for nid in elem['nodes']]
        obj = rad.ObjHexahedron(vertices, [0, 0, 0])
        orig_objs.append(obj)

    # X-mirrored elements (mx)
    mx_objs = []
    for elem in elements:
        vertices = [nodes[nid] for nid in elem['nodes']]
        mirrored = mirror_vertices(vertices, 'X')
        obj = rad.ObjHexahedron(mirrored, [0, 0, 0])
        mx_objs.append(obj)

    # Z-mirrored elements (mz)
    mz_objs = []
    for elem in elements:
        vertices = [nodes[nid] for nid in elem['nodes']]
        mirrored = mirror_vertices(vertices, 'Z')
        obj = rad.ObjHexahedron(mirrored, [0, 0, 0])
        mz_objs.append(obj)

    # XZ-mirrored elements (mxz)
    mxz_objs = []
    for elem in elements:
        vertices = [nodes[nid] for nid in elem['nodes']]
        mirrored = mirror_vertices(vertices, 'XZ')
        obj = rad.ObjHexahedron(mirrored, [0, 0, 0])
        mxz_objs.append(obj)

    # Create combined container (52 elements total)
    all_objs = orig_objs + mx_objs + mz_objs + mxz_objs
    container = rad.ObjCnt(all_objs)
    mat = rad.MatLin(1000)
    rad.MatApl(container, mat)

    print(f"\nCreated {len(all_objs)} total elements (13 x 4)")

    # Get full interaction matrix (208 x 208)
    intrc = rad.PreRelax(container, container)
    full_matrix_raw, dof = rad.GetInteractMatrix(intrc)
    full_matrix = np.array(full_matrix_raw)

    print(f"Full matrix shape: {full_matrix.shape}")

    # Extract submatrices
    # Layout: [o, mx, mz, mxz] each n_dof wide
    A_oo = full_matrix[:n_dof, :n_dof]              # orig -> orig
    A_o_mx = full_matrix[:n_dof, n_dof:2*n_dof]     # mx -> orig
    A_o_mz = full_matrix[:n_dof, 2*n_dof:3*n_dof]   # mz -> orig
    A_o_mxz = full_matrix[:n_dof, 3*n_dof:4*n_dof]  # mxz -> orig

    print(f"\nSubmatrix dimensions: {A_oo.shape}")

    # ELF^T for comparison
    elf_T = elf_matrix.T

    # Test different formulas
    print("\n" + "-"*70)
    print("Testing symmetry formulas:")
    print("-"*70)

    # Formula 1: A = A_oo + A_o_mx - A_o_mz - A_o_mxz
    A_test1 = A_oo + A_o_mx - A_o_mz - A_o_mxz
    err1 = np.max(np.abs(elf_T - A_test1))
    rel1 = err1 / (np.max(np.abs(elf_matrix)) + 1e-15)
    print(f"1. A_oo + A_o_mx - A_o_mz - A_o_mxz: error = {err1:.6e} (rel: {rel1:.6e})")

    # Formula 2: A = A_oo + A_o_mx + A_o_mz + A_o_mxz
    A_test2 = A_oo + A_o_mx + A_o_mz + A_o_mxz
    err2 = np.max(np.abs(elf_T - A_test2))
    rel2 = err2 / (np.max(np.abs(elf_matrix)) + 1e-15)
    print(f"2. A_oo + A_o_mx + A_o_mz + A_o_mxz: error = {err2:.6e} (rel: {rel2:.6e})")

    # Formula 3: A = A_oo + A_o_mx - A_o_mz + A_o_mxz
    A_test3 = A_oo + A_o_mx - A_o_mz + A_o_mxz
    err3 = np.max(np.abs(elf_T - A_test3))
    rel3 = err3 / (np.max(np.abs(elf_matrix)) + 1e-15)
    print(f"3. A_oo + A_o_mx - A_o_mz + A_o_mxz: error = {err3:.6e} (rel: {rel3:.6e})")

    # Formula 4: A = A_oo - A_o_mx - A_o_mz + A_o_mxz
    A_test4 = A_oo - A_o_mx - A_o_mz + A_o_mxz
    err4 = np.max(np.abs(elf_T - A_test4))
    rel4 = err4 / (np.max(np.abs(elf_matrix)) + 1e-15)
    print(f"4. A_oo - A_o_mx - A_o_mz + A_o_mxz: error = {err4:.6e} (rel: {rel4:.6e})")

    # Formula 5: A = A_oo - A_o_mx + A_o_mz - A_o_mxz
    A_test5 = A_oo - A_o_mx + A_o_mz - A_o_mxz
    err5 = np.max(np.abs(elf_T - A_test5))
    rel5 = err5 / (np.max(np.abs(elf_matrix)) + 1e-15)
    print(f"5. A_oo - A_o_mx + A_o_mz - A_o_mxz: error = {err5:.6e} (rel: {rel5:.6e})")

    # Find best
    errors = [err1, err2, err3, err4, err5]
    formulas = [
        "A_oo + A_o_mx - A_o_mz - A_o_mxz",
        "A_oo + A_o_mx + A_o_mz + A_o_mxz",
        "A_oo + A_o_mx - A_o_mz + A_o_mxz",
        "A_oo - A_o_mx - A_o_mz + A_o_mxz",
        "A_oo - A_o_mx + A_o_mz - A_o_mxz",
    ]
    best_idx = np.argmin(errors)
    print(f"\nBest formula: {formulas[best_idx]}")
    print(f"Error: {errors[best_idx]:.6e}")

    if errors[best_idx] < 1e-5:
        print("\n--> VERIFIED! Quarter model formula found.")
    else:
        print("\n--> Not yet verified. Need to try more formulas.")

        # Show diagonal block comparison for debugging
        print("\n" + "-"*70)
        print("Diagonal block comparison for element 0:")
        print("-"*70)
        print(f"ELF^T diagonal block:\n{elf_T[:6, :6]}")
        print(f"\nA_oo diagonal block:\n{A_oo[:6, :6]}")
        print(f"\nA_o_mx diagonal block:\n{A_o_mx[:6, :6]}")
        print(f"\nA_o_mz diagonal block:\n{A_o_mz[:6, :6]}")
        print(f"\nA_o_mxz diagonal block:\n{A_o_mxz[:6, :6]}")


if __name__ == '__main__':
    main()
