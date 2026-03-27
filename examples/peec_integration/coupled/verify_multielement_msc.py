"""
verify_multielement_msc.py

Multi-element MSC verification: MMMBuilder vs Radia for Delta_L.

Same geometry as demo_fasthenry_magnetic_core.py scenario 2:
  Wire: 100mm x 1mm x 1mm, copper
  Core: 60mm x 10mm x 10mm, centered 10mm above wire, mu_r=1000
  Core divisions: 3x1x1 (3 hexahedral elements)

Compares:
  1. Radia direct (rad.Solve + Fld('a'))
  2. MMMBuilder MSC (standalone N matrix, same kernel as Radia after fix)

This validates the multi-element case where cross-element interactions
(6x6 blocks between different elements) are tested.

Part of Radia project
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad
from mmm_core import MMMBuilder, MMMSolver
from peec_msc_schur import biot_savart_finite_filament

MU_0 = 4e-7 * np.pi


def create_core_elements(center, size, divisions):
    """Create subdivided hexahedral elements.

    Args:
        center: [cx, cy, cz] in meters
        size: [sx, sy, sz] in meters
        divisions: [nx, ny, nz]

    Returns:
        list of (8,3) vertex arrays, one per element
    """
    cx, cy, cz = center
    sx, sy, sz = size
    nx, ny, nz = divisions

    elements = []
    dx, dy, dz = sx / nx, sy / ny, sz / nz
    x0 = cx - sx / 2
    y0 = cy - sy / 2
    z0 = cz - sz / 2

    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                x = x0 + ix * dx
                y = y0 + iy * dy
                z = z0 + iz * dz
                verts = np.array([
                    [x,    y,    z],
                    [x+dx, y,    z],
                    [x+dx, y+dy, z],
                    [x,    y+dy, z],
                    [x,    y,    z+dz],
                    [x+dx, y,    z+dz],
                    [x+dx, y+dy, z+dz],
                    [x,    y+dy, z+dz],
                ])
                elements.append(verts)
    return elements


def test_radia_direct(elements, chi, seg_p1, seg_p2, seg_center, seg_dir, seg_len):
    """Radia direct: create hex objects, solve, extract Delta_L."""
    rad.UtiDelAll()

    cores = []
    for verts in elements:
        core = rad.ObjHexahedron(verts.tolist(), [0, 0, 0])
        rad.MatApl(core, rad.MatLin(chi))
        cores.append(core)

    mag_container = rad.ObjCnt(cores)

    def h_bkg(pt):
        H = biot_savart_finite_filament(seg_p1, seg_p2, pt)
        return [MU_0 * H[0], MU_0 * H[1], MU_0 * H[2]]

    bkg = rad.ObjBckg(h_bkg)
    cnt = rad.ObjCnt([mag_container, bkg])
    rad.Solve(cnt, 1e-6, 5000, 0)

    A = np.array(rad.Fld(mag_container, 'a', seg_center.tolist()))
    Delta_L = np.dot(A, seg_dir) * seg_len

    # Get M per element
    M_list = []
    for i, verts in enumerate(elements):
        ec = verts.mean(axis=0)
        M = np.array(rad.Fld(cores[i], 'm', ec.tolist()))
        M_list.append(M)

    return Delta_L, M_list


def test_mmm_builder(elements, chi, seg_p1, seg_p2, seg_center, seg_dir, seg_len):
    """MMMBuilder MSC: build N matrix, solve, extract Delta_L via Radia A-field."""
    builder = MMMBuilder()
    for verts in elements:
        builder.add_hexahedron(verts)

    N_mmm, dof_off = builder.build()
    N_mmm = np.array(N_mmm)
    dof_off = np.array(dof_off)

    n_dof = N_mmm.shape[0]
    n_elem = len(elements)
    inv_chi = np.full(n_dof, 1.0 / chi)

    # System: (1/chi + N) sigma = H_ext
    K_msc = np.diag(inv_chi) + N_mmm

    # Face normals and eval points (MMMBuilder order per element)
    normals_template = np.array([
        [0, 0, -1], [0, 0, 1], [0, -1, 0], [0, 1, 0], [-1, 0, 0], [1, 0, 0]
    ], dtype=float)

    # Build eval points per element
    eval_pts = np.zeros((n_dof, 3))
    face_normals = np.zeros((n_dof, 3))
    for e, verts in enumerate(elements):
        v = verts
        face_verts_list = [
            [v[0], v[3], v[2], v[1]],  # bottom (-Z)
            [v[4], v[5], v[6], v[7]],  # top    (+Z)
            [v[0], v[1], v[5], v[4]],  # front  (-Y)
            [v[2], v[3], v[7], v[6]],  # back   (+Y)
            [v[0], v[4], v[7], v[3]],  # left   (-X)
            [v[1], v[2], v[6], v[5]],  # right  (+X)
        ]
        elem_center = verts.mean(axis=0)
        offset = dof_off[e]
        for f, fv in enumerate(face_verts_list):
            fc = np.mean(fv, axis=0)
            eval_pts[offset + f] = (fc + elem_center) / 2.0
            face_normals[offset + f] = normals_template[f]

    # B_pm: H normal at eval points from unit PEEC current
    B_pm = np.zeros(n_dof)
    for i in range(n_dof):
        H = biot_savart_finite_filament(seg_p1, seg_p2, eval_pts[i])
        B_pm[i] = np.dot(H, face_normals[i])

    # Solve MSC
    sigma = np.linalg.solve(K_msc, B_pm)

    # Reconstruct M per element and compute Delta_L via Radia A-field
    rad.UtiDelAll()
    total_moment = np.zeros(3)
    radia_cores = []

    for e, verts in enumerate(elements):
        offset = dof_off[e]
        n_face = dof_off[e + 1] - offset
        sigma_e = sigma[offset:offset + n_face]
        n_e = face_normals[offset:offset + n_face]

        # Reconstruct M from sigma
        NtN = n_e.T @ n_e
        M_recon = np.linalg.solve(NtN, n_e.T @ sigma_e)

        core = rad.ObjHexahedron(verts.tolist(), M_recon.tolist())
        radia_cores.append(core)
        V = np.linalg.det(verts[1:4] - verts[0])  # approximate volume
        total_moment += M_recon * abs(V)

    # Compute A at segment center from all elements
    mag_cnt = rad.ObjCnt(radia_cores)
    A = np.array(rad.Fld(mag_cnt, 'a', seg_center.tolist()))
    Delta_L = np.dot(A, seg_dir) * seg_len

    # Get M per element
    M_list = []
    for e, verts in enumerate(elements):
        offset = dof_off[e]
        n_face = dof_off[e + 1] - offset
        sigma_e = sigma[offset:offset + n_face]
        n_e = face_normals[offset:offset + n_face]
        M_recon = np.linalg.solve(n_e.T @ n_e, n_e.T @ sigma_e)
        M_list.append(M_recon)

    return Delta_L, M_list, sigma


def main():
    # Geometry (same as FastHenry demo scenario 2)
    seg_p1 = np.array([0.0, 0.0, 0.0])
    seg_p2 = np.array([0.1, 0.0, 0.0])
    seg_center = (seg_p1 + seg_p2) / 2.0
    seg_dir = np.array([1.0, 0.0, 0.0])
    seg_len = 0.1

    core_center = [0.05, 0.01, 0.0]  # 10mm above wire
    core_size = [0.06, 0.01, 0.01]   # 60x10x10 mm

    chi = 999.0  # mu_r = 1000

    for divisions in [(1, 1, 1), (3, 1, 1), (3, 3, 1), (3, 3, 3)]:
        elements = create_core_elements(core_center, core_size, divisions)
        n_elem = len(elements)
        n_dof = n_elem * 6

        print(f"\n{'='*70}")
        print(f"Divisions: {divisions[0]}x{divisions[1]}x{divisions[2]} "
              f"({n_elem} elements, {n_dof} MSC DOF)")
        print(f"{'='*70}")

        # Radia direct
        DL_radia, M_radia = test_radia_direct(
            elements, chi, seg_p1, seg_p2, seg_center, seg_dir, seg_len)
        print(f"Radia direct:   Delta_L = {DL_radia*1e9:.4f} nH")

        # MMMBuilder
        DL_mmm, M_mmm, sigma_mmm = test_mmm_builder(
            elements, chi, seg_p1, seg_p2, seg_center, seg_dir, seg_len)
        print(f"MMMBuilder MSC: Delta_L = {DL_mmm*1e9:.4f} nH")

        ratio = DL_mmm / DL_radia if abs(DL_radia) > 1e-20 else float('inf')
        print(f"Ratio (MMMBuilder/Radia) = {ratio:.6f}")

        # Per-element M comparison
        for e in range(min(n_elem, 3)):
            print(f"  Element {e}: Radia M = [{M_radia[e][0]:.2f}, {M_radia[e][1]:.2f}, {M_radia[e][2]:.2f}], "
                  f"MMM M = [{M_mmm[e][0]:.2f}, {M_mmm[e][1]:.2f}, {M_mmm[e][2]:.2f}]")


if __name__ == "__main__":
    main()
