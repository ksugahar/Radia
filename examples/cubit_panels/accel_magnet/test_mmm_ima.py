#!/usr/bin/env python
"""
Diagnostic: MMM IMA for tetrahedra.

Compare full model (4 tets, no IMA) vs quarter model (1 tet, IMA '+x-z').
Uses a simple configuration to isolate the sign matrix issue.

Configuration:
  - 4 identical tetrahedra in +-x, +-z quadrants (x-z symmetry)
  - mu_r = 1000 (linear)
  - Uniform background field Bz = 0.1 T via ObjBckg callback
  - Solve LU, compare B at origin
"""

import sys
import os
import math
import numpy as np

repo_root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..')
sys.path.insert(0, os.path.join(repo_root, 'src'))
sys.path.insert(0, os.path.join(repo_root, 'src', 'radia'))

import radia as rad

MU_0 = 4e-7 * math.pi
MU_R = 1000
mm = 1e-3

# Define tetrahedra forming a plate near z=0 (in +x, +z quadrant)
# The closer to the symmetry planes, the more IMA matters
# Create a grid of tets in the quarter region x=[0..20], y=[-10..10], z=[0..5] mm
def make_tet_grid(x0, x1, y0, y1, z0, z1, nx, ny, nz):
    """Create a grid of tetrahedra by subdividing each voxel into 6 tets."""
    tets = []
    dx = (x1 - x0) / nx
    dy = (y1 - y0) / ny
    dz = (z1 - z0) / nz
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                # 8 corners of the voxel
                cx = x0 + ix * dx
                cy = y0 + iy * dy
                cz = z0 + iz * dz
                c = [
                    [cx, cy, cz],           # 0
                    [cx+dx, cy, cz],         # 1
                    [cx+dx, cy+dy, cz],      # 2
                    [cx, cy+dy, cz],         # 3
                    [cx, cy, cz+dz],         # 4
                    [cx+dx, cy, cz+dz],      # 5
                    [cx+dx, cy+dy, cz+dz],   # 6
                    [cx, cy+dy, cz+dz],      # 7
                ]
                # Subdivide hex into 6 tets (Kuhn triangulation)
                tets.append([c[0], c[1], c[3], c[4]])
                tets.append([c[1], c[2], c[3], c[6]])
                tets.append([c[1], c[3], c[4], c[6]])
                tets.append([c[3], c[4], c[6], c[7]])
                tets.append([c[1], c[4], c[5], c[6]])
                tets.append([c[4], c[5], c[6], c[7]])  # 6th tet is degenerate for some configs
    return tets

# Iron plate: 20x20x5 mm, quarter = x=[1..10], y=[-10..10], z=[1..5] mm
# Offset from axes to avoid boundary issues
tet_grid_q1 = make_tet_grid(1*mm, 10*mm, -10*mm, 10*mm, 1*mm, 5*mm, 3, 4, 1)

# Single tet near symmetry planes for matrix comparison
single_tet_q1 = [[
    [2*mm, -3*mm, 2*mm],
    [8*mm, -3*mm, 2*mm],
    [5*mm,  3*mm, 2*mm],
    [5*mm,  0*mm, 6*mm],
]]

def mirror_verts(verts, mx=False, mz=False):
    """Mirror vertices about x=0 and/or z=0."""
    result = []
    for v in verts:
        x, y, z = v
        if mx: x = -x
        if mz: z = -z
        result.append([x, y, z])
    return result


def make_objs(tet_list, mu_r=MU_R):
    """Create Radia tet objects from vertex list."""
    objs = []
    for verts in tet_list:
        obj = rad.ObjTetrahedron(verts, [0, 0, 0])
        mat = rad.MatLin(mu_r)
        rad.MatApl(obj, mat)
        objs.append(obj)
    return objs


def run_full_model():
    """Full model: all quadrant tets, no IMA."""
    rad.UtiDelAll()

    # Create full model: mirror the quarter grid to all 4 quadrants
    all_tets = []
    for mx in [False, True]:
        for mz in [False, True]:
            for verts in tet_grid_q1:
                all_tets.append(mirror_verts(verts, mx=mx, mz=mz))

    all_objs = make_objs(all_tets)
    container = rad.ObjCnt(all_objs)

    bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])
    model = rad.ObjCnt([container, bkg])

    result = rad.Solve(model, 0.0001, 100, 0)
    print(f"Full model solve: {result}")

    B = np.array(rad.Fld(model, 'b', [0, 0, 0]))
    print(f"Full model B at origin: ({B[0]*1e3:.4f}, {B[1]*1e3:.4f}, {B[2]*1e3:.4f}) mT")
    print(f"Full model |Bz| = {abs(B[2]*1e3):.4f} mT")
    print(f"Full model: {len(all_objs)} tets, {len(all_objs)*3} DOF")

    # Get full model M for q1 elements (first 72)
    full_M_q1 = []
    for obj in all_objs[:len(tet_grid_q1)]:
        info = rad.ObjM(obj)
        full_M_q1.append(info['magnetization'])
    M_full = np.array(full_M_q1)
    print(f"Full model M (q1): mean={np.mean(np.abs(M_full), axis=0)}")

    return B, full_M_q1


def run_quarter_ima():
    """Quarter model: +x,+z tets, IMA '+x-z'."""
    rad.UtiDelAll()

    objs = make_objs(tet_grid_q1)
    container = rad.ObjCnt(objs)

    bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])
    model = rad.ObjCnt([container, bkg])

    result = rad.Solve(model, 0.0001, 100, 0, image='+x-z')
    print(f"IMA quarter solve: {result}")

    B = np.array(rad.Fld(model, 'b', [0, 0, 0]))
    print(f"IMA quarter B at origin: ({B[0]*1e3:.4f}, {B[1]*1e3:.4f}, {B[2]*1e3:.4f}) mT")
    print(f"IMA quarter |Bz| = {abs(B[2]*1e3):.4f} mT")
    print(f"IMA quarter: {len(objs)} tets, {len(objs)*3} DOF")

    return B


def run_quarter_no_ima():
    """Quarter model: +x,+z tets, no IMA."""
    rad.UtiDelAll()

    objs = make_objs(tet_grid_q1)
    container = rad.ObjCnt(objs)

    bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])
    model = rad.ObjCnt([container, bkg])

    result = rad.Solve(model, 0.0001, 100, 0)
    print(f"Quarter no-IMA solve: {result}")

    B = np.array(rad.Fld(model, 'b', [0, 0, 0]))
    print(f"Quarter no-IMA B: ({B[0]*1e3:.4f}, {B[1]*1e3:.4f}, {B[2]*1e3:.4f}) mT")

    return B


def run_quarter_ima_with_image_field():
    """Quarter model + IMA, then add image field using ima_field module."""
    from ima_field import add_ima_images

    rad.UtiDelAll()

    q_objs = make_objs(tet_grid_q1)
    container = rad.ObjCnt(q_objs)

    bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])
    model = rad.ObjCnt([container, bkg])

    result = rad.Solve(model, 0.0001, 100, 0, image='+x-z')
    print(f"IMA+image solve: {result}")

    B_before = np.array(rad.Fld(model, 'b', [0, 0, 0]))
    print(f"B before images: ({B_before[0]*1e3:.4f}, {B_before[1]*1e3:.4f}, {B_before[2]*1e3:.4f}) mT")

    # Add image elements using the utility function
    img_objs = add_ima_images(q_objs, tet_grid_q1, image='+x-z', container=model)
    print(f"Created {len(img_objs)} image elements")

    B_after = np.array(rad.Fld(model, 'b', [0, 0, 0]))
    print(f"B after images:  ({B_after[0]*1e3:.4f}, {B_after[1]*1e3:.4f}, {B_after[2]*1e3:.4f}) mT")

    # Get solved M values for comparison
    solved_M = []
    for obj in q_objs:
        info = rad.ObjM(obj)
        solved_M.append(info['magnetization'])

    M_summary = np.array(solved_M)
    print(f"\nIMA solved M: mean={np.mean(np.abs(M_summary), axis=0)}")
    print(f"  max|M| = {np.max(np.linalg.norm(M_summary, axis=1)):.2f} A/m")

    return B_after, solved_M


# ================================================================
# Main
# ================================================================
print("=" * 70)
print("MMM IMA Diagnostic: Full vs Quarter+IMA")
print(f"Quarter grid: {len(tet_grid_q1)} tets in x=[1,10]mm, y=[-10,10]mm, z=[1,5]mm")
print("=" * 70)

print("\n--- Full model (all quadrant tets, no IMA) ---")
B_full, M_full_q1 = run_full_model()

print("\n--- Quarter model + IMA '+x-z' ---")
B_ima = run_quarter_ima()

print("\n--- Quarter model (no IMA, baseline) ---")
B_quarter = run_quarter_no_ima()

print("\n--- Quarter model + IMA + manual image field ---")
B_ima_fixed, M_ima = run_quarter_ima_with_image_field()

# Multi-element matrix check: compare IMA vs full model
# Sign matrices for '+x-z'
S_x_g = np.array([-1, 1, 1])
S_z_g = np.array([-1, -1, 1])
S_xz_g = np.array([1, -1, 1])

print("\n--- Multi-element matrix check ---")
rad.UtiDelAll()
n_q = len(tet_grid_q1)
full_tets_all = []
for mx in [False, True]:
    for mz in [False, True]:
        for verts in tet_grid_q1:
            full_tets_all.append(mirror_verts(verts, mx=mx, mz=mz))
full_objs_all = make_objs(full_tets_all)
full_cnt_all = rad.ObjCnt(full_objs_all)
h_multi_full = rad.BuildMatrix(full_cnt_all)
N_multi_full, ndof_mf = rad.GetInteractMatrix(h_multi_full)
print(f"Full multi: {ndof_mf} DOF")

rad.UtiDelAll()
ima_objs_all = make_objs(tet_grid_q1)
ima_cnt_all = rad.ObjCnt(ima_objs_all)
h_multi_ima = rad.BuildMatrix(ima_cnt_all, image='+x-z')
N_multi_ima, ndof_mi = rad.GetInteractMatrix(h_multi_ima)
print(f"IMA multi: {ndof_mi} DOF")

# Compare: for each tet i in q1, sum full model blocks for all 4 quadrant images of each tet j in q1
# Full model order: q1[0..71], q_zmirror[72..143], q_xmirror[144..215], q_xzmirror[216..287]
N_expected_multi = np.zeros((n_q*3, n_q*3))
for i in range(n_q):
    for j in range(n_q):
        # Direct: q1-q1 block
        r0, r1 = i*3, (i+1)*3
        c0, c1 = j*3, (j+1)*3

        block = N_multi_full[r0:r1, c0:c1].copy()

        # z-mirror block (q1 row, q_zmirror col)
        c0_z = (n_q + j) * 3
        c1_z = c0_z + 3
        b_z = N_multi_full[r0:r1, c0_z:c1_z]
        for b in range(3):
            block[:, b] += S_z_g[b] * b_z[:, b]

        # x-mirror block (q1 row, q_xmirror col)
        c0_x = (2*n_q + j) * 3
        c1_x = c0_x + 3
        b_x = N_multi_full[r0:r1, c0_x:c1_x]
        for b in range(3):
            block[:, b] += S_x_g[b] * b_x[:, b]

        # xz-mirror block (q1 row, q_xzmirror col)
        c0_xz = (3*n_q + j) * 3
        c1_xz = c0_xz + 3
        b_xz = N_multi_full[r0:r1, c0_xz:c1_xz]
        for b in range(3):
            block[:, b] += S_xz_g[b] * b_xz[:, b]

        N_expected_multi[r0:r1, c0:c1] = block

multi_err = np.max(np.abs(N_multi_ima - N_expected_multi))
multi_rel = multi_err / np.max(np.abs(N_multi_ima)) * 100
print(f"Multi-element matrix max|diff|: {multi_err:.6e}")
print(f"Multi-element matrix relative:  {multi_rel:.4f}%")
if multi_err < 1e-8:
    print("*** MULTI-ELEMENT MATRIX MATCH ***")
elif multi_rel < 0.1:
    print(f"*** MULTI-ELEMENT: small numerical difference ({multi_rel:.4f}%) ***")
else:
    print(f"*** MULTI-ELEMENT MATRIX MISMATCH: {multi_rel:.2f}% ***")

# Compare M values
M_full_arr = np.array(M_full_q1)
M_ima_arr = np.array(M_ima)
M_diff = M_ima_arr - M_full_arr
print(f"\n--- M value comparison (IMA vs Full q1) ---")
print(f"Max |M_IMA - M_full|: {np.max(np.abs(M_diff)):.4f} A/m")
print(f"Mean |M_IMA|: {np.mean(np.linalg.norm(M_ima_arr, axis=1)):.2f} A/m")
print(f"Relative M error: {np.max(np.abs(M_diff)) / np.mean(np.linalg.norm(M_ima_arr, axis=1)) * 100:.4f}%")

# Compare
print("\n" + "=" * 70)
print("COMPARISON")
print("=" * 70)
print(f"Full model Bz:        {B_full[2]*1e3:.4f} mT")
print(f"IMA quarter Bz:       {B_ima[2]*1e3:.4f} mT")
print(f"IMA + image Bz:       {B_ima_fixed[2]*1e3:.4f} mT")
print(f"Quarter no-IMA Bz:    {B_quarter[2]*1e3:.4f} mT")
print(f"Background only Bz:   100.0000 mT")

if abs(B_full[2]) > 1e-10:
    ratio = B_ima[2] / B_full[2]
    ratio_fixed = B_ima_fixed[2] / B_full[2]
    print(f"\nIMA/Full ratio (Fld only):      {ratio:.4f}")
    print(f"IMA+images/Full ratio:          {ratio_fixed:.4f} (should be ~1.0)")
    diff_pct = abs(ratio_fixed - 1.0) * 100
    if diff_pct < 5:
        print(f"*** PASS: IMA+images matches full model ({diff_pct:.1f}%) ***")
    else:
        print(f"*** FAIL: IMA+images differs by {diff_pct:.1f}% ***")

# ================================================================
# Matrix comparison: single tet
# ================================================================
print("\n" + "=" * 70)
print("MATRIX COMPARISON (single tet)")
print("=" * 70)

# Full model: 4 tets in all quadrants
rad.UtiDelAll()
full_tets = []
for mx in [False, True]:
    for mz in [False, True]:
        full_tets.append(mirror_verts(single_tet_q1[0], mx=mx, mz=mz))
full_objs = make_objs(full_tets)
full_cnt = rad.ObjCnt(full_objs)
h_full = rad.BuildMatrix(full_cnt)
N_full, ndof_full = rad.GetInteractMatrix(h_full)
print(f"Full model: {ndof_full} DOF, matrix {N_full.shape}")

# IMA quarter: 1 tet
rad.UtiDelAll()
ima_objs = make_objs(single_tet_q1)
ima_cnt = rad.ObjCnt(ima_objs)
h_ima = rad.BuildMatrix(ima_cnt, image='+x-z')
N_ima, ndof_ima = rad.GetInteractMatrix(h_ima)
print(f"IMA quarter: {ndof_ima} DOF, matrix {N_ima.shape}")

# Extract blocks from full matrix
# IMPORTANT: loop order is mx=[F,T], mz=[F,T]
# so: (F,F)=q1, (F,T)=q_zmirror, (T,F)=q_xmirror, (T,T)=q_xzmirror
# Corrected labels:
print(f"\nFull model N[q1, q1] (self):")
N_self = N_full[0:3, 0:3]
print(N_self)

print(f"\nFull model N[q1, q_zmir] (from z-mirror, block [0:3,3:6]):")
N_zmir = N_full[0:3, 3:6]
print(N_zmir)

print(f"\nFull model N[q1, q_xmir] (from x-mirror, block [0:3,6:9]):")
N_xmir = N_full[0:3, 6:9]
print(N_xmir)

print(f"\nFull model N[q1, q_xzmir] (from xz-mirror, block [0:3,9:12]):")
N_xzmir = N_full[0:3, 9:12]
print(N_xzmir)

# IMA matrix should equal: N_self + S_x*N_xmir + S_z*N_zmir + S_xz*N_xzmir
# where S_x, S_z, S_xz are the sign matrices for each mirror
print(f"\nIMA N_IMA:")
print(N_ima)

# What the full model gives if we just sum all blocks (no sign matrix):
N_sum = N_self + N_xmir + N_zmir + N_xzmir
print(f"\nSum (no sign matrix): N_self + N_xmir + N_zmir + N_xzmir:")
print(N_sum)

# Compute difference
print(f"\nDifference (N_IMA - N_sum):")
diff = N_ima - N_sum
print(diff)

# Now check: what sign matrix would make it match?
# N_IMA should = N_self + S_x @ N_xmir + S_z @ N_zmir + S_xz @ N_xzmir
# where S is diagonal sign matrix applied to columns (source M components)
# N_IMA[a,b] = N_self[a,b] + S_x[b]*N_xmir[a,b] + S_z[b]*N_zmir[a,b] + S_xz[b]*N_xzmir[a,b]
#
# But the full model doesn't use sign matrices - the physical M is self-consistent
# So the IMA matrix should match: N_self + N_xmir + N_zmir + N_xzmir
# This would only be true if S = diag(1,1,1) for all mirrors!
#
# Hmm, but that can't be right. The sign matrix accounts for the pseudovector
# transformation AND the IMA sign.
#
# The correct comparison is:
# IMA matrix = N_self + Σ_mirrors S_mirror * N_mirror_geom
# where N_mirror_geom uses MIRRORED source geometry
# But N_full[q1, q2] already uses the mirrored geometry (since q2 IS the x-mirror)
# So: N_IMA should = N_self + S_x · N_xmir + S_z · N_zmir + S_xz · N_xzmir
# where · means column-wise multiplication (S[b] * N[a,b])

# Current code signs (with CORRECT element ordering):
# +x: S = [-1, 1, 1]
# -z: S = [-1, -1, 1]
# +x-z: S = [1, -1, 1]
S_x = np.array([-1, 1, 1])
S_z = np.array([-1, -1, 1])
S_xz = np.array([1, -1, 1])

N_expected = N_self.copy()
for b in range(3):
    N_expected[:, b] += S_x[b] * N_xmir[:, b]
    N_expected[:, b] += S_z[b] * N_zmir[:, b]
    N_expected[:, b] += S_xz[b] * N_xzmir[:, b]

print(f"\nExpected (correct ordering + current signs):")
print(N_expected)
err = np.max(np.abs(N_ima - N_expected))
print(f"Match with IMA: max|diff| = {err:.6e}")
if err < 1e-8:
    print("*** MATRIX MATCH: Sign matrix is CORRECT ***")
else:
    print(f"*** MATRIX MISMATCH: error = {err:.3e} ***")
