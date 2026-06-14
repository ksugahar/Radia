#!/usr/bin/env python
"""
C-type dipole basic tests: MSC and FEM (Omega-reduced).

Both tests use:
  - Cubit hex mesh (same geometry as generate_hex_mesh.py)
  - CoilBuilder -> rad.Fld() for coil source field
  - 1/4 model (x>0, z>0 quarter)

MSC:  Radia Solve with IMA '+x-z'
FEM:  NGSolve Omega-reduced + Kelvin (1/4 sphere)

Reference (20000 AT, 6x6x6 quarter + IMA '+x-z'):
  Bz at origin ~ -961 mT

Requirements: Cubit 2025.8+, NGSolve, Radia
Usage:
    python test_occ_dipole.py msc [coarse|medium]
    python test_occ_dipole.py fem [coarse]
"""

import sys
import os
import math
import time

work_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.join(work_dir, '..', '..', '..')
sys.path.insert(0, os.path.join(repo_root, 'src'))
sys.path.insert(0, os.path.join(repo_root, 'src', 'radia'))

import numpy as np
import radia as rad

scale = 0.001  # mm -> m
MU_0 = 4e-7 * math.pi

BH_FILE = os.path.join(repo_root, 'examples', 'c_type_electromagnet',
                        'nonlinear', 'BH.txt')

MESH_CONFIGS = {
    'coarse': dict(main=2,  sec=1, third=1, factor=10),
    'medium': dict(main=6,  sec=4, third=2, factor=8),
    'fine':   dict(main=12, sec=8, third=4, factor=7),
}


# ================================================================
# Common: Cubit full model + CoilBuilder
# ================================================================

def build_full_model_cubit(resolution='coarse'):
    """Build full C-type hex mesh in Cubit (reflect + transform).

    Same as generate_hex_mesh.py. Returns cubit module.
    """
    # Auto-detect Cubit installation
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src', 'radia'))
from install_panels import find_cubit_bin as _fcb
_cubit_path = _fcb()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)
    if cubit_path not in sys.path:
        sys.path.insert(0, cubit_path)

    import contextlib, io
    with contextlib.redirect_stderr(io.StringIO()):
        import cubit
        cubit.init(['cubit', '-nojournal', '-nographics'])

    cfg = MESH_CONFIGS[resolution]
    cubit.cmd("reset")
    cubit.cmd("set journal off")

    # Quarter geometry (Trelis.jou)
    cubit.cmd("brick x 62.5 y 105 z 25")
    cubit.cmd("brick x 80 y 50 z 25")
    cubit.cmd("move Volume 2 location 71.25 -27.5 0")
    cubit.cmd("brick x 40 y 100 z 25")
    cubit.cmd("move Volume 3 location 131.25 -2.5 0")
    cubit.cmd("modify curve 35 26 36 chamfer radius 8")
    cubit.cmd("webcut volume 3 with plane yplane offset 39.5")
    cubit.cmd("webcut volume 3 with plane yplane offset 27.5")
    cubit.cmd("webcut volume 1 3 with plane yplane offset -2.5")
    cubit.cmd("imprint volume 7 3 2 1 6")
    cubit.cmd("merge volume 7 3 2 1 6")
    cubit.cmd("imprint volume 5 4")
    cubit.cmd("merge volume 5 4")

    # Hex mesh
    cubit.cmd(f"curve 62 41 42 60 47 32 interval {cfg['main']}")
    cubit.cmd("curve 62 41 42 60 47 32 scheme equal")
    cubit.cmd("mesh curve 62 41 42 60 47 32")
    cubit.cmd(f"curve 61 40 38 46 45 63 interval {cfg['sec']}")
    cubit.cmd("curve 61 40 38 46 45 63 scheme equal")
    cubit.cmd("mesh curve 61 40 38 46 45 63")
    cubit.cmd(f"curve 32 37 64 46 45 63 interval {cfg['third']}")
    cubit.cmd("curve 32 37 64 46 45 63 scheme equal")
    cubit.cmd("mesh curve 32 37 64 46 45 63")
    cubit.cmd(f"volume 5 4 size auto factor {cfg['factor']}")
    cubit.cmd("mesh volume 5 4")
    m = max(1, cfg['third'] // 2)
    cubit.cmd(f"curve 55 85 34 11 69 9 12 23 96 interval {m}")
    cubit.cmd("curve 55 85 34 11 69 9 12 23 96 scheme equal")
    cubit.cmd("mesh curve 55 85 34 11 69 9 12 23 96")
    cubit.cmd(f"curve 55 85 34 88 90 13 15 72 11 69 9 74 72 82 83 "
              f"12 23 80 81 8 2 70 4 16 14 20 interval {max(1, cfg['sec']//2)}")
    cubit.cmd("curve 55 85 34 88 90 13 15 72 11 69 9 74 72 82 83 "
              "12 23 80 81 8 2 70 4 16 14 20 scheme equal")
    cubit.cmd("mesh curve 55 85 34 88 90 13 15 72 11 69 9 74 72 82 83 "
              "12 23 80 81 8 2 70 4 16 14 20")
    cubit.cmd(f"volume 1 2 3 6 7 size auto factor {cfg['factor']}")
    cubit.cmd("mesh volume 1 2 3 6 7")
    cubit.cmd("imprint volume all")
    cubit.cmd("merge volume all")

    # Reflect + transform
    cubit.cmd("Volume all copy reflect vertex 55 53")
    cubit.cmd("Volume all copy reflect vertex 3 47")
    cubit.cmd("imprint volume all")
    cubit.cmd("merge volume all")
    cubit.cmd("move Volume all x -131.25 y -52.5 z 12.5 include_merged")
    cubit.cmd("rotate Volume all angle 90 about X include_merged")
    cubit.cmd("rotate Volume all angle -90 about Z include_merged")

    print(f"  Cubit full model: {cubit.get_hex_count()} hex")
    return cubit


def create_coil_radia(current_at=20000.0):
    """Create racetrack coil via CoilBuilder, return Radia container.

    Uses CoilBuilder (NOT ObjRecCur from coil_model.py).
    Same geometry as coil_model.py: center (0, 131.25, 0) mm.
    Coordinate system: transformed ELF (same as Cubit after transform).
    """
    from coil_builder import CoilBuilder

    # ELF coil geometry (mm -> m)
    mm = 1e-3
    Y_CENTER = 131.25 * mm     # coil center
    STRAIGHT_Y = 62.5 * mm     # Y-straight length
    R_MEAN = 22.5 * mm         # centerline arc radius
    WIDTH = 35 * mm            # cross-section radial
    HEIGHT = 105 * mm          # cross-section axial

    # Racetrack: 4 straights (2 Y-legs + 2 X-caps) + 4 quarter-circle arcs
    # Path: Y-up -> arc(90) -> X-cap(-x) -> arc(90) -> Y-down -> arc(90) -> X-cap(+x) -> arc(90)
    X_CAP = 50 * mm         # X-cap straight length
    Y_INNER = Y_CENTER - STRAIGHT_Y / 2   # y = 100 mm
    X_CENTER = 47.5 * mm    # inner Y-straight x position

    coil = (CoilBuilder(current=current_at)
        .set_start([X_CENTER, Y_INNER, 0])
        .set_cross_section(width=WIDTH, height=HEIGHT)
        .add_straight(STRAIGHT_Y)              # Y up:    (47.5, 100) -> (47.5, 162.5)
        .add_arc(radius=R_MEAN, arc_angle=90)  # top-right arc
        .add_straight(X_CAP)                   # X-cap:   top, -x direction
        .add_arc(radius=R_MEAN, arc_angle=90)  # top-left arc
        .add_straight(STRAIGHT_Y)              # Y down:  (-47.5, 162.5) -> (-47.5, 100)
        .add_arc(radius=R_MEAN, arc_angle=90)  # bottom-left arc
        .add_straight(X_CAP)                   # X-cap:   bottom, +x direction
        .add_arc(radius=R_MEAN, arc_angle=90))  # bottom-right arc (close)

    objs = coil.to_radia()
    container = rad.ObjCnt(objs)
    return container, coil


def extract_quarter_hex(cubit):
    """Extract x>0, z>0 hex elements from Cubit as (node_coords_m, hex_connectivity)."""
    all_hexes = []
    for vid in cubit.get_entities("volume"):
        for hid in cubit.get_volume_hexes(vid):
            nids = cubit.get_connectivity("hex", hid)
            coords = [cubit.get_nodal_coordinates(n) for n in nids]
            cx = sum(c[0] for c in coords) / 8
            cz = sum(c[2] for c in coords) / 8
            if cx < -0.01 or cz < -0.01:
                continue
            verts_m = [[c[0]*scale, c[1]*scale, c[2]*scale] for c in coords]
            all_hexes.append(verts_m)
    return all_hexes


# ================================================================
# MSC test: Radia Solve + IMA
# ================================================================

def run_msc(resolution='coarse', current_at=20000.0, solver=0):
    """MSC solve: Cubit hex quarter + CoilBuilder + IMA '+x-z'."""
    print(f"\n{'='*60}")
    print(f"  MSC: {resolution}, {current_at:.0f} AT")
    print(f"{'='*60}")

    cubit = build_full_model_cubit(resolution)
    quarter_hexes = extract_quarter_hex(cubit)

    rad.UtiDelAll()
    bh_data = np.loadtxt(BH_FILE, comments='#').tolist()
    mat = rad.MatSatIsoTab(bh_data)

    all_objs = []
    for verts in quarter_hexes:
        obj = rad.ObjHexahedron(verts, [0, 0, 0])
        rad.MatApl(obj, mat)
        all_objs.append(obj)
    yoke = rad.ObjCnt(all_objs)

    coil_container, _ = create_coil_radia(current_at)
    model = rad.ObjCnt([yoke, coil_container])

    n_elem = len(all_objs)
    print(f"  {n_elem} hex, {n_elem*6} DOF")

    rad.SolverConfig(newton_method=False, keep_magnetization=False)
    t0 = time.time()
    result = rad.Solve(model, 0.001, 100, solver, image='+x-z')
    t_solve = time.time() - t0

    B = rad.Fld(model, 'b', [0, 0, 0])
    print(f"  Solve: {t_solve:.2f}s, residual={result[1]:.2e}")
    print(f"  B = ({B[0]*1e3:.2f}, {B[1]*1e3:.2f}, {B[2]*1e3:.2f}) mT")
    print(f"  Bz = {B[2]*1e3:.2f} mT (ref: -961)")
    return B


# ================================================================
# FEM test: Omega-reduced + Kelvin (1/4 sphere)
# ================================================================


def _add_sym_blocks(cubit, air_vids, yoke_vids, kelvin_vids):
    """Add symmetry surface blocks for quarter C-type dipole (IMA '+x-z').

    For Omega-reduced formulation:
      sym_normal (z=0):  Dirichlet Omega=0  (B normal to face, '-z' antisymmetry)
      sym_tangential (x=0): natural BC       (B parallel to face, '+x' symmetry)

    Only sym_normal gets Dirichlet constraint. sym_tangential is for reference
    (A-formulation would use it for Dirichlet A x n = 0).
    """
    tol = 0.5  # mm tolerance
    z0_sids = []   # z=0 plane -> sym_normal
    x0_sids = []   # x=0 plane -> sym_tangential
    all_vids = list(air_vids) + list(yoke_vids) + list(kelvin_vids)
    for vid in all_vids:
        for sid in cubit.get_relatives("volume", vid, "surface"):
            adj = cubit.get_relatives("surface", sid, "volume")
            if len(adj) > 1:
                continue  # internal face, skip
            bb = cubit.get_bounding_box("surface", sid)
            # z=0 plane: zmin~0, zmax~0
            if abs(bb[6]) < tol and abs(bb[7]) < tol:
                z0_sids.append(sid)
            # x=0 plane: xmin~0, xmax~0
            elif abs(bb[0]) < tol and abs(bb[1]) < tol:
                x0_sids.append(sid)

    z0_sids = list(set(z0_sids))
    x0_sids = list(set(x0_sids))

    bid = max(cubit.get_block_id_list()) + 1
    if z0_sids:
        for sid in z0_sids:
            tris = cubit.get_surface_tris(sid)
            quads = cubit.get_surface_quads(sid)
            if tris:
                cubit.cmd(f"block {bid} add tri in surface {sid}")
            if quads:
                cubit.cmd(f"block {bid} add face in surface {sid}")
        cubit.cmd(f'block {bid} name "sym_normal"')
        bid += 1
        print(f"  sym_normal (z=0): {len(z0_sids)} surfaces")
    if x0_sids:
        for sid in x0_sids:
            tris = cubit.get_surface_tris(sid)
            quads = cubit.get_surface_quads(sid)
            if tris:
                cubit.cmd(f"block {bid} add tri in surface {sid}")
            if quads:
                cubit.cmd(f"block {bid} add face in surface {sid}")
        cubit.cmd(f'block {bid} name "sym_tangential"')
        print(f"  sym_tangential (x=0): {len(x0_sids)} surfaces")
    return z0_sids, x0_sids


def run_fem(resolution='coarse', current_at=20000.0, fes_order=1,
            nonlinear=False, max_iter=20, tol=1e-3, relax=0.3):
    """FEM Omega-reduced: Cubit hex yoke + tet air + 1/4 Kelvin sphere.

    Key: mesh air FIRST, then create kelvin and copy mesh (preserves topology).

    Formulation (unified reduced potential):
      H = Hs - grad(Omega)
      B = mu * (Hs - grad(Omega))
      Weak: mu*grad(u)*grad(v) dx = mu*Hs*grad(v) dx  (all regions)
      Kelvin weight: mu -> mu*(a/r')^2 in exterior domain.

    Args:
        fes_order: FE polynomial order (1 or 2)
        nonlinear: Use BH curve with Picard iteration
        max_iter: Max Picard iterations
        tol: Convergence tolerance for relative mu change
        relax: Under-relaxation factor (0 = no update, 1 = full update)
    """
    print(f"\n{'='*60}")
    print(f"  FEM (Omega-reduced, Cubit + Kelvin): {current_at:.0f} AT")
    print(f"{'='*60}")

    from ngsolve import (Mesh, H1, Periodic, BilinearForm, LinearForm,
                         GridFunction, grad, dx, CF, VOL, TaskManager)
    from ngsolve import x as ngx, y as ngy, z as ngz

    # 1. Coil -> Radia Biot-Savart source
    rad.UtiDelAll()
    coil_container, _ = create_coil_radia(current_at)
    Hs = rad.RadiaField(coil_container, 'h')

    # 2. Build full hex model in Cubit (then use quarter only)
    cubit = build_full_model_cubit(resolution)
    cubit.cmd("set journal off")

    # Filter yoke to quarter model (x>0, z>0) -- same as MSC extract_quarter_hex
    all_yoke_vids = sorted(v for v in cubit.get_entities("volume")
                           if len(cubit.get_volume_hexes(v)) > 0)

    yoke_vids = []
    for vid in all_yoke_vids:
        bb = cubit.get_bounding_box("volume", vid)
        cx = (bb[0] + bb[1]) / 2
        cz = (bb[6] + bb[7]) / 2
        if cx >= -0.01 and cz >= -0.01:
            yoke_vids.append(vid)
        else:
            # Delete non-quarter yoke volumes to avoid mesh export
            cubit.cmd(f"delete volume {vid}")

    # Bounding box of yoke (quarter)
    bb_min = [1e30, 1e30, 1e30]
    bb_max = [-1e30, -1e30, -1e30]
    for vid in yoke_vids:
        bb = cubit.get_bounding_box("volume", vid)
        bb_min[0] = min(bb_min[0], bb[0]); bb_max[0] = max(bb_max[0], bb[1])
        bb_min[1] = min(bb_min[1], bb[3]); bb_max[1] = max(bb_max[1], bb[4])
        bb_min[2] = min(bb_min[2], bb[6]); bb_max[2] = max(bb_max[2], bb[7])

    sx, sy, sz = 0.0, (bb_min[1] + bb_max[1]) / 2, 0.0
    corners = [(x, y, z) for x in (bb_min[0], bb_max[0])
                          for y in (bb_min[1], bb_max[1])
                          for z in (bb_min[2], bb_max[2])]
    max_dist = max(math.sqrt((x-sx)**2 + (y-sy)**2 + (z-sz)**2)
                   for x, y, z in corners)
    R_kelvin = max_dist * 1.3
    offset_x = R_kelvin * 10

    print(f"  Yoke: {len(yoke_vids)} vols, R={R_kelvin:.1f} mm, offset={offset_x:.0f} mm")

    # ================================================================
    # 3. Air sphere: create, webcut, imprint with yoke, tet mesh
    # ================================================================
    cubit.cmd(f"create sphere radius {R_kelvin}")
    int_sphere = cubit.get_last_id("volume")
    if abs(sy) > 0.01:
        cubit.cmd(f"move volume {int_sphere} y {sy}")

    # Webcut for quarter (z=0, then x=0)
    cubit.cmd(f"webcut volume {int_sphere} with plane zplane offset 0")
    new_vids = [v for v in cubit.get_entities("volume")
                if v not in yoke_vids and len(cubit.get_volume_hexes(v)) == 0]
    for vid in new_vids:
        cubit.cmd(f"webcut volume {vid} with plane xplane offset 0")

    # Delete x<0 and z<0 pieces
    to_delete = []
    for vid in cubit.get_entities("volume"):
        if vid in yoke_vids or len(cubit.get_volume_hexes(vid)) > 0:
            continue
        bb = cubit.get_bounding_box("volume", vid)
        if (bb[0]+bb[1])/2 < -0.01 or (bb[6]+bb[7])/2 < -0.01:
            to_delete.append(vid)
    for vid in to_delete:
        cubit.cmd(f"delete volume {vid}")

    # Create air_gap box around pole gap (same idea as Keiko's SAir)
    # Gap = z=[0, 5]mm (quarter), pole width x=[0, 25]mm, y=[-20, 20]mm
    # Extend slightly beyond to capture fringe field
    gap_margin = 5.0  # mm extra around gap
    cubit.cmd(f"create brick x 42 y 73 z 13")  # slightly larger than gap region
    gap_box_id = cubit.get_last_id("volume")
    cubit.cmd(f"move volume {gap_box_id} x 21 y 0 z 6.5")

    # Imprint/merge everything
    cubit.cmd("imprint volume all")
    cubit.cmd("merge volume all")

    # Find air volumes: separate air_gap (small, near yoke) from air (far)
    non_yoke = sorted(v for v in cubit.get_entities("volume")
                      if v not in yoke_vids and len(cubit.get_volume_hexes(v)) == 0)

    # Classify: volumes touching yoke AND inside gap region -> air_gap
    air_gap_vids = []
    air_vids = []
    for vid in non_yoke:
        bb = cubit.get_bounding_box("volume", vid)
        vol_center_z = (bb[6] + bb[7]) / 2
        vol_center_y = (bb[3] + bb[4]) / 2
        vol_size = max(bb[1]-bb[0], bb[4]-bb[3], bb[7]-bb[6])
        # Small volume near z=0 gap (z < 15mm) and near yoke (y near 0)
        if vol_size < 50 and vol_center_z < 15 and abs(vol_center_y) < 40:
            air_gap_vids.append(vid)
        else:
            air_vids.append(vid)

    air_mesh_size = R_kelvin / 5  # ~37mm for outer air
    gap_mesh_size = 2.0           # 2mm for air_gap (resolves 5mm gap)
    kelvin_mesh_size = R_kelvin / 3

    # Mesh air_gap first (fine)
    for av in air_gap_vids:
        cubit.cmd(f"volume {av} scheme tetmesh")
        cubit.cmd(f"volume {av} size {gap_mesh_size}")
        cubit.cmd(f"mesh volume {av}")

    # Then mesh outer air (coarse)
    for av in air_vids:
        cubit.cmd(f"volume {av} scheme tetmesh")
        cubit.cmd(f"volume {av} size {air_mesh_size}")
        cubit.cmd(f"mesh volume {av}")

    all_air = air_gap_vids + air_vids

    n_gap = sum(len(cubit.get_volume_tets(v)) for v in air_gap_vids)
    n_air = sum(len(cubit.get_volume_tets(v)) for v in air_vids)
    print(f"  Air: {len(air_vids)} vols ({n_air} tet), air_gap: {len(air_gap_vids)} vols ({n_gap} tet)")

    # ================================================================
    # 4. Kelvin sphere: create at offset, webcut, copy surface mesh
    # ================================================================
    # Track all existing volumes
    all_known = set(cubit.get_entities("volume"))

    cubit.cmd(f"create sphere radius {R_kelvin}")
    ext_sphere = cubit.get_last_id("volume")
    cubit.cmd(f"move volume {ext_sphere} x {offset_x} y {sy}")

    # Webcut exterior: z=0, then x=offset_x
    cubit.cmd(f"webcut volume {ext_sphere} with plane zplane offset 0")
    ext_new = [v for v in cubit.get_entities("volume") if v not in all_known]
    if ext_new:
        cubit.cmd(f"webcut volume {' '.join(str(v) for v in ext_new)} "
                  f"with plane xplane offset {offset_x}")

    # Delete non-quarter pieces
    to_delete = []
    for vid in cubit.get_entities("volume"):
        if vid in all_known:
            continue
        bb = cubit.get_bounding_box("volume", vid)
        nc_x = (bb[0]+bb[1])/2
        nc_z = (bb[6]+bb[7])/2
        if nc_x < offset_x - 0.01 or nc_z < -0.01:
            to_delete.append(vid)
    for vid in to_delete:
        cubit.cmd(f"delete volume {vid}")

    kelvin_vids = sorted(v for v in cubit.get_entities("volume")
                         if v not in all_known)
    print(f"  Kelvin vols: {len(kelvin_vids)}")

    # ================================================================
    # 5. Find hemisphere surfaces, copy mesh from air to kelvin
    # ================================================================
    def _find_sphere_surface(vol_ids, cx, cy, cz, R):
        """Find the quarter-sphere (curved) surface on given volumes.

        Distinguishes from flat faces by checking bounding box span:
        a curved hemisphere spans a range in ALL 3 coordinates,
        while flat faces (z=0, x=0) have zero span in one axis.
        """
        found = []
        for vid in vol_ids:
            for sid in cubit.get_relatives("volume", vid, "surface"):
                verts = cubit.get_relatives("surface", sid, "vertex")
                if not verts or len(verts) < 2:
                    continue
                # Check vertices on sphere
                on_sphere = True
                for v in verts:
                    c = cubit.vertex(v).coordinates()
                    d = math.sqrt((c[0]-cx)**2 + (c[1]-cy)**2 + (c[2]-cz)**2)
                    if abs(d - R) > R * 0.05:
                        on_sphere = False
                        break
                if not on_sphere:
                    continue
                # Reject flat surfaces (zero span in one axis)
                bb = cubit.get_bounding_box("surface", sid)
                dx = bb[1] - bb[0]  # x span
                dy = bb[4] - bb[3]  # y span
                dz = bb[7] - bb[6]  # z span
                if dx < R * 0.1 or dy < R * 0.1 or dz < R * 0.1:
                    continue  # flat face, skip
                found.append(sid)
        return found

    int_hemis = _find_sphere_surface(all_air, 0, sy, 0, R_kelvin)
    ext_hemis = _find_sphere_surface(kelvin_vids, offset_x, sy, 0, R_kelvin)
    print(f"  Sphere surfaces: {len(int_hemis)} int, {len(ext_hemis)} ext")

    # Interior hemisphere: already meshed (from air tet mesh)
    # Copy to each exterior hemisphere surface
    copy_ok = True
    for src_sid, dst_sid in zip(int_hemis, ext_hemis):
        src_curves = cubit.get_relatives("surface", src_sid, "curve")
        dst_curves = cubit.get_relatives("surface", dst_sid, "curve")
        if not src_curves or not dst_curves:
            copy_ok = False
            continue
        src_c = src_curves[0]
        dst_c = dst_curves[0]
        try:
            src_v = cubit.get_relatives("curve", src_c, "vertex")[0]
            dst_v = cubit.get_relatives("curve", dst_c, "vertex")[0]
            cubit.cmd(f"copy mesh surface {src_sid} onto surface {dst_sid} "
                      f"source curve {src_c} source vertex {src_v} "
                      f"target curve {dst_c} target vertex {dst_v}")
            if not cubit.get_surface_tris(dst_sid):
                copy_ok = False
        except Exception as e:
            print(f"  Copy mesh error: {e}")
            copy_ok = False
    print(f"  Copy mesh: {'OK' if copy_ok else 'FAILED'}")

    # Tet mesh kelvin volumes (surface mesh constrains boundary)
    for kv in kelvin_vids:
        cubit.cmd(f"volume {kv} scheme tetmesh")
        cubit.cmd(f"volume {kv} size {kelvin_mesh_size}")
        cubit.cmd(f"mesh volume {kv}")

    n_kel = sum(len(cubit.get_volume_tets(v)) for v in kelvin_vids)
    print(f"  Kelvin: {n_kel} tets")

    # ================================================================
    # 6. Blocks
    # ================================================================
    cubit.cmd("set duplicate block elements on")
    bid = 1
    for v in yoke_vids:
        cubit.cmd(f"block {bid} add volume {v}")
    cubit.cmd(f'block {bid} name "yoke"')
    bid += 1
    for v in air_gap_vids:
        cubit.cmd(f"block {bid} add volume {v}")
    cubit.cmd(f'block {bid} name "air_gap"')
    bid += 1
    for v in air_vids:
        cubit.cmd(f"block {bid} add volume {v}")
    cubit.cmd(f'block {bid} name "air"')
    bid += 1
    for v in kelvin_vids:
        cubit.cmd(f"block {bid} add volume {v}")
    cubit.cmd(f'block {bid} name "kelvin"')

    # kelvin_int / kelvin_ext blocks
    # Re-find hemisphere surfaces (may have new IDs after imprint/merge)
    int_hemis2 = _find_sphere_surface(all_air, 0, sy, 0, R_kelvin)
    ext_hemis2 = _find_sphere_surface(kelvin_vids, offset_x, sy, 0, R_kelvin)
    bid += 1
    for sid in int_hemis2:
        tris = cubit.get_surface_tris(sid)
        if tris:
            cubit.cmd(f"block {bid} add tri in surface {sid}")
    cubit.cmd(f'block {bid} name "kelvin_int"')
    bid += 1
    for sid in ext_hemis2:
        tris = cubit.get_surface_tris(sid)
        if tris:
            cubit.cmd(f"block {bid} add tri in surface {sid}")
    cubit.cmd(f'block {bid} name "kelvin_ext"')

    # Symmetry blocks
    _add_sym_blocks(cubit, all_air, yoke_vids, kelvin_vids)

    print(f"  Total: hex={cubit.get_hex_count()}, tet={cubit.get_tet_count()}")

    # ================================================================
    # 7. Export to NGSolve
    # ================================================================
    radia_src = os.path.join(repo_root, 'src', 'radia')
    if radia_src not in sys.path:
        sys.path.insert(0, radia_src)
    from cubit_mesh_export import extract_curved_mesh
    mesh = Mesh(extract_curved_mesh(cubit, order=1))

    ngmesh = mesh.ngmesh
    ngmesh.Scale(scale)
    mesh = Mesh(ngmesh)

    materials = mesh.GetMaterials()
    boundaries = mesh.GetBoundaries()
    print(f"  NGSolve: ne={mesh.ne}, mat={set(materials)}")

    # 8. Kelvin periodic BC
    panels_dir = os.path.join(repo_root, 'src', 'radia', 'panels')
    if panels_dir not in sys.path:
        sys.path.insert(0, panels_dir)
    from calc_common import add_periodic_kelvin
    # Offset = TRANSLATION from int sphere to ext sphere (NOT ext center!)
    # int sphere at (0, sy, 0), ext at (offset_x, sy, 0) -> translation = (offset_x, 0, 0)
    kelvin_center = (offset_x * scale, 0.0, 0.0)
    has_kelvin = add_periodic_kelvin(mesh, kelvin_center)
    if has_kelvin:
        mesh = Mesh(mesh.ngmesh)
        materials = mesh.GetMaterials()
        boundaries = mesh.GetBoundaries()
        print(f"  Kelvin periodic: offset={[f'{c:.4f}' for c in kelvin_center]}")
    else:
        print("  WARNING: Kelvin periodic BC failed")

    # 9. BH curve and material setup
    from ngsolve import L2, InnerProduct, sqrt
    R_m = R_kelvin * scale
    # Kelvin weight center = exterior sphere center (NOT the translation vector!)
    kx_w = offset_x * scale  # exterior sphere x-center
    ky_w = sy * scale         # exterior sphere y-center (same as interior)
    kz_w = 0.0                # exterior sphere z-center
    rp_sq = (ngx - kx_w)**2 + (ngy - ky_w)**2 + (ngz - kz_w)**2 + 1e-20
    kelvin_fac = R_m**2 / rp_sq

    # Load BH curve if nonlinear
    mu_r_init = 1000.0
    bh_interp = None
    if nonlinear and os.path.exists(BH_FILE):
        from scipy.interpolate import interp1d
        bh_data = np.loadtxt(BH_FILE, comments='#')
        H_arr, B_arr = bh_data[:, 0], bh_data[:, 1]
        mu_arr = np.where(H_arr > 1e-10, B_arr / H_arr, B_arr[1] / max(H_arr[1], 1e-10))
        bh_interp = interp1d(H_arr, mu_arr, fill_value=(mu_arr[0], mu_arr[-1]),
                             bounds_error=False)
        mu_r_init = float(mu_arr[0] / MU_0)
        print(f"  BH curve: {len(H_arr)} points, mu_r_init={mu_r_init:.0f}")

    # Per-element mu (L2 order=0), Kelvin-weighted for bilinear form
    fes_l2 = L2(mesh, order=0)
    gf_mu = GridFunction(fes_l2)

    # Element material mapping
    elem_mat = {}
    for el in mesh.Elements(VOL):
        elem_mat[el.nr] = materials[el.index]

    yoke_elems = [nr for nr, m in elem_mat.items() if "yoke" in m.lower()]

    def _set_mu(mu_r_yoke):
        """Set per-element mu values."""
        for nr, m in elem_mat.items():
            if "yoke" in m.lower():
                gf_mu.vec[nr] = mu_r_yoke * MU_0 if isinstance(mu_r_yoke, (int, float)) else MU_0
            elif "kelvin" in m.lower():
                gf_mu.vec[nr] = MU_0  # Kelvin weight applied separately
            else:
                gf_mu.vec[nr] = MU_0

    _set_mu(mu_r_init)

    # Kelvin-weighted mu CF
    kw_dict = {}
    for m in set(materials):
        kw_dict[m] = kelvin_fac if "kelvin" in m.lower() else CF(1.0)
    kelvin_weight = mesh.MaterialCF(kw_dict, default=CF(1.0))
    Mu = gf_mu * kelvin_weight

    # 10. Dirichlet BC
    dir_parts = []
    if "sym_normal" in boundaries:
        dir_parts.append("sym_normal")
    if "GND" in boundaries:
        dir_parts.append("GND")
    dirichlet_bnd = "|".join(dir_parts) if dir_parts else ""
    print(f"  Dirichlet: {dirichlet_bnd or '(none)'}")

    # 11. FE space
    fes_base = H1(mesh, order=fes_order, dirichlet=dirichlet_bnd)
    fes = Periodic(fes_base) if has_kelvin else fes_base
    u, v = fes.TnT()

    freedof_before = sum(1 for d in fes_base.FreeDofs() if d)
    freedof_after = sum(1 for d in fes.FreeDofs() if d)
    if has_kelvin:
        print(f"  Periodic: free {freedof_before} -> {freedof_after}")

    gfu = GridFunction(fes)
    print(f"  FES: order={fes_order}, ndof={fes.ndof}")

    # 12. Solve (Picard iteration if nonlinear)
    t0 = time.time()
    n_iter = 1 if not nonlinear or bh_interp is None else max_iter

    for iteration in range(n_iter):
        a = BilinearForm(fes)
        a += Mu * grad(u) * grad(v) * dx(bonus_intorder=4)
        a.Assemble()

        f_lf = LinearForm(fes)
        f_lf += Mu * Hs * grad(v) * dx(bonus_intorder=4)
        f_lf.Assemble()

        with TaskManager():
            gfu.vec.data = a.mat.Inverse(fes.FreeDofs()) * f_lf.vec

        if not nonlinear or bh_interp is None:
            break

        # Update mu from BH curve (Picard: chord permeability)
        H_total = Hs - grad(gfu)
        gf_Hmag = GridFunction(fes_l2)
        gf_Hmag.Set(sqrt(InnerProduct(H_total, H_total)))

        max_change = 0.0
        for nr in yoke_elems:
            H_val = max(float(gf_Hmag.vec[nr]), 1e-10)
            mu_new = float(bh_interp(H_val))
            mu_old = float(gf_mu.vec[nr])
            mu_upd = relax * mu_new + (1 - relax) * mu_old
            change = abs(mu_upd - mu_old) / max(abs(mu_old), 1e-30)
            max_change = max(max_change, change)
            gf_mu.vec[nr] = mu_upd

        print(f"  Picard iter {iteration}: max_change={max_change:.4e}")
        if max_change < tol and iteration > 0:
            print(f"  Converged at iter {iteration}")
            break

    t_solve = time.time() - t0
    print(f"  Solve: {t_solve:.2f}s")

    # 13. Evaluate B in physical domain (mu without Kelvin weight)
    mu_phys_dict = {}
    for m in set(materials):
        if "yoke" in m.lower():
            mu_phys_dict[m] = gf_mu  # uses current per-element mu
        else:
            mu_phys_dict[m] = MU_0
    # Simpler: just use gf_mu directly (no Kelvin weight for physical evaluation)
    Mu_phys = gf_mu

    H_total = Hs - grad(gfu)
    B_total = Mu_phys * H_total

    eval_pt = (0.0125, 0.0, 0.002)
    B_coil = rad.Fld(coil_container, 'b', list(eval_pt))
    print(f"  B_coil = ({B_coil[0]*1e3:.1f}, {B_coil[1]*1e3:.1f}, {B_coil[2]*1e3:.1f}) mT")

    B = None
    for pt in [eval_pt, (0.010, 0.0, 0.0025), (0.005, 0.0, 0.001)]:
        try:
            B = B_total(mesh(*pt))
            print(f"  B_fem at {pt} = ({B[0]*1e3:.1f}, {B[1]*1e3:.1f}, {B[2]*1e3:.1f}) mT")
            print(f"  Bz = {B[2]*1e3:.1f} mT (ref: -961)")
            if abs(B_coil[2]) > 1e-6:
                print(f"  Iron amplification: {B[2]/B_coil[2]:.1f}x")
            break
        except Exception as e:
            print(f"  Point {pt}: {e}")
            continue

    if B is None:
        print("  No valid air point found")
    return B


# ================================================================
# Main
# ================================================================

if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'msc'
    res = sys.argv[2] if len(sys.argv) > 2 else 'coarse'

    if mode == 'msc':
        run_msc(res, current_at=20000, solver=0)
    elif mode == 'fem':
        run_fem(res, current_at=20000)
    elif mode == 'fem1nl':
        # Order=1 + nonlinear BH (fast convergence test)
        run_fem(res, current_at=20000, fes_order=1, nonlinear=True,
                max_iter=80, relax=0.8, tol=0.01)
    elif mode == 'fem2':
        # Higher order + nonlinear BH
        run_fem(res, current_at=20000, fes_order=2, nonlinear=True,
                max_iter=80, relax=0.8, tol=0.01)
    elif mode == 'both':
        run_msc(res, current_at=20000, solver=0)
        run_fem(res, current_at=20000)
    elif mode == 'full':
        run_msc(res, current_at=20000, solver=0)
        run_fem(res, current_at=20000, fes_order=2, nonlinear=True)

    print(f"\n{'='*60}")
    print("DONE")
    print(f"{'='*60}")
