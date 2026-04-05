"""
Cubit → Netgen bridge: extract mesh via cubit Python API, curve via C++.

This module replaces the old radia_cubit_mesh.extract_curved_mesh() which
called Cubit C++ APIs (MeshExportInterface, GeometryQueryTool) that hang
when called from external Python.

Instead, mesh data and geometry callbacks are extracted via cubit's Python
API (which works from external Python 3.12), then passed to the C++ pybind11
module for high-order curving via Netgen's BuildCurvedElements.

Usage:
    import netgen  # MUST import before cubit to avoid DLL conflicts
    import cubit
    from cubit_mesh_export import extract_curved_mesh

    cubit.init(['cubit', '-nojournal', '-batch', '-nographics'])
    cubit.cmd('open "model.cub5"')
    ng_mesh = extract_curved_mesh(cubit, order=3)
    mesh = ngsolve.Mesh(ng_mesh)
"""

import numpy as np


def extract_curved_mesh(cubit_mod, order=2, surface_only=False, split_quads=False):
    """Extract a curved Netgen mesh from Cubit via Python API + C++ curving.

    Args:
        cubit_mod: The cubit module (already initialized with cubit.init()).
        order: Polynomial order for curving (>= 2).
        surface_only: Remove volume elements (BEM workflow).
        split_quads: Split quad faces into triangles.

    Returns:
        netgen.meshing.Mesh: Curved mesh ready for ngsolve.Mesh().
    """
    import netgen  # noqa: F401 — must import before radia_cubit_mesh (DLL setup)
    try:
        # v0.2.0+: .pyd bundled inside cubit_mesh_export package
        from cubit_mesh_export import radia_cubit_mesh
    except ImportError:
        try:
            # Fallback: .pyd installed separately (e.g., radia package)
            import radia_cubit_mesh
        except ImportError:
            # Fallback: .pyd in same directory
            import os, sys
            pkg_dir = os.path.dirname(os.path.abspath(__file__))
            if pkg_dir not in sys.path:
                sys.path.insert(0, pkg_dir)
            import radia_cubit_mesh

    # --- Phase A: Extract nodes ---
    node_id_list = list(cubit_mod.parse_cubit_list("node", "all"))
    N = len(node_id_list)
    node_coords = np.empty((N, 3), dtype=np.float64)
    node_ids = np.array(node_id_list, dtype=np.int32)
    for i, nid in enumerate(node_id_list):
        node_coords[i] = cubit_mod.get_nodal_coordinates(nid)

    # --- Phase B: Extract volume elements (with per-volume domain index) ---
    vol_elements = []
    volume_ids = list(cubit_mod.get_entities("volume"))

    tet_ids = cubit_mod.parse_cubit_list("tet", "all")
    hex_ids = cubit_mod.parse_cubit_list("hex", "all")

    for domain_idx, vid in enumerate(volume_ids, start=1):
        for tid in cubit_mod.get_volume_tets(vid):
            conn = list(cubit_mod.get_connectivity("tet", tid))
            vol_elements.append((0, conn, domain_idx))  # TET
        for hid in cubit_mod.get_volume_hexes(vid):
            conn = list(cubit_mod.get_connectivity("hex", hid))
            vol_elements.append((1, conn, domain_idx))  # HEX
        # Wedges and pyramids (hex-tet transition)
        try:
            for wid in cubit_mod.get_volume_wedges(vid):
                conn = list(cubit_mod.get_connectivity("wedge", wid))
                vol_elements.append((2, conn, domain_idx))  # PRISM/WEDGE
        except Exception:
            pass
        try:
            for pid in cubit_mod.get_volume_pyramids(vid):
                conn = list(cubit_mod.get_connectivity("pyramid", pid))
                vol_elements.append((3, conn, domain_idx))  # PYRAMID
        except Exception:
            pass

    # --- Phase C: Extract surface mesh + UVs ---
    surface_ids = list(cubit_mod.parse_cubit_list("surface", "all"))
    surface_data = []
    surface_objects = {}  # cache cubit.surface() objects

    for sid in surface_ids:
        surf = cubit_mod.surface(sid)
        surface_objects[sid] = surf

        tris = []
        quads = []

        # Always check both tri and face — mixed meshes (tet+wedge, hex+tet)
        # have both element types on surfaces.
        tri_ids_on_surf = cubit_mod.parse_cubit_list("tri", f"in surface {sid}")
        for tid in tri_ids_on_surf:
            tris.extend(cubit_mod.get_connectivity("tri", tid))

        face_ids_on_surf = cubit_mod.parse_cubit_list("face", f"in surface {sid}")
        for fid in face_ids_on_surf:
            quads.extend(cubit_mod.get_connectivity("face", fid))

        # UV coordinates for surface nodes
        surf_node_set = set(tris) | set(quads)
        uvs = {}
        for nid in surf_node_set:
            coords = cubit_mod.get_nodal_coordinates(nid)
            u, v = surf.u_v_from_position(coords)
            uvs[nid] = (u, v)

        surface_data.append({
            "id": sid,
            "tris": list(tris),
            "quads": list(quads),
            "uvs": uvs,
        })

    # --- Phase D: Build projection/normal callbacks ---
    # surfnr is 1-based FaceDescriptor index; map to Cubit surface ID
    def project_func(surfnr, x, y, z, u_hint, v_hint):
        sid = surface_ids[surfnr - 1]
        surf = surface_objects[sid]
        xp, yp, zp = surf.closest_point_trimmed((x, y, z))
        u, v = surf.u_v_from_position((xp, yp, zp))
        return (xp, yp, zp, u, v)

    def normal_func(surfnr, x, y, z):
        sid = surface_ids[surfnr - 1]
        surf = surface_objects[sid]
        nx, ny, nz = surf.normal_at((x, y, z))
        return (nx, ny, nz)

    # --- Phase D2: Build edge (curve) projection callback ---
    # Map (surfnr1, surfnr2) -> Cubit curve ID for edge projection
    curve_ids = list(cubit_mod.parse_cubit_list("curve", "all"))
    edge_map = {}  # (surfnr1, surfnr2) -> cubit curve object
    for cid in curve_ids:
        parent_surfs = list(cubit_mod.parse_cubit_list("surface", f"in curve {cid}"))
        if len(parent_surfs) >= 2:
            # Map surface pair (as FD indices) to curve
            for i in range(len(parent_surfs)):
                for j in range(i + 1, len(parent_surfs)):
                    si = parent_surfs[i]
                    sj = parent_surfs[j]
                    # FD index = 1-based position in surface_ids list
                    if si in surface_ids and sj in surface_ids:
                        fi = surface_ids.index(si) + 1
                        fj = surface_ids.index(sj) + 1
                        crv = cubit_mod.curve(cid)
                        edge_map[(fi, fj)] = crv
                        edge_map[(fj, fi)] = crv

    # --- Phase D3: Extract edge (segment) elements on curves ---
    # Each segment lies on a curve shared by two surfaces.
    # Netgen needs these as 1D elements for BuildCurvedElements edge snapping.
    edge_elements = []  # list of dicts: {surfnr1, surfnr2, segments: [[n0,n1], ...]}
    edge_curve_ids = []  # Cubit curve ID for each edge_element (edgenr = index + 1)
    for cid in curve_ids:
        parent_surfs = list(cubit_mod.parse_cubit_list("surface", f"in curve {cid}"))
        if len(parent_surfs) < 2:
            continue
        si, sj = parent_surfs[0], parent_surfs[1]
        if si not in surface_ids or sj not in surface_ids:
            continue
        fi = surface_ids.index(si) + 1  # 1-based FD index
        fj = surface_ids.index(sj) + 1

        crv_obj = cubit_mod.curve(cid)
        crv_length = crv_obj.length()
        edge_ids_on_curve = cubit_mod.parse_cubit_list("edge", f"in curve {cid}")
        segs = []
        dists = []
        for eid in edge_ids_on_curve:
            conn = list(cubit_mod.get_connectivity("edge", eid))
            if len(conn) >= 2:
                segs.append(conn[:2])
                # Normalized curve parameter (0~1) for each endpoint
                c0 = cubit_mod.get_nodal_coordinates(conn[0])
                c1 = cubit_mod.get_nodal_coordinates(conn[1])
                u0 = crv_obj.u_from_position(c0)
                u1 = crv_obj.u_from_position(c1)
                # Normalize to [0, 1] using arc length fraction
                # length_from_u(u_start, u) returns arc length between parameters
                u_start = crv_obj.u_from_position(
                    crv_obj.position_from_fraction(0))
                arc0 = crv_obj.length_from_u(u_start, u0)
                arc1 = crv_obj.length_from_u(u_start, u1)
                d0 = arc0 / crv_length if crv_length > 0 else 0
                d1 = arc1 / crv_length if crv_length > 0 else 0
                dists.append([d0, d1])

        if segs:
            edge_elements.append({
                "surfnr1": fi,
                "surfnr2": fj,
                "segments": segs,
                "dists": dists,
            })
            edge_curve_ids.append(cid)

    def edge_project_func(surfnr1, surfnr2, x, y, z):
        key = (surfnr1, surfnr2)
        crv = edge_map.get(key)
        if crv is not None:
            xp, yp, zp = crv.closest_point_trimmed((x, y, z))
            return (xp, yp, zp)
        # Fallback: project onto surfnr1
        sid = surface_ids[surfnr1 - 1]
        surf = surface_objects[sid]
        xp, yp, zp = surf.closest_point_trimmed((x, y, z))
        return (xp, yp, zp)

    # --- Phase E: Call C++ for high-order curving ---
    # For surface_only: pass empty volume elements to avoid
    # ClearVolumeElements() segfault in NGSolve Mesh() constructor.
    ve = [] if surface_only else vol_elements
    ng_mesh = radia_cubit_mesh.build_curved_mesh(
        node_coords, node_ids, ve, surface_data,
        project_func, normal_func,
        order=order, surface_only=False, split_quads=split_quads,
        edge_project_func=edge_project_func if edge_map else None,
        edge_elements=edge_elements if edge_elements else [],
    )

    # --- Phase F: Set labels from Cubit blocks ---
    _set_labels(ng_mesh, cubit_mod, surface_ids, edge_curve_ids)

    return ng_mesh


def _set_labels(ng_mesh, cubit_mod, surface_ids, edge_curve_ids=None):
    """Set material and boundary labels from Cubit blocks/sidesets/entity names.

    Cubit allows blocks/sidesets/nodesets to contain any combination of
    entities, so we detect material vs boundary usage by content:

    Material (volume domain) — priority:
      1. Block that contains the volume (parse_cubit_list "volume in block")
      2. Cubit volume entity name (if user renamed it)
      3. Fallback: "volume_N"

    Boundary (surface bcname) — priority:
      1. Sideset containing the surface (parse_cubit_list "surface in sideset")
      2. Block whose tris/faces overlap the surface's tris/faces
      3. Cubit surface entity name (if user renamed it)
      4. Fallback: "surface_N"

    A single block may serve BOTH roles (e.g., block contains volume AND
    surface tris). Material and boundary detection are independent.
    """
    block_ids = list(cubit_mod.get_block_id_list())
    vol_ids = list(cubit_mod.get_entities("volume"))

    # --- Material labels (per-domain, multi-volume) ---
    # Detect by volume membership, NOT element type
    vol_block_names = {}  # vid -> name
    for bid in block_ids:
        bvols = list(cubit_mod.parse_cubit_list("volume", f"in block {bid}"))
        if bvols:
            bname = cubit_mod.get_exodus_entity_name("block", bid)
            for vid in bvols:
                if vid not in vol_block_names:
                    vol_block_names[vid] = bname

    for domain_idx, vid in enumerate(vol_ids, start=1):
        if vid in vol_block_names:
            name = vol_block_names[vid]
        else:
            ename = cubit_mod.get_entity_name("volume", vid)
            if ename and not ename.startswith("Volume"):
                name = ename
            else:
                name = f"volume_{vid}"
        ng_mesh.SetMaterial(domain_idx, name)

    # --- Boundary labels ---
    surf_names = {}  # sid -> name

    # Priority 1: Sideset names (surface membership)
    try:
        for ssid in cubit_mod.get_sideset_id_list():
            ssname = cubit_mod.get_exodus_entity_name("sideset", ssid)
            for sid in cubit_mod.parse_cubit_list(
                    "surface", f"in sideset {ssid}"):
                if sid in surface_ids and sid not in surf_names:
                    surf_names[sid] = ssname
    except Exception:
        pass

    # Priority 2: Block names (tri/face overlap with surface)
    # Independent of material detection — same block can provide both
    for bid in block_ids:
        btris = set(cubit_mod.get_block_tris(bid))
        try:
            bfaces = set(cubit_mod.get_block_faces(bid))
        except Exception:
            bfaces = set()
        if not btris and not bfaces:
            continue
        bname = cubit_mod.get_exodus_entity_name("block", bid)
        for sid in surface_ids:
            if sid in surf_names:
                continue  # sideset takes priority
            stris = set(cubit_mod.parse_cubit_list(
                "tri", f"in surface {sid}"))
            sfaces = set(cubit_mod.parse_cubit_list(
                "face", f"in surface {sid}"))
            if (stris and stris & btris) or (sfaces and sfaces & bfaces):
                surf_names[sid] = bname

    # Priority 3/4: Entity name or fallback
    for i, sid in enumerate(surface_ids):
        if sid in surf_names:
            name = surf_names[sid]
        else:
            ename = cubit_mod.get_entity_name("surface", sid)
            if ename and not ename.startswith("Surface"):
                name = ename
            else:
                name = f"surface_{sid}"
        ng_mesh.SetBCName(i, name)

    # --- Edge labels (BBND / codim-2) ---
    # edgenr in C++ is 1-based, one per edge_curve_ids entry.
    # SetCD2Name uses 1-based index (matching edgenr and read_gmsh convention).
    if edge_curve_ids:
        # Priority 1: Sideset on curves
        curve_sideset_names = {}  # cid -> name
        try:
            for ssid in cubit_mod.get_sideset_id_list():
                ssname = cubit_mod.get_exodus_entity_name("sideset", ssid)
                for cid in cubit_mod.parse_cubit_list(
                        "curve", f"in sideset {ssid}"):
                    if cid not in curve_sideset_names:
                        curve_sideset_names[cid] = ssname
        except Exception:
            pass

        for i, cid in enumerate(edge_curve_ids):
            if cid in curve_sideset_names:
                name = curve_sideset_names[cid]
            else:
                ename = cubit_mod.get_entity_name("curve", cid)
                if ename and not ename.startswith("Curve"):
                    name = ename
                else:
                    name = f"curve_{cid}"
            # edgenr is 1-based (matching segment edgenr in C++)
            cd2nr = i + 1
            try:
                ng_mesh.SetCD2Name(cd2nr, name)
            except Exception as e:
                print(f"Warning: SetCD2Name({cd2nr}, '{name}') failed: {e}")
