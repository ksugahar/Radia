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
    from radia.cubit_netgen_bridge import extract_curved_mesh

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
        import radia_cubit_mesh
    except ImportError:
        # PyPI install: .pyd is inside radia package directory
        import os, sys
        radia_dir = os.path.dirname(os.path.abspath(__file__))
        if radia_dir not in sys.path:
            sys.path.insert(0, radia_dir)
        import radia_cubit_mesh

    # --- Phase A: Extract nodes ---
    node_id_list = list(cubit_mod.parse_cubit_list("node", "all"))
    N = len(node_id_list)
    node_coords = np.empty((N, 3), dtype=np.float64)
    node_ids = np.array(node_id_list, dtype=np.int32)
    for i, nid in enumerate(node_id_list):
        node_coords[i] = cubit_mod.get_nodal_coordinates(nid)

    # --- Phase B: Extract volume elements ---
    vol_elements = []
    TYPE_MAP = {
        "TETRA4": 0, "TETRA": 0, "TET4": 0, "TET": 0,
        "HEX8": 1, "HEX": 1,
        "WEDGE6": 2, "WEDGE": 2,
        "PYRAMID5": 3, "PYRAMID": 3,
    }
    tet_ids = cubit_mod.parse_cubit_list("tet", "all")
    for tid in tet_ids:
        conn = list(cubit_mod.get_connectivity("tet", tid))
        vol_elements.append((0, conn))  # TET

    hex_ids = cubit_mod.parse_cubit_list("hex", "all")
    for hid in hex_ids:
        conn = list(cubit_mod.get_connectivity("hex", hid))
        vol_elements.append((1, conn))  # HEX

    # --- Phase C: Extract surface mesh + UVs ---
    surface_ids = list(cubit_mod.parse_cubit_list("surface", "all"))
    surface_data = []
    surface_objects = {}  # cache cubit.surface() objects

    # Determine surface element type: Cubit uses "tri" for tet meshes,
    # "face" for hex meshes (returns 4-node quads)
    has_tets = len(tet_ids) > 0
    has_hexes = len(hex_ids) > 0

    for sid in surface_ids:
        surf = cubit_mod.surface(sid)
        surface_objects[sid] = surf

        tris = []
        quads = []

        if has_tets:
            # Tet mesh: surface tris via "tri" entity type
            tri_ids_on_surf = cubit_mod.parse_cubit_list("tri", f"in surface {sid}")
            for tid in tri_ids_on_surf:
                tris.extend(cubit_mod.get_connectivity("tri", tid))

        if has_hexes:
            # Hex mesh: surface quads via "face" entity type
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
    )
    return ng_mesh
