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

    # --- Phase E: Call C++ for high-order curving ---
    # For surface_only: pass empty volume elements to avoid
    # ClearVolumeElements() segfault in NGSolve Mesh() constructor.
    ve = [] if surface_only else vol_elements
    ng_mesh = radia_cubit_mesh.build_curved_mesh(
        node_coords, node_ids, ve, surface_data,
        project_func, normal_func,
        order=order, surface_only=False, split_quads=split_quads,
    )
    return ng_mesh
