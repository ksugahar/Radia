"""Export the active Cubit mesh as a static solver input (.jcf).

The analysis settings are a JSON object that names Cubit blocks and sidesets:

    {"unit": "mm",
     "material_library": {"pm":   {"magnet": {"br": 1.4, "hcb": 1061000, "direction": [0, 0, 1]}},
                          "core": {"bh": [[0, 0], [100, 1.0], ...], "lamination_factor": 98},
                          "cu":   {"mu_r": 1, "resistivity": 1.673e-8},
                          "air":  {"mu_r": 1}},
     "block_materials": {"magnet": "pm", "yoke": "core", "coil": "cu", "air": "air"},
     "external_field": {"axis": [0, 0, 1], "value": 0.5},
     "coils": [{"blocks": ["coil"], "current": 10, "turns": 100,
                "flowin": {"sideset": "coil_cut", "normal": [0, 0, 1]}},
               {"blocks": ["coil2"], "current": 5, "direction": [0, 1, 0], "origin": [0, 0, 0]}],
     "boundaries": {"symmetry_xy": "symmetry"}}

``material_library`` may also be a file name relative to the settings file (resolved by the
caller). Every block needs a material. Lengths (``origin``) are in ``unit``. A flow-in sideset
lists the faces the current crosses and ``normal`` is the current direction across them. Each
face goes to the coil element whose outward normal points along ``normal``, the element the
current leaves; this is the convention checked against the solver. Only static analysis is
written.

``boundaries`` replaced the former ``symmetry_sidesets`` list, which is now rejected:
``"symmetry_sidesets": ["sym"]`` becomes ``"boundaries": {"sym": "symmetry"}``.

Nothing else is imported: the writer is :mod:`cubit_mesh_export.jcf_writer`.
"""
from __future__ import annotations

from . import jcf_writer as jw

SETTINGS_KEYS = {"unit", "material_library", "block_materials", "external_field", "coils",
                 "boundaries", "control"}
_FREE_MESH = ("Cubit deletes mesh that no geometry owns (free mesh, e.g. from MeshImport) when "
              "control returns to the GUI; mesh the geometry instead")


def _vector(v, what):
    if (not isinstance(v, (list, tuple)) or len(v) != 3
            or any(type(x) not in (int, float) for x in v)):
        raise ValueError(f"{what} must be a list of three numbers")
    return [float(x) for x in v]


def _named(cubit, kind, name):
    ids = cubit.get_block_id_list() if kind == "block" else cubit.get_sideset_id_list()
    hits = [int(i) for i in ids if cubit.get_exodus_entity_name(kind, i) == name]
    if len(hits) != 1:
        raise ValueError(f"expected one Cubit {kind} named {name!r}, got {len(hits)}")
    return hits[0]


def _sideset_faces(cubit, name):
    """Sorted node-id tuples of every mesh face in the named sideset (surface members and
    explicit tris/quads)."""
    sid = _named(cubit, "sideset", name)
    if cubit.get_sideset_curves(sid) or cubit.get_sideset_edges(sid):
        raise ValueError(f"sideset {name!r} contains non-face members")
    faces = set()
    surfaces = list(cubit.get_sideset_surfaces(sid))
    for surface in surfaces:
        for kind, ids in (("tri", cubit.get_surface_tris(surface)), ("quad", cubit.get_surface_quads(surface))):
            faces.update(tuple(sorted(cubit.get_expanded_connectivity(kind, f))) for f in ids)
    for kind, ids in (("tri", cubit.get_sideset_tris(sid)), ("quad", cubit.get_sideset_quads(sid))):
        faces.update(tuple(sorted(cubit.get_expanded_connectivity(kind, f))) for f in ids)
    if not faces:
        why = (f"its surfaces {surfaces} are not meshed" if surfaces else
               f"it holds no surfaces and no mesh faces; if its faces were free mesh: {_FREE_MESH}")
        raise ValueError(f"sideset {name!r} has no mesh faces: {why}")
    return sorted(faces)


def settings_to_spec(cubit, settings: dict) -> dict:
    """Resolve block and sideset names into the writer's numeric specification."""
    if not isinstance(settings, dict):
        raise ValueError("analysis settings must be a JSON object")
    unknown = set(settings) - SETTINGS_KEYS
    if unknown:
        raise ValueError(f"unknown analysis settings {sorted(unknown)}; allowed {sorted(SETTINGS_KEYS)}")
    library = settings.get("material_library")
    assignments = settings.get("block_materials")
    if not isinstance(library, dict) or not isinstance(assignments, dict) or not assignments:
        raise ValueError("material_library and a nonempty block_materials must be JSON objects")
    parts, materials, magnets = {}, [], []
    for name, key in assignments.items():
        if not isinstance(key, str) or not isinstance(library.get(key), dict):
            raise ValueError(f"block {name!r} has unknown material {key!r}")
        part = _named(cubit, "block", name)
        material = {**library[key], "part": part, "name": name}
        if "magnet" in material:
            magnet = dict(material["magnet"]) if isinstance(material["magnet"], dict) else None
            if magnet is None or "direction" not in magnet:
                raise ValueError(f"magnet material {key!r} requires a direction")
            magnets.append({"part": part, "direction": _vector(magnet.pop("direction"), f"direction of {key!r}")})
            material["magnet"] = magnet
        parts[name] = part
        materials.append(material)
    blocks = {int(b) for b in cubit.get_block_id_list()}
    unnamed = sorted(blocks - set(parts.values()))
    if unnamed:
        raise ValueError(f"blocks {unnamed} have no material in block_materials")
    spec = {"unit": settings.get("unit", "mm"), "materials": materials, "magnets": magnets}
    if settings.get("external_field") is not None:
        spec["external_field"] = settings["external_field"]
    if settings.get("control"):
        spec["control"] = settings["control"]
    coils = []
    for c in settings.get("coils") or []:
        names = c.get("blocks") or []
        missing = [n for n in names if n not in parts]
        if not names or missing:
            raise ValueError(f"coil blocks {missing or names} must be blocks with materials")
        coil = {"parts": [parts[n] for n in names], "current": c.get("current"), "turns": c.get("turns", 1)}
        if coil["current"] is None:
            raise ValueError("a coil needs 'current' [A]")
        fi = c.get("flowin")
        if fi is not None:
            if set(fi) != {"sideset", "normal"}:
                raise ValueError("coil flowin needs exactly 'sideset' and 'normal'")
            coil["flowin"] = {"face_nodes": _sideset_faces(cubit, fi["sideset"]),
                              "normal": _vector(fi["normal"], "flowin normal"), "label": fi["sideset"]}
        else:
            coil["direction"] = _vector(c.get("direction"), "coil direction")
            coil["origin"] = _vector(c.get("origin"), "coil origin")
        coils.append(coil)
    if coils:
        spec["coils"] = coils
    boundaries = settings.get("boundaries") or {}
    if not isinstance(boundaries, dict):
        raise ValueError("boundaries must map sideset names to a boundary type")
    spec["boundaries"] = [{"type": kind, "face_nodes": _sideset_faces(cubit, name), "label": name}
                          for name, kind in boundaries.items()]
    return spec


def export_jcf(path, settings: dict) -> dict:
    """Write the active Cubit mesh and the named analysis settings to ``path`` (.jcf)."""
    import cubit

    spec = settings_to_spec(cubit, settings)
    mesh = read_active_mesh(cubit, spec["unit"])
    return jw.write_jcf(mesh, path, spec)


def read_active_mesh(cubit, unit):
    kinds = (
        ("hex", jw.HEX, 8, cubit.get_block_hexes, cubit.get_volume_hexes),
        ("tet", jw.TET, 4, cubit.get_block_tets, cubit.get_volume_tets),
        ("wedge", jw.PRISM, 6, cubit.get_block_wedges, cubit.get_volume_wedges),
        ("pyramid", jw.PYRAMID, 5, cubit.get_block_pyramids, cubit.get_volume_pyramids),
    )
    block_volumes = {block: cubit.get_block_volumes(block) for block in cubit.get_block_id_list()}
    elements = []
    for kind, code, arity, block_members, volume_members in kinds:
        assigned = set()
        for block, volumes in block_volumes.items():
            members = set(block_members(block))
            for volume in volumes:
                members.update(volume_members(volume))
            for local_id in members:
                if local_id in assigned:
                    raise ValueError(f"Cubit {kind} {local_id} belongs to multiple blocks")
                connectivity = tuple(cubit.get_expanded_connectivity(kind, local_id))
                if len(connectivity) != arity:
                    raise ValueError(f"Cubit {kind} {local_id} must have {arity} nodes")
                elements.append((cubit.get_global_element_id(kind, local_id), code, block, connectivity))
                assigned.add(local_id)
        unassigned = set(cubit.get_entities(kind)) - assigned
        if unassigned:
            raise ValueError(f"Cubit {kind} {min(unassigned)} has no block assignment")
    if not elements:
        raise ValueError(f"the Cubit model has no volume elements; {_FREE_MESH}")
    node_ids = {node for _, _, _, conn in elements for node in conn}
    nodes = {node: cubit.get_nodal_coordinates(node) for node in node_ids}
    return jw.mesh_from_data(nodes, sorted(elements), unit)
