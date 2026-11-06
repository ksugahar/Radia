"""Export the active Cubit mesh as a JMAG static solver input."""
import math


def _resolve_block_materials(cubit, spec):
    """Resolve named assignments into the numeric JCF builder specification."""
    if "block_materials" not in spec and "material_library" not in spec:
        return spec
    if "materials" in spec or "magnets" in spec:
        raise ValueError("use named block assignments or numeric materials and magnets")
    library = spec.get("material_library")
    assignments = spec.get("block_materials")
    if not isinstance(library, dict) or not isinstance(assignments, dict) or not assignments:
        raise ValueError("material_library and nonempty block_materials must be JSON objects")

    by_name = {}
    for block in cubit.get_block_id_list():
        name = cubit.get_exodus_entity_name("block", block)
        by_name.setdefault(name, []).append(int(block))

    materials, magnets = [], []
    for name, key in assignments.items():
        matches = by_name.get(name, [])
        if len(matches) != 1:
            raise ValueError(f"expected one Cubit block named {name!r}, got {len(matches)}")
        if not isinstance(key, str) or not isinstance(library.get(key), dict):
            raise ValueError(f"block {name!r} has unknown material {key!r}")
        part = matches[0]
        material = {**library[key], "part": part, "name": name}
        if "magnet" in material:
            magnet = material["magnet"]
            if not isinstance(magnet, dict) or "direction" not in magnet:
                raise ValueError(f"magnet material {key!r} requires a direction")
            magnet = dict(magnet)
            direction = magnet.pop("direction")
            if (not isinstance(direction, list) or len(direction) != 3
                or any(type(v) not in (float, int) or not math.isfinite(v) for v in direction)
                or not any(direction)):
                raise ValueError(f"magnet material {key!r} requires a finite nonzero 3-vector")
            material["magnet"] = magnet
            magnets.append({"part": part, "direction": direction})
        materials.append(material)
    return {**spec, "materials": materials, "magnets": magnets}


def export_jcf(path, spec, *, symmetry_sidesets=()):
    """Write a JCF using named exterior sidesets, independent of surface sense."""
    import cubit
    from mcp_server_jmag.jcf_build import build_jcf

    spec = _resolve_block_materials(cubit, spec)
    boundaries = list(spec.get("boundaries") or [])
    for name in symmetry_sidesets:
        matches = [sid for sid in cubit.get_sideset_id_list()
                   if cubit.get_exodus_entity_name("sideset", sid) == name]
        if len(matches) != 1:
            raise ValueError(f"expected one Cubit sideset named {name!r}, got {len(matches)}")
        faces = set()
        for surface in cubit.get_sideset_surfaces(matches[0]):
            for kind, ids in (("tri", cubit.get_surface_tris(surface)),
                              ("quad", cubit.get_surface_quads(surface))):
                for face in ids:
                    faces.add(tuple(sorted(cubit.get_expanded_connectivity(kind, face))))
        for kind, ids in (("tri", cubit.get_sideset_tris(matches[0])),
                          ("quad", cubit.get_sideset_quads(matches[0]))):
            for face in ids:
                faces.add(tuple(sorted(cubit.get_expanded_connectivity(kind, face))))
        if cubit.get_sideset_curves(matches[0]) or cubit.get_sideset_edges(matches[0]):
            raise ValueError(f"sideset {name!r} contains non-face members")
        if not faces:
            raise ValueError(f"sideset {name!r} has no mesh faces")
        boundaries.append({"type": "symmetry", "face_nodes": sorted(faces), "label": name})

    mesh = _read_active_mesh(cubit, spec.get("unit", "mm"))
    return build_jcf(mesh, path, {**spec, "boundaries": boundaries})


def _read_active_mesh(cubit, unit):
    from mcp_server_jmag.jcf_build import HEX, PRISM, PYRAMID, TET, mesh_from_cubit_data

    kinds = (("hex", HEX, 8, cubit.get_block_hexes), ("tet", TET, 4, cubit.get_block_tets),
             ("wedge", PRISM, 6, cubit.get_block_wedges), ("pyramid", PYRAMID, 5, cubit.get_block_pyramids))
    elements = []
    for kind, code, arity, block_members in kinds:
        assigned = set()
        for block in cubit.get_block_id_list():
            for local_id in block_members(block):
                if local_id in assigned:
                    raise ValueError(f"Cubit {kind} {local_id} belongs to multiple blocks")
                connectivity = tuple(cubit.get_expanded_connectivity(kind, local_id))
                if len(connectivity) != arity:
                    raise ValueError(f"Cubit {kind} {local_id} must have {arity} nodes")
                global_id = cubit.get_global_element_id(kind, local_id)
                elements.append((global_id, code, block, connectivity))
                assigned.add(local_id)
        unassigned = set(cubit.get_entities(kind)) - assigned
        if unassigned:
            raise ValueError(f"Cubit {kind} {min(unassigned)} has no block assignment")

    node_ids = {node for _, _, _, conn in elements for node in conn}
    nodes = {node: cubit.get_nodal_coordinates(node) for node in node_ids}
    return mesh_from_cubit_data(nodes, sorted(elements), unit)
