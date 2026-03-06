#!/usr/bin/env python
"""
nastran_mesh_import.py - Import Nastran mesh files to Radia

Supports fixed-width and long-format Nastran (.bdf, .nas, .dat)
for hexahedral (CHEXA) and tetrahedral (CTETRA) elements.

Handles continuation lines marked with '+' for multi-line cards.

Author: Radia Development Team
Created: 2025-11-22
Version: 0.2.0
"""


def parse_nastran_fixed_width(line, card_type):
    """
    Parse Nastran fixed-width format card.

    Fixed-width format:
    - Field width: 8 characters
    - Fields: 0-7, 8-15, 16-23, 24-31, ...
    """
    fields = []
    for i in range(0, len(line), 8):
        field = line[i:i+8].strip()
        if field and field != '+':  # Ignore continuation markers
            fields.append(field)
    return fields


def import_nastran_mesh(nas_file, units='m', verbose=True):
    """
    Import Nastran mesh file to Radia.

    Parameters
    ----------
    nas_file : str
        Path to Nastran file (.bdf, .nas, .dat)
    units : str, default='m'
        Unit system: 'm' (meters) or 'mm' (millimeters)
        Must match rad.FldUnits() setting
    verbose : bool, default=True
        Print progress information

    Returns
    -------
    dict
        {
            'vertices': list of [x, y, z],
            'hex_elements': list of element vertex indices,
            'tet_elements': list of element vertex indices,
            'wedge_elements': list of element vertex indices,
            'pyramid_elements': list of element vertex indices
        }

    Notes
    -----
    Supports fixed-width and long-format Nastran:
    - GRID/GRID* cards: Node definitions (with continuation)
    - CHEXA cards: 8-node hexahedral elements (with continuation)
    - CTETRA cards: 4-node tetrahedral elements (with continuation)
    - CPENTA cards: 6-node wedge/pentahedron elements (with continuation)
    - CPYRAM cards: 5-node pyramid elements (with continuation)
    """
    if verbose:
        print(f"[Nastran Import] Reading file: {nas_file}")

    vertices = []
    hex_elements = []
    tet_elements = []
    wedge_elements = []
    pyramid_elements = []
    tria_groups = {}  # material_id -> {'faces': [[n1,n2,n3], ...], 'node_ids': set()}
    vertex_map = {}  # node_id -> vertex_index

    with open(nas_file, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i]

        # Skip comments and empty lines
        if line.startswith('$') or not line.strip():
            i += 1
            continue

        # End of bulk data
        if line.strip().startswith('ENDDATA'):
            break

        # Parse card type (first 8 characters)
        card_type = line[:8].strip()

        # GRID/GRID* card: Node definition with continuation support
        if card_type == 'GRID' or card_type.startswith('GRID'):
            # Check for long format (GRID*)
            is_long = '*' in card_type

            if is_long:
                # GRID* long format: 16-character fields
                # First line: GRID* (0-7), node_id (8-23), coord_sys (24-39), X (40-55), Y (56-71)
                # Second line: * (0-7), Z (8-23)
                try:
                    node_id = int(line[8:24].strip())
                    coord_sys = int(line[24:40].strip()) if len(line) > 24 else 0
                    x = float(line[40:56].strip())
                    y = float(line[56:72].strip()) if len(line) > 56 else 0.0

                    # Get continuation line for Z coordinate
                    z = 0.0
                    if i + 1 < len(lines) and lines[i + 1].strip().startswith('*'):
                        i += 1
                        cont_line = lines[i]
                        z = float(cont_line[8:24].strip())

                    vertex_idx = len(vertices)
                    vertices.append([x, y, z])
                    vertex_map[node_id] = vertex_idx
                except (ValueError, IndexError) as e:
                    if verbose:
                        print(f"  Warning: Failed to parse GRID* at line {i+1}: {e}")
            else:
                # Standard GRID format: 8-character fields
                # GRID (0-7), node_id (8-15), coord_sys (16-23), X (24-31), Y (32-39), Z (40-47)
                try:
                    node_id = int(line[8:16].strip())
                    coord_sys = int(line[16:24].strip()) if len(line) > 16 and line[16:24].strip() else 0
                    x = float(line[24:32].strip())
                    y = float(line[32:40].strip()) if len(line) > 32 else 0.0
                    z = float(line[40:48].strip()) if len(line) > 40 else 0.0

                    vertex_idx = len(vertices)
                    vertices.append([x, y, z])
                    vertex_map[node_id] = vertex_idx
                except (ValueError, IndexError) as e:
                    if verbose:
                        print(f"  Warning: Failed to parse GRID at line {i+1}: {e}")

        # CHEXA card: 8-node hexahedral element with continuation support
        elif card_type == 'CHEXA' or card_type.startswith('CHEXA'):
            # Collect continuation lines (marked with + at end)
            full_line = line
            while i + 1 < len(lines) and full_line.rstrip().endswith('+'):
                i += 1
                full_line = full_line.rstrip()[:-1]  # Remove trailing '+'
                cont_line = lines[i]
                if cont_line.startswith('+'):
                    cont_line = cont_line[1:]  # Remove leading '+'
                full_line += cont_line

            fields = parse_nastran_fixed_width(full_line, card_type)
            if len(fields) >= 10:
                try:
                    # Element ID (field 1), Property ID (field 2), then 8 nodes
                    node_ids = [int(fields[j]) for j in range(3, 11)]
                    element_verts = [vertex_map[nid] for nid in node_ids]
                    hex_elements.append(element_verts)
                except (ValueError, KeyError, IndexError) as e:
                    if verbose:
                        print(f"  Warning: Failed to parse CHEXA at line {i+1}: {e}")

        # CTETRA card: 4-node tetrahedral element with continuation support
        elif card_type == 'CTETRA' or card_type.startswith('CTETRA'):
            # Collect continuation lines
            full_line = line
            while i + 1 < len(lines) and full_line.rstrip().endswith('+'):
                i += 1
                full_line = full_line.rstrip()[:-1]  # Remove trailing '+'
                cont_line = lines[i]
                if cont_line.startswith('+'):
                    cont_line = cont_line[1:]  # Remove leading '+'
                full_line += cont_line

            fields = parse_nastran_fixed_width(full_line, card_type)
            if len(fields) >= 6:
                try:
                    node_ids = [int(fields[j]) for j in range(3, 7)]
                    element_verts = [vertex_map[nid] for nid in node_ids]
                    tet_elements.append(element_verts)
                except (ValueError, KeyError, IndexError) as e:
                    if verbose:
                        print(f"  Warning: Failed to parse CTETRA at line {i+1}: {e}")

        # CPENTA card: 6-node wedge/pentahedron element with continuation support
        elif card_type == 'CPENTA' or card_type.startswith('CPENTA'):
            # Collect continuation lines
            full_line = line
            while i + 1 < len(lines) and full_line.rstrip().endswith('+'):
                i += 1
                full_line = full_line.rstrip()[:-1]  # Remove trailing '+'
                cont_line = lines[i]
                if cont_line.startswith('+'):
                    cont_line = cont_line[1:]  # Remove leading '+'
                full_line += cont_line

            fields = parse_nastran_fixed_width(full_line, card_type)
            if len(fields) >= 8:
                try:
                    # Element ID (field 1), Property ID (field 2), then 6 nodes
                    node_ids = [int(fields[j]) for j in range(3, 9)]
                    element_verts = [vertex_map[nid] for nid in node_ids]
                    wedge_elements.append(element_verts)
                except (ValueError, KeyError, IndexError) as e:
                    if verbose:
                        print(f"  Warning: Failed to parse CPENTA at line {i+1}: {e}")

        # CPYRAM card: 5-node pyramid element with continuation support
        elif card_type == 'CPYRAM' or card_type.startswith('CPYRAM'):
            # Collect continuation lines
            full_line = line
            while i + 1 < len(lines) and full_line.rstrip().endswith('+'):
                i += 1
                full_line = full_line.rstrip()[:-1]  # Remove trailing '+'
                cont_line = lines[i]
                if cont_line.startswith('+'):
                    cont_line = cont_line[1:]  # Remove leading '+'
                full_line += cont_line

            fields = parse_nastran_fixed_width(full_line, card_type)
            if len(fields) >= 7:
                try:
                    # Element ID (field 1), Property ID (field 2), then 5 nodes
                    node_ids = [int(fields[j]) for j in range(3, 8)]
                    element_verts = [vertex_map[nid] for nid in node_ids]
                    pyramid_elements.append(element_verts)
                except (ValueError, KeyError, IndexError) as e:
                    if verbose:
                        print(f"  Warning: Failed to parse CPYRAM at line {i+1}: {e}")

        # CTRIA3 card: 3-node triangle surface element
        elif card_type == 'CTRIA3' or card_type.startswith('CTRIA3'):
            # CTRIA3 format: CTRIA3, elem_id, prop_id, n1, n2, n3
            # Fixed format: 8-character fields
            try:
                elem_id = int(line[8:16].strip())
                prop_id = int(line[16:24].strip())  # Material/Property ID
                n1 = int(line[24:32].strip())
                n2 = int(line[32:40].strip())
                n3 = int(line[40:48].strip())

                # Group by material ID (prop_id)
                if prop_id not in tria_groups:
                    tria_groups[prop_id] = {'faces': [], 'node_ids': set()}

                # Store face (using original node IDs, will convert later)
                tria_groups[prop_id]['faces'].append([n1, n2, n3])
                tria_groups[prop_id]['node_ids'].update([n1, n2, n3])

            except (ValueError, IndexError) as e:
                if verbose:
                    print(f"  Warning: Failed to parse CTRIA3 at line {i+1}: {e}")

        i += 1

    # Calculate total triangle faces across all material groups
    total_tria_faces = sum(len(group['faces']) for group in tria_groups.values())

    if verbose:
        print(f"[Nastran Import] Parsing complete")
        print(f"                 Vertices: {len(vertices)}")
        print(f"                 Hexahedral elements: {len(hex_elements)}")
        print(f"                 Wedge elements: {len(wedge_elements)}")
        print(f"                 Pyramid elements: {len(pyramid_elements)}")
        print(f"                 Tetrahedral elements: {len(tet_elements)}")
        print(f"                 Triangle surface groups: {len(tria_groups)} (total {total_tria_faces} faces)")

    return {
        'vertices': vertices,
        'hex_elements': hex_elements,
        'wedge_elements': wedge_elements,
        'pyramid_elements': pyramid_elements,
        'tet_elements': tet_elements,
        'tria_groups': tria_groups,
        'vertex_map': vertex_map
    }


def create_radia_from_nastran(nas_file, material=None, units='m',
                                combine=True, verbose=True):
    """
    Create Radia geometry from Nastran mesh file.

    Parameters
    ----------
    nas_file : str
        Path to Nastran file
    material : dict or None
        Material specification: {'magnetization': [Mx, My, Mz]}
    units : str, default='m'
        'm' (meters) or 'mm' (millimeters)
    combine : bool, default=True
        If True, return rad.ObjCnt() container
        If False, return list of object IDs
    verbose : bool, default=True
        Print progress

    Returns
    -------
    int or list
        Radia object ID (if combine=True) or list of IDs (if combine=False)
    """
    import radia as rad

    # Import mesh data
    mesh_data = import_nastran_mesh(nas_file, units=units, verbose=verbose)

    vertices = mesh_data['vertices']
    hex_elements = mesh_data['hex_elements']
    wedge_elements = mesh_data['wedge_elements']
    pyramid_elements = mesh_data['pyramid_elements']
    tet_elements = mesh_data['tet_elements']
    tria_groups = mesh_data['tria_groups']
    vertex_map = mesh_data['vertex_map']

    total_elements = len(hex_elements) + len(wedge_elements) + len(pyramid_elements) + len(tet_elements)
    total_tria_groups = len(tria_groups)
    if total_elements == 0 and total_tria_groups == 0:
        raise ValueError("No valid elements found in Nastran file")

    # Material specification
    if material is None:
        magnetization = [0, 0, 0]
    elif isinstance(material, dict) and 'magnetization' in material:
        magnetization = material['magnetization']
    else:
        raise ValueError("material must be dict with 'magnetization' key")

    # Import face topologies
    from netgen_mesh_import import HEX_FACES, WEDGE_FACES, PYRAMID_FACES, TETRA_FACES

    # Create Radia objects
    radia_objects = []

    # Create hexahedral elements
    if len(hex_elements) > 0 and verbose:
        print(f"[Nastran Import] Creating {len(hex_elements)} hexahedral elements...")

    for i, elem_verts in enumerate(hex_elements):
        elem_coords = [vertices[vi] for vi in elem_verts]
        obj_id = rad.ObjPolyhdr(elem_coords, HEX_FACES, magnetization)
        radia_objects.append(obj_id)

        if verbose and (i + 1) % 25 == 0:
            print(f"                 Progress: {i+1}/{len(hex_elements)}", end='\r')

    if len(hex_elements) > 0 and verbose:
        print(f"                 Progress: {len(hex_elements)}/{len(hex_elements)}")

    # Create wedge elements
    if len(wedge_elements) > 0 and verbose:
        print(f"[Nastran Import] Creating {len(wedge_elements)} wedge elements...")

    for i, elem_verts in enumerate(wedge_elements):
        elem_coords = [vertices[vi] for vi in elem_verts]
        obj_id = rad.ObjPolyhdr(elem_coords, WEDGE_FACES, magnetization)
        radia_objects.append(obj_id)

        if verbose and (i + 1) % 50 == 0:
            print(f"                 Progress: {i+1}/{len(wedge_elements)}", end='\r')

    if len(wedge_elements) > 0 and verbose:
        print(f"                 Progress: {len(wedge_elements)}/{len(wedge_elements)}")

    # Create pyramid elements
    if len(pyramid_elements) > 0 and verbose:
        print(f"[Nastran Import] Creating {len(pyramid_elements)} pyramid elements...")

    for i, elem_verts in enumerate(pyramid_elements):
        elem_coords = [vertices[vi] for vi in elem_verts]
        obj_id = rad.ObjPolyhdr(elem_coords, PYRAMID_FACES, magnetization)
        radia_objects.append(obj_id)

        if verbose and (i + 1) % 50 == 0:
            print(f"                 Progress: {i+1}/{len(pyramid_elements)}", end='\r')

    if len(pyramid_elements) > 0 and verbose:
        print(f"                 Progress: {len(pyramid_elements)}/{len(pyramid_elements)}")

    # Create tetrahedral elements
    if len(tet_elements) > 0 and verbose:
        print(f"[Nastran Import] Creating {len(tet_elements)} tetrahedral elements...")

    for i, elem_verts in enumerate(tet_elements):
        elem_coords = [vertices[vi] for vi in elem_verts]
        obj_id = rad.ObjPolyhdr(elem_coords, TETRA_FACES, magnetization)
        radia_objects.append(obj_id)

        if verbose and (i + 1) % 100 == 0:
            print(f"                 Progress: {i+1}/{len(tet_elements)}", end='\r')

    if len(tet_elements) > 0 and verbose:
        print(f"                 Progress: {len(tet_elements)}/{len(tet_elements)}")

    # Create surface polyhedra from triangle groups
    if len(tria_groups) > 0 and verbose:
        print(f"[Nastran Import] Creating {len(tria_groups)} surface polyhedra from triangle groups...")

    for mat_id, group in tria_groups.items():
        faces_global = group['faces']  # Faces with global node IDs
        node_ids = sorted(group['node_ids'])  # All unique node IDs for this material

        # Create local node ID mapping (global node ID -> local index 1-based for Radia)
        node_id_to_local = {nid: idx + 1 for idx, nid in enumerate(node_ids)}

        # Get vertex coordinates for this group
        group_vertices = [vertices[vertex_map[nid]] for nid in node_ids]

        # Convert faces from global node IDs to local 1-based indices
        faces_local = []
        for face in faces_global:
            local_face = [node_id_to_local[nid] for nid in face]
            faces_local.append(local_face)

        # Create polyhedron from surface mesh
        obj_id = rad.ObjPolyhdr(group_vertices, faces_local, magnetization)
        radia_objects.append(obj_id)

        if verbose:
            print(f"                 Material {mat_id}: {len(faces_local)} faces, {len(group_vertices)} vertices")

    # Return result
    if combine:
        if verbose:
            print(f"[Nastran Import] Combining into container...")

        container = rad.ObjCnt(radia_objects)

        if verbose:
            print(f"                 [OK] Container object ID: {container}")

        return container
    else:
        return radia_objects


__all__ = [
    'import_nastran_mesh',
    'create_radia_from_nastran',
]
