"""
gmsh_mesh_import.py

Direct GMSH .msh file reader for Radia MMM/MSC solver.

Reads GMSH 2.x ASCII format and creates Radia hexahedral/tetrahedral objects.
Does NOT require NGSolve as a dependency.

Supported element types:
  - Tet4  (GMSH type 4): 4-node tetrahedron -> ObjTetrahedron (3 DOF, MMM)
  - Hex8  (GMSH type 5): 8-node hexahedron  -> ObjHexahedron  (6 DOF, MSC)
  - Wedge6 (GMSH type 6): 6-node prism      -> ObjWedge       (5 DOF, MSC)

Surface elements (Tri3, Quad4) are skipped for volume mesh import.

Usage:
    import radia as rad
    from gmsh_mesh_import import gmsh_to_radia

    rad.FldUnits('m')
    core = gmsh_to_radia('ferrite_core.msh', mu_r=1000)
    rad.Solve(core, 0.0001, 1000, 0)

Part of Radia project
"""

import os
import numpy as np


# GMSH element type codes
GMSH_TRI3 = 2
GMSH_QUAD4 = 3
GMSH_TET4 = 4
GMSH_HEX8 = 5
GMSH_WEDGE6 = 6

# Number of nodes per element type
GMSH_NODES_PER_TYPE = {
    1: 2,   # Line2
    2: 3,   # Tri3
    3: 4,   # Quad4
    4: 4,   # Tet4
    5: 8,   # Hex8
    6: 6,   # Wedge6 (Prism)
    7: 5,   # Pyramid5
    8: 3,   # Line3
    9: 6,   # Tri6
    10: 9,  # Quad9
    11: 10, # Tet10
    15: 1,  # Point
}


def read_gmsh(filename, physical_group=None):
    """
    Read a GMSH .msh file and extract volume elements.

    Supports GMSH 2.x ASCII format.

    Args:
        filename: Path to .msh file
        physical_group: Optional physical group ID to filter elements.
                       If None, all volume elements are returned.

    Returns:
        dict with:
            'nodes': dict of {node_id: [x, y, z]}
            'elements': list of {'type': int, 'nodes': [id, ...], 'tags': [int, ...]}
            'format_version': str
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"GMSH file not found: {filename}")

    with open(filename, 'r') as f:
        lines = f.readlines()

    nodes = {}
    elements = []
    format_version = '2.0'

    idx = 0
    while idx < len(lines):
        line = lines[idx].strip()

        if line == '$MeshFormat':
            idx += 1
            parts = lines[idx].strip().split()
            format_version = parts[0]
            if not format_version.startswith('2'):
                raise ValueError(
                    f"Only GMSH 2.x format supported, got {format_version}. "
                    f"Use gmsh.option.setNumber('Mesh.MshFileVersion', 2.2) "
                    f"before saving.")
            idx += 1  # Skip to $EndMeshFormat

        elif line == '$Nodes':
            idx += 1
            n_nodes = int(lines[idx].strip())
            idx += 1
            for _ in range(n_nodes):
                parts = lines[idx].strip().split()
                nid = int(parts[0])
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                nodes[nid] = [x, y, z]
                idx += 1

        elif line == '$Elements':
            idx += 1
            n_elements = int(lines[idx].strip())
            idx += 1
            for _ in range(n_elements):
                parts = lines[idx].strip().split()
                eid = int(parts[0])
                etype = int(parts[1])
                n_tags = int(parts[2])
                tags = [int(parts[3 + t]) for t in range(n_tags)]
                node_start = 3 + n_tags
                n_elem_nodes = GMSH_NODES_PER_TYPE.get(etype, 0)
                elem_nodes = [int(parts[node_start + n])
                              for n in range(n_elem_nodes)]

                # Filter: keep only volume elements (Tet4, Hex8, Wedge6)
                if etype in (GMSH_TET4, GMSH_HEX8, GMSH_WEDGE6):
                    # Filter by physical group if specified
                    if physical_group is not None:
                        # First tag is physical group ID
                        if tags and tags[0] != physical_group:
                            idx += 1
                            continue

                    elements.append({
                        'type': etype,
                        'nodes': elem_nodes,
                        'tags': tags,
                    })
                idx += 1

        else:
            idx += 1

    return {
        'nodes': nodes,
        'elements': elements,
        'format_version': format_version,
    }


def gmsh_to_radia(filename, mu_r=1000, magnetization=None,
                   physical_group=None, unit_scale=1.0):
    """
    Read GMSH .msh file and create Radia magnetic objects.

    Args:
        filename: Path to GMSH .msh file
        mu_r: Relative permeability (for soft magnetic material)
              Set to None to skip material assignment (e.g., for permanent magnets)
        magnetization: Initial magnetization [Mx, My, Mz] in A/m.
                      Default [0, 0, 0] for soft magnetic material.
        physical_group: Optional GMSH physical group ID to import
        unit_scale: Coordinate scale factor (e.g., 0.001 for mm -> m)

    Returns:
        Radia container object ID

    Raises:
        FileNotFoundError: If file does not exist
        ValueError: If no volume elements found
        ImportError: If radia module not available
    """
    import sys
    try:
        sys.path.insert(0, os.path.dirname(__file__))
        from _radia_pybind import (ObjHexahedron, ObjTetrahedron, ObjWedge,
                                    ObjCnt, MatLin, MatApl)
    except ImportError:
        import radia as rad
        ObjHexahedron = rad.ObjHexahedron
        ObjTetrahedron = rad.ObjTetrahedron
        ObjWedge = rad.ObjWedge
        ObjCnt = rad.ObjCnt
        MatLin = rad.MatLin
        MatApl = rad.MatApl

    if magnetization is None:
        magnetization = [0, 0, 0]

    mesh_data = read_gmsh(filename, physical_group=physical_group)
    nodes = mesh_data['nodes']
    elements = mesh_data['elements']

    if not elements:
        raise ValueError(
            f"No volume elements found in {filename}. "
            f"Ensure the mesh contains Tet4 (type 4) or Hex8 (type 5) elements.")

    radia_objects = []

    for elem in elements:
        # Get vertex coordinates with unit scaling
        verts = [[nodes[nid][c] * unit_scale for c in range(3)]
                 for nid in elem['nodes']]

        if elem['type'] == GMSH_TET4:
            obj = ObjTetrahedron(verts, magnetization)
        elif elem['type'] == GMSH_HEX8:
            obj = ObjHexahedron(verts, magnetization)
        elif elem['type'] == GMSH_WEDGE6:
            obj = ObjWedge(verts, magnetization)
        else:
            continue

        # Apply material
        if mu_r is not None and mu_r > 1:
            mat = MatLin(mu_r)
            MatApl(obj, mat)

        radia_objects.append(obj)

    if len(radia_objects) == 1:
        return radia_objects[0]

    return ObjCnt(radia_objects)


def get_mesh_info(filename):
    """
    Get summary information about a GMSH mesh file.

    Args:
        filename: Path to .msh file

    Returns:
        dict with mesh statistics
    """
    mesh_data = read_gmsh(filename)

    n_tet = sum(1 for e in mesh_data['elements'] if e['type'] == GMSH_TET4)
    n_hex = sum(1 for e in mesh_data['elements'] if e['type'] == GMSH_HEX8)
    n_wedge = sum(1 for e in mesh_data['elements'] if e['type'] == GMSH_WEDGE6)

    # Bounding box
    if mesh_data['nodes']:
        coords = np.array(list(mesh_data['nodes'].values()))
        bbox_min = coords.min(axis=0)
        bbox_max = coords.max(axis=0)
        bbox_size = bbox_max - bbox_min
    else:
        bbox_min = bbox_max = bbox_size = np.zeros(3)

    # Physical groups
    phys_groups = set()
    for elem in mesh_data['elements']:
        if elem['tags']:
            phys_groups.add(elem['tags'][0])

    return {
        'format': mesh_data['format_version'],
        'n_nodes': len(mesh_data['nodes']),
        'n_tet4': n_tet,
        'n_hex8': n_hex,
        'n_wedge6': n_wedge,
        'n_total_volume': n_tet + n_hex + n_wedge,
        'bbox_min': bbox_min.tolist(),
        'bbox_max': bbox_max.tolist(),
        'bbox_size': bbox_size.tolist(),
        'physical_groups': sorted(phys_groups),
    }


# ===================================================================
# Surface mesh support (for ngbem / PEEC conductor BEM)
# ===================================================================

def read_gmsh_surface(filename, physical_group=None):
    """Read GMSH .msh file and extract SURFACE elements (Tri3, Quad4).

    Unlike read_gmsh() which extracts volume elements, this function
    extracts surface elements for PEEC/BEM conductor analysis.

    Args:
        filename: Path to .msh file
        physical_group: Optional physical group ID to filter elements

    Returns:
        dict with:
            'nodes': dict of {node_id: [x, y, z]}
            'elements': list of {'type': int, 'nodes': [...], 'tags': [...]}
            'format_version': str
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"GMSH file not found: {filename}")

    with open(filename, 'r') as f:
        lines = f.readlines()

    nodes = {}
    elements = []
    format_version = '2.0'

    idx = 0
    while idx < len(lines):
        line = lines[idx].strip()

        if line == '$MeshFormat':
            idx += 1
            parts = lines[idx].strip().split()
            format_version = parts[0]
            if not format_version.startswith('2'):
                raise ValueError(
                    f"Only GMSH 2.x format supported, got {format_version}.")
            idx += 1

        elif line == '$Nodes':
            idx += 1
            n_nodes = int(lines[idx].strip())
            idx += 1
            for _ in range(n_nodes):
                parts = lines[idx].strip().split()
                nid = int(parts[0])
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                nodes[nid] = [x, y, z]
                idx += 1

        elif line == '$Elements':
            idx += 1
            n_elements = int(lines[idx].strip())
            idx += 1
            for _ in range(n_elements):
                parts = lines[idx].strip().split()
                etype = int(parts[1])
                n_tags = int(parts[2])
                tags = [int(parts[3 + t]) for t in range(n_tags)]
                node_start = 3 + n_tags
                n_elem_nodes = GMSH_NODES_PER_TYPE.get(etype, 0)
                elem_nodes = [int(parts[node_start + n])
                              for n in range(n_elem_nodes)]

                # Filter: keep only surface elements (Tri3, Quad4)
                if etype in (GMSH_TRI3, GMSH_QUAD4):
                    if physical_group is not None:
                        if tags and tags[0] != physical_group:
                            idx += 1
                            continue

                    elements.append({
                        'type': etype,
                        'nodes': elem_nodes,
                        'tags': tags,
                    })
                idx += 1

        else:
            idx += 1

    return {
        'nodes': nodes,
        'elements': elements,
        'format_version': format_version,
    }


def gmsh_surface_to_ngsolve(filename, label="conductor",
                             physical_group=None, unit_scale=1.0):
    """Read GMSH surface mesh and create NGSolve Mesh for ngbem.

    Creates a 2D surface mesh embedded in 3D, suitable for
    HDivSurface / SurfaceL2 function spaces.

    Quad4 elements are split into 2 triangles.

    Args:
        filename: Path to GMSH .msh file
        label: Surface label for conductor boundary
        physical_group: Optional GMSH physical group to filter
        unit_scale: Coordinate scale factor (e.g., 0.001 for mm -> m)

    Returns:
        ngsolve.Mesh: Surface mesh in 3D space
    """
    from netgen.meshing import Mesh as NetgenMesh, MeshPoint, Pnt
    from netgen.meshing import Element2D, FaceDescriptor
    from ngsolve import Mesh

    mesh_data = read_gmsh_surface(filename, physical_group)
    nodes = mesh_data['nodes']
    elements = mesh_data['elements']

    if not elements:
        raise ValueError(
            f"No surface elements found in {filename}. "
            f"Check physical_group={physical_group}.")

    # Build Netgen mesh
    ngmesh = NetgenMesh(dim=3)

    # Add face descriptor
    fd = FaceDescriptor(surfnr=1, domin=0, domout=0, bc=1)
    fd.bcname = label
    ngmesh.Add(fd)

    # Add nodes (map GMSH IDs to sequential Netgen IDs, 1-indexed)
    node_map = {}
    for nid in sorted(nodes.keys()):
        coords = nodes[nid]
        scaled = [c * unit_scale for c in coords]
        pnt_id = ngmesh.Add(MeshPoint(Pnt(*scaled)))
        node_map[nid] = pnt_id

    # Add surface elements
    for elem in elements:
        verts = [node_map[nid] for nid in elem['nodes']]
        if elem['type'] == GMSH_TRI3:
            ngmesh.Add(Element2D(index=1, vertices=verts))
        elif elem['type'] == GMSH_QUAD4:
            # Split quad into 2 triangles
            ngmesh.Add(Element2D(index=1,
                                 vertices=[verts[0], verts[1], verts[2]]))
            ngmesh.Add(Element2D(index=1,
                                 vertices=[verts[0], verts[2], verts[3]]))

    ngmesh.SetBCName(0, label)

    return Mesh(ngmesh)


def extract_surface_edges(filename, physical_group=None, unit_scale=1.0):
    """Extract unique edges from GMSH surface mesh.

    Each edge becomes a potential PEEC filament (Loop DOF).
    Useful for comparing ngbem edge DOFs with PEEC filament discretization.

    Args:
        filename: Path to GMSH .msh file
        physical_group: Optional physical group filter
        unit_scale: Coordinate scale factor

    Returns:
        dict with:
            'edges': list of (node_id_1, node_id_2) tuples
            'edge_centers': (n_edges, 3) array
            'edge_directions': (n_edges, 3) array
            'edge_lengths': (n_edges,) array
            'n_edges': int
            'n_faces': int
    """
    mesh_data = read_gmsh_surface(filename, physical_group)
    nodes = mesh_data['nodes']
    elements = mesh_data['elements']

    # Collect unique edges (sorted vertex pairs)
    edge_set = set()
    for elem in elements:
        verts = elem['nodes']
        n_v = len(verts)
        for k in range(n_v):
            v0 = verts[k]
            v1 = verts[(k + 1) % n_v]
            edge = (min(v0, v1), max(v0, v1))
            edge_set.add(edge)

    edges = sorted(edge_set)

    # Compute geometry
    centers = []
    directions = []
    lengths = []
    for v0, v1 in edges:
        p0 = np.array(nodes[v0]) * unit_scale
        p1 = np.array(nodes[v1]) * unit_scale
        center = 0.5 * (p0 + p1)
        diff = p1 - p0
        length = np.linalg.norm(diff)
        direction = diff / length if length > 1e-20 else np.zeros(3)

        centers.append(center)
        directions.append(direction)
        lengths.append(length)

    return {
        'edges': edges,
        'edge_centers': np.array(centers),
        'edge_directions': np.array(directions),
        'edge_lengths': np.array(lengths),
        'n_edges': len(edges),
        'n_faces': len(elements),
    }
