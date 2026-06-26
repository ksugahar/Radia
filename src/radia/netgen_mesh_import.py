#!/usr/bin/env python
"""
netgen_mesh_import.py - Convert Netgen meshes to Radia geometry

This module provides functionality to import NGSolve/Netgen tetrahedral meshes
into Radia's magnetic field computation framework.

Author: Radia Development Team
Created: 2025-11-21
Version: 0.1.0

Functions
---------
netgen_mesh_to_radia : Convert NGSolve mesh to Radia geometry
extract_tetrahedra : Extract tetrahedral elements from mesh
create_radia_tetrahedron : Create single tetrahedron in Radia

Example
-------
>>> from ngsolve import Mesh
>>> from netgen.occ import Box, OCCGeometry
>>> import radia as rad
>>> from netgen_mesh_import import netgen_mesh_to_radia
>>>
>>> # Radia always uses meters
>>> geo = OCCGeometry(Box((0, 0, 0), (0.01, 0.01, 0.01)))
>>> mesh = Mesh(geo.GenerateMesh(maxh=0.003))
>>> mag_obj = netgen_mesh_to_radia(mesh,
...                                 material={'magnetization': [0, 0, 1.2]},
...                                 units='m')
"""

import sys
import radia as rad
from ngsolve import VOL, ET


# Standard tetrahedral face topology (1-indexed for Radia)
# Face winding from ELF/MAGIC KK4T definition for compatibility
# Netgen vertices: v0, v1, v2, v3 (0-indexed) -> 1, 2, 3, 4 (Radia 1-indexed)
# Each face is a triangle opposite to one vertex
TETRA_FACES = [
    [1, 2, 4],  # Face 0: opposite to v2 (triangle v0-v1-v3)
    [2, 3, 4],  # Face 1: opposite to v0 (triangle v1-v2-v3)
    [3, 1, 4],  # Face 2: opposite to v1 (triangle v2-v0-v3)
    [1, 3, 2],  # Face 3: opposite to v3 (triangle v0-v2-v1)
]

# Hexahedral face topology (1-indexed for Radia)
# Standard brick/hexahedron with 8 vertices
# Vertices numbered as: v0-v7 (0-indexed) -> 1-8 (Radia 1-indexed)
#
# WARNING: Hexahedral elements have known issues in Radia MMM
# ============================================================
# 1. Non-convex (concave) meshes: Hexahedral meshes from Netgen may produce
#    non-convex (concave) elements, which cause errors in Radia MMM.
#    Radia requires convex polyhedra for correct field computation.
#
# 2. Netgen 3D hex meshing limitations: Netgen's 3D hexahedral meshing
#    functionality is very limited and often produces poor quality meshes.
#
# EXCEPTION: Regular cubic hexahedral meshes (structured grids) are safe
#            and do not produce concave elements. These work correctly in Radia.
#
# Recommended: Use tetrahedral meshes for general geometries
HEX_FACES = [
    [1, 4, 3, 2],  # Bottom face (z=0)
    [5, 6, 7, 8],  # Top face (z=1)
    [1, 2, 6, 5],  # Front face (y=0)
    [3, 4, 8, 7],  # Back face (y=1)
    [1, 5, 8, 4],  # Left face (x=0)
    [2, 3, 7, 6]   # Right face (x=1)
]

# Wedge/Pentahedron face topology (1-indexed for Radia)
# Standard wedge/prism with 6 vertices (triangular prism)
# Vertices: v0-v5 (0-indexed) -> 1-6 (Radia 1-indexed)
# Bottom triangle: v0, v1, v2
# Top triangle: v3, v4, v5
WEDGE_FACES = [
    [1, 3, 2],     # Bottom triangle face (v0-v2-v1)
    [4, 5, 6],     # Top triangle face (v3-v4-v5)
    [1, 2, 5, 4],  # Quad face (v0-v1-v4-v3)
    [2, 3, 6, 5],  # Quad face (v1-v2-v5-v4)
    [3, 1, 4, 6]   # Quad face (v2-v0-v3-v5)
]

# Pyramid face topology (1-indexed for Radia)
# Standard pyramid with 5 vertices
# Vertices: v0-v4 (0-indexed) -> 1-5 (Radia 1-indexed)
# Base quad: v0, v1, v2, v3
# Apex: v4
PYRAMID_FACES = [
    [1, 4, 3, 2],  # Base quad face (v0-v3-v2-v1)
    [1, 2, 5],     # Triangle face (v0-v1-v4)
    [2, 3, 5],     # Triangle face (v1-v2-v4)
    [3, 4, 5],     # Triangle face (v2-v3-v4)
    [4, 1, 5]      # Triangle face (v3-v0-v4)
]


def create_radia_hexahedron(vertices, magnetization=None):
    """
    Create a single hexahedral (brick) polyhedron in Radia.

    WARNING: Hexahedral elements may cause numerical issues in Radia MMM.
    Tetrahedral meshes are recommended for better stability.

    Parameters
    ----------
    vertices : list of list
        8 vertices: [[x1,y1,z1], ..., [x8,y8,z8]]
        Vertices should follow standard hexahedron numbering
    magnetization : list, optional
        Magnetization vector [Mx, My, Mz] in Tesla
        Default: [0, 0, 0]

    Returns
    -------
    int
        Radia object ID

    Raises
    ------
    RuntimeError
        If Radia polyhedron creation fails

    Notes
    -----
    Uses rad.ObjHexahedron() which auto-generates face topology internally.
    """
    if magnetization is None:
        magnetization = [0, 0, 0]

    try:
        # Create hexahedron using new API (auto-generates faces)
        poly_id = rad.ObjHexahedron(vertices, magnetization)
        return poly_id
    except Exception as e:
        raise RuntimeError(
            f"Failed to create Radia hexahedron: {e}\n"
            f"Vertices: {vertices}"
        )


def create_radia_tetrahedron(vertices, magnetization=None):
    """
    Create a single tetrahedral polyhedron in Radia.

    Parameters
    ----------
    vertices : list of list
        4 vertices: [[x1,y1,z1], [x2,y2,z2], [x3,y3,z3], [x4,y4,z4]]
        Coordinates should be in meters.
    magnetization : list, optional
        Magnetization vector [Mx, My, Mz] in Tesla.
        Default: [0, 0, 0] (no magnetization)

    Returns
    -------
    int
        Radia object ID

    Raises
    ------
    RuntimeError
        If Radia polyhedron creation fails

    Notes
    -----
    Uses rad.ObjTetrahedron() which auto-generates face topology internally.
    Tetrahedra are always convex, making them ideal for Radia polyhedra.
    """
    if magnetization is None:
        magnetization = [0, 0, 0]

    if len(vertices) != 4:
        raise ValueError(f"Tetrahedron must have exactly 4 vertices, got {len(vertices)}")

    try:
        obj_id = rad.ObjTetrahedron(vertices, magnetization)
        return obj_id
    except Exception as e:
        raise RuntimeError(f"Failed to create Radia tetrahedron: {e}")


def create_radia_wedge(vertices, magnetization=None):
    """
    Create a single wedge (triangular prism, 6 vertices) polyhedron in Radia.

    Parameters
    ----------
    vertices : list of list
        6 vertices: [[x1,y1,z1], ..., [x6,y6,z6]] in NGSolve ET.PRISM order
        (bottom triangle v0,v1,v2 then top triangle v3,v4,v5, with v[i+3] above v[i]).
        Coordinates in meters.
    magnetization : list, optional
        Magnetization vector [Mx, My, Mz] in A/m.  Default: [0, 0, 0].

    Returns
    -------
    int
        Radia object ID

    Notes
    -----
    Uses rad.ObjWedge() which auto-generates face topology internally.  The NGSolve
    ET.PRISM vertex order is the rad.ObjWedge order with no reordering (validated: a
    uniform-M cube assembled from prisms matches a single ObjHexahedron to ~1e-15).
    """
    if magnetization is None:
        magnetization = [0, 0, 0]

    if len(vertices) != 6:
        raise ValueError(f"Wedge must have exactly 6 vertices, got {len(vertices)}")

    try:
        return rad.ObjWedge(vertices, magnetization)
    except Exception as e:
        raise RuntimeError(f"Failed to create Radia wedge: {e}\nVertices: {vertices}")


def compute_element_centroid(vertices):
    """
    Compute the centroid of an element from its vertices.

    Parameters
    ----------
    vertices : list of list
        List of vertex coordinates [[x1,y1,z1], [x2,y2,z2], ...]

    Returns
    -------
    list
        Centroid coordinates [cx, cy, cz]
    """
    n = len(vertices)
    cx = sum(v[0] for v in vertices) / n
    cy = sum(v[1] for v in vertices) / n
    cz = sum(v[2] for v in vertices) / n
    return [cx, cy, cz]


def extract_elements(mesh, material_filter=None, allow_hex=False, allow_wedge=False):
    """
    Extract volume elements (tetrahedra and optionally hexahedra) from NGSolve mesh.

    Parameters
    ----------
    mesh : ngsolve.Mesh
        Input mesh from Netgen/NGSolve
    material_filter : str, list of str, or None, optional
        Filter elements by material name(s):
        - str: Import only elements with this material name
        - list of str: Import only elements with these material names
        - None: Import all elements (default)
    allow_hex : bool, optional
        If True, allow hexahedral elements (ET.HEX)
        If False, raise error on non-tetrahedral elements
        Default: False
        WARNING: Hexahedral elements may cause issues in Radia MMM

    Returns
    -------
    tuple of (list of dict, int)
        (elements, skipped_count) where elements is list of dicts:
        {
            'vertices': [[x1,y1,z1], ...],  # 4 for TET, 8 for HEX
            'element_index': int,
            'element_type': str,  # 'TET' or 'HEX'
            'material': str  # Material name
        }

    Raises
    ------
    ValueError
        If mesh contains unsupported element types

    Notes
    -----
    - Vertex coordinates extracted from mesh.vertices
    - Tetrahedral elements (ET.TET): 4 vertices
    - Hexahedral elements (ET.HEX): 8 vertices (if allow_hex=True)
    - Element indices are 0-based
    - Material filtering reduces import time for multi-material meshes
    """
    # Normalize material_filter to set
    if material_filter is None:
        allowed_materials = None
    elif isinstance(material_filter, str):
        allowed_materials = {material_filter}
    elif isinstance(material_filter, (list, tuple)):
        allowed_materials = set(material_filter)
    else:
        raise ValueError(
            f"material_filter must be str, list, or None, got {type(material_filter)}"
        )

    elements = []
    skipped_count = 0
    hex_count = 0
    tet_count = 0
    wedge_count = 0

    for el_idx, el in enumerate(mesh.Elements(VOL)):
        # Check material filter
        if allowed_materials is not None:
            if el.mat not in allowed_materials:
                skipped_count += 1
                continue

        # Check element type
        if el.type == ET.TET:
            element_type = 'TET'
            expected_vertices = 4
            tet_count += 1
        elif el.type == ET.HEX and allow_hex:
            element_type = 'HEX'
            expected_vertices = 8
            hex_count += 1
        elif el.type == ET.PRISM and allow_wedge:
            element_type = 'WEDGE'
            expected_vertices = 6
            wedge_count += 1
        else:
            if el.type == ET.HEX:
                raise ValueError(
                    f"Element {el_idx} is hexahedral (ET.HEX). "
                    f"Hexahedral elements are not allowed by default. "
                    f"Set allow_hex=True to enable (WARNING: may cause MMM issues)."
                )
            elif el.type == ET.PRISM:
                raise ValueError(
                    f"Element {el_idx} is a wedge/prism (ET.PRISM). "
                    f"Wedge elements are not allowed by default. "
                    f"Set allow_wedge=True to enable (WARNING: may cause MMM issues)."
                )
            else:
                raise ValueError(
                    f"Element {el_idx} has unsupported type {el.type}. "
                    f"Only ET.TET (tetrahedra) and optionally ET.HEX (hexahedra) / "
                    f"ET.PRISM (wedges) supported."
                )

        # Extract vertex NodeId objects (NGSolve format: V0, V1, etc.)
        vert_node_ids = el.vertices

        if len(vert_node_ids) != expected_vertices:
            raise ValueError(
                f"Element {el_idx}: Expected {expected_vertices} vertices for {element_type}, "
                f"got {len(vert_node_ids)}"
            )

        # Get vertex coordinates in original NGSolve order
        vertices = []
        for v_node in vert_node_ids:
            v_idx = v_node.nr  # Get integer index from NodeId
            v = mesh.vertices[v_idx]
            coord = v.point
            vertices.append([coord[0], coord[1], coord[2]])

        elements.append({
            'vertices': vertices,
            'element_index': el_idx,
            'element_type': element_type,
            'material': el.mat
        })

    if hex_count > 0:
        print(f"[WARNING] Imported {hex_count} hexahedral elements. "
              f"Hexahedra may cause numerical issues in Radia MMM.")
    if wedge_count > 0:
        print(f"[WARNING] Imported {wedge_count} wedge/prism elements. "
              f"Wedges may cause numerical issues in Radia MMM.")

    return elements, skipped_count


def netgen_mesh_to_radia(mesh, material=None, units='m', combine=True, verbose=True,
                          material_filter=None, allow_hex=False, allow_wedge=False):
    """
    Convert NGSolve/Netgen mesh to Radia geometry.

    All coordinates must be in meters. Radia always uses meters.

    Parameters
    ----------
    mesh : ngsolve.Mesh
        Input mesh from Netgen/NGSolve
    material : dict or callable, optional
        Material specification:

        - dict: {'magnetization': [Mx, My, Mz]}
          Applies uniform magnetization to all elements

        - callable: function(element_index) -> {'magnetization': [Mx, My, Mz]}
          Per-element material specification

        - None: Use default [0, 0, 0] (no magnetization)

    units : str, default='m'
        Must be 'm'. Kept for API compatibility.

    combine : bool, default=True
        If True, return rad.ObjCnt() container of all elements.
        If False, return list of individual polyhedra object IDs.

    verbose : bool, default=True
        If True, print progress information during conversion.

    material_filter : str, list of str, or None, optional
        Filter elements by mesh material name(s):
        - str: Import only elements with this material name
        - list of str: Import only elements with these material names
        - None: Import all elements (default)

        Example: material_filter='magnetic' imports only 'magnetic' material elements

    allow_hex : bool, optional
        If True, allow hexahedral elements (ET.HEX) in addition to tetrahedra
        If False, raise error if mesh contains hexahedral elements
        Default: False

        WARNING: Hexahedral meshes may produce concave elements causing MMM errors.
        Exception: Regular cubic grids (structured hex meshes) are safe.

    Returns
    -------
    int or list
        - If combine=True: Radia container object ID (int)
        - If combine=False: List of individual polyhedron object IDs (list of int)

    Raises
    ------
    ValueError
        If mesh contains hexahedral elements and allow_hex=False
        If mesh contains unsupported element types
        If material specification is invalid
        If units parameter is not 'm' or 'mm'
    RuntimeError
        If Radia polyhedron creation fails

    Examples
    --------
    Basic usage with uniform magnetization:

    >>> from ngsolve import Mesh
    >>> from netgen.occ import Box, OCCGeometry
    >>> import radia as rad
    >>> from netgen_mesh_import import netgen_mesh_to_radia
    >>>
    >>> # Radia always uses meters
    >>>
    >>> # Create Netgen mesh
    >>> geo = OCCGeometry(Box((0, 0, 0), (0.01, 0.01, 0.01)))
    >>> mesh = Mesh(geo.GenerateMesh(maxh=0.003))
    >>>
    >>> # Convert to Radia with uniform magnetization
    >>> mag_obj = netgen_mesh_to_radia(mesh,
    ...                                 material={'magnetization': [0, 0, 1.2]},
    ...                                 units='m')
    >>> print(f"Created Radia object: {mag_obj}")

    Per-element material specification:

    >>> def material_func(el_idx):
    ...     # Left half: magnetized, right half: air
    ...     if el_idx < 100:
    ...         return {'magnetization': [0, 0, 1.2]}
    ...     else:
    ...         return {'magnetization': [0, 0, 0]}
    >>>
    >>> mag_obj = netgen_mesh_to_radia(mesh, material=material_func)

    Using millimeters (Radia default):

    >>> # All coordinates in meters
    >>> mag_obj = netgen_mesh_to_radia(mesh, units='mm')  # Auto-scales coordinates

    Notes
    -----
    - Supports tetrahedral (ET.TET) and optionally hexahedral (ET.HEX) elements
    - Vertex coordinates are extracted in Netgen's native units (meters)
    - Scaling applied automatically if units='mm'
    - All tetrahedra are convex, suitable for rad.ObjTetrahedron()
    - Hexahedra may be concave (avoid except for structured cubic grids)
    - Progress printed every 100 elements if verbose=True
    """
    # Validate units parameter
    if units not in ['m', 'mm']:
        raise ValueError(f"units must be 'm' or 'mm', got '{units}'")

    # Extract elements (tetrahedra and optionally hexahedra)
    if verbose:
        print(f"[Netgen Import] Extracting elements from mesh...")
        print(f"                Mesh: {mesh.ne} elements")
        if material_filter is not None:
            filter_str = material_filter if isinstance(material_filter, str) else ', '.join(material_filter)
            print(f"                Material filter: {filter_str}")
        if allow_hex:
            print(f"                [WARNING] Hexahedral elements enabled (may cause MMM issues)")

    try:
        elements, skipped_count = extract_elements(mesh, material_filter=material_filter,
                                                    allow_hex=allow_hex, allow_wedge=allow_wedge)
    except ValueError as e:
        print(f"[ERROR] {e}")
        raise

    num_elements = len(elements)
    if verbose:
        # Count element types
        tet_count = sum(1 for el in elements if el['element_type'] == 'TET')
        hex_count = sum(1 for el in elements if el['element_type'] == 'HEX')
        wedge_count = sum(1 for el in elements if el['element_type'] == 'WEDGE')

        if hex_count > 0 or wedge_count > 0:
            print(f"                Extracted: {num_elements} elements "
                  f"({tet_count} TET, {hex_count} HEX, {wedge_count} WEDGE)")
        else:
            print(f"                Extracted: {num_elements} tetrahedra")

        if skipped_count > 0:
            print(f"                Skipped: {skipped_count} elements (filtered by material)")

    # Radia always uses meters - no coordinate scaling needed
    coord_scale = 1.0

    if verbose:
        print(f"                Units: {units} (Radia always uses meters)")

    # Process material specification
    if material is None:
        # Default: no magnetization
        def get_material(el_idx):
            return {'magnetization': [0, 0, 0]}
    elif isinstance(material, dict):
        # Uniform material
        def get_material(el_idx):
            return material
    elif callable(material):
        # Per-element material function
        get_material = material
    else:
        raise ValueError(
            f"material must be dict, callable, or None. Got {type(material)}"
        )

    # Create Radia polyhedra
    polyhedra = []

    if verbose:
        print(f"[Netgen Import] Creating Radia polyhedra...")

    for i, element in enumerate(elements):
        # Scale coordinates if needed
        vertices = element['vertices']
        if coord_scale != 1.0:
            vertices = [[x*coord_scale, y*coord_scale, z*coord_scale]
                       for x, y, z in vertices]

        # Get material for this element
        el_idx = element['element_index']
        try:
            mat = get_material(el_idx)
        except Exception as e:
            raise RuntimeError(
                f"Material function failed for element {el_idx}: {e}"
            )

        # Validate material specification
        if not isinstance(mat, dict) or 'magnetization' not in mat:
            raise ValueError(
                f"Material function for element {el_idx} must return "
                f"dict with 'magnetization' key. Got: {mat}"
            )

        magnetization = mat['magnetization']

        # Create polyhedron based on element type
        element_type = element['element_type']
        try:
            if element_type == 'TET':
                obj_id = create_radia_tetrahedron(vertices, magnetization)
            elif element_type == 'HEX':
                obj_id = create_radia_hexahedron(vertices, magnetization)
            elif element_type == 'WEDGE':
                obj_id = create_radia_wedge(vertices, magnetization)
            else:
                raise ValueError(f"Unknown element type: {element_type}")

            polyhedra.append(obj_id)
        except RuntimeError as e:
            raise RuntimeError(
                f"Failed to create {element_type} polyhedron for element {el_idx}: {e}"
            )

        # Progress reporting
        if verbose and (i + 1) % 100 == 0:
            print(f"                Progress: {i+1}/{num_elements}", end='\r')

    if verbose:
        print(f"                Progress: {num_elements}/{num_elements}")
        print(f"                [OK] Created {len(polyhedra)} polyhedra")

    # Return result
    if combine:
        if verbose:
            print(f"[Netgen Import] Combining into container...")

        container = rad.ObjCnt(polyhedra)

        if verbose:
            print(f"                [OK] Container object ID: {container}")

        return container
    else:
        if verbose:
            print(f"[Netgen Import] Returning list of {len(polyhedra)} object IDs")

        return polyhedra


def cubit_hex_to_radia(hex_elements, magnetization=None, mu_r=None, combine=True, verbose=True):
    """
    Convert Cubit hexahedral mesh data to Radia geometry for CplMag solver.

    This function takes hex element data (list of vertex lists) from Cubit
    and creates Radia ObjHexahedron objects suitable for CplMag (PEEC-MMM coupling).

    Parameters
    ----------
    hex_elements : list of list
        List of hexahedral elements, each element is a list of 8 vertices:
        [[[x1,y1,z1], [x2,y2,z2], ..., [x8,y8,z8]], ...]
        Vertices should be in meters.

    magnetization : list, optional
        Initial magnetization vector [Mx, My, Mz] in A/m.
        Default: [0, 0, 0] (soft magnetic material, no remanent magnetization)

    mu_r : float, optional
        Relative permeability for soft magnetic material.
        If provided, creates and applies MatLin(mu_r) to all elements.
        Default: None (no material applied, user must apply material separately)

    combine : bool, default=True
        If True, return rad.ObjCnt() container of all elements.
        If False, return list of individual hexahedron object IDs.

    verbose : bool, default=True
        If True, print progress information during conversion.

    Returns
    -------
    int or list
        - If combine=True: Radia container object ID (int)
        - If combine=False: List of individual hexahedron object IDs (list of int)

    Examples
    --------
    Basic usage with Cubit hex mesh:

    >>> import radia as rad
    >>> from netgen_mesh_import import cubit_hex_to_radia
    >>>
    >>> # Radia always uses meters
    >>>
    >>> # Hex elements from Cubit (or manual creation)
    >>> hex_elements = [
    ...     [[-0.01,-0.01,-0.01], [0.01,-0.01,-0.01], [0.01,0.01,-0.01], [-0.01,0.01,-0.01],
    ...      [-0.01,-0.01,0.01], [0.01,-0.01,0.01], [0.01,0.01,0.01], [-0.01,0.01,0.01]],
    ... ]
    >>>
    >>> # Create Radia objects with mu_r=1000
    >>> core = cubit_hex_to_radia(hex_elements, mu_r=1000)

    Notes
    -----
    - All coordinates must be in meters
    - Vertices should follow standard hexahedron numbering convention
    - Multi-element meshes enable proper MMM demagnetization coupling
    - Legacy CplMag APIs (CndLoop, CplMagCreate, etc.) have been removed.
      Use PEECBuilder + CoupledPEECSolver instead.

    See Also
    --------
    netgen_mesh_to_radia : For Netgen mesh import
    create_radia_hexahedron : Create single hexahedron
    """
    if magnetization is None:
        magnetization = [0, 0, 0]

    if verbose:
        print(f"[Cubit Hex Import] Processing {len(hex_elements)} hexahedral elements...")

    polyhedra = []

    for i, vertices in enumerate(hex_elements):
        if len(vertices) != 8:
            raise ValueError(
                f"Element {i}: Expected 8 vertices for hexahedron, got {len(vertices)}"
            )

        try:
            obj_id = create_radia_hexahedron(vertices, magnetization)
            polyhedra.append(obj_id)
        except RuntimeError as e:
            raise RuntimeError(f"Failed to create hexahedron for element {i}: {e}")

        if verbose and (i + 1) % 100 == 0:
            print(f"                  Progress: {i+1}/{len(hex_elements)}", end='\r')

    if verbose:
        print(f"                  Progress: {len(hex_elements)}/{len(hex_elements)}")
        print(f"                  [OK] Created {len(polyhedra)} hexahedra")

    # Apply material if mu_r is specified
    if mu_r is not None:
        if verbose:
            print(f"[Cubit Hex Import] Applying material (mu_r = {mu_r})...")

        mat = rad.MatLin(mu_r)
        for obj_id in polyhedra:
            rad.MatApl(obj_id, mat)

        if verbose:
            print(f"                  [OK] Material applied to all elements")

    # Return result
    if combine:
        if verbose:
            print(f"[Cubit Hex Import] Combining into container...")

        container = rad.ObjCnt(polyhedra)

        if verbose:
            print(f"                  [OK] Container object ID: {container}")

        return container
    else:
        if verbose:
            print(f"[Cubit Hex Import] Returning list of {len(polyhedra)} object IDs")

        return polyhedra


# Module-level constants for external use
__version__ = '0.3.0'
__all__ = [
    'netgen_mesh_to_radia',
    'extract_elements',
    'compute_element_centroid',
    'create_radia_tetrahedron',
    'create_radia_hexahedron',
    'create_radia_wedge',
    'cubit_hex_to_radia',
    'TETRA_FACES',
    'HEX_FACES',
    'WEDGE_FACES',
    'PYRAMID_FACES'
]
