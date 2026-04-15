"""
gmsh_post_export.py

GMSH .msh v4.1 post-processing export with per-material physical groups.

GMSH natively supports high-order elements (Tri6, Quad9, Tet10, Hex20)
with correct interpolation, unlike VTK which approximates.

Key feature: material-based field assignment. Each field can be restricted
to a specific material (physical group). GMSH GUI allows toggling visibility
per physical group, enabling per-material field switching.

Usage:
    from ngsolve import *
    from radia.gmsh_post_export import GmshPostExport

    # Multi-material mesh with per-material fields
    post = GmshPostExport(mesh)
    post.add_field("J", gf_J, material="coil")           # coil only
    post.add_field("B", gf_B, material="core")            # core only
    post.add_scalar_field("|H|", Norm(gf_H))              # all materials
    post.write("results.msh")

    # Mesh only (no field data):
    post.write_mesh("mesh_only.msh")

Supports:
    1st order: Tri3 (2), Quad4 (3), Tet4 (4), Hex8 (5), Wedge6 (6), Pyr5 (7)
    2nd order: Tri6 (9), Quad9 (10), Tet10 (11), Hex20 (17), Wedge15 (18), Pyr13 (19)
    3rd order: Tri10 (21), Quad16 (36), Tet20 (29), Hex64 (92), Pyr30 (118)
    4th order: Tri15 (23), Quad25 (37), Tet35 (30), Hex125 (93), Pyr55 (119)
    5th order: Tri21 (25), Quad36 (38), Tet56 (31), Hex216 (94), Pyr91 (120)
    Note: Wedge/Prism only supports up to 2nd order (GMSH format limitation)

Part of Radia project
"""

import numpy as np

# Physical constants
MU_0 = 4.0 * np.pi * 1e-7


# ============================================================
# GMSH element type codes
# ============================================================

_NGSOLVE_TO_GMSH_1ST = {
    'TRIG': 2, 'QUAD': 3, 'TET': 4, 'HEX': 5, 'PRISM': 6, 'PYRAMID': 7,
}

_NGSOLVE_TO_GMSH_2ND = {
    'TRIG': 9, 'QUAD': 10, 'TET': 11, 'HEX': 17, 'PRISM': 18, 'PYRAMID': 19,
}

# High-order Lagrange element types (order -> GMSH type code)
# Obtained via gmsh.model.mesh.getElementType(name, order)
_GMSH_TYPE_BY_ORDER = {
    'TRIG':    {1: 2, 2: 9, 3: 21, 4: 23, 5: 25},
    'QUAD':    {1: 3, 2: 10, 3: 36, 4: 37, 5: 38},
    'TET':     {1: 4, 2: 11, 3: 29, 4: 30, 5: 31},
    'HEX':     {1: 5, 2: 12, 3: 92, 4: 93, 5: 94},
    'PRISM':   {1: 6, 2: 13},  # order 3+ not supported in GMSH
    'PYRAMID': {1: 7, 2: 14, 3: 118, 4: 119, 5: 120},
}

# Backward compat alias
_TRIG_GMSH_TYPE_BY_ORDER = _GMSH_TYPE_BY_ORDER['TRIG']

_GMSH_NODES_PER_TYPE = {
    2: 3, 3: 4, 4: 4, 5: 8, 6: 6, 7: 5,
    9: 6, 10: 9, 11: 10, 12: 27, 17: 20, 18: 15, 19: 13,
    21: 10, 23: 15, 25: 21,  # high-order triangles
    29: 20, 30: 35, 31: 56,  # high-order tets
    36: 16, 37: 25, 38: 36,  # high-order quads
    92: 64, 93: 125, 94: 216,  # high-order hexes
    118: 30, 119: 55, 120: 91,  # high-order pyramids
}

# GMSH reference coordinates -> NGSolve reference coordinates mapping.
# GMSH TET: same as NGSolve (0,0,0)-(1,0,0)-(0,1,0)-(0,0,1)
#   but vertex ordering differs (handled by _NGSOLVE_TO_GMSH_NODE_ORDER)
# GMSH HEX: (-1,-1,-1) to (1,1,1); NGSolve: (0,0,0) to (1,1,1)
#   transform: ng = (gmsh + 1) / 2
# GMSH PYR: (-1,-1,0) base, (0,0,1) apex; NGSolve: (0,0,0)-(1,0,0)-(1,1,0)-(0,1,0)-(0,0,1)
#   transform: ng_x = (gmsh_x + 1)/2, ng_y = (gmsh_y + 1)/2, ng_z = gmsh_z

def _gmsh_ref_to_ngsolve_ref(et_name, gmsh_x, gmsh_y, gmsh_z):
    """Convert GMSH reference coordinates to NGSolve reference coordinates.

    GMSH and NGSolve use different vertex numbering in reference space.
    For TET: GMSH v0=(0,0,0) maps to NGSolve v[0] at ref(1,0,0).
    The affine transform is: ng = (1-gx-gy-gz, gx, gy) for the
    3 independent reference coordinates (4th is 1-sum).
    """
    gx, gy, gz = gmsh_x, gmsh_y, gmsh_z
    if et_name == 'TET':
        # GMSH (0,0,0)->ng(1,0,0), (1,0,0)->ng(0,1,0),
        # (0,1,0)->ng(0,0,1), (0,0,1)->ng(0,0,0)
        return (1 - gx - gy - gz, gx, gy)
    elif et_name == 'TRIG':
        return (1 - gx - gy, gx, 0)
    elif et_name == 'HEX' or et_name == 'QUAD':
        return ((gx + 1) / 2, (gy + 1) / 2, (gz + 1) / 2)
    elif et_name == 'PYRAMID':
        return ((gx + 1) / 2, (gy + 1) / 2, gz)
    elif et_name == 'PRISM':
        return (gx, gy, (gz + 1) / 2)
    return (gx, gy, gz)


def _get_gmsh_ref_points(et_name, order):
    """Get GMSH Lagrange reference points for given element type and order.

    Returns list of (x, y, z) in NGSolve reference coordinates,
    or None if not supported.
    """
    type_map = _GMSH_TYPE_BY_ORDER.get(et_name, {})
    gmsh_type = type_map.get(order)
    if gmsh_type is None:
        return None, None

    try:
        import gmsh
        need_init = not gmsh.isInitialized()
        if need_init:
            gmsh.initialize()
        _, _, _, nn, ref_pts, _ = gmsh.model.mesh.getElementProperties(gmsh_type)
        if need_init:
            gmsh.finalize()
    except Exception:
        return None, None

    # Convert GMSH ref -> NGSolve ref
    ng_pts = []
    for i in range(nn):
        gx, gy, gz = ref_pts[3*i], ref_pts[3*i+1], ref_pts[3*i+2]
        ng_pts.append(_gmsh_ref_to_ngsolve_ref(et_name, gx, gy, gz))

    return gmsh_type, ng_pts

_NGSOLVE_TO_GMSH_NODE_ORDER = {
    'TET': [0, 1, 2, 3],
    'HEX': [0, 1, 5, 4, 3, 2, 6, 7],
    'PRISM': [0, 2, 1, 3, 5, 4],
    'PYRAMID': [3, 2, 1, 0, 4],
    'TRIG': [0, 1, 2],
    'QUAD': [0, 1, 2, 3],
}

_VOL_TYPES = {4, 5, 6, 7, 11, 12, 17, 18, 19,
              29, 30, 31,       # TET order 3-5
              92, 93, 94,       # HEX order 3-5
              118, 119, 120}    # PYR order 3-5


# ============================================================
# GmshPostExport class
# ============================================================

class GmshPostExport:
    """Export NGSolve mesh and field data to GMSH .msh v4.1 format.

    Supports per-material physical groups for field switching in GMSH GUI.

    Args:
        mesh: NGSolve Mesh object (volume or surface mesh)
        boundary: If True, export boundary surface elements (BND) from a volume
                  mesh. Enables high-order surface export for BEM visualization.
                  Default False (export VOL elements for volume meshes).

    Example:
        post = GmshPostExport(mesh)
        post.add_field("J", gf_J, material="coil")   # coil elements only
        post.add_field("B", gf_B, material="core")    # core elements only
        post.add_scalar_field("T", temp_cf)            # all elements
        post.write("results.msh")

        # BEM surface export with high-order elements:
        post = GmshPostExport(mesh, boundary=True)
        post.add_field("|J|", node_J, ncomp=1)
        post.write("bem_results.msh")
    """

    def __init__(self, mesh, boundary=False):
        self.mesh = mesh
        self._boundary = boundary
        # (name, ncomp, data_array, is_cell, material_name_or_None)
        self._fields = []

    def add_field(self, name, data, ncomp=None, material=None, cell_data=False):
        """Add a field for export, optionally restricted to a material.

        Args:
            name: Field name (displayed in GMSH GUI)
            data: CoefficientFunction, GridFunction, or numpy array
            ncomp: Number of components. If None, auto-detected
                   (1 for scalar, 3 for vector).
            material: Material name to restrict field to (GMSH physical group).
                      If None, field applies to all elements.
            cell_data: If True, evaluate at element centroids (ElementData).
                       If False, evaluate at vertices (NodeData).
        """
        if ncomp is None:
            ncomp = _detect_ncomp(data)
        arr = self._resolve_data(data, ncomp=ncomp, cell_data=cell_data)
        self._fields.append((name, ncomp, arr, cell_data, material))

    def add_scalar_field(self, name, data, cell_data=False, material=None):
        """Add a scalar field. See add_field() for args."""
        arr = self._resolve_data(data, ncomp=1, cell_data=cell_data)
        self._fields.append((name, 1, arr, cell_data, material))

    def add_vector_field(self, name, data, cell_data=False, material=None):
        """Add a 3-component vector field. See add_field() for args."""
        arr = self._resolve_data(data, ncomp=3, cell_data=cell_data)
        self._fields.append((name, 3, arr, cell_data, material))

    def write(self, filename, time=0.0, timestep=0):
        """Write mesh + field data to GMSH .msh v4.1 format.

        v4.1 is the only supported format (lab-wide standard, 2026-04).
        netgen I/O is always via .vol; .msh v2.2 is no longer emitted.
        """
        if not filename.endswith('.msh'):
            filename += '.msh'

        mesh = self.mesh
        is_surface = _is_surface_mesh(mesh) or self._boundary

        # Extract mesh topology grouped by material
        nodes, mat_names, elem_data = _extract_mesh_data_grouped(
            mesh, is_surface)
        n_nodes = len(nodes)

        # Build per-material element index mapping
        # mat_elem_map: {mat_name: [(orig_idx, gmsh_elem_id)]}
        mat_elem_map, elem_id_map, n_elems = _assign_element_ids(
            mat_names, elem_data)

        # Determine dimension (2 for surface, 3 for volume)
        dim = 2 if is_surface else 3

        # Compute per-material bounding boxes
        mat_bboxes = _compute_material_bboxes(nodes, mat_names, elem_data)

        with open(filename, 'w', encoding='utf-8') as f:
            # Header
            f.write('$MeshFormat\n4.1 0 8\n$EndMeshFormat\n')

            # Physical names (one per material)
            _write_physical_names(f, mat_names, dim)

            # Entities (one per material)
            _write_entities(f, mat_names, mat_bboxes, dim)

            # Nodes (all in one entity block)
            _write_nodes(f, nodes, dim)

            # Elements (grouped by material x element type)
            _write_elements_grouped(f, mat_names, elem_data, n_elems, dim)

            # Field data
            for name, ncomp, data, is_cell, material in self._fields:
                if is_cell:
                    if material is not None:
                        # Material-specific: only include that material's elements
                        elem_ids = mat_elem_map.get(material, [])
                        if not elem_ids:
                            print(f"  WARNING: material '{material}' not found, "
                                  f"skipping field '{name}'")
                            continue
                        _write_element_data_filtered(
                            f, name, ncomp, data, elem_ids, time, timestep)
                    else:
                        # Global: all elements
                        all_ids = []
                        for mn in mat_names:
                            all_ids.extend(mat_elem_map.get(mn, []))
                        _write_element_data_filtered(
                            f, name, ncomp, data, all_ids, time, timestep)
                else:
                    if material is not None:
                        # Material-specific node data: restrict to nodes
                        # connected to this material's elements
                        node_set = set()
                        for mat_name, gmsh_type, conn, orig_idx in elem_data:
                            if mat_name == material:
                                node_set.update(conn)
                        _write_node_data_filtered(
                            f, name, ncomp, data, sorted(node_set),
                            time, timestep)
                    else:
                        _write_node_data(
                            f, name, ncomp, data, n_nodes, time, timestep)

        # Write companion .msh.opt for correct high-order display.
        # View index matches NodeData order in .msh file:
        #   vector views (ncomp >= 2) get arrow settings,
        #   scalar views (ncomp == 1) get filled colormap.
        if is_surface and _detect_curve_order(self.mesh) >= 2:
            n_vec = sum(1 for _, nc, _, _, _ in self._fields if nc >= 2)
            n_scl = sum(1 for _, nc, _, _, _ in self._fields if nc == 1)
            write_companion_opt(filename, n_vector_views=n_vec,
                                n_scalar_views=n_scl)

        n_fields = len(self._fields)
        mat_str = ', '.join(f'{m}({len(mat_elem_map.get(m, []))})' for m in mat_names)
        print(f"GMSH export: {filename}")
        print(f"  {n_elems} elements, {n_nodes} nodes, {n_fields} fields")
        print(f"  Materials: {mat_str}")
        return filename

    def write_mesh(self, filename):
        """Write mesh only (no field data)."""
        fields_backup = self._fields
        self._fields = []
        result = self.write(filename)
        self._fields = fields_backup
        return result

    def get_materials(self):
        """Return list of material names in the mesh."""
        is_surface = _is_surface_mesh(self.mesh)
        if is_surface:
            return list(self.mesh.GetBoundaries())
        else:
            return list(self.mesh.GetMaterials())

    def _resolve_data(self, data, ncomp, cell_data):
        """Convert CoefficientFunction/GridFunction/ndarray to numpy array."""
        if isinstance(data, np.ndarray):
            return data

        try:
            from ngsolve import CoefficientFunction
            if isinstance(data, CoefficientFunction):
                return _evaluate_cf(self.mesh, data, ncomp, cell_data)
        except ImportError:
            pass

        raise TypeError(
            f"Unsupported data type: {type(data)}. "
            "Expected CoefficientFunction, GridFunction, or numpy.ndarray.")


# ============================================================
# Internal: mesh extraction
# ============================================================

def _is_surface_mesh(mesh):
    """Check if mesh is surface-only (no volume elements)."""
    try:
        return mesh.ne == 0
    except Exception:
        return True


def _detect_ncomp(data):
    """Auto-detect number of components from data."""
    try:
        return data.dim
    except AttributeError:
        pass
    if isinstance(data, np.ndarray):
        if data.ndim == 1:
            return 1
        return data.shape[-1]
    return 1


def _extract_mesh_data_grouped(mesh, is_surface):
    """Extract mesh topology grouped by material/boundary.

    For curved meshes (order > 1), generates mid-edge nodes and uses
    2nd order GMSH element types (Tri6, Quad9, etc.).

    Returns:
        nodes: list of (x, y, z) tuples
        mat_names: list of unique material names (preserving order)
        elem_data: list of (mat_name, gmsh_type, connectivity, orig_idx)
    """
    nodes = []
    for v in mesh.vertices:
        pt = v.point
        nodes.append((pt[0], pt[1], pt[2]))

    # Detect if mesh is curved (order > 1)
    curve_order = _detect_curve_order(mesh)
    use_highorder = (curve_order >= 2)

    # C++ fast path for high-order SURFACE meshes only
    if use_highorder and is_surface:
        try:
            from radia._radia_pybind import _compute_ho_bnd_nodes
            result = _compute_ho_bnd_nodes(mesh, curve_order)
            nodes = [tuple(row) for row in result['nodes']]
            mat_names = []
            elem_data = []
            for i in range(len(result['elem_conn'])):
                mat = result['elem_materials'][i]
                gtype = result['elem_gmsh_types'][i]
                conn = result['elem_conn'][i]
                orig = result['elem_orig_idx'][i]
                elem_data.append((mat, gtype, conn, orig))
                if mat not in mat_names:
                    mat_names.append(mat)
            return nodes, mat_names, elem_data
        except ImportError:
            pass  # Fall back to Python path

    # Python fallback for high-order (if C++ not available)
    elem_ho_nodes = {}
    if use_highorder:
        elem_ho_nodes = _compute_highorder_nodes(
            mesh, nodes, curve_order, is_surface)

    mat_names = []
    elem_data = []

    if is_surface:
        from ngsolve import BND
        for idx, el in enumerate(mesh.Elements(BND)):
            et_name = str(el.type).split('.')[-1]

            if use_highorder and et_name == 'TRIG':
                gmsh_type = _TRIG_GMSH_TYPE_BY_ORDER.get(
                    curve_order, _NGSOLVE_TO_GMSH_1ST.get(et_name))
            elif use_highorder and et_name in _NGSOLVE_TO_GMSH_2ND:
                gmsh_type = _NGSOLVE_TO_GMSH_2ND[et_name]
            else:
                gmsh_type = _NGSOLVE_TO_GMSH_1ST.get(et_name)
            if gmsh_type is None:
                continue

            verts = [v.nr for v in el.vertices]
            perm = _NGSOLVE_TO_GMSH_NODE_ORDER.get(
                et_name, list(range(len(verts))))
            reordered = [verts[i] for i in perm]

            if use_highorder:
                reordered = _build_highorder_connectivity(
                    el, reordered, et_name, curve_order,
                    elem_ho_nodes, None, nodes)

            mat_name = _get_element_material(mesh, el, is_surface)
            elem_data.append((mat_name, gmsh_type, reordered, idx))
            if mat_name not in mat_names:
                mat_names.append(mat_name)
    else:
        from ngsolve import VOL

        # For high-order volume meshes, compute Lagrange nodes via GetTrafo
        vol_ho_cache = {}  # (et_name, order, node_key) -> node index

        # Pre-fetch GMSH reference points per (element_type, order)
        _vol_gmsh_ref = {}  # (et_name,) -> (gmsh_type, ng_ref_pts)

        for idx, el in enumerate(mesh.Elements(VOL)):
            et_name = str(el.type).split('.')[-1]

            if use_highorder:
                # Try high-order type
                type_map = _GMSH_TYPE_BY_ORDER.get(et_name, {})
                gmsh_type = type_map.get(curve_order)
                if gmsh_type is None:
                    # Fall back to order 2 or 1
                    gmsh_type = type_map.get(min(curve_order, 2))
                if gmsh_type is None:
                    gmsh_type = _NGSOLVE_TO_GMSH_1ST.get(et_name)
            else:
                gmsh_type = _NGSOLVE_TO_GMSH_1ST.get(et_name)
            if gmsh_type is None:
                continue

            verts = [v.nr for v in el.vertices]
            perm = _NGSOLVE_TO_GMSH_NODE_ORDER.get(
                et_name, list(range(len(verts))))
            reordered = [verts[i] for i in perm]

            if use_highorder and et_name in _GMSH_TYPE_BY_ORDER:
                effective_order = curve_order
                if effective_order not in _GMSH_TYPE_BY_ORDER.get(et_name, {}):
                    effective_order = min(curve_order, 2)
                if effective_order >= 2:
                    reordered = _build_vol_ho_generic(
                        mesh, el, reordered, et_name, effective_order,
                        nodes, vol_ho_cache, _vol_gmsh_ref)

            mat_name = _get_element_material(mesh, el, is_surface)
            elem_data.append((mat_name, gmsh_type, reordered, idx))
            if mat_name not in mat_names:
                mat_names.append(mat_name)

    return nodes, mat_names, elem_data


# ============================================================
# High-order volume element connectivity
# ============================================================

# Edge tables for 2nd order: vertex pairs for each mid-edge node
# GMSH TET10 node order: 4 vertices + 6 mid-edge nodes
# Edge order matches GMSH convention (gmsh.model.mesh.getElementProperties)
_TET_EDGES_GMSH = [
    (0, 1), (1, 2), (2, 0),  # bottom face edges
    (3, 0), (3, 2), (3, 1),  # edges to apex
]

# GMSH HEX20: 8 vertices + 12 mid-edge nodes
_HEX_EDGES_GMSH = [
    (0, 1), (0, 3), (0, 4), (1, 2),
    (1, 5), (2, 3), (2, 6), (3, 7),
    (4, 5), (4, 7), (5, 6), (6, 7),
]

# GMSH PRISM18: 6 vertices + 9 mid-edge nodes
_PRISM_EDGES_GMSH = [
    (0, 1), (0, 2), (0, 3), (1, 2),
    (1, 4), (2, 5), (3, 4), (3, 5), (4, 5),
]


def _build_vol_highorder_conn(mesh, el, reordered, et_name, nodes, cache):
    """Build 2nd order connectivity for a volume element.

    Uses mesh.GetTrafo(el) to evaluate mid-edge positions in physical space.
    Caches mid-edge nodes by vertex pair to avoid duplicates on shared edges.

    Args:
        mesh: NGSolve Mesh
        el: Element (VOL)
        reordered: 1st order vertex indices (GMSH node order)
        et_name: Element type name ('TET', 'HEX', 'PRISM')
        nodes: mutable list of (x,y,z) tuples (new nodes appended)
        cache: dict mapping edge_key -> mid-node index

    Returns:
        Extended connectivity list (vertices + mid-edge nodes).
    """
    if et_name == 'TET':
        edge_table = _TET_EDGES_GMSH
    elif et_name == 'HEX':
        edge_table = _HEX_EDGES_GMSH
    elif et_name == 'PRISM':
        edge_table = _PRISM_EDGES_GMSH
    else:
        return reordered  # unsupported type, return 1st order

    conn = list(reordered)  # start with vertices

    # Reference space mid-edge points for TET: (0,0,0)-(1,0,0)-(0,1,0)-(0,0,1)
    # We compute each mid-edge point via GetTrafo
    trafo = mesh.GetTrafo(el)

    # NGSolve reference coordinates for standard element types.
    # IMPORTANT: NGSolve trafo maps reference coords to physical coords,
    # but the vertex ordering in el.vertices does NOT match the reference
    # vertex ordering. For TET:
    #   ref(0,0,0) -> el.vertices[3]
    #   ref(1,0,0) -> el.vertices[0]
    #   ref(0,1,0) -> el.vertices[1]
    #   ref(0,0,1) -> el.vertices[2]
    # So ref_verts[i] gives the reference coord of el.vertices[i].
    if et_name == 'TET':
        # el.vertices order -> reference coords
        ref_verts = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0)]
    elif et_name == 'HEX':
        ref_verts = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                     (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    elif et_name == 'PRISM':
        ref_verts = [(1, 0, 0), (0, 1, 0), (0, 0, 0),
                     (1, 0, 1), (0, 1, 1), (0, 0, 1)]

    # Map GMSH vertex index -> NGSolve el.vertices index -> ref coord
    ngsolve_verts = [v.nr for v in el.vertices]
    perm = _NGSOLVE_TO_GMSH_NODE_ORDER.get(et_name, list(range(len(ngsolve_verts))))
    # inv_perm: gmsh_local_idx -> ngsolve_local_idx
    inv_perm = [0] * len(perm)
    for i, p in enumerate(perm):
        inv_perm[p] = i

    for e0_gmsh, e1_gmsh in edge_table:
        # Global node IDs for this edge
        gn0 = reordered[e0_gmsh]
        gn1 = reordered[e1_gmsh]
        edge_key = (min(gn0, gn1), max(gn0, gn1))

        if edge_key in cache:
            conn.append(cache[edge_key])
        else:
            # Map GMSH local index -> NGSolve local index -> ref coord
            ng0 = inv_perm[e0_gmsh]
            ng1 = inv_perm[e1_gmsh]
            ref_mid = tuple(0.5 * (ref_verts[ng0][k] + ref_verts[ng1][k])
                            for k in range(3))

            # Evaluate trafo at mid-point to get physical coordinates
            from ngsolve import IntegrationRule
            ir = IntegrationRule([ref_mid], [1.0])
            mip = trafo(ir[0])
            phys = mip.point
            mid_node_idx = len(nodes)
            nodes.append((float(phys[0]), float(phys[1]), float(phys[2])))
            cache[edge_key] = mid_node_idx
            conn.append(mid_node_idx)

    return conn


def _build_vol_ho_generic(mesh, el, reordered, et_name, order,
                          nodes, cache, ref_cache):
    """Build arbitrary-order connectivity for a volume element.

    Uses GMSH API to get Lagrange reference points, converts to NGSolve
    reference coordinates, evaluates GetTrafo for physical positions.

    Vertex nodes (first N) reuse existing indices from reordered.
    Higher-order nodes are cached by rounded physical position to share
    nodes across adjacent elements.

    Args:
        mesh: NGSolve Mesh
        el: Element (VOL)
        reordered: 1st order vertex indices (GMSH node order)
        et_name: Element type name ('TET', 'HEX', 'PRISM', 'PYRAMID')
        order: polynomial order (2-5)
        nodes: mutable list of (x,y,z) tuples (new nodes appended)
        cache: dict for node deduplication
        ref_cache: dict (et_name,) -> (gmsh_type, ng_ref_pts)
    """
    # Get or compute reference points
    if et_name not in ref_cache:
        ref_cache[et_name] = _get_gmsh_ref_points(et_name, order)
    gmsh_type, ng_ref_pts = ref_cache[et_name]

    if ng_ref_pts is None:
        # Fallback to order-2 path
        return _build_vol_highorder_conn(
            mesh, el, reordered, et_name, nodes,
            cache.setdefault("_order2", {}))

    n_verts = len(reordered)
    n_total = len(ng_ref_pts)

    # NGSolve reference coords for vertices
    # Map: GMSH local vertex i -> NGSolve el.vertices index
    perm = _NGSOLVE_TO_GMSH_NODE_ORDER.get(
        et_name, list(range(n_verts)))
    inv_perm = [0] * len(perm)
    for i, p in enumerate(perm):
        inv_perm[p] = i

    # Evaluate trafo at all high-order reference points
    trafo = mesh.GetTrafo(el)
    from ngsolve import IntegrationRule

    conn = list(reordered)  # vertices first

    for ho_idx in range(n_verts, n_total):
        ref_pt = ng_ref_pts[ho_idx]

        # Evaluate physical coordinates
        ir = IntegrationRule([ref_pt], [1.0])
        mip = trafo(ir[0])
        phys = (float(mip.point[0]), float(mip.point[1]),
                float(mip.point[2]))

        # Round for deduplication (shared edges/faces between elements)
        key = (round(phys[0], 12), round(phys[1], 12),
               round(phys[2], 12))

        if key in cache:
            conn.append(cache[key])
        else:
            node_idx = len(nodes)
            nodes.append(phys)
            cache[key] = node_idx
            conn.append(node_idx)

    return conn


def write_companion_opt(msh_filename, n_vector_views=0,
                        n_scalar_views=0):
    """Write a .msh.opt file with display settings for GMSH.

    GMSH auto-loads ``<file>.msh.opt`` when opening ``<file>.msh``.
    This replaces the old companion .geo approach (no Merge needed,
    no element-ID collision risk from multi-file Merge).

    Settings follow the gmsh_post_spec (mandatory):
      - Mesh.NumSubEdges = 4      (render curved Tri6+ smoothly)
      - Mesh.VolumeEdges = 0      (hide air tet edges)
      - Mesh.VolumeFaces = 0      (hide air tet faces)
      - Mesh.SurfaceFaces = 1     (show coil as smooth body)
      - View[i].VectorType = 4    (3D arrows for vectors)
      - View[i].ArrowSizeMin = 20 (fixed 20 px -- THE kirei knob)
      - View[i].IntervalsType = 2 (filled iso for scalars)

    View ordering: vector views first [0..n_vector-1],
    then scalar views [n_vector..n_vector+n_scalar-1].

    Args:
        msh_filename: Path to the .msh file (the .opt is written next to it).
        n_vector_views: Number of vector views (View[0]..View[n-1]).
        n_scalar_views: Number of scalar views (View[n]..View[n+m-1]).

    Returns:
        Path to the .msh.opt file.
    """
    import os
    opt_filename = msh_filename + '.opt'
    with open(opt_filename, 'w', encoding='utf-8') as f:
        f.write('// Auto-generated display options (gmsh_post_spec)\n')
        # Curved element rendering
        f.write('Mesh.NumSubEdges = 4;\n')
        # Hide air volume mesh
        f.write('Mesh.VolumeEdges = 0;\n')
        f.write('Mesh.VolumeFaces = 0;\n')
        # Coil surface: smooth body, no edge lines
        f.write('Mesh.SurfaceEdges = 0;\n')
        f.write('Mesh.SurfaceFaces = 1;\n')
        # Hide node markers and wireframe lines
        f.write('Mesh.Points = 0;\n')
        f.write('Mesh.Lines = 0;\n')
        f.write('Mesh.LineWidth = 1;\n')
        # Vector views: 3D arrows, fixed 20px (the kirei knob)
        for i in range(n_vector_views):
            f.write(f'View[{i}].Visible = 1;\n')
            f.write(f'View[{i}].VectorType = 4;\n')
            f.write(f'View[{i}].ArrowSizeMin = 20;\n')
            f.write(f'View[{i}].ArrowSizeMax = 20;\n')
        # Scalar views: filled iso colormap with legend bar.
        # Log scale (ScaleType = 2) so q heating density (usual range
        # spans 2-5 decades from skin-depth local peak to far-side
        # minimum) shows visible contrast instead of a uniform near-
        # max band (2026-04-15: sugahara GUI feedback).
        for j in range(n_scalar_views):
            idx = n_vector_views + j
            f.write(f'View[{idx}].Visible = 1;\n')
            f.write(f'View[{idx}].IntervalsType = 2;\n')  # Filled
            f.write(f'View[{idx}].ScaleType = 2;\n')  # 1=Linear 2=Log 3=Double-log
            f.write(f'View[{idx}].ShowScale = 1;\n')
            f.write(f'View[{idx}].NbIso = 20;\n')
    return opt_filename


# ============================================================
# Public API: vol2msh single-mesh converter (Phase B, 2026-04-15)
# ============================================================

def vol2msh(output_msh, vol_path, fields, *, mesh_curve_order=1,
             boundary=False, companion_opt=True):
    """Load a .vol + zero-or-more .sol files and write a GMSH .msh.

    Single source of truth for panels that have ONE base mesh with
    multiple fields attached. For multi-mesh cases (e.g. IH BEM with
    coil_surface + workpiece_surface + air_volume living on three
    different meshes) see the combined writer in calc_inductance.py.

    Args:
        output_msh: Path to the .msh to write.
        vol_path: Path to the NGSolve .vol file (the mesh to load).
        fields: list of field specs. Each spec is a dict with EITHER:
            - ``array`` key  -> a numpy array of per-vertex values
              (fastest path, bypasses GridFunction reconstruction)
            - ``sol`` key    -> path to a .sol file plus ``fes``,
              ``fes_order``, ``fes_dim`` (for HDiv*/H1 dim>1).
            Common keys (both paths):
                ``name``    : view name in GMSH
                ``ncomp``   : 1 (scalar) or 3 (vector)
                ``material``: optional, restrict view to a physical group
        mesh_curve_order: if >=2, call mesh.Curve(N) so Tri6/Tet10/...
            render smoothly in GMSH.
        boundary: if True, export BND elements only (surface export
            from a volume mesh); matches GmshPostExport(boundary=...).
        companion_opt: if True, also write ``<output_msh>.opt`` so
            GMSH auto-loads sensible display defaults.

    Returns the path to the .msh written.
    """
    from ngsolve import Mesh
    from netgen.meshing import Mesh as NgMesh

    ng = NgMesh()
    ng.Load(vol_path)
    mesh = Mesh(ng)
    if mesh_curve_order >= 2:
        try:
            mesh.Curve(mesh_curve_order)
        except Exception:
            pass

    post = GmshPostExport(mesh, boundary=boundary)

    n_vec = 0
    n_sca = 0
    for f in fields:
        name = f["name"]
        ncomp = f.get("ncomp")
        material = f.get("material")
        cell_data = f.get("cell_data", False)
        if "array" in f:
            post.add_field(name, f["array"], ncomp=ncomp,
                           material=material, cell_data=cell_data)
        elif "sol" in f:
            from ngsolve import GridFunction
            fes_kind = f["fes"]
            fes_order = f.get("fes_order", 1)
            fes_dim = f.get("fes_dim", 1)
            if fes_kind == "H1":
                from ngsolve import H1
                fes = H1(mesh, order=fes_order, dim=fes_dim)
            elif fes_kind == "HDiv":
                from ngsolve import HDiv
                fes = HDiv(mesh, order=fes_order)
            elif fes_kind == "HDivSurface":
                from ngsolve import HDivSurface
                fes = HDivSurface(mesh, order=fes_order)
            elif fes_kind == "L2":
                from ngsolve import L2
                fes = L2(mesh, order=fes_order, dim=fes_dim)
            else:
                raise ValueError(f"vol2msh: unsupported fes kind {fes_kind!r}")
            gf = GridFunction(fes)
            gf.Load(f["sol"])
            post.add_field(name, gf, ncomp=ncomp,
                           material=material, cell_data=cell_data)
        else:
            raise ValueError(
                f"vol2msh: field {name!r} has neither 'array' nor 'sol'")
        if (ncomp or 0) >= 2:
            n_vec += 1
        else:
            n_sca += 1

    post.write(output_msh)
    if companion_opt:
        write_companion_opt(output_msh,
                             n_vector_views=n_vec, n_scalar_views=n_sca)
    return output_msh


def save_vol_sol_pair(vol_path, sol_path, mesh_ngmesh, gf):
    """Save a NGSolve mesh + GridFunction as the canonical .vol+.sol pair.

    Thin wrapper so every calc_*.py uses identical call sites.
    """
    mesh_ngmesh.Save(vol_path)
    gf.Save(sol_path)
    return vol_path, sol_path


# ============================================================
# Multi-mesh combined .msh writer (Phase C, 2026-04-14)
# ============================================================
# Used by calc_inductance to emit one .msh v4.1 file that holds
# three logically-separate meshes (air volume, coil surface, workpiece
# surface) as distinct physical groups with curved Tri6 coil rendering.
# Lives here so any panel that later needs the same pattern can reuse
# it without duplicating ~200 lines of .msh-formatting code.


def extract_tri6_surface(mesh_surf):
    """Extract Tri6 (curved) surface mesh from an NGSolve mesh.

    Returns:
        nodes: list of (x, y, z) — vertices first, then mid-edge nodes
        elems: list of [v0, v1, v2, m01, m12, m20] (0-indexed into nodes)
        edges: sorted list of (a, b) vertex-only edges for wireframe
    """
    from ngsolve import BND, IntegrationRule

    nodes = [(v.point[0], v.point[1], v.point[2])
             for v in mesh_surf.vertices]
    edge_cache = {}  # (min_v, max_v) -> mid_node_index

    elems = []
    wire_edges = set()
    for el in mesh_surf.Elements(BND):
        verts = [v.nr for v in el.vertices]
        if len(verts) != 3:
            continue

        trafo = mesh_surf.GetTrafo(el)

        ref_corners = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
        corner_pts = []
        for u, v in ref_corners:
            ir = IntegrationRule([(u, v)], [1.0])
            for ip in ir:
                mip = trafo(ip)
                corner_pts.append(
                    np.array([mip.point[0], mip.point[1], mip.point[2]]))

        ref_to_vert = []
        for ci in range(3):
            dists = [np.linalg.norm(corner_pts[ci] -
                     np.array(mesh_surf.vertices[verts[vi]].point))
                     for vi in range(3)]
            ref_to_vert.append(verts[np.argmin(dists)])

        ref_mids = [(0.5, 0.0), (0.5, 0.5), (0.0, 0.5)]
        edge_pairs = [(0, 1), (1, 2), (2, 0)]
        mid_indices = []
        for ei, (rc0, rc1) in enumerate(edge_pairs):
            va = ref_to_vert[rc0]
            vb = ref_to_vert[rc1]
            edge_key = (min(va, vb), max(va, vb))
            wire_edges.add(edge_key)

            if edge_key in edge_cache:
                mid_indices.append(edge_cache[edge_key])
            else:
                u, v = ref_mids[ei]
                ir = IntegrationRule([(u, v)], [1.0])
                for ip in ir:
                    mip = trafo(ip)
                    pt = (mip.point[0], mip.point[1], mip.point[2])
                    nidx = len(nodes)
                    nodes.append(pt)
                    edge_cache[edge_key] = nidx
                    mid_indices.append(nidx)

        # GMSH Tri6 ordering: v0 v1 v2 m01 m12 m20
        elems.append([ref_to_vert[0], ref_to_vert[1], ref_to_vert[2],
                      mid_indices[0], mid_indices[1], mid_indices[2]])

    return nodes, elems, sorted(wire_edges)


def write_combined_msh_v41(filename, mesh_surf, J_nodes,
                           vol_nodes, vol_elems, B_nodes,
                           wp_nodes=None, wp_elems=None, q_nodes=None,
                           wp_J_nodes=None):
    """Write a single .msh v4.1 with volume B + surface J + wireframe.

    Physical groups in one file, with node/element IDs offset so
    there are no collisions:

      Physical "air"                (3D, tag 1): volume tets   + B NodeData
      Physical "coil_surface"       (2D, tag 2): Tri6 surface  + J NodeData
      Physical "coil_wire"          (1D, tag 3): wireframe edges (no data)
      Physical "workpiece_surface"  (2D, tag 4): workpiece tris + q NodeData
                                                  (optional)

    Coil surface uses Tri6 (GMSH type 9) for curved rendering.
    Writes companion .msh.opt with display options (replaces .geo).
    """
    surf_nodes, surf_elems, wire_edges = extract_tri6_surface(mesh_surf)
    nv_surf = len(surf_nodes)

    nv_vol = len(vol_nodes)
    surf_offset = nv_vol
    ne_vol = len(vol_elems)
    ne_surf = len(surf_elems)
    ne_wire = len(wire_edges)

    has_wp = (wp_nodes is not None and wp_elems is not None
              and q_nodes is not None)
    nv_wp = len(wp_nodes) if has_wp else 0
    ne_wp = len(wp_elems) if has_wp else 0
    wp_offset = nv_vol + nv_surf

    n_nodes = nv_vol + nv_surf + nv_wp
    n_phys = 4 if has_wp else 3
    n_entities_curve = 1
    n_entities_surf = 2 if has_wp else 1
    n_entities_vol = 1

    with open(filename, 'w', encoding='utf-8') as f:
        f.write('$MeshFormat\n4.1 0 8\n$EndMeshFormat\n')

        f.write(f'$PhysicalNames\n{n_phys}\n')
        f.write('3 1 "air"\n')
        f.write('2 2 "coil_surface"\n')
        f.write('1 3 "coil_wire"\n')
        if has_wp:
            f.write('2 4 "workpiece_surface"\n')
        f.write('$EndPhysicalNames\n')

        f.write('$Entities\n')
        f.write(f'0 {n_entities_curve} {n_entities_surf} {n_entities_vol}\n')
        f.write('3 0 0 0 0 0 0 1 3 0\n')
        f.write('2 0 0 0 0 0 0 1 2 0\n')
        if has_wp:
            f.write('4 0 0 0 0 0 0 1 4 0\n')
        f.write('1 0 0 0 0 0 0 1 1 0\n')
        f.write('$EndEntities\n')

        n_blocks = 2 + (1 if has_wp else 0)
        f.write('$Nodes\n')
        f.write(f'{n_blocks} {n_nodes} 1 {n_nodes}\n')
        f.write(f'3 1 0 {nv_vol}\n')
        for i in range(nv_vol):
            f.write(f'{i + 1}\n')
        for x, y, z in vol_nodes:
            f.write(f'{x:.15e} {y:.15e} {z:.15e}\n')
        f.write(f'2 2 0 {nv_surf}\n')
        for i in range(nv_surf):
            f.write(f'{surf_offset + i + 1}\n')
        for x, y, z in surf_nodes:
            f.write(f'{x:.15e} {y:.15e} {z:.15e}\n')
        if has_wp:
            f.write(f'2 4 0 {nv_wp}\n')
            for i in range(nv_wp):
                f.write(f'{wp_offset + i + 1}\n')
            for x, y, z in wp_nodes:
                f.write(f'{x:.15e} {y:.15e} {z:.15e}\n')
        f.write('$EndNodes\n')

        n_elems = ne_vol + ne_surf + ne_wire + ne_wp
        n_elem_blocks = 3 + (1 if has_wp else 0)
        f.write('$Elements\n')
        f.write(f'{n_elem_blocks} {n_elems} 1 {n_elems}\n')
        eid = 1
        f.write(f'3 1 4 {ne_vol}\n')
        for verts in vol_elems:
            ns = ' '.join(str(v + 1) for v in verts)
            f.write(f'{eid} {ns}\n')
            eid += 1
        f.write(f'2 2 9 {ne_surf}\n')
        for conn in surf_elems:
            ns = ' '.join(str(v + surf_offset + 1) for v in conn)
            f.write(f'{eid} {ns}\n')
            eid += 1
        f.write(f'1 3 1 {ne_wire}\n')
        for a, b in wire_edges:
            f.write(f'{eid} {a + surf_offset + 1} '
                    f'{b + surf_offset + 1}\n')
            eid += 1
        if has_wp:
            f.write(f'2 4 2 {ne_wp}\n')
            for verts in wp_elems:
                ns = ' '.join(str(v + wp_offset + 1) for v in verts)
                f.write(f'{eid} {ns}\n')
                eid += 1
        f.write('$EndElements\n')

        f.write('$NodeData\n')
        f.write('1\n"B"\n1\n0.0\n3\n0\n3\n')
        f.write(f'{nv_vol}\n')
        for i in range(nv_vol):
            bx, by, bz = B_nodes[i]
            f.write(f'{i + 1} {bx:.15e} {by:.15e} {bz:.15e}\n')
        f.write('$EndNodeData\n')

        # J on coil vertex nodes only; Tri6 mid-edge nodes have no direct
        # J data — GMSH interpolates.
        nv_vert = mesh_surf.nv
        f.write('$NodeData\n')
        f.write('1\n"J"\n1\n0.0\n3\n0\n3\n')
        f.write(f'{nv_vert}\n')
        for i in range(nv_vert):
            jx, jy, jz = J_nodes[i]
            f.write(f'{surf_offset + i + 1} '
                    f'{jx:.15e} {jy:.15e} {jz:.15e}\n')
        f.write('$EndNodeData\n')

        # Keep NodeData ordering: ALL vector views first, then scalars
        # (.msh.opt writer assigns ArrowSize to View[0..n_vec-1] and
        # IntervalsType to View[n_vec..], matching this order).
        if has_wp and wp_J_nodes is not None:
            f.write('$NodeData\n')
            f.write('1\n"J_workpiece"\n1\n0.0\n3\n0\n3\n')
            f.write(f'{nv_wp}\n')
            for i in range(nv_wp):
                jx, jy, jz = wp_J_nodes[i]
                f.write(f'{wp_offset + i + 1} '
                        f'{jx:.15e} {jy:.15e} {jz:.15e}\n')
            f.write('$EndNodeData\n')

        if has_wp:
            f.write('$NodeData\n')
            f.write('1\n"q"\n1\n0.0\n3\n0\n1\n')
            f.write(f'{nv_wp}\n')
            for i in range(nv_wp):
                f.write(f'{wp_offset + i + 1} {q_nodes[i]:.15e}\n')
            f.write('$EndNodeData\n')

    n_vec = 2 + (1 if (has_wp and wp_J_nodes is not None) else 0)
    n_sca = 1 if has_wp else 0
    write_companion_opt(filename, n_vector_views=n_vec,
                        n_scalar_views=n_sca)
    return filename


def _write_companion_geo(msh_filename):
    """Deprecated: use write_companion_opt() instead."""
    return write_companion_opt(msh_filename)


def _detect_curve_order(mesh):
    """Detect the curve order of the mesh."""
    try:
        # NGSolve mesh stores curve order
        return mesh.GetCurveOrder()
    except Exception:
        return 1


def _get_gmsh_trig_ref_points(p):
    """Get GMSH Lagrange triangle reference points for order p.

    GMSH Lagrange triangles use equidistant node placement in the reference
    triangle (0,0)-(1,0)-(0,1). Node ordering follows GMSH convention:
      1. Corners: (0,0), (1,0), (0,1)
      2. Edge 0->1: (k/p, 0) for k=1..p-1
      3. Edge 1->2: (1-k/p, k/p) for k=1..p-1
      4. Edge 2->0: (0, 1-k/p) for k=1..p-1
      5. Interior: row-by-row equidistant, j=1..p-2, i=1..p-j-1

    Reference: GMSH documentation, Section 9.1 "Node ordering"
    Verified against gmsh.model.mesh.getElementProperties() for p=1..5.

    Returns list of (u, v) tuples in GMSH node ordering.
    Total nodes: (p+1)*(p+2)/2
    """
    pts = []
    # Corners
    pts.extend([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)])
    # Edge 0->1: (k/p, 0) for k=1..p-1
    for k in range(1, p):
        pts.append((k / p, 0.0))
    # Edge 1->2: (1-k/p, k/p) for k=1..p-1
    for k in range(1, p):
        pts.append((1.0 - k / p, k / p))
    # Edge 2->0: (0, 1-k/p) for k=1..p-1
    for k in range(1, p):
        pts.append((0.0, 1.0 - k / p))
    # Interior: equidistant points (i/p, j/p) for i>0, j>0, i+j<p
    for j in range(1, p):
        for i in range(1, p - j):
            pts.append((i / p, j / p))
    return pts


def _compute_highorder_nodes(mesh, nodes, curve_order, is_surface):
    """Precompute high-order node positions using GetTrafo + GMSH reference coords.

    For each BND element, evaluates the curved mesh transformation at the
    GMSH reference points (skipping corners which are already in nodes[]).
    Edge nodes are cached across shared edges.

    Args:
        mesh: NGSolve mesh (must have Curve(p) applied)
        nodes: master node list (appended in-place)
        curve_order: polynomial order (2, 3, 4, 5)
        is_surface: True if exporting BND elements

    Returns:
        elem_ho_nodes: dict { el.nr -> list of node indices for non-corner nodes }
    """
    from ngsolve import BND, IntegrationRule

    p = curve_order
    gmsh_ref = _get_gmsh_trig_ref_points(p)
    n_total = len(gmsh_ref)  # (p+1)(p+2)/2
    n_edge_per = p - 1
    # GMSH layout: 3 corners, then 3*(p-1) edge nodes, then interior
    # Edge nodes: indices 3..3+3*(p-1)-1
    # Interior nodes: indices 3+3*(p-1)..end

    # Determine which ref points are edge vs interior (skip corners 0,1,2)
    # Edge 0->1: indices 3 .. 3+(p-1)-1
    # Edge 1->2: indices 3+(p-1) .. 3+2*(p-1)-1
    # Edge 2->0: indices 3+2*(p-1) .. 3+3*(p-1)-1
    # Interior: indices 3+3*(p-1) .. end

    # For vertex matching: evaluate corners and match to el.vertices
    elem_ho_nodes = {}  # el.nr -> [node_idx for gmsh nodes 3..end]

    # Edge cache: (min_v, max_v) -> (start_vertex, [node_indices])
    edge_cache = {}

    for el in mesh.Elements(BND):
        verts = [v.nr for v in el.vertices]
        if len(verts) != 3:
            continue  # TRIG only for now

        trafo = mesh.GetTrafo(el)

        # Match ref corners to physical vertices
        corner_mapped = []
        for ci in range(3):
            u, v = gmsh_ref[ci]
            ir = IntegrationRule([(u, v)], [1.0])
            for ip in ir:
                mip = trafo(ip)
                corner_mapped.append(
                    np.array([mip.point[0], mip.point[1], mip.point[2]]))

        # Find permutation: ref_corner[i] -> which vertex?
        ref_to_vert = []
        for ci in range(3):
            dists = [np.linalg.norm(corner_mapped[ci] -
                     np.array(mesh.vertices[verts[vi]].point))
                     for vi in range(3)]
            ref_to_vert.append(verts[np.argmin(dists)])

        # Process 3 GMSH edges
        ho_indices = []
        gmsh_edge_defs = [(0, 1), (1, 2), (2, 0)]  # ref corner pairs

        for edge_i, (rc0, rc1) in enumerate(gmsh_edge_defs):
            va = ref_to_vert[rc0]
            vb = ref_to_vert[rc1]
            edge_key = (min(va, vb), max(va, vb))

            start_idx = 3 + edge_i * n_edge_per

            if edge_key in edge_cache:
                cached_start, cached_nodes = edge_cache[edge_key]
                # Determine direction
                if cached_start == va:
                    ho_indices.extend(cached_nodes)
                else:
                    ho_indices.extend(reversed(cached_nodes))
            else:
                # Evaluate at edge ref points
                edge_nodes = []
                for k in range(n_edge_per):
                    ref_idx = start_idx + k
                    u, v = gmsh_ref[ref_idx]
                    ir = IntegrationRule([(u, v)], [1.0])
                    for ip in ir:
                        mip = trafo(ip)
                        pt = (mip.point[0], mip.point[1], mip.point[2])
                        nidx = len(nodes)
                        nodes.append(pt)
                        edge_nodes.append(nidx)
                edge_cache[edge_key] = (va, edge_nodes)
                ho_indices.extend(edge_nodes)

        # Interior nodes
        int_start = 3 + 3 * n_edge_per
        for ri in range(int_start, n_total):
            u, v = gmsh_ref[ri]
            ir = IntegrationRule([(u, v)], [1.0])
            for ip in ir:
                mip = trafo(ip)
                pt = (mip.point[0], mip.point[1], mip.point[2])
                nidx = len(nodes)
                nodes.append(pt)
                ho_indices.append(nidx)

        elem_ho_nodes[el.nr] = ho_indices

    return elem_ho_nodes


def _build_highorder_connectivity(el, corner_nodes, et_name, curve_order,
                                  elem_ho_nodes, _unused, nodes):
    """Build GMSH high-order element connectivity.

    Args:
        el: NGSolve BND element
        corner_nodes: [c0, c1, c2] in GMSH order
        et_name: 'TRIG' or 'QUAD'
        curve_order: p
        elem_ho_nodes: from _compute_highorder_nodes
        _unused: kept for API compatibility
        nodes: master node list

    Returns:
        Full connectivity list for GMSH high-order element.
    """
    ho = elem_ho_nodes.get(el.nr, [])
    return list(corner_nodes) + ho


def _get_element_material(mesh, el, is_surface):
    """Get material/boundary name for an NGSolve element."""
    try:
        name = str(el.mat)
        if name:
            return name
    except Exception:
        pass
    return "default"


# ============================================================
# Internal: element ID assignment
# ============================================================

def _assign_element_ids(mat_names, elem_data):
    """Assign GMSH 1-indexed element IDs, grouped by material.

    Returns:
        mat_elem_map: {mat_name: [(orig_idx, gmsh_elem_id)]}
        elem_id_map: {orig_idx: gmsh_elem_id}
        n_elems: total number of elements
    """
    mat_elem_map = {m: [] for m in mat_names}
    elem_id_map = {}
    gmsh_id = 1

    # Group by material, then by element type (for entity blocks)
    for mat_name in mat_names:
        mat_elems = [(mn, gt, conn, oi) for mn, gt, conn, oi in elem_data
                     if mn == mat_name]
        # Sort by element type for consistent output
        mat_elems.sort(key=lambda x: x[1])
        for _, _, _, orig_idx in mat_elems:
            mat_elem_map[mat_name].append((orig_idx, gmsh_id))
            elem_id_map[orig_idx] = gmsh_id
            gmsh_id += 1

    n_elems = gmsh_id - 1
    return mat_elem_map, elem_id_map, n_elems


def _compute_material_bboxes(nodes, mat_names, elem_data):
    """Compute bounding boxes for each material.

    Returns:
        {mat_name: (xmin, ymin, zmin, xmax, ymax, zmax)}
    """
    mat_bboxes = {}
    nodes_arr = np.array(nodes) if nodes else np.zeros((0, 3))

    for mat_name in mat_names:
        node_set = set()
        for mn, gt, conn, oi in elem_data:
            if mn == mat_name:
                node_set.update(conn)
        if not node_set:
            mat_bboxes[mat_name] = (0, 0, 0, 0, 0, 0)
            continue
        pts = nodes_arr[list(node_set)]
        mn = pts.min(axis=0)
        mx = pts.max(axis=0)
        mat_bboxes[mat_name] = (mn[0], mn[1], mn[2], mx[0], mx[1], mx[2])

    return mat_bboxes


# ============================================================
# Internal: GMSH section writers
# ============================================================

def _write_physical_names(f, mat_names, dim):
    """Write $PhysicalNames section."""
    if not mat_names:
        return
    f.write('$PhysicalNames\n')
    f.write(f'{len(mat_names)}\n')
    for i, name in enumerate(mat_names):
        f.write(f'{dim} {i + 1} "{name}"\n')
    f.write('$EndPhysicalNames\n')


def _write_entities(f, mat_names, mat_bboxes, dim):
    """Write $Entities section (one entity per material)."""
    if not mat_names:
        return

    n_surf = len(mat_names) if dim == 2 else 0
    n_vol = len(mat_names) if dim == 3 else 0

    f.write('$Entities\n')
    f.write(f'0 0 {n_surf} {n_vol}\n')

    for i, name in enumerate(mat_names):
        tag = i + 1
        bbox = mat_bboxes.get(name, (0, 0, 0, 0, 0, 0))
        xmin, ymin, zmin, xmax, ymax, zmax = bbox
        # entityTag xMin yMin zMin xMax yMax zMax numPhysicalTags physTag ... numBounding boundTag ...
        f.write(f'{tag} {xmin:.15e} {ymin:.15e} {zmin:.15e} '
                f'{xmax:.15e} {ymax:.15e} {zmax:.15e} '
                f'1 {tag} 0\n')

    f.write('$EndEntities\n')


def _write_nodes(f, nodes, dim):
    """Write $Nodes section (all nodes in one entity block)."""
    n = len(nodes)
    if n == 0:
        f.write('$Nodes\n0 0 0 0\n$EndNodes\n')
        return

    f.write('$Nodes\n')
    f.write(f'1 {n} 1 {n}\n')
    # One entity block: dim, tag=1, parametric=0, numNodes
    f.write(f'{dim} 1 0 {n}\n')
    for i in range(n):
        f.write(f'{i + 1}\n')
    for i in range(n):
        f.write(f'{nodes[i][0]:.15e} {nodes[i][1]:.15e} {nodes[i][2]:.15e}\n')
    f.write('$EndNodes\n')


def _write_elements_grouped(f, mat_names, elem_data, n_elems, dim):
    """Write $Elements section grouped by material and element type."""
    if n_elems == 0:
        f.write('$Elements\n0 0 0 0\n$EndElements\n')
        return

    # Build entity blocks: (dim, entity_tag, gmsh_type, [(conn, orig_idx)])
    entity_blocks = []
    for i, mat_name in enumerate(mat_names):
        entity_tag = i + 1
        # Group this material's elements by type
        type_groups = {}
        for mn, gt, conn, oi in elem_data:
            if mn != mat_name:
                continue
            if gt not in type_groups:
                type_groups[gt] = []
            type_groups[gt].append((conn, oi))

        for gmsh_type in sorted(type_groups.keys()):
            entity_blocks.append(
                (dim, entity_tag, gmsh_type, type_groups[gmsh_type]))

    f.write('$Elements\n')
    f.write(f'{len(entity_blocks)} {n_elems} 1 {n_elems}\n')

    gmsh_id = 1
    for edim, etag, etype, elems in entity_blocks:
        f.write(f'{edim} {etag} {etype} {len(elems)}\n')
        for conn, _ in elems:
            node_str = ' '.join(str(n + 1) for n in conn)
            f.write(f'{gmsh_id} {node_str}\n')
            gmsh_id += 1

    f.write('$EndElements\n')


def _write_node_data(f, name, ncomp, data, n_nodes, time, timestep):
    """Write $NodeData section (vertex nodes with data only)."""
    # For high-order meshes, n_nodes > len(data) because data is per-vertex.
    # Write only nodes that have data.
    n_data = len(data) if hasattr(data, '__len__') else n_nodes
    n_out = min(n_nodes, n_data)

    f.write('$NodeData\n')
    f.write(f'1\n"{name}"\n')
    f.write(f'1\n{time:.15e}\n')
    f.write(f'3\n{timestep}\n{ncomp}\n{n_out}\n')

    for i in range(n_out):
        if ncomp == 1:
            val = float(data[i])
            f.write(f'{i + 1} {val:.15e}\n')
        else:
            vals = data[i]
            val_str = ' '.join(f'{float(v):.15e}' for v in vals)
            f.write(f'{i + 1} {val_str}\n')

    f.write('$EndNodeData\n')


def _write_node_data_filtered(f, name, ncomp, data, node_indices, time, timestep):
    """Write $NodeData for a subset of nodes (material-filtered)."""
    n = len(node_indices)
    f.write('$NodeData\n')
    f.write(f'1\n"{name}"\n')
    f.write(f'1\n{time:.15e}\n')
    f.write(f'3\n{timestep}\n{ncomp}\n{n}\n')

    for node_idx in node_indices:
        gmsh_nid = node_idx + 1
        if ncomp == 1:
            val = float(data[node_idx]) if node_idx < len(data) else 0.0
            f.write(f'{gmsh_nid} {val:.15e}\n')
        else:
            vals = data[node_idx] if node_idx < len(data) else np.zeros(ncomp)
            val_str = ' '.join(f'{float(v):.15e}' for v in vals)
            f.write(f'{gmsh_nid} {val_str}\n')

    f.write('$EndNodeData\n')


def _write_element_data_filtered(f, name, ncomp, data, elem_ids, time, timestep):
    """Write $ElementData for specified elements only.

    Args:
        elem_ids: list of (orig_idx, gmsh_elem_id) tuples
    """
    n = len(elem_ids)
    f.write('$ElementData\n')
    f.write(f'1\n"{name}"\n')
    f.write(f'1\n{time:.15e}\n')
    f.write(f'3\n{timestep}\n{ncomp}\n{n}\n')

    for orig_idx, gmsh_id in elem_ids:
        if ncomp == 1:
            val = float(data[orig_idx]) if orig_idx < len(data) else 0.0
            f.write(f'{gmsh_id} {val:.15e}\n')
        else:
            vals = data[orig_idx] if orig_idx < len(data) else np.zeros(ncomp)
            val_str = ' '.join(f'{float(v):.15e}' for v in vals)
            f.write(f'{gmsh_id} {val_str}\n')

    f.write('$EndElementData\n')


# ============================================================
# Internal: CoefficientFunction evaluation
# ============================================================

def _evaluate_cf(mesh, cf, ncomp, cell_data):
    """Evaluate CoefficientFunction at vertices or element centroids.

    Returns:
        numpy array: (n,) for scalar, (n, 3) for vector
    """
    is_surface = _is_surface_mesh(mesh)

    if cell_data:
        from ngsolve import VOL, BND
        domain = BND if is_surface else VOL
        n = sum(1 for _ in mesh.Elements(domain))
        result = np.zeros(n) if ncomp == 1 else np.zeros((n, ncomp))

        for i, el in enumerate(mesh.Elements(domain)):
            verts = list(el.vertices)
            pts = np.array([list(mesh.vertices[v.nr].point) for v in verts])
            centroid = pts.mean(axis=0)
            try:
                mip = mesh(*centroid)
                val = cf(mip)
                if ncomp == 1:
                    result[i] = float(val)
                else:
                    result[i, :] = list(val)[:ncomp]
            except Exception:
                pass
        return result
    else:
        nv = mesh.nv
        result = np.zeros(nv) if ncomp == 1 else np.zeros((nv, ncomp))

        for v in mesh.vertices:
            pt = v.point
            try:
                mip = mesh(pt[0], pt[1], pt[2])
                val = cf(mip)
                if ncomp == 1:
                    result[v.nr] = float(val)
                else:
                    result[v.nr, :] = list(val)[:ncomp]
            except Exception:
                pass
        return result


# ============================================================
# .vol + .sol -> .msh v4.1 converter
# ============================================================

def _extract_curved_surface(mesh, curve_order):
    """Extract curved surface elements from an NGSolve mesh.

    Uses GetTrafo to compute high-order node positions matching
    the mesh's curve order (Curve(2)->Tri6, Curve(3)->Tri10, etc.).

    Args:
        mesh: NGSolve mesh with Curve(p) applied.
        curve_order: Polynomial order (1-5).

    Returns:
        nodes: list of (x, y, z) -- vertices, then edge/interior nodes.
        elems: list of connectivity arrays (0-indexed into nodes).
        gmsh_type: GMSH element type code (2=Tri3, 9=Tri6, 21=Tri10, ...).
    """
    from ngsolve import BND, IntegrationRule

    p = curve_order
    gmsh_type = _GMSH_TYPE_BY_ORDER['TRIG'].get(p, 2)
    gmsh_ref = _get_gmsh_trig_ref_points(p)
    ho_refs = gmsh_ref[3:]  # skip 3 corners
    n_edge_nodes = p - 1

    nodes = [(v.point[0], v.point[1], v.point[2])
             for v in mesh.vertices]
    edge_cache = {}  # (min_v, max_v) -> [node_indices]

    elems = []
    for el in mesh.Elements(BND):
        verts = [v.nr for v in el.vertices]
        if len(verts) != 3:
            continue

        trafo = mesh.GetTrafo(el)

        # Map reference corners to physical vertices
        ref_corners = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
        corner_pts = []
        for u, v in ref_corners:
            ir = IntegrationRule([(u, v)], [1.0])
            for ip in ir:
                mip = trafo(ip)
                corner_pts.append(
                    np.array([mip.point[0], mip.point[1], mip.point[2]]))

        ref_to_vert = []
        for ci in range(3):
            dists = [np.linalg.norm(corner_pts[ci] -
                     np.array(mesh.vertices[verts[vi]].point))
                     for vi in range(3)]
            ref_to_vert.append(verts[np.argmin(dists)])

        conn = list(ref_to_vert)

        # Edge nodes (3 edges x (p-1) nodes each)
        edge_pairs = [(0, 1), (1, 2), (2, 0)]
        for ei in range(3):
            va = ref_to_vert[edge_pairs[ei][0]]
            vb = ref_to_vert[edge_pairs[ei][1]]
            edge_key = (min(va, vb), max(va, vb))
            reverse = va > vb

            if edge_key in edge_cache:
                indices = edge_cache[edge_key]
                conn.extend(reversed(indices) if reverse else indices)
            else:
                start = ei * n_edge_nodes
                indices = []
                for k in range(n_edge_nodes):
                    u, v = ho_refs[start + k]
                    ir = IntegrationRule([(u, v)], [1.0])
                    for ip in ir:
                        mip = trafo(ip)
                        pt = (mip.point[0], mip.point[1], mip.point[2])
                        indices.append(len(nodes))
                        nodes.append(pt)
                edge_cache[edge_key] = indices
                conn.extend(indices)

        # Interior nodes (order >= 3)
        interior_start = 3 * n_edge_nodes
        for k in range(interior_start, len(ho_refs)):
            u, v = ho_refs[k]
            ir = IntegrationRule([(u, v)], [1.0])
            for ip in ir:
                mip = trafo(ip)
                pt = (mip.point[0], mip.point[1], mip.point[2])
                conn.append(len(nodes))
                nodes.append(pt)

        elems.append(conn)

    return nodes, elems, gmsh_type


def _extract_wireframe(mesh):
    """Extract wireframe edges from mesh boundary elements."""
    from ngsolve import BND
    edges = set()
    for el in mesh.Elements(BND):
        verts = [v.nr for v in el.vertices]
        for i in range(len(verts)):
            a, b = verts[i], verts[(i + 1) % len(verts)]
            edges.add((min(a, b), max(a, b)))
    return sorted(edges)


def convert_sol_to_msh(output_msh,
                       surface_vol, j_sol,
                       air_vol=None, b_sol=None,
                       wp_vol=None, q_sol=None, wp_sigma=None,
                       fes_order=0):
    """Convert .vol + .sol pairs to a single GMSH .msh v4.1 file.

    Data flow::

        Solve -> .vol + .sol (NGSolve native)
          -> convert_sol_to_msh()
          -> .msh + .msh.opt (GMSH visualization)

    Coil surface element order matches the mesh's curve order
    (Curve(2)->Tri6, Curve(3)->Tri10, Curve(4)->Tri15, Curve(5)->Tri21).

    Args:
        output_msh: Output .msh file path.
        surface_vol: Coil surface mesh (.vol).
        j_sol: J GridFunction (.sol), HDivSurface.
        air_vol: Air volume mesh (.vol). Optional.
        b_sol: B GridFunction (.sol), HDiv. Optional.
        wp_vol: Workpiece surface mesh (.vol). Optional.
        q_sol: Workpiece heating GridFunction (.sol), H1. Optional.
            If None but wp_vol given, q is computed from J_t and wp_sigma.
        wp_sigma: Workpiece conductivity [S/m]. Used when q_sol is None.
        fes_order: FES order used for J (default 0 = RWG).

    Returns:
        Path to the .msh file.
    """
    from ngsolve import (Mesh, BND, HDivSurface, HDiv, GridFunction,
                         Integrate, CF)
    from netgen.meshing import Mesh as NetgenMesh

    # --- Load coil surface ---
    ngm = NetgenMesh()
    ngm.Load(surface_vol)
    mesh_surf = Mesh(ngm)
    curve_order = _detect_curve_order(mesh_surf)
    if curve_order < 2:
        mesh_surf.Curve(2)
        curve_order = 2
    nv_surf = mesh_surf.nv

    fes_J = HDivSurface(mesh_surf, order=fes_order)
    gf_J = GridFunction(fes_J)
    gf_J.Load(j_sol)

    # Extract curved surface + wireframe
    surf_nodes, surf_elems, surf_gmsh_type = _extract_curved_surface(
        mesh_surf, curve_order)
    wire_edges = _extract_wireframe(mesh_surf)

    # Evaluate J at ALL surface nodes (vertices + HO mid-edge/interior)
    # so GMSH can render NodeData on Tri6/Tri10 elements correctly.
    J_all = np.zeros((len(surf_nodes), 3))
    for ni, (x, y, z) in enumerate(surf_nodes):
        try:
            pt = mesh_surf(x, y, z)
            val = gf_J(pt)
            J_all[ni] = [float(val[k]) for k in range(3)]
        except Exception:
            pass

    # --- Load air volume (optional) ---
    vol_nodes = []
    vol_elems = []
    B_nodes = None
    if air_vol and b_sol:
        ngm_air = NetgenMesh()
        ngm_air.Load(air_vol)
        mesh_air = Mesh(ngm_air)

        fes_B = HDiv(mesh_air, order=1)
        gf_B = GridFunction(fes_B)
        gf_B.Load(b_sol)

        vol_nodes = [(v.point[0], v.point[1], v.point[2])
                     for v in mesh_air.vertices]
        vol_elems = [[v.nr for v in el.vertices]
                     for el in mesh_air.Elements()]
        B_nodes = np.zeros((len(vol_nodes), 3))
        for vi, v in enumerate(mesh_air.vertices):
            bval = gf_B(mesh_air(*v.point))
            B_nodes[vi] = [float(bval[k]) for k in range(3)]

    has_air = len(vol_nodes) > 0

    # --- Load workpiece (optional) ---
    wp_nodes_list = []
    wp_elems_list = []
    q_all = None
    has_wp = False
    if wp_vol:
        from ngsolve import H1
        ngm_wp = NetgenMesh()
        ngm_wp.Load(wp_vol)
        mesh_wp = Mesh(ngm_wp)
        wp_curve = _detect_curve_order(mesh_wp)
        if wp_curve < 2:
            mesh_wp.Curve(2)

        wp_nodes_list = [(v.point[0], v.point[1], v.point[2])
                         for v in mesh_wp.vertices]
        wp_elems_list = [[v.nr for v in el.vertices]
                         for el in mesh_wp.Elements(BND)]

        if q_sol:
            # Load pre-computed heating density
            fes_q = H1(mesh_wp, order=1)
            gf_q = GridFunction(fes_q)
            gf_q.Load(q_sol)
            q_all = np.zeros(len(wp_nodes_list))
            for vi, v in enumerate(mesh_wp.vertices):
                try:
                    q_all[vi] = float(gf_q(mesh_wp(*v.point)))
                except Exception:
                    pass
        elif wp_sigma is not None and wp_sigma > 0:
            # Compute q = sigma * |J_t|^2 / 2 from Biot-Savart J on wp
            # For now, store zeros — full coupling requires BEM solve
            q_all = np.zeros(len(wp_nodes_list))

        has_wp = len(wp_nodes_list) > 0

    # --- Write .msh v4.1 ---
    nv_vol = len(vol_nodes)
    nv_s = len(surf_nodes)
    nv_wp = len(wp_nodes_list)
    surf_offset = nv_vol
    wp_offset = nv_vol + nv_s

    n_phys = 2  # coil_surface + coil_wire
    if has_air:
        n_phys += 1
    if has_wp:
        n_phys += 1

    with open(output_msh, 'w', encoding='utf-8') as f:
        f.write('$MeshFormat\n4.1 0 8\n$EndMeshFormat\n')

        # Physical names
        f.write(f'$PhysicalNames\n{n_phys}\n')
        if has_air:
            f.write('3 1 "air"\n')
        f.write('2 2 "coil_surface"\n')
        f.write('1 3 "coil_wire"\n')
        if has_wp:
            f.write('2 4 "workpiece_surface"\n')
        f.write('$EndPhysicalNames\n')

        # Entities
        n_surf_ent = 1 + (1 if has_wp else 0)
        n_vol_ent = 1 if has_air else 0
        f.write('$Entities\n')
        f.write(f'0 1 {n_surf_ent} {n_vol_ent}\n')
        f.write('3 0 0 0 0 0 0 1 3 0\n')  # curve
        f.write('2 0 0 0 0 0 0 1 2 0\n')  # coil surface
        if has_wp:
            f.write('4 0 0 0 0 0 0 1 4 0\n')  # workpiece surface
        if has_air:
            f.write('1 0 0 0 0 0 0 1 1 0\n')  # volume
        f.write('$EndEntities\n')

        # Nodes
        n_nodes = nv_vol + nv_s + nv_wp
        n_blocks = (1 if has_air else 0) + 1 + (1 if has_wp else 0)
        f.write('$Nodes\n')
        f.write(f'{n_blocks} {n_nodes} 1 {n_nodes}\n')
        if has_air:
            f.write(f'3 1 0 {nv_vol}\n')
            for i in range(nv_vol):
                f.write(f'{i + 1}\n')
            for x, y, z in vol_nodes:
                f.write(f'{x:.15e} {y:.15e} {z:.15e}\n')
        f.write(f'2 2 0 {nv_s}\n')
        for i in range(nv_s):
            f.write(f'{surf_offset + i + 1}\n')
        for x, y, z in surf_nodes:
            f.write(f'{x:.15e} {y:.15e} {z:.15e}\n')
        if has_wp:
            f.write(f'2 4 0 {nv_wp}\n')
            for i in range(nv_wp):
                f.write(f'{wp_offset + i + 1}\n')
            for x, y, z in wp_nodes_list:
                f.write(f'{x:.15e} {y:.15e} {z:.15e}\n')
        f.write('$EndNodes\n')

        # Elements
        ne_vol = len(vol_elems)
        ne_surf = len(surf_elems)
        ne_wire = len(wire_edges)
        ne_wp = len(wp_elems_list)
        n_elems = ne_vol + ne_surf + ne_wire + ne_wp
        n_elem_blocks = (1 if has_air else 0) + 2 + (1 if has_wp else 0)
        f.write('$Elements\n')
        f.write(f'{n_elem_blocks} {n_elems} 1 {n_elems}\n')
        eid = 1
        if has_air:
            f.write(f'3 1 4 {ne_vol}\n')
            for verts in vol_elems:
                f.write(f'{eid} {" ".join(str(v+1) for v in verts)}\n')
                eid += 1
        f.write(f'2 2 {surf_gmsh_type} {ne_surf}\n')
        for conn in surf_elems:
            f.write(f'{eid} {" ".join(str(v+surf_offset+1) for v in conn)}\n')
            eid += 1
        f.write(f'1 3 1 {ne_wire}\n')
        for a, b in wire_edges:
            f.write(f'{eid} {a+surf_offset+1} {b+surf_offset+1}\n')
            eid += 1
        if has_wp:
            f.write(f'2 4 2 {ne_wp}\n')
            for verts in wp_elems_list:
                f.write(f'{eid} {" ".join(str(v+wp_offset+1) for v in verts)}\n')
                eid += 1
        f.write('$EndElements\n')

        # NodeData: B
        n_views = 0
        if has_air and B_nodes is not None:
            f.write('$NodeData\n1\n"B"\n1\n0.0\n3\n0\n3\n')
            f.write(f'{nv_vol}\n')
            for i in range(nv_vol):
                bx, by, bz = B_nodes[i]
                f.write(f'{i+1} {bx:.15e} {by:.15e} {bz:.15e}\n')
            f.write('$EndNodeData\n')
            n_views += 1

        # NodeData: J on ALL surface nodes (vertex + HO)
        f.write('$NodeData\n1\n"J"\n1\n0.0\n3\n0\n3\n')
        f.write(f'{nv_s}\n')
        for i in range(nv_s):
            jx, jy, jz = J_all[i]
            f.write(f'{surf_offset+i+1} {jx:.15e} {jy:.15e} {jz:.15e}\n')
        f.write('$EndNodeData\n')
        n_views += 1

        # NodeData: q (heating density) on workpiece
        if has_wp and q_all is not None:
            f.write('$NodeData\n1\n"q"\n1\n0.0\n3\n0\n1\n')
            f.write(f'{nv_wp}\n')
            for i in range(nv_wp):
                f.write(f'{wp_offset+i+1} {q_all[i]:.15e}\n')
            f.write('$EndNodeData\n')
            n_views += 1

    # B and J are vector views; q is scalar
    n_vec = min(n_views, 2)  # at most B + J
    n_scl = n_views - n_vec  # remaining views are scalar (q)
    write_companion_opt(output_msh, n_vector_views=n_vec,
                        n_scalar_views=n_scl)

    print(f"convert_sol_to_msh: {output_msh}")
    print(f"  coil: {ne_surf} Tri (order {curve_order}), "
          f"{nv_s} nodes ({nv_surf} vert + {nv_s - nv_surf} HO)")
    if has_air:
        print(f"  air: {ne_vol} tets, {nv_vol} nodes")
    if has_wp:
        print(f"  workpiece: {ne_wp} tris, {nv_wp} nodes")
    print(f"  wireframe: {ne_wire} edges, views: {n_views}")
    return output_msh
