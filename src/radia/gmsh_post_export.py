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

# High-order Lagrange triangle types (order -> GMSH type code)
# Verified via gmsh.model.mesh.getElementType('Triangle', order)
_TRIG_GMSH_TYPE_BY_ORDER = {1: 2, 2: 9, 3: 21, 4: 23, 5: 25}

_GMSH_NODES_PER_TYPE = {
    2: 3, 3: 4, 4: 4, 5: 8, 6: 6, 7: 5,
    9: 6, 10: 9, 11: 10, 17: 20, 18: 15, 19: 13,
    21: 10, 23: 15, 25: 21,  # high-order triangles
}

_NGSOLVE_TO_GMSH_NODE_ORDER = {
    'TET': [0, 1, 2, 3],
    'HEX': [0, 1, 5, 4, 3, 2, 6, 7],
    'PRISM': [0, 2, 1, 3, 5, 4],
    'PYRAMID': [3, 2, 1, 0, 4],
    'TRIG': [0, 1, 2],
    'QUAD': [0, 1, 2, 3],
}

_VOL_TYPES = {4, 5, 6, 7, 11, 17, 18, 19}


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

    def write(self, filename, time=0.0, timestep=0, version="4.1"):
        """Write mesh + field data to GMSH .msh format.

        Args:
            filename: Output .msh file path
            time: Time value for time series
            timestep: Time step index
            version: GMSH format version ("4.1" or "2.2")
                     v4.1 (default): structured physical groups, NodeData, transient support
                     v2.2: legacy format
        """
        if version == "2.2":
            return self.write_v22(filename)

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

        # Write companion .geo file for correct high-order display
        if is_surface and _detect_curve_order(self.mesh) >= 2:
            _write_companion_geo(filename)

        n_fields = len(self._fields)
        mat_str = ', '.join(f'{m}({len(mat_elem_map.get(m, []))})' for m in mat_names)
        print(f"GMSH export: {filename}")
        print(f"  {n_elems} elements, {n_nodes} nodes, {n_fields} fields")
        print(f"  Materials: {mat_str}")
        return filename

    def write_v22(self, filename):
        """Write mesh + field data to GMSH .msh v2.2 format.

        Supports high-order elements (p >= 2) via the same internal
        node computation as write() (v4.1).

        v2.2 is the legacy format (use .vol for NGSolve mesh input).

        Args:
            filename: Output .msh file path
        """
        if not filename.endswith('.msh'):
            filename += '.msh'

        mesh = self.mesh
        is_surface = _is_surface_mesh(mesh) or self._boundary

        nodes, mat_names, elem_data = _extract_mesh_data_grouped(
            mesh, is_surface)
        n_nodes = len(nodes)
        n_elems = len(elem_data)

        # Build material -> physical tag mapping (1-indexed)
        mat_to_tag = {name: i + 1 for i, name in enumerate(mat_names)}

        # Determine dimension
        dim = 2 if is_surface else 3

        with open(filename, 'w', encoding='utf-8') as f:
            # Header
            f.write('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n')

            # Physical names
            f.write('$PhysicalNames\n')
            f.write(f'{len(mat_names)}\n')
            for name in mat_names:
                tag = mat_to_tag[name]
                f.write(f'{dim} {tag} "{name}"\n')
            f.write('$EndPhysicalNames\n')

            # Nodes
            f.write('$Nodes\n')
            f.write(f'{n_nodes}\n')
            for i, (x, y, z) in enumerate(nodes):
                f.write(f'{i + 1} {x:.15e} {y:.15e} {z:.15e}\n')
            f.write('$EndNodes\n')

            # Elements
            f.write('$Elements\n')
            f.write(f'{n_elems}\n')
            for elem_idx, (mat_name, gmsh_type, conn, _) in enumerate(elem_data):
                elem_id = elem_idx + 1
                phys_tag = mat_to_tag[mat_name]
                # v2.2 format: id type n_tags tag1 tag2 node1 node2 ...
                # tag1 = physical group, tag2 = elementary entity (= same)
                node_str = ' '.join(str(c + 1) for c in conn)
                f.write(f'{elem_id} {gmsh_type} 2 {phys_tag} {phys_tag} '
                        f'{node_str}\n')
            f.write('$EndElements\n')

            # NodeData (field data)
            for name, ncomp, data, is_cell, material in self._fields:
                if is_cell:
                    continue  # v2.2 NodeData only for simplicity
                arr = self._resolve_data(data, ncomp, False)
                if ncomp is None:
                    ncomp = 1 if arr.ndim == 1 else arr.shape[1]

                # Limit to nodes that have data (vertex nodes only for
                # high-order meshes where arr has nv entries, not n_nodes)
                n_data = arr.shape[0]
                if material is not None:
                    node_set = set()
                    for mn, _, conn, _ in elem_data:
                        if mn == material:
                            node_set.update(c for c in conn if c < n_data)
                    out_nodes = sorted(node_set)
                else:
                    out_nodes = list(range(min(n_nodes, n_data)))

                f.write('$NodeData\n')
                f.write(f'1\n"{name}"\n')
                f.write(f'1\n0.0\n')
                f.write(f'3\n0\n{ncomp}\n{len(out_nodes)}\n')
                for ni in out_nodes:
                    if ncomp == 1:
                        val = float(arr[ni]) if arr.ndim == 1 else float(arr[ni, 0])
                        f.write(f'{ni + 1} {val:.15e}\n')
                    else:
                        vals = ' '.join(f'{float(arr[ni, c]):.15e}'
                                        for c in range(ncomp))
                        f.write(f'{ni + 1} {vals}\n')
                f.write('$EndNodeData\n')

        # Write companion .geo file for correct high-order display
        if is_surface and _detect_curve_order(mesh) >= 2:
            _write_companion_geo(filename)

        print(f"GMSH v2.2 export: {filename}")
        print(f"  {n_elems} elements, {n_nodes} nodes, "
              f"{len(self._fields)} fields")
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

        # For high-order volume meshes, compute mid-edge nodes via GetTrafo
        vol_ho_cache = {}  # edge_key -> [mid-node indices]

        for idx, el in enumerate(mesh.Elements(VOL)):
            et_name = str(el.type).split('.')[-1]

            if use_highorder and et_name in _NGSOLVE_TO_GMSH_2ND:
                gmsh_type = _NGSOLVE_TO_GMSH_2ND[et_name]
            else:
                gmsh_type = _NGSOLVE_TO_GMSH_1ST.get(et_name)
            if gmsh_type is None:
                continue

            verts = [v.nr for v in el.vertices]
            perm = _NGSOLVE_TO_GMSH_NODE_ORDER.get(
                et_name, list(range(len(verts))))
            reordered = [verts[i] for i in perm]

            if use_highorder and et_name in _NGSOLVE_TO_GMSH_2ND:
                reordered = _build_vol_highorder_conn(
                    mesh, el, reordered, et_name, nodes, vol_ho_cache)

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


def _write_companion_geo(msh_filename):
    """Write a companion .geo file that merges .msh with correct display settings.

    GMSH default Mesh.NumSubEdges=1 draws curved (Tri6) elements as flat.
    The .geo file sets NumSubEdges=4 for proper curved surface rendering.
    """
    import os
    geo_filename = os.path.splitext(msh_filename)[0] + '.geo'
    msh_basename = os.path.basename(msh_filename)
    with open(geo_filename, 'w', encoding='utf-8') as f:
        f.write(f'// Auto-generated companion for {msh_basename}\n')
        f.write(f'Merge "{msh_basename}";\n')
        f.write('Mesh.NumSubEdges = 4;\n')
    return geo_filename


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
