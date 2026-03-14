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

_GMSH_NODES_PER_TYPE = {
    2: 3, 3: 4, 4: 4, 5: 8, 6: 6, 7: 5,
    9: 6, 10: 9, 11: 10, 17: 20, 18: 15, 19: 13,
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

    Example:
        post = GmshPostExport(mesh)
        post.add_field("J", gf_J, material="coil")   # coil elements only
        post.add_field("B", gf_B, material="core")    # core elements only
        post.add_scalar_field("T", temp_cf)            # all elements
        post.write("results.msh")
    """

    def __init__(self, mesh):
        self.mesh = mesh
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
        """Write mesh + field data to GMSH .msh v4.1 with physical groups.

        Args:
            filename: Output .msh file path
            time: Time value for time series
            timestep: Time step index
        """
        if not filename.endswith('.msh'):
            filename += '.msh'

        mesh = self.mesh
        is_surface = _is_surface_mesh(mesh)

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

    Returns:
        nodes: list of (x, y, z) tuples
        mat_names: list of unique material names (preserving order)
        elem_data: list of (mat_name, gmsh_type, connectivity, orig_idx)
    """
    nodes = []
    for v in mesh.vertices:
        pt = v.point
        nodes.append((pt[0], pt[1], pt[2]))

    mat_names = []
    elem_data = []

    if is_surface:
        from ngsolve import BND
        for idx, el in enumerate(mesh.Elements(BND)):
            et_name = str(el.type).split('.')[-1]
            gmsh_type = _NGSOLVE_TO_GMSH_1ST.get(et_name)
            if gmsh_type is None:
                continue
            verts = [v.nr for v in el.vertices]
            perm = _NGSOLVE_TO_GMSH_NODE_ORDER.get(
                et_name, list(range(len(verts))))
            reordered = [verts[i] for i in perm]
            mat_name = _get_element_material(mesh, el, is_surface)
            elem_data.append((mat_name, gmsh_type, reordered, idx))
            if mat_name not in mat_names:
                mat_names.append(mat_name)
    else:
        from ngsolve import VOL
        for idx, el in enumerate(mesh.Elements(VOL)):
            et_name = str(el.type).split('.')[-1]
            gmsh_type = _NGSOLVE_TO_GMSH_1ST.get(et_name)
            if gmsh_type is None:
                continue
            verts = [v.nr for v in el.vertices]
            perm = _NGSOLVE_TO_GMSH_NODE_ORDER.get(
                et_name, list(range(len(verts))))
            reordered = [verts[i] for i in perm]
            mat_name = _get_element_material(mesh, el, is_surface)
            elem_data.append((mat_name, gmsh_type, reordered, idx))
            if mat_name not in mat_names:
                mat_names.append(mat_name)

    return nodes, mat_names, elem_data


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
    """Write $NodeData section (all nodes)."""
    f.write('$NodeData\n')
    f.write(f'1\n"{name}"\n')
    f.write(f'1\n{time:.15e}\n')
    f.write(f'3\n{timestep}\n{ncomp}\n{n_nodes}\n')

    for i in range(n_nodes):
        if ncomp == 1:
            val = float(data[i]) if i < len(data) else 0.0
            f.write(f'{i + 1} {val:.15e}\n')
        else:
            vals = data[i] if i < len(data) else np.zeros(ncomp)
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
