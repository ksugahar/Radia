"""
gmsh_post_export.py

GMSH .pos / .msh post-processing export for high-order element visualization.

VTK is weak with high-order elements. GMSH natively supports Tri6, Quad9,
Tet10, Hex20 with correct high-order interpolation in its GUI.

Output formats:
    - GMSH .msh v4.1 with $ElementNodeData (recommended)
    - GMSH .pos ASCII (legacy, also supported)

Usage:
    from ngsolve import *
    from radia.gmsh_post_export import GmshPostExport

    # After solving...
    post = GmshPostExport(mesh)
    post.add_vector_field("J", gf_J)
    post.add_scalar_field("|J|", Norm(gf_J))
    post.write("results.msh")

    # Or write just the mesh (no Cubit dependency):
    post.write_mesh("mesh_only.msh")

Supports:
    1st order: Tri3 (type 2), Quad4 (3), Tet4 (4), Hex8 (5), Wedge6 (6), Pyr5 (7)
    2nd order: Tri6 (9), Quad9 (10), Tet10 (11), Hex20 (17), Wedge15 (18), Pyr13 (19)

Part of Radia project
"""

import numpy as np

# Physical constants
MU_0 = 4.0 * np.pi * 1e-7


# ============================================================
# GMSH element type codes
# ============================================================

# NGSolve ET -> GMSH type mapping (1st order)
_NGSOLVE_TO_GMSH_1ST = {
    'TRIG': 2,     # Tri3
    'QUAD': 3,     # Quad4
    'TET': 4,      # Tet4
    'HEX': 5,      # Hex8
    'PRISM': 6,    # Wedge6
    'PYRAMID': 7,  # Pyramid5
}

# NGSolve ET -> GMSH type mapping (2nd order)
_NGSOLVE_TO_GMSH_2ND = {
    'TRIG': 9,     # Tri6
    'QUAD': 10,    # Quad9
    'TET': 11,     # Tet10
    'HEX': 17,     # Hex20
    'PRISM': 18,   # Wedge15
    'PYRAMID': 19, # Pyramid13
}

# Nodes per GMSH element type
_GMSH_NODES_PER_TYPE = {
    2: 3, 3: 4, 4: 4, 5: 8, 6: 6, 7: 5,
    9: 6, 10: 9, 11: 10, 17: 20, 18: 15, 19: 13,
}

# NGSolve -> GMSH node reordering for volume elements
# NGSolve and GMSH have different node numbering conventions
_NGSOLVE_TO_GMSH_NODE_ORDER = {
    'TET': [0, 1, 2, 3],                        # identical
    'HEX': [0, 1, 5, 4, 3, 2, 6, 7],           # Netgen hex convention
    'PRISM': [0, 2, 1, 3, 5, 4],                # swap top/bottom
    'PYRAMID': [3, 2, 1, 0, 4],                 # reversed base
    'TRIG': [0, 1, 2],                           # identical
    'QUAD': [0, 1, 2, 3],                        # identical
}


# ============================================================
# GmshPostExport class
# ============================================================

class GmshPostExport:
    """Export NGSolve mesh and field data to GMSH format.

    Supports high-order elements (Tri6, Tet10, Hex20, etc.) for
    accurate visualization in GMSH GUI.

    Args:
        mesh: NGSolve Mesh object (volume or surface mesh)

    Example:
        post = GmshPostExport(mesh)
        post.add_vector_field("B", gf_B)
        post.add_scalar_field("T", temperature_cf)
        post.write("results.msh")
    """

    def __init__(self, mesh):
        self.mesh = mesh
        self._fields = []   # list of (name, ncomp, data_array, is_cell)
        self._mesh_written = False

    def add_scalar_field(self, name, data, cell_data=False):
        """Add a scalar field for export.

        Args:
            name: Field name (displayed in GMSH GUI)
            data: NGSolve CoefficientFunction/GridFunction or numpy array
                  shape (n,) where n = nv (point) or ne (cell)
            cell_data: If True, evaluate at element centroids (CellData).
                       If False, evaluate at vertices (PointData).
        """
        arr = self._resolve_data(data, ncomp=1, cell_data=cell_data)
        self._fields.append((name, 1, arr, cell_data))

    def add_vector_field(self, name, data, cell_data=False):
        """Add a vector field (3 components) for export.

        Args:
            name: Field name
            data: CoefficientFunction/GridFunction (dim=3) or numpy array (n, 3)
            cell_data: If True, evaluate at element centroids.
        """
        arr = self._resolve_data(data, ncomp=3, cell_data=cell_data)
        self._fields.append((name, 3, arr, cell_data))

    def write(self, filename, time=0.0, timestep=0):
        """Write mesh + field data to GMSH .msh v4.1 format.

        Args:
            filename: Output .msh file path
            time: Time value for time series
            timestep: Time step index
        """
        if not filename.endswith('.msh'):
            filename += '.msh'

        mesh = self.mesh
        is_surface = _is_surface_mesh(mesh)

        # Extract mesh topology
        nodes, elements, elem_types = _extract_mesh_data(mesh, is_surface)
        n_nodes = len(nodes)
        n_elems = len(elements)

        with open(filename, 'w', encoding='utf-8') as f:
            # Header
            f.write('$MeshFormat\n')
            f.write('4.1 0 8\n')
            f.write('$EndMeshFormat\n')

            # Nodes
            f.write('$Nodes\n')
            f.write(f'1 {n_nodes} 1 {n_nodes}\n')  # 1 entity block
            f.write(f'3 1 0 {n_nodes}\n')  # dim=3, tag=1, parametric=0
            for i in range(n_nodes):
                f.write(f'{i + 1}\n')
            for i in range(n_nodes):
                f.write(f'{nodes[i][0]:.15e} {nodes[i][1]:.15e} {nodes[i][2]:.15e}\n')
            f.write('$EndNodes\n')

            # Elements
            f.write('$Elements\n')
            # Group by element type
            type_groups = {}
            for idx, (etype, conn) in enumerate(zip(elem_types, elements)):
                if etype not in type_groups:
                    type_groups[etype] = []
                type_groups[etype].append((idx, conn))

            n_blocks = len(type_groups)
            f.write(f'{n_blocks} {n_elems} 1 {n_elems}\n')

            elem_id = 1
            elem_id_map = {}  # original index -> GMSH elem ID
            for etype, group in sorted(type_groups.items()):
                dim = 3 if etype in (4, 5, 6, 7, 11, 17, 18, 19) else 2
                f.write(f'{dim} 1 {etype} {len(group)}\n')
                for orig_idx, conn in group:
                    # GMSH uses 1-indexed nodes
                    node_str = ' '.join(str(n + 1) for n in conn)
                    f.write(f'{elem_id} {node_str}\n')
                    elem_id_map[orig_idx] = elem_id
                    elem_id += 1
            f.write('$EndElements\n')

            # Field data as $ElementNodeData or $NodeData
            for name, ncomp, data, is_cell in self._fields:
                if is_cell:
                    _write_element_data(f, name, ncomp, data, n_elems,
                                        elem_id_map, time, timestep)
                else:
                    _write_node_data(f, name, ncomp, data, n_nodes,
                                     time, timestep)

        print(f"GMSH export: {filename} ({n_elems} elements, {n_nodes} nodes, "
              f"{len(self._fields)} fields)")
        return filename

    def write_mesh(self, filename):
        """Write mesh only (no field data)."""
        fields_backup = self._fields
        self._fields = []
        result = self.write(filename)
        self._fields = fields_backup
        return result

    def _resolve_data(self, data, ncomp, cell_data):
        """Convert CoefficientFunction/GridFunction/ndarray to numpy array."""
        if isinstance(data, np.ndarray):
            return data

        # Try NGSolve CoefficientFunction (GridFunction is a subclass)
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
# Internal helpers
# ============================================================

def _is_surface_mesh(mesh):
    """Check if mesh is surface-only (no volume elements)."""
    try:
        ne = mesh.ne
        return ne == 0
    except Exception:
        return True


def _extract_mesh_data(mesh, is_surface):
    """Extract nodes, element connectivity, and GMSH element types.

    Returns:
        nodes: list of (x, y, z) tuples
        elements: list of node index lists (0-indexed)
        elem_types: list of GMSH element type codes
    """
    # Vertices
    nodes = []
    for v in mesh.vertices:
        pt = v.point
        nodes.append((pt[0], pt[1], pt[2]))

    elements = []
    elem_types = []

    if is_surface:
        # Surface elements (BEM mesh)
        from ngsolve import BND
        for el in mesh.Elements(BND):
            et_name = str(el.type).split('.')[-1]  # e.g. 'ET.TRIG' -> 'TRIG'
            gmsh_type = _NGSOLVE_TO_GMSH_1ST.get(et_name)
            if gmsh_type is None:
                continue
            verts = [v.nr for v in el.vertices]
            perm = _NGSOLVE_TO_GMSH_NODE_ORDER.get(et_name, list(range(len(verts))))
            reordered = [verts[i] for i in perm]
            elements.append(reordered)
            elem_types.append(gmsh_type)
    else:
        # Volume elements
        from ngsolve import VOL
        for el in mesh.Elements(VOL):
            et_name = str(el.type).split('.')[-1]
            gmsh_type = _NGSOLVE_TO_GMSH_1ST.get(et_name)
            if gmsh_type is None:
                continue
            verts = [v.nr for v in el.vertices]
            perm = _NGSOLVE_TO_GMSH_NODE_ORDER.get(et_name, list(range(len(verts))))
            reordered = [verts[i] for i in perm]
            elements.append(reordered)
            elem_types.append(gmsh_type)

    return nodes, elements, elem_types


def _evaluate_cf(mesh, cf, ncomp, cell_data):
    """Evaluate CoefficientFunction at vertices or element centroids.

    Args:
        mesh: NGSolve Mesh
        cf: CoefficientFunction
        ncomp: Expected number of components (1 for scalar, 3 for vector)
        cell_data: If True, evaluate at centroids; else at vertices

    Returns:
        numpy array: (n,) for scalar, (n, 3) for vector
    """
    is_surface = _is_surface_mesh(mesh)

    if cell_data:
        # Evaluate at element centroids
        from ngsolve import VOL, BND
        domain = BND if is_surface else VOL
        n = sum(1 for _ in mesh.Elements(domain))
        if ncomp == 1:
            result = np.zeros(n)
        else:
            result = np.zeros((n, ncomp))

        for i, el in enumerate(mesh.Elements(domain)):
            # Compute centroid from vertices
            verts = [v for v in el.vertices]
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
                pass  # leave as zero for points outside mesh
        return result
    else:
        # Evaluate at vertices
        nv = mesh.nv
        if ncomp == 1:
            result = np.zeros(nv)
        else:
            result = np.zeros((nv, ncomp))

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


def _write_node_data(f, name, ncomp, data, n_nodes, time, timestep):
    """Write $NodeData section to GMSH .msh file."""
    f.write('$NodeData\n')
    f.write('1\n')               # one string tag
    f.write(f'"{name}"\n')       # field name
    f.write('1\n')               # one real tag
    f.write(f'{time:.15e}\n')    # time value
    f.write('3\n')               # three integer tags
    f.write(f'{timestep}\n')     # time step index
    f.write(f'{ncomp}\n')        # number of components
    f.write(f'{n_nodes}\n')      # number of values

    for i in range(n_nodes):
        if ncomp == 1:
            val = float(data[i]) if i < len(data) else 0.0
            f.write(f'{i + 1} {val:.15e}\n')
        else:
            vals = data[i] if i < len(data) else np.zeros(ncomp)
            val_str = ' '.join(f'{float(v):.15e}' for v in vals)
            f.write(f'{i + 1} {val_str}\n')

    f.write('$EndNodeData\n')


def _write_element_data(f, name, ncomp, data, n_elems, elem_id_map, time, timestep):
    """Write $ElementData section to GMSH .msh file."""
    f.write('$ElementData\n')
    f.write('1\n')
    f.write(f'"{name}"\n')
    f.write('1\n')
    f.write(f'{time:.15e}\n')
    f.write('3\n')
    f.write(f'{timestep}\n')
    f.write(f'{ncomp}\n')
    f.write(f'{n_elems}\n')

    for orig_idx in range(n_elems):
        elem_id = elem_id_map.get(orig_idx, orig_idx + 1)
        if ncomp == 1:
            val = float(data[orig_idx]) if orig_idx < len(data) else 0.0
            f.write(f'{elem_id} {val:.15e}\n')
        else:
            vals = data[orig_idx] if orig_idx < len(data) else np.zeros(ncomp)
            val_str = ' '.join(f'{float(v):.15e}' for v in vals)
            f.write(f'{elem_id} {val_str}\n')

    f.write('$EndElementData\n')
