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

from dataclasses import dataclass, field
from typing import List, Tuple

import numpy as np

# Physical constants
MU_0 = 4.0 * np.pi * 1e-7


def _point3(point):
    """Return an NGSolve 1-D/2-D/3-D point as an explicit 3-D tuple."""
    dimension = len(point)
    return (
        float(point[0]),
        float(point[1]) if dimension >= 2 else 0.0,
        float(point[2]) if dimension >= 3 else 0.0,
    )


# ============================================================
# Unified GMSH v4.1 payload model (Phase E.1, 2026-04-15)
# ============================================================
# Both the single-mesh writer (GmshPostExport.write) and the multi-mesh
# writer (write_combined_msh_v41) build PayloadV41 instances and emit
# through _emit_v41_msh. Every on-disk field of the .msh format lives
# here once, so format changes only need touching one place.


@dataclass
class PhysicalGroup:
    dim: int      # 0/1/2/3
    tag: int      # global physical tag
    name: str


@dataclass
class Entity:
    dim: int
    tag: int
    bbox: Tuple[float, float, float, float, float, float]
    physical_tags: List[int] = field(default_factory=list)
    bounding_tags: List[int] = field(default_factory=list)


@dataclass
class NodeBlock:
    entity_dim: int
    entity_tag: int
    start_id: int   # 1-indexed first node ID in this block
    coords: List[Tuple[float, float, float]]


@dataclass
class ElementBlock:
    entity_dim: int
    entity_tag: int
    gmsh_type: int
    # connectivity entries are lists of GLOBAL 1-indexed node IDs
    connectivity: List[List[int]]


@dataclass
class DataField:
    kind: str        # "$NodeData" or "$ElementData"
    name: str
    ncomp: int       # 1 (scalar) or 3 (vector)
    # items[i] = (gmsh_id, value_or_tuple). For scalar ncomp=1 value is
    # float; for ncomp>=2 value is iterable of floats.
    items: List[Tuple[int, object]]
    time: float = 0.0
    timestep: int = 0


@dataclass
class PayloadV41:
    groups: List[PhysicalGroup] = field(default_factory=list)
    entities: List[Entity] = field(default_factory=list)
    node_blocks: List[NodeBlock] = field(default_factory=list)
    element_blocks: List[ElementBlock] = field(default_factory=list)
    data_fields: List[DataField] = field(default_factory=list)


def _emit_v41_msh(f, payload: PayloadV41) -> None:
    """Single v4.1 emitter. Writes a full .msh from a PayloadV41."""
    f.write('$MeshFormat\n4.1 0 8\n$EndMeshFormat\n')

    # $PhysicalNames
    if payload.groups:
        f.write(f'$PhysicalNames\n{len(payload.groups)}\n')
        for g in payload.groups:
            f.write(f'{g.dim} {g.tag} "{g.name}"\n')
        f.write('$EndPhysicalNames\n')

    # $Entities — format per dim: point=3 coords, curve/surf/vol=6 bbox
    n_pt = sum(1 for e in payload.entities if e.dim == 0)
    n_cu = sum(1 for e in payload.entities if e.dim == 1)
    n_su = sum(1 for e in payload.entities if e.dim == 2)
    n_vl = sum(1 for e in payload.entities if e.dim == 3)
    f.write(f'$Entities\n{n_pt} {n_cu} {n_su} {n_vl}\n')
    for e in payload.entities:
        xmin, ymin, zmin, xmax, ymax, zmax = e.bbox
        parts = [str(e.tag)]
        if e.dim == 0:
            parts += [f'{xmin:.15e}', f'{ymin:.15e}', f'{zmin:.15e}']
        else:
            parts += [f'{xmin:.15e}', f'{ymin:.15e}', f'{zmin:.15e}',
                      f'{xmax:.15e}', f'{ymax:.15e}', f'{zmax:.15e}']
        parts.append(str(len(e.physical_tags)))
        parts += [str(t) for t in e.physical_tags]
        parts.append(str(len(e.bounding_tags)))
        parts += [str(t) for t in e.bounding_tags]
        f.write(' '.join(parts) + '\n')
    f.write('$EndEntities\n')

    # $Nodes
    total_nodes = sum(len(b.coords) for b in payload.node_blocks)
    if total_nodes == 0:
        f.write('$Nodes\n0 0 0 0\n$EndNodes\n')
    else:
        f.write(f'$Nodes\n{len(payload.node_blocks)} {total_nodes} '
                f'1 {total_nodes}\n')
        for b in payload.node_blocks:
            f.write(f'{b.entity_dim} {b.entity_tag} 0 {len(b.coords)}\n')
            for i in range(len(b.coords)):
                f.write(f'{b.start_id + i}\n')
            for x, y, z in b.coords:
                f.write(f'{x:.15e} {y:.15e} {z:.15e}\n')
        f.write('$EndNodes\n')

    # $Elements — IDs auto-assigned sequentially from 1
    total_elems = sum(len(b.connectivity) for b in payload.element_blocks)
    if total_elems == 0:
        f.write('$Elements\n0 0 0 0\n$EndElements\n')
    else:
        f.write(f'$Elements\n{len(payload.element_blocks)} {total_elems} '
                f'1 {total_elems}\n')
        eid = 1
        for b in payload.element_blocks:
            f.write(f'{b.entity_dim} {b.entity_tag} {b.gmsh_type} '
                    f'{len(b.connectivity)}\n')
            for conn in b.connectivity:
                f.write(f'{eid} ' + ' '.join(str(n) for n in conn) + '\n')
                eid += 1
        f.write('$EndElements\n')

    # $NodeData / $ElementData
    for df in payload.data_fields:
        _emit_data_section(f, df)


def _emit_data_section(f, df: DataField) -> None:
    """Emit a $NodeData or $ElementData block from a DataField."""
    end_kind = df.kind.replace('$', '$End')
    f.write(f'{df.kind}\n')
    f.write(f'1\n"{df.name}"\n')
    f.write(f'1\n{df.time:.15e}\n')
    f.write(f'3\n{df.timestep}\n{df.ncomp}\n')
    f.write(f'{len(df.items)}\n')
    if df.ncomp == 1:
        for tag, val in df.items:
            f.write(f'{tag} {float(val):.15e}\n')
    else:
        for tag, vec in df.items:
            vals = ' '.join(f'{float(c):.15e}' for c in vec)
            f.write(f'{tag} {vals}\n')
    f.write(f'{end_kind}\n')


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
# GMSH TET: identical reference simplex (0,0,0)-(1,0,0)-(0,1,0)-(0,0,1);
#   with _NGSOLVE_TO_GMSH_NODE_ORDER['TET'] = [3,0,1,2] the GMSH corner i
#   sits exactly at NGSolve reference vertex i, so the map is identity.
#   (NGSolve lists tet vertices so that ref(1,0,0)->v[0], ref(0,1,0)->v[1],
#   ref(0,0,1)->v[2], ref(0,0,0)->v[3]; writing corners in el.vertices
#   order produced a NEGATIVE Jacobian for every element -- the 2026-08
#   docs/gmsh_animation inversion incident.)
# GMSH HEX: (-1,-1,-1) to (1,1,1); NGSolve: (0,0,0) to (1,1,1)
#   transform: ng = (gmsh + 1) / 2
# GMSH PYR: (-1,-1,0) base, (0,0,1) apex; NGSolve: (0,0,0)-(1,0,0)-(1,1,0)-(0,1,0)-(0,0,1)
#   transform: ng_x = (gmsh_x + 1)/2, ng_y = (gmsh_y + 1)/2, ng_z = gmsh_z

def _gmsh_ref_to_ngsolve_ref(et_name, gmsh_x, gmsh_y, gmsh_z,
                             orient="neg"):
    """Convert GMSH reference coordinates to NGSolve reference coordinates.

    GMSH and NGSolve use different vertex numbering in reference space.
    For TET the two reference simplices coincide point-by-point under the
    [3,0,1,2] corner permutation, so no coordinate transform is needed.

    orient: "neg" for the netgen-native vertex listing (the mapping
    below, as probed 2026-08); "pos" for the MIRRORED listing that
    Cubit-plugin .vol files carry (probed 2026-08-06: same ref-corner
    correspondence, opposite geometric orientation).  The mirror is one
    gx<->gy swap of the GMSH reference, uniform across element types,
    paired with _NGSOLVE_TO_GMSH_NODE_ORDER_POS.
    """
    gx, gy, gz = gmsh_x, gmsh_y, gmsh_z
    if orient == "pos":
        gx, gy = gy, gx
    if et_name == 'TET':
        # Identity: GMSH corner i == NGSolve reference vertex i under
        # _NGSOLVE_TO_GMSH_NODE_ORDER['TET'] = [3,0,1,2].
        return (gx, gy, gz)
    elif et_name == 'TRIG':
        return (1 - gx - gy, gx, 0)
    elif et_name == 'HEX':
        # Must match _NGSOLVE_TO_GMSH_NODE_ORDER['HEX'] = [0,1,5,4,3,2,6,7]:
        # the trilinear map with T(gmsh corner g) == NGSolve ref of
        # el.vertices[perm[g]] swaps the y and z axes (probed 2026-08 on
        # MakeStructured3DMesh; the previous identity-axes transform put
        # every Hex27 high-order node on the wrong edge/face -> all
        # Jacobian points negative, volume -7/9 for the unit cube).
        return ((gx + 1) / 2, (gz + 1) / 2, (gy + 1) / 2)
    elif et_name == 'QUAD':
        return ((gx + 1) / 2, (gy + 1) / 2, (gz + 1) / 2)
    elif et_name == 'PYRAMID':
        # Bijection between the GMSH frustum-style reference pyramid
        # (|gx|,|gy| <= 1-gz) and NGSolve's (0 <= x,y <= 1-z): shrink the
        # xy square towards the apex.  Corners AND apex map exactly onto
        # NGSolve reference vertices with the identity corner permutation
        # (probed 2026-08 on a hand-built netgen pyramid: el.vertices[i]
        # sits at NGSolve ref vertex i).
        return ((gx + 1 - gz) / 2, (gy + 1 - gz) / 2, gz)
    elif et_name == 'PRISM':
        # Must match _NGSOLVE_TO_GMSH_NODE_ORDER['PRISM'] = [0,2,1,3,5,4]:
        # triangle barycentric cycle (1-gx-gy, gy) plus the (-1,1) -> (0,1)
        # z rescale (probed 2026-08; the previous identity triangle map
        # was inconsistent with the corner permutation).
        return (1 - gx - gy, gy, (gz + 1) / 2)
    return (gx, gy, gz)


def _get_gmsh_ref_points(et_name, order, orient="neg"):
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
        _, dimension, _, nn, ref_pts, _ = gmsh.model.mesh.getElementProperties(gmsh_type)
        if need_init:
            gmsh.finalize()
    except Exception:
        return None, None

    # Convert GMSH ref -> NGSolve ref
    ng_pts = []
    for i in range(nn):
        offset = dimension * i
        gx = ref_pts[offset]
        gy = ref_pts[offset + 1] if dimension >= 2 else 0.0
        gz = ref_pts[offset + 2] if dimension >= 3 else 0.0
        ng_pts.append(_gmsh_ref_to_ngsolve_ref(et_name, gx, gy, gz,
                                               orient=orient))

    return gmsh_type, ng_pts

# GMSH corner position i holds NGSolve el.vertices[perm[i]].
# TET is [3,0,1,2] (NOT identity): NGSolve's el.vertices listing order is
# NEGATIVELY oriented as a GMSH tet (det[v1-v0,v2-v0,v3-v0] = -det(trafo)),
# so writing corners in listing order inverts every element.  [3,0,1,2]
# places corners at GMSH==NGSolve reference positions and keeps det > 0.
# Verified via gmsh getJacobians (0 negative points, volume vs CAD).
_NGSOLVE_TO_GMSH_NODE_ORDER = {
    'TET': [3, 0, 1, 2],
    'HEX': [0, 1, 5, 4, 3, 2, 6, 7],
    'PRISM': [0, 2, 1, 3, 5, 4],
    # PYRAMID is identity: NGSolve pyramid vertices are listed in
    # reference order (base counter-clockwise + apex), matching GMSH.
    # The previous [3,2,1,0,4] reversed the base cycle -> negative
    # Jacobian on every PYR5 (2026-08 audit).
    'PYRAMID': [0, 1, 2, 3, 4],
    'TRIG': [0, 1, 2],
    'QUAD': [0, 1, 2, 3],
}

# MIRRORED corner permutations for positively-listed meshes (Cubit
# plugin .vol): the gx<->gy swap of the GMSH reference induces these
# corner transpositions (TET (1 2); HEX (1 3)(5 7); PRISM (1 2)(4 5);
# PYRAMID (1 3)) composed with the "neg" tables above.  Selected PER
# ELEMENT from the measured corner-frame orientation; every element is
# verified after reordering (fail loud on a non-positive frame).
_NGSOLVE_TO_GMSH_NODE_ORDER_POS = {
    'TET': [3, 1, 0, 2],
    'HEX': [0, 4, 5, 1, 3, 7, 6, 2],
    'PRISM': [0, 1, 2, 3, 4, 5],
    'PYRAMID': [0, 3, 2, 1, 4],
    'TRIG': [0, 1, 2],
    'QUAD': [0, 1, 2, 3],
}

# GMSH-convention corner frame whose determinant must be POSITIVE for a
# valid (non-inverted) linear element, per type: (origin, a, b, c).
_GMSH_FRAME_CORNERS = {
    'TET': (0, 1, 2, 3),
    'HEX': (0, 1, 3, 4),
    'PRISM': (0, 1, 2, 3),
    'PYRAMID': (0, 1, 3, 4),
}


def _gmsh_corner_frame_det(et_name, pts):
    """det of the GMSH corner frame; > 0 for a valid element."""
    frame = _GMSH_FRAME_CORNERS.get(et_name)
    if frame is None:
        return 1.0
    o, ia, ib, ic = frame
    ax = [pts[ia][k] - pts[o][k] for k in range(3)]
    bx = [pts[ib][k] - pts[o][k] for k in range(3)]
    cx = [pts[ic][k] - pts[o][k] for k in range(3)]
    return (ax[0] * (bx[1] * cx[2] - bx[2] * cx[1])
            - ax[1] * (bx[0] * cx[2] - bx[2] * cx[0])
            + ax[2] * (bx[0] * cx[1] - bx[1] * cx[0]))


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

        Phase E.1 (2026-04-15): now goes through the unified
        _emit_v41_msh emitter via a PayloadV41 built here.
        """
        if not filename.endswith('.msh'):
            filename += '.msh'

        mesh = self.mesh
        is_surface = _is_surface_mesh(mesh) or self._boundary

        nodes, mat_names, elem_data = _extract_mesh_data_grouped(
            mesh, is_surface)
        n_nodes = len(nodes)
        mat_elem_map, _elem_id_map, n_elems = _assign_element_ids(
            mat_names, elem_data)
        # Top-element dimension: 2 for a BEM surface export OR a native 2-D
        # domain mesh (its VOL elements ARE triangles); 3 only for a real
        # 3-D volume mesh.  Without the mesh.dim==2 check a 2-D mesh wrote
        # Tri3 elements into a dim-3 "volume" entity and GMSH rejected the
        # file ("Trying to add unsupported element in volume 1").
        dim = 2 if (is_surface or mesh.dim == 2) else 3
        mat_bboxes = _compute_material_bboxes(nodes, mat_names, elem_data)

        payload = self._build_single_mesh_payload(
            nodes, mat_names, elem_data, mat_elem_map, n_elems,
            dim, mat_bboxes, time, timestep)

        with open(filename, 'w', encoding='utf-8') as f:
            _emit_v41_msh(f, payload)

        # Companion launch artifacts:
        #   case.geo      -- user-facing Gmsh launch recipe
        #   case.geo.opt  -- exact sidecar auto-loaded when case.geo opens
        #   case.msh.opt  -- raw mesh/data inspection sidecar
        # Per-view ncomp follows the order GmshPostExport.add_field() was
        # called, avoiding the old "vectors come first" assumption.
        view_ncomps = [nc for _, nc, _, _, _ in self._fields]
        write_companion_opt(filename, view_ncomps=view_ncomps)

        n_fields = len(self._fields)
        mat_str = ', '.join(f'{m}({len(mat_elem_map.get(m, []))})' for m in mat_names)
        print(f"GMSH export: {filename}")
        print(f"  {n_elems} elements, {n_nodes} nodes, {n_fields} fields")
        print(f"  Materials: {mat_str}")
        return filename

    def _build_single_mesh_payload(self, nodes, mat_names, elem_data,
                                     mat_elem_map, n_elems, dim, mat_bboxes,
                                     time, timestep) -> PayloadV41:
        """Translate single-mesh state into a PayloadV41.

        One PhysicalGroup + one Entity per material. Nodes share a
        single block (entity tag = 1, covering the whole mesh). Elements
        are grouped by (material, gmsh_type) — each group becomes one
        ElementBlock under the corresponding material's entity tag.
        """
        # Physical groups + entities, one per material
        groups = []
        entities = []
        for i, name in enumerate(mat_names):
            tag = i + 1
            groups.append(PhysicalGroup(dim=dim, tag=tag, name=name))
            bbox = mat_bboxes.get(name, (0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
            entities.append(Entity(
                dim=dim, tag=tag, bbox=bbox,
                physical_tags=[tag], bounding_tags=[]))

        # Nodes: one block covering all vertices (entity_tag=1)
        node_blocks = [NodeBlock(
            entity_dim=dim, entity_tag=1, start_id=1, coords=list(nodes))]

        # Elements: one block per (material, gmsh_type), entity_tag =
        # that material's index+1. Connectivity uses global 1-indexed
        # node IDs (equal to local index + 1 since nodes are flat).
        element_blocks = []
        for i, mat_name in enumerate(mat_names):
            entity_tag = i + 1
            type_groups = {}
            for mn, gt, conn, oi in elem_data:
                if mn != mat_name:
                    continue
                type_groups.setdefault(gt, []).append(conn)
            for gmsh_type in sorted(type_groups.keys()):
                conns_1indexed = [[n + 1 for n in conn]
                                   for conn in type_groups[gmsh_type]]
                element_blocks.append(ElementBlock(
                    entity_dim=dim, entity_tag=entity_tag,
                    gmsh_type=gmsh_type, connectivity=conns_1indexed))

        # Data fields (preserving the material-filter semantics of the
        # original GmshPostExport.write path)
        data_fields = []
        n_nodes = len(nodes)
        for name, ncomp, data, is_cell, material in self._fields:
            if not isinstance(data, np.ndarray):
                # CoefficientFunction, evaluated lazily now that the
                # emitted node list is known.
                if is_cell:
                    data = _evaluate_cf(self.mesh, data, ncomp, True)
                else:
                    data = _evaluate_cf_at_points(self.mesh, data, ncomp,
                                                  nodes)
            elif not is_cell and len(data) < n_nodes:
                # Vertex-length array on a mesh that emits high-order
                # nodes: embed P1 data exactly (edge mid = corner mean,
                # face/center = corner means). Leaving those nodes
                # without data corrupts the quadratic view interpolation.
                data = _expand_vertex_data_to_nodes(
                    np.asarray(data), ncomp, n_nodes, elem_data,
                    field_name=name)
            if is_cell:
                if material is not None:
                    elem_ids = mat_elem_map.get(material, [])
                    if not elem_ids:
                        print(f"  WARNING: material '{material}' not found, "
                              f"skipping field '{name}'")
                        continue
                else:
                    elem_ids = []
                    for mn in mat_names:
                        elem_ids.extend(mat_elem_map.get(mn, []))
                items = []
                for orig_idx, gmsh_id in elem_ids:
                    if ncomp == 1:
                        v = (float(data[orig_idx]) if orig_idx < len(data)
                             else 0.0)
                    else:
                        v = (data[orig_idx] if orig_idx < len(data)
                             else np.zeros(ncomp))
                    items.append((gmsh_id, v))
                data_fields.append(DataField(
                    kind='$ElementData', name=name, ncomp=ncomp,
                    items=items, time=time, timestep=timestep))
            else:
                if material is not None:
                    node_set = set()
                    for mat_name, _gt, conn, _oi in elem_data:
                        if mat_name == material:
                            node_set.update(conn)
                    indices = sorted(node_set)
                else:
                    # vertex-only data: only emit nodes that have entries
                    n_data = (len(data) if hasattr(data, '__len__')
                              else n_nodes)
                    indices = list(range(min(n_nodes, n_data)))
                items = []
                for node_idx in indices:
                    if ncomp == 1:
                        v = (float(data[node_idx]) if node_idx < len(data)
                             else 0.0)
                    else:
                        v = (data[node_idx] if node_idx < len(data)
                             else np.zeros(ncomp))
                    items.append((node_idx + 1, v))
                data_fields.append(DataField(
                    kind='$NodeData', name=name, ncomp=ncomp,
                    items=items, time=time, timestep=timestep))

        return PayloadV41(
            groups=groups, entities=entities,
            node_blocks=node_blocks, element_blocks=element_blocks,
            data_fields=data_fields)

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
        """Validate field data; keep CoefficientFunctions for lazy eval.

        Arrays pass through.  A CoefficientFunction/GridFunction is kept
        as-is and evaluated in write(): nodal fields are evaluated at
        EVERY emitted node (including the high-order nodes of a curved
        mesh -- vertex-only values corrupt the quadratic interpolation),
        so the node list must exist first.
        """
        if isinstance(data, np.ndarray):
            return data

        try:
            from ngsolve import CoefficientFunction
        except ImportError:
            CoefficientFunction = ()
        if isinstance(data, CoefficientFunction):
            return data

        raise TypeError(
            f"Unsupported data type: {type(data)}. "
            "Expected CoefficientFunction, GridFunction, or numpy.ndarray.")


def export_element_activation_animation(
        mesh, active_steps, filename, *, field_name="Iron shape",
        times=None, require_monotone=True):
    """Write a Gmsh v4.1 time-series view of element activation.

    The design-superset mesh is written once and every accepted topology is a
    separate ``ElementData`` time step.  Only active volume elements receive
    data, so inactive candidate cells remain hidden in the Gmsh view.
    """
    import os

    filename = os.fspath(filename)
    if not filename.lower().endswith(".msh"):
        filename += ".msh"
    if getattr(mesh, "dim", None) != 3:
        raise ValueError("element activation animation requires a 3-D volume mesh")

    raw_steps = list(active_steps)
    if not raw_steps:
        raise ValueError("active_steps must contain at least one design state")
    masks = []
    for step in raw_steps:
        raw = np.asarray(step)
        if raw.dtype == bool:
            mask = raw.reshape(-1).copy()
        else:
            try:
                numeric = np.asarray(step, dtype=float).reshape(-1)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    "active steps must contain only boolean or binary values") from exc
            if (not np.all(np.isfinite(numeric))
                    or np.any((numeric != 0.0) & (numeric != 1.0))):
                raise ValueError(
                    "active steps must contain only boolean or binary values")
            mask = numeric.astype(bool)
        masks.append(mask)
    if any(mask.size != int(mesh.ne) for mask in masks):
        sizes = [int(mask.size) for mask in masks]
        raise ValueError(
            f"each active step must have mesh.ne={mesh.ne} entries; got {sizes}")
    if any(not np.any(mask) for mask in masks):
        raise ValueError("every animation step must contain at least one active element")
    if require_monotone:
        for index, (previous, current) in enumerate(zip(masks, masks[1:]), 1):
            removed = np.flatnonzero(previous & ~current)
            if removed.size:
                raise ValueError(
                    f"animation step {index} removes active elements "
                    f"{removed[:8].tolist()}")

    if times is None:
        step_times = np.arange(len(masks), dtype=float)
    else:
        step_times = np.asarray(times, dtype=float).reshape(-1)
        if step_times.size != len(masks):
            raise ValueError(
                f"times has {step_times.size} entries for {len(masks)} steps")
        if not np.all(np.isfinite(step_times)):
            raise ValueError("times must contain only finite values")

    nodes, mat_names, elem_data = _extract_mesh_data_grouped(mesh, False)
    mat_elem_map, elem_id_map, n_elems = _assign_element_ids(
        mat_names, elem_data)
    if n_elems != int(mesh.ne) or len(elem_id_map) != int(mesh.ne):
        raise RuntimeError(
            "Gmsh export omitted volume elements; activation indices would "
            f"be ambiguous ({n_elems} exported, mesh.ne={mesh.ne})")
    bboxes = _compute_material_bboxes(nodes, mat_names, elem_data)
    payload = GmshPostExport(mesh)._build_single_mesh_payload(
        nodes, mat_names, elem_data, mat_elem_map, n_elems, 3, bboxes,
        float(step_times[0]), 0)
    payload.data_fields = []
    for step, (mask, time_value) in enumerate(zip(masks, step_times)):
        items = [
            (elem_id_map[int(element)], 1.0)
            for element in np.flatnonzero(mask)
        ]
        payload.data_fields.append(DataField(
            kind="$ElementData", name=str(field_name), ncomp=1,
            items=items, time=float(time_value), timestep=step))

    parent = os.path.dirname(os.path.abspath(filename))
    os.makedirs(parent, exist_ok=True)
    with open(filename, "w", encoding="utf-8") as stream:
        _emit_v41_msh(stream, payload)
    opt_path = write_companion_opt(filename, view_ncomps=[1])
    geo_path = _geo_filename_for_msh(filename)
    return {
        "msh": os.path.abspath(filename),
        "geo": os.path.abspath(geo_path),
        "msh_opt": os.path.abspath(opt_path),
        "geo_opt": os.path.abspath(geo_path + ".opt"),
        "steps": len(masks),
        "active_counts": [int(np.count_nonzero(mask)) for mask in masks],
        "field_name": str(field_name),
    }


def export_nodal_deformation_animation(
        mesh, displacement_steps, filename, *, field_name="Shape displacement",
        times=None, displacement_factor=1.0):
    """Write a moving-mesh Gmsh animation as vector ``NodeData`` steps.

    Each step may be an array in exported Gmsh-node order or an NGSolve
    callable vector field.  The callable form is the stable contract for
    curved/high-order meshes because the exporter samples every generated
    Lagrange node through the NGSolve mesh map.
    """
    import os

    filename = os.fspath(filename)
    if not filename.lower().endswith(".msh"):
        filename += ".msh"
    if getattr(mesh, "dim", None) != 3:
        raise ValueError("nodal deformation animation requires a 3-D volume mesh")
    raw_steps = list(displacement_steps)
    if not raw_steps:
        raise ValueError("displacement_steps must not be empty")
    if times is None:
        step_times = np.arange(len(raw_steps), dtype=float)
    else:
        step_times = np.asarray(times, dtype=float).reshape(-1)
        if (step_times.size != len(raw_steps)
                or not np.all(np.isfinite(step_times))):
            raise ValueError(
                "times must provide one finite value per displacement step")
    factor = float(displacement_factor)
    if not np.isfinite(factor) or factor <= 0.0:
        raise ValueError("displacement_factor must be positive and finite")

    nodes, mat_names, elem_data = _extract_mesh_data_grouped(mesh, False)
    expected = (len(nodes), 3)
    steps = []
    for value in raw_steps:
        if callable(value):
            try:
                sampled = np.asarray(
                    [value(mesh(*point)) for point in nodes], dtype=float)
            except Exception as exc:
                raise ValueError(
                    "callable displacement steps must evaluate to finite "
                    "three-vectors at every exported node") from exc
        else:
            sampled = np.asarray(value, dtype=float)
        if sampled.shape != expected or not np.all(np.isfinite(sampled)):
            raise ValueError(
                f"each displacement step must have shape {expected}, or be "
                "an NGSolve-callable vector field evaluated at exported "
                f"nodes; got {sampled.shape}")
        steps.append(np.ascontiguousarray(sampled))
    mat_elem_map, _, n_elems = _assign_element_ids(mat_names, elem_data)
    bboxes = _compute_material_bboxes(nodes, mat_names, elem_data)
    payload = GmshPostExport(mesh)._build_single_mesh_payload(
        nodes, mat_names, elem_data, mat_elem_map, n_elems, 3, bboxes,
        float(step_times[0]), 0)
    payload.data_fields = []
    for step, (values, time_value) in enumerate(zip(steps, step_times)):
        payload.data_fields.append(DataField(
            kind="$NodeData", name=str(field_name), ncomp=3,
            items=[(index + 1, tuple(vector))
                   for index, vector in enumerate(values)],
            time=float(time_value), timestep=step))
    parent = os.path.dirname(os.path.abspath(filename))
    os.makedirs(parent, exist_ok=True)
    with open(filename, "w", encoding="utf-8") as stream:
        _emit_v41_msh(stream, payload)
    opt_path = write_companion_opt(filename, view_ncomps=[3])
    geo_path = _geo_filename_for_msh(filename)
    displacement_options = (
        "View[0].VectorType = 5;\n"
        f"View[0].DisplacementFactor = {factor:.15g};\n")
    for option_path in (opt_path, geo_path, geo_path + ".opt"):
        with open(option_path, "a", encoding="utf-8") as stream:
            stream.write("// Animate the vector view as mesh displacement\n")
            stream.write(displacement_options)
    return {
        "msh": os.path.abspath(filename),
        "geo": os.path.abspath(geo_path),
        "msh_opt": os.path.abspath(opt_path),
        "geo_opt": os.path.abspath(geo_path + ".opt"),
        "steps": len(steps),
        "field_name": str(field_name),
        "maximum_displacement_m": [
            float(np.max(np.linalg.norm(value, axis=1))) for value in steps],
    }


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
        # 2D meshes (axisym, planar) carry only (x, y); pad z=0 for
        # the 3D-anchored .msh format (GMSH itself is happy with a
        # flat mesh in z=0).
        z = pt[2] if len(pt) >= 3 else 0.0
        nodes.append((pt[0], pt[1], z))

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

    # Python fallback for high-order (if C++ not available).
    # elem_ho_nodes is consumed ONLY by the surface (BND) branch below;
    # the volume branch generates its own high-order nodes via
    # _build_vol_ho_generic. Calling this for volume exports appended
    # orphan surface mid-edge nodes that no element referenced.
    elem_ho_nodes = {}
    if use_highorder and is_surface:
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
            # Vertex-listing orientation is MESH-ORIGIN dependent
            # (measured 2026-08-06: netgen-generated meshes list tet
            # corners negatively, Cubit-plugin .vol files positively --
            # same reference correspondence, mirrored geometry), so the
            # corner permutation is selected PER ELEMENT: try the
            # netgen-native table, and if the GMSH corner frame comes
            # out negative use the mirrored table.  A frame that is
            # non-positive under BOTH tables is a genuinely degenerate/
            # inverted element -> fail loud.
            perm = _NGSOLVE_TO_GMSH_NODE_ORDER.get(
                et_name, list(range(len(verts))))
            reordered = [verts[i] for i in perm]
            orient = "neg"
            if et_name in _GMSH_FRAME_CORNERS:
                det = _gmsh_corner_frame_det(
                    et_name, [nodes[v] for v in reordered])
                if det <= 0:
                    perm = _NGSOLVE_TO_GMSH_NODE_ORDER_POS.get(
                        et_name, perm)
                    reordered = [verts[i] for i in perm]
                    orient = "pos"
                    det = _gmsh_corner_frame_det(
                        et_name, [nodes[v] for v in reordered])
                if det <= 0:
                    raise RuntimeError(
                        f"GmshPostExport: volume element {idx} "
                        f"({et_name}) has a non-positive GMSH corner "
                        f"frame under both orientation tables "
                        f"(det={det:.3e}) -- degenerate or inverted "
                        "mesh element")

            if use_highorder and et_name in _GMSH_TYPE_BY_ORDER:
                effective_order = curve_order
                if effective_order not in _GMSH_TYPE_BY_ORDER.get(et_name, {}):
                    effective_order = min(curve_order, 2)
                if effective_order >= 2:
                    reordered = _build_vol_ho_generic(
                        mesh, el, reordered, et_name, effective_order,
                        nodes, vol_ho_cache, _vol_gmsh_ref,
                        orient=orient)

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

    # Map GMSH corner index -> NGSolve el.vertices index -> ref coord.
    # reordered[g] = verts[perm[g]], so GMSH local g corresponds to
    # NGSolve local perm[g].  (The pre-2026-08 code indexed the INVERSE
    # permutation here; every permutation in the table happened to be an
    # involution back then, which made the mix-up unobservable.  The TET
    # fix [3,0,1,2] is not an involution, so the direction matters.)
    perm = _NGSOLVE_TO_GMSH_NODE_ORDER.get(
        et_name, list(range(len(reordered))))

    for e0_gmsh, e1_gmsh in edge_table:
        # Global node IDs for this edge
        gn0 = reordered[e0_gmsh]
        gn1 = reordered[e1_gmsh]
        edge_key = (min(gn0, gn1), max(gn0, gn1))

        if edge_key in cache:
            conn.append(cache[edge_key])
        else:
            # Map GMSH local index -> NGSolve local index -> ref coord
            ng0 = perm[e0_gmsh]
            ng1 = perm[e1_gmsh]
            ref_mid = tuple(0.5 * (ref_verts[ng0][k] + ref_verts[ng1][k])
                            for k in range(3))

            # Evaluate trafo at mid-point to get physical coordinates
            from ngsolve import IntegrationRule
            ir = IntegrationRule([ref_mid], [1.0])
            mip = trafo(ir[0])
            phys = _point3(mip.point)
            mid_node_idx = len(nodes)
            nodes.append(phys)
            cache[edge_key] = mid_node_idx
            conn.append(mid_node_idx)

    return conn


def _build_vol_ho_generic(mesh, el, reordered, et_name, order,
                          nodes, cache, ref_cache, orient="neg"):
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
    # Get or compute reference points (cached per type AND orientation)
    ref_key = (et_name, orient)
    if ref_key not in ref_cache:
        ref_cache[ref_key] = _get_gmsh_ref_points(et_name, order,
                                                  orient=orient)
    gmsh_type, ng_ref_pts = ref_cache[ref_key]

    if ng_ref_pts is None:
        # Fallback to order-2 path
        return _build_vol_highorder_conn(
            mesh, el, reordered, et_name, nodes,
            cache.setdefault("_order2", {}))

    n_verts = len(reordered)
    n_total = len(ng_ref_pts)

    # Vertex nodes come straight from `reordered` (already in GMSH corner
    # order); ng_ref_pts[i] below is likewise already expressed in NGSolve
    # reference coordinates, so no further permutation is needed here.

    # Evaluate trafo at all high-order reference points
    trafo = mesh.GetTrafo(el)
    from ngsolve import IntegrationRule

    conn = list(reordered)  # vertices first

    for ho_idx in range(n_verts, n_total):
        ref_pt = ng_ref_pts[ho_idx]

        # Evaluate physical coordinates
        ir = IntegrationRule([ref_pt], [1.0])
        mip = trafo(ir[0])
        phys = _point3(mip.point)

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


def _gmsh_display_options_text(vector_indices, scalar_indices):
    """Return the standard Gmsh display options shared by .geo and .opt."""
    lines = [
        "General.Orthographic = 1;",
        "General.Trackball = 0;",
        "General.RotationX = -68;",
        "General.RotationY = 0;",
        "General.RotationZ = 0;",
        "General.RotationCenterGravity = 1;",
        "Mesh.NumSubEdges = 4;",
        "Mesh.VolumeEdges = 0;",
        "Mesh.VolumeFaces = 0;",
        "Mesh.SurfaceEdges = 0;",
        "Mesh.SurfaceFaces = 1;",
        "Mesh.Points = 0;",
        "Mesh.Lines = 0;",
        "Mesh.LineWidth = 1;",
        "Mesh.ColorCarousel = 2;",
    ]
    for i in vector_indices:
        lines.extend([
            f"View[{i}].Visible = 1;",
            f"View[{i}].VectorType = 4;",
            f"View[{i}].ArrowSizeMin = 20;",
            f"View[{i}].ArrowSizeMax = 20;",
        ])
    for idx in scalar_indices:
        lines.extend([
            f"View[{idx}].Visible = 1;",
            f"View[{idx}].IntervalsType = 2;",
            f"View[{idx}].ScaleType = 2;",
            f"View[{idx}].ShowScale = 1;",
            f"View[{idx}].NbIso = 20;",
        ])
    return "\n".join(lines) + "\n"


def _geo_filename_for_msh(msh_filename):
    import os
    base, ext = os.path.splitext(msh_filename)
    if ext.lower() == ".msh":
        return base + ".geo"
    return msh_filename + ".geo"


def _relative_gmsh_merge_path(msh_filename, geo_filename):
    import os
    rel = os.path.relpath(
        os.path.abspath(msh_filename),
        os.path.dirname(os.path.abspath(geo_filename)) or ".",
    )
    return rel.replace("\\", "/")


def write_companion_opt(msh_filename, n_vector_views=0,
                        n_scalar_views=0, view_ncomps=None):
    """Write Gmsh launch companions for a .msh post-processing artifact.

    The Radia user-facing launch target is ``case.geo``. GMSH auto-loads
    ``case.geo.opt`` when opening that file. ``case.msh.opt`` is kept as the
    raw mesh/data inspection sidecar for users who intentionally open the
    .msh directly.

    Settings follow the gmsh_post_spec (mandatory):
      - Mesh.NumSubEdges = 4      (render curved Tri6+ smoothly)
      - Mesh.VolumeEdges = 0      (hide air tet edges)
      - Mesh.VolumeFaces = 0      (hide air tet faces)
      - Mesh.SurfaceFaces = 1     (show coil as smooth body)
      - View[i].VectorType = 4    (3D arrows for vectors)
      - View[i].ArrowSizeMin = 20 (fixed 20 px -- THE kirei knob)
      - View[i].IntervalsType = 2 (filled iso for scalars)

    Two ways to specify which views are vectors vs scalars:

    1. ``view_ncomps``: list of ncomp per view IN ORDER, e.g.
       ``[1, 1, 3]`` means View[0]=scalar, View[1]=scalar, View[2]=vec.
       PREFERRED -- avoids assumptions about ordering.
    2. ``n_vector_views`` + ``n_scalar_views`` (legacy): assumes vectors
       come FIRST.  Wrong when the caller appends views in any other
       order (this caused q_surf to be misclassified as a vector and
       displayed with linear scale, hiding its log-distributed hot
       spot -- 2026-04-29 calc_peec_bem GUI bug).

    Args:
        msh_filename: Path to the .msh file. The companions are written next
            to it as ``case.geo``, ``case.geo.opt``, and ``case.msh.opt``.
        n_vector_views: legacy; vectors are View[0..n-1].
        n_scalar_views: legacy; scalars are View[n..n+m-1].
        view_ncomps: per-view component counts (preferred over the
            legacy split args).  Indices in this list match the view
            indices in the .msh.

    Returns:
        Path to the .msh.opt file for backward compatibility. New callers
        should open the generated ``case.geo``.
    """
    import os
    msh_filename = os.fspath(msh_filename)
    # Resolve which view indices are vectors vs scalars.
    if view_ncomps is not None:
        vector_indices = [i for i, n in enumerate(view_ncomps)
                          if (n or 0) >= 2]
        scalar_indices = [i for i, n in enumerate(view_ncomps)
                          if (n or 0) < 2]
    else:
        # Legacy split (vectors first).
        vector_indices = list(range(n_vector_views))
        scalar_indices = list(range(n_vector_views,
                                     n_vector_views + n_scalar_views))
    options_text = _gmsh_display_options_text(vector_indices, scalar_indices)

    opt_filename = msh_filename + '.opt'
    with open(opt_filename, 'w', encoding='utf-8') as f:
        f.write('// Auto-generated display options (gmsh_post_spec)\n')
        f.write(options_text)

    geo_filename = _geo_filename_for_msh(msh_filename)
    merge_path = _relative_gmsh_merge_path(msh_filename, geo_filename)
    with open(geo_filename, 'w', encoding='utf-8') as f:
        f.write('// Auto-generated Radia Gmsh launch companion\n')
        f.write('// Open this .geo for normal review; open .msh only for raw inspection.\n')
        f.write(f'Merge "{merge_path}";\n\n')
        f.write(options_text)

    geo_opt_filename = geo_filename + '.opt'
    with open(geo_opt_filename, 'w', encoding='utf-8') as f:
        f.write('// Auto-generated display options for the .geo launch artifact\n')
        f.write(options_text)

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
    companion_opt: if True, also write ``case.geo``, ``case.geo.opt``, and
        ``case.msh.opt``. Open ``case.geo`` for normal review.

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
    view_ncomps = []
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
        view_ncomps.append(ncomp or 1)
        if (ncomp or 0) >= 2:
            n_vec += 1
        else:
            n_sca += 1

    post.write(output_msh)
    if companion_opt:
        # Pass view_ncomps so write_companion_opt classifies each view
        # by its actual component count (vector vs scalar) rather than
        # assuming vectors-first ordering.  This was the cause of the
        # 2026-04-29 calc_peec_bem GUI bug: q_surf (scalar appended at
        # idx 0) was misclassified as a vector and rendered with
        # default linear color scale -- hot helical bands invisible.
        write_companion_opt(output_msh,
                             view_ncomps=view_ncomps)
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

    Physical groups (contiguous node ID space across three meshes):

      Physical "air"                (3D, tag 1): volume tets   + B NodeData
      Physical "coil_surface"       (2D, tag 2): Tri6 surface  + J NodeData
      Physical "coil_wire"          (1D, tag 3): wireframe edges (no data)
      Physical "workpiece_surface"  (2D, tag 4): workpiece tris + q NodeData
                                                  (optional)

    Coil surface uses Tri6 (GMSH type 9) for curved rendering.
    Writes companion .geo/.geo.opt/.msh.opt display artifacts.

    Phase E.1 (2026-04-15): builds a PayloadV41 and emits via
    _emit_v41_msh (shared with GmshPostExport.write).
    """
    surf_nodes, surf_elems, wire_edges = extract_tri6_surface(mesh_surf)
    nv_surf = len(surf_nodes)
    nv_vol = len(vol_nodes)
    surf_offset = nv_vol
    has_wp = (wp_nodes is not None and wp_elems is not None
              and q_nodes is not None)
    nv_wp = len(wp_nodes) if has_wp else 0
    wp_offset = nv_vol + nv_surf

    groups = [
        PhysicalGroup(dim=3, tag=1, name='air'),
        PhysicalGroup(dim=2, tag=2, name='coil_surface'),
        PhysicalGroup(dim=1, tag=3, name='coil_wire'),
    ]
    if has_wp:
        groups.append(PhysicalGroup(dim=2, tag=4, name='workpiece_surface'))

    # Preserve historical entity-ordering (curve before surfaces before
    # volume) by building the list in that order.
    entities = [
        Entity(dim=1, tag=3, bbox=(0, 0, 0, 0, 0, 0),
               physical_tags=[3], bounding_tags=[]),
        Entity(dim=2, tag=2, bbox=(0, 0, 0, 0, 0, 0),
               physical_tags=[2], bounding_tags=[]),
    ]
    if has_wp:
        entities.append(Entity(dim=2, tag=4, bbox=(0, 0, 0, 0, 0, 0),
                                physical_tags=[4], bounding_tags=[]))
    entities.append(Entity(dim=3, tag=1, bbox=(0, 0, 0, 0, 0, 0),
                            physical_tags=[1], bounding_tags=[]))

    node_blocks = [
        NodeBlock(entity_dim=3, entity_tag=1, start_id=1,
                  coords=list(vol_nodes)),
        NodeBlock(entity_dim=2, entity_tag=2, start_id=surf_offset + 1,
                  coords=list(surf_nodes)),
    ]
    if has_wp:
        node_blocks.append(NodeBlock(
            entity_dim=2, entity_tag=4, start_id=wp_offset + 1,
            coords=list(wp_nodes)))

    # Elements (all connectivity in GLOBAL 1-indexed node IDs)
    vol_conn = [[v + 1 for v in verts] for verts in vol_elems]
    surf_conn = [[v + surf_offset + 1 for v in conn] for conn in surf_elems]
    wire_conn = [[a + surf_offset + 1, b + surf_offset + 1]
                  for a, b in wire_edges]
    element_blocks = [
        ElementBlock(entity_dim=3, entity_tag=1, gmsh_type=4,
                     connectivity=vol_conn),
        ElementBlock(entity_dim=2, entity_tag=2, gmsh_type=9,
                     connectivity=surf_conn),
        ElementBlock(entity_dim=1, entity_tag=3, gmsh_type=1,
                     connectivity=wire_conn),
    ]
    if has_wp:
        wp_conn = [[v + wp_offset + 1 for v in verts] for verts in wp_elems]
        element_blocks.append(ElementBlock(
            entity_dim=2, entity_tag=4, gmsh_type=2, connectivity=wp_conn))

    # Data fields. B on all volume nodes, J only on coil vertex IDs
    # (Tri6 mid-edge nodes are interpolated by GMSH). Scalar q last so
    # companion .geo/.geo.opt/.msh.opt assigns IntervalsType correctly.
    data_fields = [
        DataField(
            kind='$NodeData', name='B', ncomp=3,
            items=[(i + 1, B_nodes[i]) for i in range(nv_vol)]),
        DataField(
            kind='$NodeData', name='J', ncomp=3,
            items=[(surf_offset + i + 1, J_nodes[i])
                    for i in range(mesh_surf.nv)]),
    ]
    if has_wp and wp_J_nodes is not None:
        data_fields.append(DataField(
            kind='$NodeData', name='J_workpiece', ncomp=3,
            items=[(wp_offset + i + 1, wp_J_nodes[i])
                    for i in range(nv_wp)]))
    if has_wp:
        data_fields.append(DataField(
            kind='$NodeData', name='q', ncomp=1,
            items=[(wp_offset + i + 1, q_nodes[i]) for i in range(nv_wp)]))

    payload = PayloadV41(
        groups=groups, entities=entities,
        node_blocks=node_blocks, element_blocks=element_blocks,
        data_fields=data_fields)

    with open(filename, 'w', encoding='utf-8') as f:
        _emit_v41_msh(f, payload)

    n_vec = 2 + (1 if (has_wp and wp_J_nodes is not None) else 0)
    n_sca = 1 if has_wp else 0
    write_companion_opt(filename, n_vector_views=n_vec,
                        n_scalar_views=n_sca)
    return filename


def _write_companion_geo(msh_filename):
    """Backward-compatible alias that now writes .geo and exact .opt twins."""
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
# Phase E.1 (2026-04-15): the per-section writers (_write_physical_names,
# _write_entities, _write_nodes, _write_elements_grouped,
# _write_node_data, _write_node_data_filtered,
# _write_element_data_filtered) were collapsed into _emit_v41_msh /
# _emit_data_section at the top of this file. Both GmshPostExport.write
# and write_combined_msh_v41 now build a PayloadV41 and call the unified
# emitter, so changing the .msh on-disk layout is a one-touch edit.


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
            mip = mesh(*centroid)
            val = cf(mip)
            if ncomp == 1:
                result[i] = float(val)
            else:
                result[i, :] = list(val)[:ncomp]
        return result
    else:
        nv = mesh.nv
        points = [list(v.point) for v in mesh.vertices]
        return _evaluate_cf_at_points(mesh, cf, ncomp, points)


def _evaluate_cf_at_points(mesh, cf, ncomp, points):
    """Evaluate a CoefficientFunction at explicit physical points.

    Used to sample nodal fields at EVERY emitted .msh node (vertices
    plus the high-order nodes of a curved export). Batch evaluation
    through mesh(x_arr, y_arr, z_arr); a point that cannot be located
    raises -- silently writing zeros produced corrupted quadratic
    views (the 2026-08 vertex-only NodeData bug).
    """
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    mips = mesh(pts[:, 0], pts[:, 1], pts[:, 2])
    vals = np.asarray(cf(mips), dtype=float)
    if ncomp == 1:
        return vals.reshape(-1)
    return vals.reshape(-1, ncomp)


# P1 -> high-order node fill rules per gmsh element type:
# (n_corners, edge pairs, quad faces, has volume center).
# Values on the extra Lagrange nodes of an order-2 element are the
# corner means (exact for P1 data). Order >= 3 layouts are not listed:
# pass a CoefficientFunction instead.
_HO_FILL_RULES = {
    9: (3, [(0, 1), (1, 2), (2, 0)], [], False),                # tri6
    16: (4, [(0, 1), (1, 2), (2, 3), (3, 0)], [], False),       # quad8
    10: (4, [(0, 1), (1, 2), (2, 3), (3, 0)],
         [(0, 1, 2, 3)], False),                                # quad9
    11: (4, _TET_EDGES_GMSH, [], False),                        # tet10
    17: (8, _HEX_EDGES_GMSH, [], False),                        # hex20
    12: (8, _HEX_EDGES_GMSH,
         [(0, 3, 2, 1), (0, 1, 5, 4), (0, 4, 7, 3),
          (1, 2, 6, 5), (2, 3, 7, 6), (4, 5, 6, 7)], True),     # hex27
    18: (6, _PRISM_EDGES_GMSH, [], False),                      # prism15
    13: (6, _PRISM_EDGES_GMSH,
         [(0, 1, 4, 3), (0, 3, 5, 2), (1, 2, 5, 4)], False),    # prism18
}


def _expand_vertex_data_to_nodes(data, ncomp, n_nodes, elem_data,
                                 field_name=""):
    """Embed vertex (P1) field data onto every emitted node.

    Exact for P1 data: each extra order-2 Lagrange node value is the
    mean of its parent corner values (edge midpoint = 2-corner mean,
    quad face node = 4-corner mean, hex center = 8-corner mean).
    Fails loudly for element layouts without a fill rule (order >= 3):
    those exports must pass a CoefficientFunction instead.
    """
    n_vertex = len(data)
    shape = (n_nodes,) if ncomp == 1 else (n_nodes, ncomp)
    out = np.full(shape, np.nan)
    out[:n_vertex] = data

    for _mat, gmsh_type, conn, _oi in elem_data:
        if not any(idx >= n_vertex for idx in conn):
            continue  # linear element: corners only
        rule = _HO_FILL_RULES.get(gmsh_type)
        if rule is None:
            raise ValueError(
                f"field '{field_name}': cannot expand vertex-length data "
                f"onto gmsh element type {gmsh_type}; pass a "
                f"CoefficientFunction (or a full {n_nodes}-node array)")
        n_corners, edges, faces, has_center = rule
        pos = n_corners
        for a, b in edges:
            idx = conn[pos]
            pos += 1
            if idx >= n_vertex:
                out[idx] = 0.5 * (out[conn[a]] + out[conn[b]])
        for face in faces:
            idx = conn[pos]
            pos += 1
            if idx >= n_vertex:
                out[idx] = sum(out[conn[c]] for c in face) / len(face)
        if has_center:
            idx = conn[pos]
            if idx >= n_vertex:
                out[idx] = sum(out[conn[c]]
                               for c in range(n_corners)) / n_corners

    if np.isnan(out).any():
        missing = int(np.isnan(np.atleast_2d(out.T)[0]).sum())
        raise ValueError(
            f"field '{field_name}': {missing} emitted node(s) not covered "
            f"by the vertex-data expansion; pass a CoefficientFunction")
    return out


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
          -> .msh + .geo/.geo.opt/.msh.opt (GMSH visualization)

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


# =====================================================================
# Filament visualization (PEEC sub-filament paths as GMSH line elements)
# =====================================================================

# Line-element .msh tag offsets.  Element tags start high so filament
# lines do not collide with the volume-mesh element tags when GMSH
# merges this file as a sibling overlay (open_gmsh.py auto-merge);
# without the offset GMSH warns "Skipping duplicate element N" for every
# line and silently DROPS the filament geometry.
#
# ENTITY tags need the same offset: merging a STEP file creates OCC
# curves with tags 1..N in the same (dim, tag) space, and a discrete
# curve sharing a tag gets BOUND to the OCC curve as its "mesh" -- gmsh
# then hard-crashes (0xC000041D) as soon as shaded geometry faces are
# enabled (Geometry.Surfaces=1 + Geometry.SurfaceType>0).  Measured
# gmsh 4.15.2; the tag offset resolves it.
_LINE_ELEM_TAG_OFFSET = 10_000_000
_LINE_ENT_TAG_OFFSET = 10_000_000


def _write_line_element_msh(output_msh, segments, groups, group_names,
                            element_views=()):
    """Write 2-node line elements as a GMSH .msh v4.1 post file.

    Shared low-level writer behind ``export_filaments_msh`` (per-filament
    polylines) and ``export_peec_topology_msh`` (per-segment PEEC
    branches).  Coincident endpoints are merged into one node so a
    filament chain renders as a continuous polyline.

    Args:
        output_msh: output .msh path.
        segments: list of ((x1,y1,z1), (x2,y2,z2)) endpoint pairs, in
            the order the elements should be tagged.
        groups: per-segment 0-based physical-group index.
        group_names: physical group names, one per group.
        element_views: iterable of view dicts, each
            ``{"name": str, "ncomp": 1 or 3, "step": int, "time": float,
            "values": array of shape (n_seg,) or (n_seg, ncomp)}``.
            Values are PER ELEMENT (callers expand per-filament data).

    Returns:
        dict with n_nodes / n_elements / element_tags.
    """
    import numpy as np

    n_groups = len(group_names)
    if len(groups) != len(segments):
        raise ValueError(
            f"groups length {len(groups)} != segments length {len(segments)}")

    node_map = {}
    nodes = []

    def _get_node(pt):
        key = (round(pt[0], 12), round(pt[1], 12), round(pt[2], 12))
        if key not in node_map:
            node_map[key] = len(nodes) + 1  # 1-based
            nodes.append(key)
        return node_map[key]

    elements = []          # (tag, group, n1, n2) in `segments` order
    elem_tag = _LINE_ELEM_TAG_OFFSET
    for (p1, p2), g in zip(segments, groups):
        n1 = _get_node(p1)
        n2 = _get_node(p2)
        elem_tag += 1
        elements.append((elem_tag, int(g), n1, n2))

    n_nodes = len(nodes)
    n_elems = len(elements)

    elems_by_group = [[] for _ in range(n_groups)]
    for (tag, g, n1, n2) in elements:
        elems_by_group[g].append((tag, n1, n2))

    with open(output_msh, "w") as f:
        f.write("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")

        f.write("$PhysicalNames\n")
        f.write(f"{n_groups}\n")
        for k, name in enumerate(group_names):
            f.write(f'1 {k+1} "{name}"\n')
        f.write("$EndPhysicalNames\n")

        # One curve entity per group, REAL bounding boxes, tags offset
        # out of the OCC range (see _LINE_ENT_TAG_OFFSET note).
        f.write("$Entities\n")
        f.write(f"0 {n_groups} 0 0\n")   # 0 points, n curves, 0 surf, 0 vol
        for k in range(n_groups):
            gnodes = {n for (_t, a, b) in elems_by_group[k] for n in (a, b)}
            pts = [nodes[n - 1] for n in gnodes]
            if pts:
                lo = [min(p[i] for p in pts) for i in range(3)]
                hi = [max(p[i] for p in pts) for i in range(3)]
            else:
                lo = hi = [0.0, 0.0, 0.0]
            f.write(f"{k + 1 + _LINE_ENT_TAG_OFFSET} "
                    f"{lo[0]:.9e} {lo[1]:.9e} {lo[2]:.9e} "
                    f"{hi[0]:.9e} {hi[1]:.9e} {hi[2]:.9e} 1 {k+1} 0\n")
        f.write("$EndEntities\n")

        # Nodes (one block, all nodes, on the first curve entity)
        f.write("$Nodes\n")
        f.write(f"1 {n_nodes} 1 {n_nodes}\n")
        f.write(f"1 {1 + _LINE_ENT_TAG_OFFSET} 0 {n_nodes}\n")
        for i in range(n_nodes):
            f.write(f"{i+1}\n")
        for (x, y, z) in nodes:
            f.write(f"{x:.15e} {y:.15e} {z:.15e}\n")
        f.write("$EndNodes\n")

        # Elements (one block per group)
        f.write("$Elements\n")
        min_tag = _LINE_ELEM_TAG_OFFSET + 1
        f.write(f"{n_groups} {n_elems} {min_tag} {elem_tag}\n")
        for k in range(n_groups):
            es = elems_by_group[k]
            # dim=1, entityTag (offset), elemType=1 (2-node line), count
            f.write(f"1 {k + 1 + _LINE_ENT_TAG_OFFSET} 1 {len(es)}\n")
            for (tag, n1, n2) in es:
                f.write(f"{tag} {n1} {n2}\n")
        f.write("$EndElements\n")

        for view in element_views:
            ncomp = int(view["ncomp"])
            vals = np.asarray(view["values"], dtype=np.float64)
            if ncomp == 1:
                vals = vals.reshape(-1, 1)
            if vals.shape != (n_elems, ncomp):
                raise ValueError(
                    f"view {view['name']!r}: values shape {vals.shape} "
                    f"!= expected ({n_elems}, {ncomp})")
            step = int(view.get("step", 0))
            f.write("$ElementData\n")
            f.write(f'1\n"{view["name"]}"\n')
            f.write(f"1\n{float(view.get('time', step)):g}\n")
            f.write(f"3\n{step}\n{ncomp}\n{n_elems}\n")
            for (tag, _g, _n1, _n2), row in zip(elements, vals):
                f.write(f"{tag} " + " ".join(f"{v:.6e}" for v in row) + "\n")
            f.write("$EndElementData\n")

    return {"n_nodes": n_nodes, "n_elements": n_elems,
            "element_tags": [e[0] for e in elements]}


def export_filaments_msh(filament_paths, output_msh, *,
                         currents=None, label="filament",
                         direction_view=False,
                         complex_steps=False,
                         viz_subdivide_long=True,
                         viz_subdivide_threshold=5.0,
                         viz_subdivide_n=8):
    """Export PEEC filament paths as GMSH .msh v4.1 line elements.

    Each filament is a polyline of ((x1,y1,z1),(x2,y2,z2)) segments.
    Writes 2-node line elements (GMSH type 1) grouped per filament
    as separate physical groups.  Optionally writes |I| as ElementData
    for current magnitude visualization.

    Usage:
        from radia.gmsh_post_export import export_filaments_msh

        result = filaments_from_step("coil.step", nwinc=3, nhinc=3)
        export_filaments_msh(result["filament_paths"], "filaments.msh")

        # In GMSH GUI:
        #   Merge "coil.step"       <- solid coil overlay
        #   Merge "filaments.msh"   <- filament lines + current
        # or headless via radia-mcp:
        #   gmsh_render("filaments.msh", merge_files=["coil.step"])

    Args:
        filament_paths: list of N filaments.  Each filament is a list
            of ((x1,y1,z1), (x2,y2,z2)) segment endpoint pairs.
        output_msh: Output .msh file path.
        currents: optional complex ndarray (N,) of per-filament currents.
            Written as |I| ElementData.
        direction_view: if True (and currents given), additionally
            write an "I direction" vector ElementData view: the unit
            tangent of every line element scaled by the SIGNED current
            Re(I) of its filament.  Rendered as arrows by gmsh
            (View.VectorType), this shows the CURRENT FLOW direction
            along the winding, and a filament carrying reverse current
            points backwards (the |I| view cannot show that).
        complex_steps: if True (and currents are complex), additionally
            write an "I_complex" ElementData view with TWO time steps
            (step 0 = Re I, step 1 = Im I).  This is the
            harmonic-to-time contract: expanding it (radia-mcp
            gmsh_harmonic_to_time) animates the AC current
            distribution over one period -- skin/proximity phase lead
            of outer filaments becomes visible; gmsh_modulus_phase
            turns the same data into amplitude/phase maps.
        label: Physical group base name (default "filament").
        viz_subdivide_long: if True (default), subdivide segments
            whose length exceeds ``viz_subdivide_threshold`` x the
            per-filament median.  This is a VISUAL-ONLY refinement
            (the PEEC solver kept the original ``filament_paths``);
            without it, a coil whose chain has a lead-to-helix
            endpoint jump (e.g. kubota 3turn coil's 62mm lead
            segment) renders in GMSH as a single line that looks
            like the filament leaves the conductor body.  Set False
            to write raw chain segments.
        viz_subdivide_threshold: factor of per-filament median above
            which a segment is treated as a long jump (default 5.0).
        viz_subdivide_n: number of sub-segments to split a long
            segment into (default 8 -- coarse enough to keep the .msh
            small, fine enough that GMSH renders a continuous line).
    """
    import numpy as np

    try:
        filament_paths = [list(path) for path in filament_paths]
    except TypeError as exc:
        raise ValueError("filament_paths must be an iterable of paths") from exc
    if not filament_paths:
        raise ValueError("filament_paths must contain at least one filament")

    normalized_paths = []
    for fil_idx, path in enumerate(filament_paths):
        if not path:
            raise ValueError(f"filament {fil_idx} has no segments")
        normalized = []
        for seg_idx, segment in enumerate(path):
            try:
                if len(segment) != 2:
                    raise ValueError
                p0 = np.asarray(segment[0], dtype=np.float64)
                p1 = np.asarray(segment[1], dtype=np.float64)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"filament {fil_idx} segment {seg_idx} must contain "
                    "two numeric 3D endpoints") from exc
            if p0.shape != (3,) or p1.shape != (3,):
                raise ValueError(
                    f"filament {fil_idx} segment {seg_idx} endpoints "
                    "must each have exactly three coordinates")
            if not np.all(np.isfinite(p0)) or not np.all(np.isfinite(p1)):
                raise ValueError(
                    f"filament {fil_idx} segment {seg_idx} contains "
                    "non-finite coordinates")
            if float(np.linalg.norm(p1 - p0)) <= 0.0:
                raise ValueError(
                    f"filament {fil_idx} segment {seg_idx} has zero length")
            normalized.append((tuple(p0), tuple(p1)))
        normalized_paths.append(normalized)
    filament_paths = normalized_paths

    try:
        viz_subdivide_threshold = float(viz_subdivide_threshold)
        viz_subdivide_n_int = int(viz_subdivide_n)
    except (TypeError, ValueError) as exc:
        raise ValueError("invalid visualization subdivision controls") from exc
    if (not np.isfinite(viz_subdivide_threshold)
            or viz_subdivide_threshold <= 0.0):
        raise ValueError("viz_subdivide_threshold must be finite and positive")
    if (viz_subdivide_n_int < 1
            or float(viz_subdivide_n) != float(viz_subdivide_n_int)):
        raise ValueError("viz_subdivide_n must be a positive integer")
    viz_subdivide_n = viz_subdivide_n_int

    n_fil = len(filament_paths)
    I_arr = None
    if currents is not None:
        try:
            I_arr = np.asarray(currents)
            dtype = np.complex128 if np.iscomplexobj(I_arr) else np.float64
            I_arr = I_arr.astype(dtype, copy=False)
        except (TypeError, ValueError) as exc:
            raise ValueError("currents must be a numeric 1D array") from exc
        if I_arr.ndim != 1 or len(I_arr) != n_fil:
            raise ValueError(
                f"currents must have shape ({n_fil},), got {I_arr.shape}")
        if not np.all(np.isfinite(I_arr)):
            raise ValueError("currents contains non-finite values")
    elif direction_view or complex_steps:
        raise ValueError(
            "currents is required when direction_view or complex_steps is set")

    if not isinstance(label, str) or not label or any(
            c in label for c in ('"', '\n', '\r')):
        raise ValueError("label must be a non-empty single-line string without quotes")

    if viz_subdivide_long:
        filament_paths_viz = []
        for path in filament_paths:
            if not path:
                filament_paths_viz.append(path)
                continue
            seglens = np.array(
                [np.linalg.norm(np.subtract(s[1], s[0])) for s in path])
            if len(seglens) == 0:
                filament_paths_viz.append(path)
                continue
            med = float(np.median(seglens))
            if med <= 0:
                filament_paths_viz.append(path)
                continue
            thresh = viz_subdivide_threshold * med
            new_path = []
            for s in path:
                p0 = np.asarray(s[0], dtype=np.float64)
                p1 = np.asarray(s[1], dtype=np.float64)
                d = float(np.linalg.norm(p1 - p0))
                if d > thresh:
                    n = viz_subdivide_n
                    for k in range(n):
                        a = p0 + (p1 - p0) * (k / n)
                        b = p0 + (p1 - p0) * ((k + 1) / n)
                        new_path.append((tuple(a), tuple(b)))
                else:
                    new_path.append(s)
            filament_paths_viz.append(new_path)
        filament_paths = filament_paths_viz

    # Flatten the polylines into per-element segments; the physical
    # group of an element is its filament index.
    n_fil = len(filament_paths)
    segments = []
    groups = []
    for k, path in enumerate(filament_paths):
        for seg in path:
            segments.append(seg)
            groups.append(k)
    group_names = [f"{label}_{k}" for k in range(n_fil)]
    n_elems = len(segments)

    views = []
    if currents is not None:
        I_arr = np.asarray(currents)
        g_idx = np.asarray(groups, dtype=int)
        I_abs_el = np.abs(I_arr)[g_idx]
        views.append({"name": "|I| [A]", "ncomp": 1, "step": 0,
                      "time": 0.0, "values": I_abs_el})

        if direction_view:
            # unit tangent x Re(I): SIGNED arrows, so a filament that
            # carries reverse current points backwards (|I| cannot show
            # that -- see the export_peec_topology_msh note).
            I_signed_el = np.real(I_arr)[g_idx]
            tang = np.zeros((n_elems, 3), dtype=np.float64)
            for i, (p1, p2) in enumerate(segments):
                t = np.asarray(p2, dtype=np.float64) - np.asarray(
                    p1, dtype=np.float64)
                norm = float(np.linalg.norm(t))
                if norm >= 1e-30:
                    tang[i] = (t / norm) * float(I_signed_el[i])
            views.append({"name": "I direction [A]", "ncomp": 3, "step": 0,
                          "time": 0.0, "values": tang})

        if complex_steps:
            # two-step re/im view: feeds gmsh harmonic-to-time /
            # modulus-phase post (AC current animation)
            for step, part in ((0, I_arr.real), (1, I_arr.imag)):
                views.append({"name": "I_complex [A]", "ncomp": 1,
                              "step": step, "time": float(step),
                              "values": np.asarray(part)[g_idx]})

    info = _write_line_element_msh(output_msh, segments, groups,
                                   group_names, views)

    print(f"export_filaments_msh: {output_msh}")
    print(f"  {n_fil} filaments, {info['n_elements']} line elements, "
          f"{info['n_nodes']} nodes")
    if currents is not None:
        print(f"  |I| range: [{float(np.min(np.abs(currents))):.3e}, "
              f"{float(np.max(np.abs(currents))):.3e}] A")


def export_peec_topology_msh(topology, output_msh, *,
                             currents=None, Zs=None,
                             current_convention="amplitude",
                             groups=None, group_names=None,
                             label="peec",
                             direction_view=False, complex_steps=False):
    """Export a solved PEEC topology as a GMSH .msh v4.1 post file.

    This is the PEEC-solution counterpart of ``export_filaments_msh``:
    it takes the ``PEECBuilder.build_topology()`` dict directly (so the
    per-segment geometry is the one the solver actually used) plus the
    branch currents from
    ``PEECCircuitSolver.compute_branch_currents(freq, port_currents)``,
    and writes one line element per PEEC segment carrying the electrical
    quantities as ElementData views.

    A multi-filament segment (``nwinc``/``nhinc`` > 1) is expanded by the
    builder into offset sub-filaments, so each strand becomes its own
    line element: the skin/proximity current redistribution is directly
    visible.  Sub-filaments meet only at the topology NODE, which is a
    circuit connection and not a geometric one -- the small visual gap
    between strands at a node is the model, not an export defect.

    REVERSE CURRENTS -- ``|I|`` hides them.  In the inductance-limited
    (high-frequency) regime the interior strands of a multi-filament bar
    carry current in the OPPOSITE direction: measured on a 4 mm square
    bar at 100 MHz, the uniform partition gives 1 reverse strand at
    nwinc=nhinc=3 (Re I = -0.091 A of 1 A total), 4 at 4x4 and 8 at 5x5,
    with the high-frequency limit reproducing L^-1 1 exactly.  A ``|I|``
    colormap shows those strands as ordinary positive current.  Use the
    ``I direction`` arrows (scaled by the SIGNED Re(I), so a reverse
    strand points backwards) or step 0 of ``I_complex`` (= Re I, read
    with a diverging colormap) to see the reversal;
    ``n_reverse_segments`` in the returned dict counts them.

    Views written (all ElementData, one value per segment):
        ``|I| [A]``        branch current magnitude (needs ``currents``)
        ``|J| [A/m^2]``    |I| / (width * height) -- the current-DENSITY
                           map, which is what actually reveals skin and
                           proximity effect when sub-filaments have
                           unequal cross-sections
        ``P [W]``          per-segment time-average ohmic dissipation
        ``I direction [A]`` (opt-in) unit tangent x SIGNED Re(I) arrows
        ``I_complex [A]``  (opt-in) two steps Re/Im -> feeds the
                           radia-mcp ``gmsh_harmonic_to_time`` /
                           ``gmsh_modulus_phase`` post-processing

    Args:
        topology: dict from ``PEECBuilder.build_topology()``.  Requires
            ``segment_centers`` / ``segment_directions`` /
            ``segment_lengths``; ``segment_widths`` /
            ``segment_heights`` enable ``|J|`` and ``R`` enables ``P``.
        output_msh: output .msh path.
        currents: optional complex (n_seg,) branch currents.  Without
            them a geometry-only file is written (useful to inspect the
            discretization before solving).
        Zs: optional complex (n_seg,) surface impedance added to the
            loss as ``Re(Zs)`` (the SIBC/proximity part of the
            dissipation); ``None`` uses the DC resistance alone.
        current_convention: ``"amplitude"`` (default) treats the
            currents as phasor AMPLITUDES, so the time-average loss is
            ``0.5 |I|^2 R``; ``"rms"`` treats them as RMS values, giving
            ``|I|^2 R``.  The choice is echoed in the returned dict so
            the reported power is never ambiguous.
        groups: optional per-segment 0-based group index (e.g. one group
            per turn or per conductor).  Default: a single group.
        group_names: names for those groups (default ``label`` /
            ``label_0`` ...).
        label: base physical-group name.
        direction_view: also write the ``I direction`` arrow view.
        complex_steps: also write the two-step ``I_complex`` view.

    Returns:
        dict with the output path, element/node counts, the view names,
        the current convention and the |I| / |J| ranges and total loss.
    """
    import numpy as np

    required = ("segment_centers", "segment_directions", "segment_lengths")
    missing = [k for k in required if topology.get(k) is None]
    if missing:
        raise KeyError(
            "topology is missing per-segment geometry %s.  Rebuild it via "
            "PEECBuilder.build_topology() from a recent Radia (the keys "
            "are written by the C++ builder)." % missing)

    centers = np.asarray(topology["segment_centers"], dtype=np.float64)
    dirs = np.asarray(topology["segment_directions"], dtype=np.float64)
    lengths = np.asarray(topology["segment_lengths"], dtype=np.float64)
    n_seg = centers.shape[0]
    if dirs.shape != (n_seg, 3) or lengths.shape != (n_seg,):
        raise ValueError(
            f"inconsistent topology geometry: centers {centers.shape}, "
            f"directions {dirs.shape}, lengths {lengths.shape}")

    # A PEEC segment is a straight stick: endpoints are exactly
    # center -/+ direction * length / 2.
    half = dirs * (lengths[:, None] / 2.0)
    p0 = centers - half
    p1 = centers + half
    segments = [(tuple(p0[i]), tuple(p1[i])) for i in range(n_seg)]

    if groups is None:
        groups = [0] * n_seg
        names = list(group_names) if group_names else [label]
        if len(names) != 1:
            raise ValueError(
                "group_names must hold exactly one name when groups is "
                f"None, got {len(names)}")
    else:
        groups = [int(g) for g in groups]
        if len(groups) != n_seg:
            raise ValueError(
                f"groups length {len(groups)} != {n_seg} segments")
        n_groups = max(groups) + 1
        names = (list(group_names) if group_names
                 else [f"{label}_{k}" for k in range(n_groups)])
        if len(names) != n_groups:
            raise ValueError(
                f"group_names length {len(names)} != {n_groups} groups")

    conv = str(current_convention).lower()
    if conv not in ("amplitude", "rms"):
        raise ValueError("current_convention must be 'amplitude' or 'rms', "
                         f"got {current_convention!r}")
    p_factor = 0.5 if conv == "amplitude" else 1.0

    views = []
    summary = {"msh": output_msh, "n_segments": n_seg,
               "current_convention": conv}

    if currents is not None:
        I_arr = np.asarray(currents).reshape(-1)
        if I_arr.shape[0] != n_seg:
            raise ValueError(
                f"currents length {I_arr.shape[0]} != {n_seg} segments")
        I_abs = np.abs(I_arr)
        views.append({"name": "|I| [A]", "ncomp": 1, "step": 0,
                      "time": 0.0, "values": I_abs})
        summary["I_abs_range"] = [float(I_abs.min()), float(I_abs.max())]

        widths = topology.get("segment_widths")
        heights = topology.get("segment_heights")
        if widths is not None and heights is not None:
            area = (np.asarray(widths, dtype=np.float64)
                    * np.asarray(heights, dtype=np.float64))
            if np.any(area <= 0.0):
                raise ValueError(
                    "non-positive segment cross-section in topology "
                    "(widths * heights); cannot form the current density")
            J_abs = I_abs / area
            views.append({"name": "|J| [A/m^2]", "ncomp": 1, "step": 0,
                          "time": 0.0, "values": J_abs})
            summary["J_abs_range"] = [float(J_abs.min()),
                                      float(J_abs.max())]

        R = topology.get("R")
        if R is not None:
            R = np.asarray(R, dtype=np.float64)
            if R.ndim == 2:            # dense diagonal form
                R = np.diag(R)
            if R.shape[0] != n_seg:
                raise ValueError(
                    f"R length {R.shape[0]} != {n_seg} segments")
            R_eff = R.copy()
            if Zs is not None:
                Zs_arr = np.asarray(Zs).reshape(-1)
                if Zs_arr.shape[0] != n_seg:
                    raise ValueError(
                        f"Zs length {Zs_arr.shape[0]} != {n_seg} segments")
                R_eff = R_eff + np.real(Zs_arr)
            P = p_factor * (I_abs ** 2) * R_eff
            views.append({"name": "P [W]", "ncomp": 1, "step": 0,
                          "time": 0.0, "values": P})
            summary["P_total_W"] = float(P.sum())

        if direction_view:
            # SIGNED arrows: Re(I) * unit tangent, so a reverse-carrying
            # strand points backwards (see the reverse-current note).
            I_signed = np.real(I_arr)
            tang = np.zeros((n_seg, 3), dtype=np.float64)
            norms = np.linalg.norm(dirs, axis=1)
            ok = norms >= 1e-30
            tang[ok] = (dirs[ok] / norms[ok, None]) * I_signed[ok, None]
            views.append({"name": "I direction [A]", "ncomp": 3, "step": 0,
                          "time": 0.0, "values": tang})
            summary["n_reverse_segments"] = int((I_signed < 0.0).sum())

        if complex_steps:
            for step, part in ((0, I_arr.real), (1, I_arr.imag)):
                views.append({"name": "I_complex [A]", "ncomp": 1,
                              "step": step, "time": float(step),
                              "values": part})

    info = _write_line_element_msh(output_msh, segments, groups, names, views)
    summary["n_nodes"] = info["n_nodes"]
    summary["n_elements"] = info["n_elements"]
    summary["views"] = list(dict.fromkeys(v["name"] for v in views))

    print(f"export_peec_topology_msh: {output_msh}")
    print(f"  {n_seg} segments, {info['n_nodes']} nodes, "
          f"views: {summary['views']}")
    if "P_total_W" in summary:
        print(f"  P_total = {summary['P_total_W']:.6e} W "
              f"({conv} current convention)")
    return summary


_TRACK_SCHEMA_PREFIX = "radia-time-domain-lorentz-track/"
_JOULE_PER_EV = 1.602176634e-19


def _beam_animation_views(values, t_mid, *, n_frames, mode, comet_window,
                          view_name):
    """Multi-step "beam" view whose visible extent grows with time.

    The trick (MEASURED on gmsh 4.15.2, locked by tests): rendered with
    a FIXED colour range and ``SaturateValues = 0``, gmsh does NOT DRAW
    an element whose value falls outside that range.  So step ``k``
    carries the real quantity on the segments the beam has already
    reached and a far out-of-range SENTINEL everywhere else -- the
    trail grows while the colour scale stays put, unlike sweeping the
    range itself (which rescales the colours every frame).

    NaN is NOT usable as the sentinel: gmsh writes the literal ``nan``
    into the .pos/.msh and its own parser then rejects the file
    ("Unknown variable 'nan'").  The sentinel is finite by necessity.
    """
    import numpy as np

    lo = float(values.min())
    hi = float(values.max())
    # The sentinel must be far outside on the scale of the VALUES, not
    # of their spread: a monoenergetic beam coloured by energy has a
    # spread of ~1e-13 eV, and a sentinel one spread above the maximum
    # would be the same colour as the beam -- i.e. not hidden at all.
    # The floor is RELATIVE (1e-3 of the magnitude), never an absolute
    # 1.0: an absolute floor poisons every downstream stats verb on a
    # small-valued quantity.  MEASURED on the time-coloured beam
    # (lo=1.68e-10, hi=6.72e-6): absolute floor -> sentinel 10.0, i.e.
    # 1.5e6x the data (field_stats then reported max 10.0 for a view
    # whose physics tops out at 6.72e-6); relative floor -> sentinel
    # 7.39e-5, 11x the data.  Hiding is unaffected -- monoenergetic
    # 1e4 eV still gives hi=1e4+0.01 and sentinel=1.01e4, far outside
    # the colour range.
    span = hi - lo
    scale = max(span, 1e-3 * max(abs(lo), abs(hi)))
    if scale <= 0.0:
        scale = 1.0                     # degenerate all-zero data, last resort
    if span < 1e-12 * scale:
        hi = lo + 1e-3 * scale          # a colour bar needs a real span
    sentinel = hi + 10.0 * scale

    t0 = float(t_mid.min())
    t1 = float(t_mid.max())
    if t1 <= t0:
        raise ValueError(
            "animation needs a nonzero time span across the tracks")
    window = comet_window * (t1 - t0) if mode == "comet" else None

    views = []
    for k in range(n_frames):
        head = t0 + (t1 - t0) * (k + 1) / n_frames
        visible = t_mid <= head
        if window is not None:
            visible &= t_mid >= head - window
        views.append({"name": view_name, "ncomp": 1, "step": k,
                      "time": head,
                      "values": np.where(visible, values, sentinel)})
    return views, lo, hi, sentinel


def export_particle_tracks_msh(tracks, output_msh, *, label="track",
                               stride=1, animation_frames=0,
                               animation_color="energy",
                               animation_mode="trail",
                               comet_window=0.15):
    """Export time-domain Lorentz particle tracks as GMSH line elements.

    Bridges ``radia.particle_tracking.track_lorentz_ivp`` (and the
    momentum-scan helpers, via their ``tracks`` lists) to the GMSH
    post-processing lane: each track becomes one physical group of
    2-node line elements carrying per-segment ElementData views

        ``time [s]``             midpoint sample time,
        ``speed [m/s]``          midpoint |v|,
        ``kinetic energy [eV]``  midpoint kinetic energy,
        ``v [m/s]``              velocity arrows (unit tangent scaled
                                 by the midpoint speed).

    Color by ``time [s]`` to read the direction of travel, by
    ``kinetic energy [eV]`` to see electrostatic acceleration, and
    overlay the source geometry with
    ``gmsh_render(..., merge_files=["coil.step"])``.

    ``animation_frames > 0`` adds a multi-step ``beam`` view: the
    particles FLY along their trajectories, one time step per frame,
    ready for ``gmsh_export_animation``.  It only renders correctly
    with the fixed colour range returned in ``animation`` (see
    ``_beam_animation_views`` for the measured mechanism) -- the
    returned ``render_hint`` is exactly the ``color=`` argument to
    pass.

    Args:
        tracks: one track dict, or a sequence of them.  A track is the
            return value of ``track_lorentz_ivp`` (schema
            ``radia-time-domain-lorentz-track/v1``); for the
            five-momentum optics helpers pass ``result["tracks"]``.
        output_msh: output .msh v4.1 path.
        label: physical-group base name; track ``k`` becomes
            ``{label}_{k}``.
        stride: keep every ``stride``-th time sample (the final sample
            is always kept).  The DOP853 ``t_eval`` grid is often much
            denser than a rendering needs.
        animation_frames: number of animation steps (0 = no beam view).
        animation_color: quantity carried by the beam view --
            ``"energy"`` (default; separates a momentum scan into
            distinct constant colours), ``"speed"``, or ``"time"``.
        animation_mode: ``"trail"`` (default -- the track accumulates
            behind the particle) or ``"comet"`` (only a moving window
            around the particle stays lit).
        comet_window: comet window as a fraction of the total flight
            time (ignored in ``"trail"`` mode).

    Returns:
        dict with ``n_tracks`` / ``n_nodes`` / ``n_elements`` /
        ``views`` / per-track segment counts / the count of skipped
        zero-length segments (repeated positions never become
        degenerate line elements), plus ``animation`` +
        ``render_hint`` when a beam view was written.
    """
    import numpy as np

    if isinstance(tracks, dict):
        tracks = [tracks]
    try:
        tracks = list(tracks)
    except TypeError as exc:
        raise ValueError(
            "tracks must be a track dict or a sequence of them") from exc
    if not tracks:
        raise ValueError("tracks must contain at least one track")
    stride = int(stride)
    if stride < 1:
        raise ValueError(f"stride must be >= 1, got {stride}")
    if not str(label).strip():
        raise ValueError("label must not be empty")
    animation_frames = int(animation_frames)
    if animation_frames < 0:
        raise ValueError(
            f"animation_frames must be >= 0, got {animation_frames}")
    if animation_frames and animation_mode not in ("trail", "comet"):
        raise ValueError(
            f"animation_mode must be 'trail' or 'comet', got "
            f"{animation_mode!r}")
    if animation_frames and animation_color not in ("energy", "speed",
                                                    "time"):
        raise ValueError(
            f"animation_color must be 'energy', 'speed' or 'time', got "
            f"{animation_color!r}")
    comet_window = float(comet_window)
    if animation_frames and animation_mode == "comet" and not (
            0.0 < comet_window <= 1.0):
        raise ValueError(
            f"comet_window must be in (0, 1], got {comet_window}")

    segments = []
    groups = []
    names = []
    t_mid_all = []
    speed_mid_all = []
    ke_mid_ev_all = []
    tang_all = []
    per_track_segments = []
    n_degenerate = 0

    for k, track in enumerate(tracks):
        if not isinstance(track, dict):
            raise ValueError(f"track {k} is not a dict")
        schema = str(track.get("schema", ""))
        if not schema.startswith(_TRACK_SCHEMA_PREFIX):
            raise ValueError(
                f"track {k} schema {schema!r} is not a "
                f"{_TRACK_SCHEMA_PREFIX}* track (pass the "
                "track_lorentz_ivp result dict)")
        missing = [key for key in ("time_s", "position_m", "velocity_m_s",
                                   "kinetic_energy_j") if key not in track]
        if missing:
            raise ValueError(f"track {k} is missing keys {missing}")
        times = np.asarray(track["time_s"], dtype=np.float64).reshape(-1)
        pos = np.asarray(track["position_m"], dtype=np.float64)
        vel = np.asarray(track["velocity_m_s"], dtype=np.float64)
        ke_j = np.asarray(track["kinetic_energy_j"],
                          dtype=np.float64).reshape(-1)
        n = times.shape[0]
        if pos.shape != (n, 3) or vel.shape != (n, 3) or ke_j.shape != (n,):
            raise ValueError(
                f"track {k}: inconsistent sample counts "
                f"(time {n}, position {pos.shape}, velocity {vel.shape}, "
                f"kinetic energy {ke_j.shape})")
        keep = list(range(0, n, stride))
        if keep[-1] != n - 1:
            keep.append(n - 1)
        if len(keep) < 2:
            raise ValueError(
                f"track {k} has {n} samples; need at least 2 after "
                f"stride {stride}")
        times, pos, vel, ke_j = times[keep], pos[keep], vel[keep], ke_j[keep]

        speed = np.linalg.norm(vel, axis=1)
        chords = pos[1:] - pos[:-1]
        lengths = np.linalg.norm(chords, axis=1)
        ok = lengths >= 1e-30
        n_degenerate += int((~ok).sum())
        if not ok.any():
            raise ValueError(
                f"track {k}: every segment has zero length "
                "(the particle never moved)")
        names.append(f"{label}_{k}")
        idx = np.flatnonzero(ok)
        for i in idx:
            segments.append((tuple(pos[i]), tuple(pos[i + 1])))
            groups.append(k)
        t_mid_all.append(0.5 * (times[idx] + times[idx + 1]))
        speed_mid_all.append(0.5 * (speed[idx] + speed[idx + 1]))
        ke_mid_ev_all.append(
            0.5 * (ke_j[idx] + ke_j[idx + 1]) / _JOULE_PER_EV)
        tang_all.append(chords[idx] / lengths[idx, None]
                        * (0.5 * (speed[idx] + speed[idx + 1]))[:, None])
        per_track_segments.append(int(idx.size))

    views = [
        {"name": "time [s]", "ncomp": 1, "step": 0, "time": 0.0,
         "values": np.concatenate(t_mid_all)},
        {"name": "speed [m/s]", "ncomp": 1, "step": 0, "time": 0.0,
         "values": np.concatenate(speed_mid_all)},
        {"name": "kinetic energy [eV]", "ncomp": 1, "step": 0, "time": 0.0,
         "values": np.concatenate(ke_mid_ev_all)},
        {"name": "v [m/s]", "ncomp": 3, "step": 0, "time": 0.0,
         "values": np.concatenate(tang_all, axis=0)},
    ]

    animation = None
    if animation_frames:
        source = {"energy": views[2], "speed": views[1],
                  "time": views[0]}[animation_color]
        beam_views, lo, hi, sentinel = _beam_animation_views(
            source["values"], views[0]["values"],
            n_frames=animation_frames, mode=animation_mode,
            comet_window=comet_window,
            view_name=f"beam ({source['name']})")
        views.extend(beam_views)
        animation = {
            "view": beam_views[0]["name"],
            "n_steps": animation_frames,
            "mode": animation_mode,
            "color_range": [lo, hi],
            "sentinel": sentinel,
            "render_hint": {"color": {"range": [lo, hi],
                                      "saturate": False}},
        }

    info = _write_line_element_msh(output_msh, segments, groups, names,
                                   views)
    summary = {
        "n_tracks": len(tracks),
        "n_nodes": info["n_nodes"],
        "n_elements": info["n_elements"],
        "per_track_segments": per_track_segments,
        "n_degenerate_skipped": n_degenerate,
        "views": list(dict.fromkeys(v["name"] for v in views)),
    }
    if animation is not None:
        summary["animation"] = animation
        summary["render_hint"] = animation["render_hint"]
    print(f"export_particle_tracks_msh: {output_msh}")
    print(f"  {len(tracks)} tracks, {info['n_elements']} segments, "
          f"{info['n_nodes']} nodes, views: {summary['views']}")
    if animation is not None:
        print(f"  beam animation: {animation_frames} steps "
              f"({animation_mode}), render with color range "
              f"[{animation['color_range'][0]:.6g}, "
              f"{animation['color_range'][1]:.6g}] and saturate=False")
    return summary
