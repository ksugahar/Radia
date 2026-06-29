"""
peec_mesh_import.py

Surface Mesh to PEEC Panel/Filament Extraction

Converts a 2D surface mesh (triangles/quads) into PEEC topology:
  - Faces -> Panels (Star DOFs, capacitive P matrix)
  - Edges -> Filaments (Loop DOFs, inductive L matrix)
  - Nodes -> Junction nodes (for MNA)

Supports:
  - GMSH .msh files (direct import)
  - NGSolve/Netgen surface meshes
  - Manual node/face specification

Usage:
    from peec_mesh_import import surface_mesh_to_peec

    # From GMSH file
    builder = surface_mesh_to_peec('conductor.msh',
                                    sigma=5.8e7,
                                    thickness=1e-3)

    # Build and solve
    topo = builder.build_topology(include_star=True)
    solver = PEECCircuitSolver(topo)

Part of Radia project
"""

import numpy as np


class GMSHCenterlineReader:
    """Read 1D line elements from a GMSH v2.2 ASCII centerline mesh."""

    def __init__(self, filename, unit_scale=1.0):
        self.filename = str(filename)
        self.unit_scale = float(unit_scale)
        self.nodes = {}
        self.edge_elements = []
        self.physical_groups = {}

    def read(self):
        """Read nodes and two-node line elements from the mesh file."""
        self.nodes = {}
        self.edge_elements = []
        self.physical_groups = {}

        with open(self.filename, "r", encoding="utf-8") as f:
            lines = f.readlines()

        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line == "$PhysicalNames":
                i = self._read_physical_names(lines, i + 1)
            elif line == "$Nodes":
                i = self._read_nodes(lines, i + 1)
            elif line == "$Elements":
                i = self._read_elements(lines, i + 1)
            i += 1

        return self.nodes, self.edge_elements

    def _read_physical_names(self, lines, start_idx):
        num_names = int(lines[start_idx].strip())
        idx = start_idx + 1
        for _ in range(num_names):
            parts = lines[idx].strip().split(maxsplit=2)
            if len(parts) >= 3:
                group_id = int(parts[1])
                self.physical_groups[group_id] = parts[2].strip('"')
            idx += 1

        while idx < len(lines) and lines[idx].strip() != "$EndPhysicalNames":
            idx += 1
        return idx

    def _read_nodes(self, lines, start_idx):
        num_nodes = int(lines[start_idx].strip())
        idx = start_idx + 1
        for _ in range(num_nodes):
            parts = lines[idx].strip().split()
            node_id = int(parts[0])
            x, y, z = (float(parts[1]) * self.unit_scale,
                       float(parts[2]) * self.unit_scale,
                       float(parts[3]) * self.unit_scale)
            self.nodes[node_id] = [x, y, z]
            idx += 1

        while idx < len(lines) and lines[idx].strip() != "$EndNodes":
            idx += 1
        return idx

    def _read_elements(self, lines, start_idx):
        num_elements = int(lines[start_idx].strip())
        idx = start_idx + 1
        for _ in range(num_elements):
            parts = lines[idx].strip().split()
            element_type = int(parts[1])
            if element_type == 1:
                num_tags = int(parts[2])
                node_start_idx = 3 + num_tags
                node1 = int(parts[node_start_idx])
                node2 = int(parts[node_start_idx + 1])
                self.edge_elements.append([node1, node2])
            idx += 1

        while idx < len(lines) and lines[idx].strip() != "$EndElements":
            idx += 1
        return idx

    def get_segments(self):
        """Return centerline segments as ``[[start_xyz, end_xyz], ...]``."""
        segments = []
        for node1_id, node2_id in self.edge_elements:
            segments.append([self.nodes[node1_id], self.nodes[node2_id]])
        return segments

    def get_total_length(self):
        """Compute the total polyline length of all 1D elements."""
        total_length = 0.0
        for start, end in self.get_segments():
            total_length += np.linalg.norm(np.array(end) - np.array(start))
        return float(total_length)

    def summary(self):
        """Return a compact, JSON-friendly mesh summary."""
        return {
            "filename": self.filename,
            "nodes": len(self.nodes),
            "edge_elements": len(self.edge_elements),
            "total_length": self.get_total_length(),
            "physical_groups": dict(self.physical_groups),
        }

    def print_info(self):
        """Print a short centerline mesh summary."""
        info = self.summary()
        print("GMSH Centerline Mesh Info:")
        print(f"  Nodes: {info['nodes']}")
        print(f"  Edge elements: {info['edge_elements']}")
        print(f"  Total length: {info['total_length']:.6f} m")
        if self.physical_groups:
            print(f"  Physical groups: {self.physical_groups}")


def read_gmsh_centerline(filename, unit_scale=1.0):
    """Read a GMSH v2.2 ASCII centerline mesh and return ``(nodes, edges)``."""
    reader = GMSHCenterlineReader(filename, unit_scale=unit_scale)
    return reader.read()


def surface_mesh_to_peec(mesh_source, sigma=5.8e7, thickness=1e-3,
                         ports=None, unit_scale=1.0, include_panels=True):
    """Convert a surface mesh to PEEC topology with panels and filaments.

    Args:
        mesh_source: One of:
            - str/Path: GMSH .msh file path
            - dict with keys 'nodes' (Nx3), 'faces' (list of index lists)
            - NGSolve Mesh object (boundary faces extracted)
        sigma: Conductivity [S/m] (default: copper 5.8e7)
        thickness: Conductor thickness [m] for filament cross-section
        ports: List of (node_id_pos, node_id_neg) tuples for port definition.
               If None, tries to auto-detect from boundary edges.
        unit_scale: Multiply coordinates by this factor (e.g., 1e-3 for mm->m)
        include_panels: If True, add mesh faces as panels for capacitive effects

    Returns:
        PEECBuilder: Configured builder ready for build_topology()
    """
    from radia.peec_matrices import PEECBuilder

    # Parse mesh source
    if isinstance(mesh_source, (str,)):
        nodes, faces = _load_gmsh(mesh_source, unit_scale)
    elif hasattr(mesh_source, '__fspath__') or hasattr(mesh_source, 'stem'):
        nodes, faces = _load_gmsh(str(mesh_source), unit_scale)
    elif isinstance(mesh_source, dict):
        nodes = np.array(mesh_source['nodes']) * unit_scale
        faces = mesh_source['faces']
    else:
        # Try NGSolve mesh
        nodes, faces = _extract_ngsolve_surface(mesh_source, unit_scale)

    n_nodes = len(nodes)
    n_faces = len(faces)

    # Extract unique edges from faces
    edges, edge_to_faces = _extract_edges(faces)
    n_edges = len(edges)

    # Find boundary edges (adjacent to only 1 face)
    boundary_edges = [e for e, adj in edge_to_faces.items() if len(adj) == 1]

    # Build PEEC model
    builder = PEECBuilder()

    # Add nodes
    node_id_map = {}
    for i, pos in enumerate(nodes):
        nid = builder.add_node_at(float(pos[0]), float(pos[1]), float(pos[2]))
        node_id_map[i] = nid

    # Add filaments (edges -> connected segments)
    for n0, n1 in edges:
        p0 = nodes[n0]
        p1 = nodes[n1]
        edge_len = np.linalg.norm(p1 - p0)
        # Filament width from dual mesh: approximate as edge_len / sqrt(3)
        # For thin conductors, use thickness as both w and h
        w = thickness
        h = thickness
        builder.add_connected_segment(
            node_id_map[n0], node_id_map[n1],
            w, h, sigma
        )

    # Add panels (faces -> triangles/quads)
    if include_panels:
        for face in faces:
            verts = [nodes[idx].tolist() for idx in face]
            builder.add_panel(verts)

    # Add ports
    if ports is not None:
        for port_def in ports:
            if len(port_def) == 2:
                n_pos, n_neg = port_def
                builder.add_port(node_id_map[n_pos], node_id_map[n_neg])
    elif len(boundary_edges) >= 2:
        # Auto-detect: use first and last boundary node as port
        boundary_nodes = set()
        for e in boundary_edges:
            boundary_nodes.add(e[0])
            boundary_nodes.add(e[1])
        boundary_list = sorted(boundary_nodes)
        if len(boundary_list) >= 2:
            builder.add_port(
                node_id_map[boundary_list[0]],
                node_id_map[boundary_list[-1]]
            )

    return builder


def _load_gmsh(filename, unit_scale=1.0):
    """Load surface mesh from GMSH .msh file.

    Returns:
        nodes: (N, 3) array of vertex positions
        faces: list of tuples (0-indexed vertex indices per face)
    """
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(filename)

    # Get nodes
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    coords = node_coords.reshape(-1, 3) * unit_scale

    # Build tag -> 0-index mapping
    tag_to_idx = {}
    for i, tag in enumerate(node_tags):
        tag_to_idx[int(tag)] = i

    # Get 2D elements (triangles and quads)
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements()

    faces = []
    for i, elem_type in enumerate(elem_types):
        if elem_type == 2:  # Triangle (3 nodes)
            ntags = elem_node_tags[i]
            n_elems = len(ntags) // 3
            for j in range(n_elems):
                tri = ntags[j * 3:(j + 1) * 3]
                faces.append(tuple(tag_to_idx[int(t)] for t in tri))
        elif elem_type == 3:  # Quadrilateral (4 nodes)
            ntags = elem_node_tags[i]
            n_elems = len(ntags) // 4
            for j in range(n_elems):
                quad = ntags[j * 4:(j + 1) * 4]
                faces.append(tuple(tag_to_idx[int(t)] for t in quad))

    gmsh.finalize()
    return coords, faces


def _extract_ngsolve_surface(mesh, unit_scale=1.0):
    """Extract surface faces from NGSolve mesh.

    Returns:
        nodes: (N, 3) array of vertex positions
        faces: list of tuples (0-indexed vertex indices per face)
    """
    # Get vertex coordinates
    n_verts = mesh.nv
    nodes = np.zeros((n_verts, 3))
    for i, v in enumerate(mesh.vertices):
        pt = v.point
        nodes[i] = [pt[0] * unit_scale, pt[1] * unit_scale, pt[2] * unit_scale]

    # Extract boundary faces (surface elements)
    faces = []
    for el in mesh.Elements(mesh.dim - 1):  # Boundary elements
        face_verts = tuple(v.nr for v in el.vertices)
        faces.append(face_verts)

    return nodes, faces


def _extract_edges(faces):
    """Extract unique edges from face list and build edge-face adjacency.

    Args:
        faces: list of tuples (vertex indices per face)

    Returns:
        edges: list of (n0, n1) tuples (sorted, unique)
        edge_to_faces: dict mapping edge tuple -> list of face indices
    """
    edge_to_faces = {}

    for face_idx, face in enumerate(faces):
        n_verts = len(face)
        for k in range(n_verts):
            n0 = face[k]
            n1 = face[(k + 1) % n_verts]
            edge_key = (min(n0, n1), max(n0, n1))
            if edge_key not in edge_to_faces:
                edge_to_faces[edge_key] = []
            edge_to_faces[edge_key].append(face_idx)

    edges = list(edge_to_faces.keys())
    return edges, edge_to_faces


def create_plate_mesh(center, size, divisions, sigma=5.8e7, thickness=1e-3,
                      ports=None, include_panels=True):
    """Create a rectangular plate surface mesh and convert to PEEC.

    Generates structured triangular mesh for a flat plate conductor.

    Args:
        center: [x, y, z] plate center position [m]
        size: [width, height] plate dimensions [m]
        divisions: [nx, ny] mesh subdivisions
        sigma: Conductivity [S/m]
        thickness: Conductor thickness [m]
        ports: List of (node_pos, node_neg) tuples
        include_panels: If True, add panels for capacitive effects

    Returns:
        PEECBuilder: Configured builder
    """
    cx, cy, cz = center
    sx, sy = size
    nx, ny = divisions

    # Generate structured grid nodes
    x = np.linspace(cx - sx / 2, cx + sx / 2, nx + 1)
    y = np.linspace(cy - sy / 2, cy + sy / 2, ny + 1)

    nodes = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            nodes.append([x[i], y[j], cz])
    nodes = np.array(nodes)

    # Generate triangular faces (2 triangles per quad cell)
    faces = []
    for j in range(ny):
        for i in range(nx):
            n00 = j * (nx + 1) + i
            n10 = j * (nx + 1) + (i + 1)
            n01 = (j + 1) * (nx + 1) + i
            n11 = (j + 1) * (nx + 1) + (i + 1)
            faces.append((n00, n10, n11))
            faces.append((n00, n11, n01))

    mesh_data = {'nodes': nodes, 'faces': faces}
    return surface_mesh_to_peec(mesh_data, sigma=sigma, thickness=thickness,
                                ports=ports, include_panels=include_panels)


def create_box_surface_mesh(center, size, divisions, sigma=5.8e7,
                            thickness=1e-3, ports=None, include_panels=True):
    """Create a box (6-face) surface mesh and convert to PEEC.

    Generates structured triangular surface mesh on all 6 faces of a box.
    Shared edges/nodes between faces are merged.

    Args:
        center: [x, y, z] box center [m]
        size: [sx, sy, sz] box dimensions [m]
        divisions: [nx, ny, nz] mesh subdivisions per face
        sigma: Conductivity [S/m]
        thickness: Conductor thickness [m]
        ports: List of (node_pos, node_neg) tuples
        include_panels: If True, add panels

    Returns:
        PEECBuilder: Configured builder
    """
    cx, cy, cz = center
    sx, sy, sz = size
    nx, ny, nz = divisions

    # Collect all raw nodes and faces, then merge duplicates
    raw_nodes = []
    raw_faces = []

    def add_face_mesh(p00, p10, p01, nu, nv):
        p00 = np.array(p00)
        p10 = np.array(p10)
        p01 = np.array(p01)
        u_vec = p10 - p00
        v_vec = p01 - p00

        offset = len(raw_nodes)
        for j in range(nv + 1):
            for i in range(nu + 1):
                pt = p00 + u_vec * (i / nu) + v_vec * (j / nv)
                raw_nodes.append(pt.tolist())

        for j in range(nv):
            for i in range(nu):
                n00 = j * (nu + 1) + i + offset
                n10 = j * (nu + 1) + (i + 1) + offset
                n01 = (j + 1) * (nu + 1) + i + offset
                n11 = (j + 1) * (nu + 1) + (i + 1) + offset
                raw_faces.append((n00, n10, n11))
                raw_faces.append((n00, n11, n01))

    x0, x1 = cx - sx / 2, cx + sx / 2
    y0, y1 = cy - sy / 2, cy + sy / 2
    z0, z1 = cz - sz / 2, cz + sz / 2

    add_face_mesh([x0, y0, z0], [x1, y0, z0], [x0, y1, z0], nx, ny)
    add_face_mesh([x0, y0, z1], [x1, y0, z1], [x0, y1, z1], nx, ny)
    add_face_mesh([x0, y0, z0], [x1, y0, z0], [x0, y0, z1], nx, nz)
    add_face_mesh([x0, y1, z0], [x1, y1, z0], [x0, y1, z1], nx, nz)
    add_face_mesh([x0, y0, z0], [x0, y1, z0], [x0, y0, z1], ny, nz)
    add_face_mesh([x1, y0, z0], [x1, y1, z0], [x1, y0, z1], ny, nz)

    # Merge duplicate nodes (shared edges between faces)
    raw_nodes_arr = np.array(raw_nodes)
    unique_nodes, remap = _merge_coincident_nodes(raw_nodes_arr, tol=1e-12)

    # Remap face indices
    merged_faces = []
    for face in raw_faces:
        new_face = tuple(remap[idx] for idx in face)
        # Skip degenerate faces (all same node)
        if len(set(new_face)) >= 3:
            merged_faces.append(new_face)

    mesh_data = {
        'nodes': unique_nodes,
        'faces': merged_faces,
    }
    return surface_mesh_to_peec(mesh_data, sigma=sigma, thickness=thickness,
                                ports=ports, include_panels=include_panels)


def _merge_coincident_nodes(nodes, tol=1e-12):
    """Merge nodes that are within tolerance distance.

    Args:
        nodes: (N, 3) array of node positions
        tol: Distance tolerance for merging

    Returns:
        unique_nodes: (M, 3) array of unique positions
        remap: (N,) array mapping old index -> new index
    """
    n = len(nodes)
    remap = np.arange(n, dtype=int)
    unique_mask = np.ones(n, dtype=bool)

    for i in range(n):
        if not unique_mask[i]:
            continue
        for j in range(i + 1, n):
            if not unique_mask[j]:
                continue
            if np.linalg.norm(nodes[i] - nodes[j]) < tol:
                unique_mask[j] = False
                remap[j] = remap[i]

    # Compact: renumber unique nodes
    unique_indices = np.where(unique_mask)[0]
    unique_nodes = nodes[unique_indices]
    new_idx = np.zeros(n, dtype=int)
    for new_i, old_i in enumerate(unique_indices):
        new_idx[old_i] = new_i

    # Update remap to use new indices
    final_remap = np.array([new_idx[remap[i]] for i in range(n)])
    return unique_nodes, final_remap
