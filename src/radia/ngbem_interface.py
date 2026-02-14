"""
ngbem_interface.py

Bridge between ngbem BEM matrices and Radia PEEC circuit extraction.

Converts NGBEMPEECSolver matrices (L, P, M_LS from HDivSurface/SurfaceL2)
into topology_dict format consumed by PEECCircuitSolver.

Architecture:
    NGBEMPEECSolver.assemble() -> L, P, M_LS, R
                |
        NGBEMBridge.to_topology_dict()
                |
        PEECCircuitSolver(topology_dict)
                |
            Z(f), SPICE, PRIMA MOR

The key mapping:
  - HDivSurface DOF i  <->  mesh edge i  <->  virtual segment i
  - SurfaceL2 DOF j    <->  mesh triangle j  <->  virtual star node j

Usage:
    from ngbem_peec import NGBEMPEECSolver, create_plate_mesh
    from ngbem_interface import NGBEMBridge

    mesh = create_plate_mesh(0.01, 0.001, 0.003)
    solver = NGBEMPEECSolver(mesh, sigma=5.8e7, thickness=1e-3)
    solver.assemble()

    bridge = NGBEMBridge(solver, port_spec=([0,0,0], [0.1,0,0]))
    topo = bridge.to_topology_dict()

    from peec_topology import PEECCircuitSolver
    circuit = PEECCircuitSolver(topo)
    Z = circuit.compute_port_impedance(1e6)

Part of Radia project
"""

import numpy as np


def extract_edge_geometry(mesh):
    """Extract edge midpoints, directions, lengths from NGSolve mesh.

    For order=0 HDivSurface, each DOF corresponds to one mesh edge.
    This function extracts the geometric data for each edge.

    Args:
        mesh: NGSolve Mesh (surface mesh)

    Returns:
        dict with:
            'centers': (n_edges, 3) edge midpoints
            'directions': (n_edges, 3) unit direction vectors
            'lengths': (n_edges,) edge lengths
            'p1': (n_edges, 3) first vertex position
            'p2': (n_edges, 3) second vertex position
            'edge_vertices': (n_edges, 2) vertex indices per edge
            'vertex_positions': dict {vertex_nr: [x, y, z]}
    """
    # Collect all vertex positions
    vertex_positions = {}
    for v in mesh.vertices:
        pt = v.point
        vertex_positions[v.nr] = np.array([pt[0], pt[1], pt[2]])

    # Extract edge data
    centers_list = []
    directions_list = []
    lengths_list = []
    p1_list = []
    p2_list = []
    edge_verts_list = []

    for edge in mesh.edges:
        verts = list(edge.vertices)
        v0_nr = verts[0].nr
        v1_nr = verts[1].nr
        p0 = vertex_positions[v0_nr]
        p1 = vertex_positions[v1_nr]

        center = 0.5 * (p0 + p1)
        diff = p1 - p0
        length = np.linalg.norm(diff)
        direction = diff / length if length > 1e-20 else np.zeros(3)

        centers_list.append(center)
        directions_list.append(direction)
        lengths_list.append(length)
        p1_list.append(p0)
        p2_list.append(p1)
        edge_verts_list.append([v0_nr, v1_nr])

    return {
        'centers': np.array(centers_list),
        'directions': np.array(directions_list),
        'lengths': np.array(lengths_list),
        'p1': np.array(p1_list),
        'p2': np.array(p2_list),
        'edge_vertices': np.array(edge_verts_list, dtype=int),
        'vertex_positions': vertex_positions,
    }


class NGBEMBridge:
    """Bridge between ngbem BEM matrices and PEEC circuit solver.

    Converts ngbem HDivSurface/SurfaceL2 matrices into topology_dict
    compatible with PEECCircuitSolver.

    The key mapping:
      - HDivSurface DOF i <-> mesh edge i <-> virtual segment i
      - SurfaceL2 DOF j <-> mesh triangle j <-> virtual star node j
      - Port: specified by user via coordinates or boundary labels
    """

    def __init__(self, ngbem_solver, port_spec=None):
        """
        Initialize bridge.

        Args:
            ngbem_solver: NGBEMPEECSolver instance (must be assembled)
            port_spec: Port specification. Options:
                - None or 'auto': auto-detect endpoints (min/max x-coord)
                - (pos_coords, neg_coords): tuple of two [x,y,z] arrays,
                  finds nearest mesh vertices
                - {'positive_label': str, 'negative_label': str}:
                  NGSolve boundary labels
        """
        if ngbem_solver.L is None:
            raise RuntimeError(
                "NGBEMPEECSolver must be assembled before creating bridge. "
                "Call solver.assemble() first.")

        self.solver = ngbem_solver
        self.mesh = ngbem_solver.mesh
        self.port_spec = port_spec

        # Extract geometry on construction
        self._edge_geom = extract_edge_geometry(self.mesh)
        self._virtual_topo = None
        self._ports = None

    def _build_virtual_topology(self):
        """Build virtual node-segment topology from mesh edges.

        Creates virtual nodes at mesh vertices and virtual segments at
        mesh edges. Each HDivSurface DOF becomes one virtual segment.

        Returns:
            dict with 'segment_nodes', 'n_nodes', 'node_positions'
        """
        edge_verts = self._edge_geom['edge_vertices']
        vert_pos = self._edge_geom['vertex_positions']

        # Get unique vertex IDs and create sequential mapping
        unique_verts = sorted(vert_pos.keys())
        vert_to_node = {v: i for i, v in enumerate(unique_verts)}
        n_nodes = len(unique_verts)

        # Build node positions array
        node_positions = np.zeros((n_nodes, 3))
        for v_nr, node_id in vert_to_node.items():
            node_positions[node_id] = vert_pos[v_nr]

        # Build segment_nodes array
        n_edges = len(edge_verts)
        segment_nodes = np.zeros((n_edges, 2), dtype=int)
        for i in range(n_edges):
            segment_nodes[i, 0] = vert_to_node[edge_verts[i, 0]]
            segment_nodes[i, 1] = vert_to_node[edge_verts[i, 1]]

        self._vert_to_node = vert_to_node
        self._virtual_topo = {
            'segment_nodes': segment_nodes,
            'n_nodes': n_nodes,
            'node_positions': node_positions,
        }
        return self._virtual_topo

    def _assign_ports(self):
        """Assign ports to the virtual topology.

        Returns:
            list of (node_pos, node_neg, port_id) tuples
        """
        if self._virtual_topo is None:
            self._build_virtual_topology()

        node_pos_array = self._virtual_topo['node_positions']
        spec = self.port_spec

        if spec is None or spec == 'auto':
            # Auto-detect: find vertices with min and max x-coordinate
            x_coords = node_pos_array[:, 0]
            node_neg = int(np.argmin(x_coords))
            node_pos = int(np.argmax(x_coords))

            # If x-coordinates are the same, try y, then z
            if abs(x_coords[node_pos] - x_coords[node_neg]) < 1e-12:
                y_coords = node_pos_array[:, 1]
                node_neg = int(np.argmin(y_coords))
                node_pos = int(np.argmax(y_coords))

            if node_pos == node_neg:
                z_coords = node_pos_array[:, 2]
                node_neg = int(np.argmin(z_coords))
                node_pos = int(np.argmax(z_coords))

        elif isinstance(spec, (tuple, list)) and len(spec) == 2:
            # Coordinate-based: find nearest vertices
            pos_target = np.asarray(spec[0], dtype=float)
            neg_target = np.asarray(spec[1], dtype=float)

            dist_pos = np.linalg.norm(node_pos_array - pos_target, axis=1)
            dist_neg = np.linalg.norm(node_pos_array - neg_target, axis=1)

            node_pos = int(np.argmin(dist_pos))
            node_neg = int(np.argmin(dist_neg))

        elif isinstance(spec, dict):
            # Boundary label-based
            pos_label = spec.get('positive_label', '')
            neg_label = spec.get('negative_label', '')
            raise NotImplementedError(
                "Boundary label-based port assignment not yet implemented. "
                "Use coordinate-based port_spec=([x1,y1,z1], [x2,y2,z2]).")
        else:
            raise ValueError(f"Invalid port_spec: {spec}")

        if node_pos == node_neg:
            raise ValueError(
                f"Port positive and negative nodes are the same ({node_pos}). "
                "Check port_spec or mesh geometry.")

        self._ports = [(node_pos, node_neg, 0)]
        return self._ports

    def to_topology_dict(self):
        """Convert ngbem matrices to topology_dict format.

        This is the main entry point. Produces a dict compatible
        with PEECCircuitSolver.__init__().

        Returns:
            topology_dict: dict with L, R, segment_nodes, n_nodes,
                          n_loop, ports, and optional P, M_LS, n_star
        """
        mats = self.solver.get_matrices()

        if self._virtual_topo is None:
            self._build_virtual_topology()
        if self._ports is None:
            self._assign_ports()

        topo = {
            # Required fields
            'L': mats['L'],
            'R': mats['R_loop'],
            'segment_nodes': self._virtual_topo['segment_nodes'],
            'n_nodes': self._virtual_topo['n_nodes'],
            'n_loop': mats['n_loop'],
            'ports': self._ports,

            # Optional: geometry for coupling
            'segment_centers': self._edge_geom['centers'],
            'segment_directions': self._edge_geom['directions'],
            'segment_lengths': self._edge_geom['lengths'],
            'node_positions': self._virtual_topo['node_positions'],

            # Optional: full Loop-Star
            'P': mats['P'],
            'M_LS': mats['M_LS'],
            'n_star': mats['n_star'],

            # Metadata
            'backend': 'ngbem',
            'order': self.solver.order,
        }

        # Build incidence data for MNA solver
        topo.update(self._build_incidence_data())

        return topo

    def _build_incidence_data(self):
        """Build CSR incidence matrix data for MNA solver.

        The MNA solver in PEECCircuitSolver uses the full node incidence
        matrix A_full where A[node, seg] = +1/-1 for outgoing/incoming.

        Returns:
            dict with incidence_data, incidence_indices, incidence_indptr
        """
        seg_nodes = self._virtual_topo['segment_nodes']
        n_nodes = self._virtual_topo['n_nodes']
        n_seg = len(seg_nodes)

        # Build dense incidence first, then convert
        # A_full[node, seg] = +1 if seg leaves node, -1 if seg enters node
        from scipy import sparse

        rows = []
        cols = []
        data = []
        for seg_idx in range(n_seg):
            node_from = seg_nodes[seg_idx, 0]
            node_to = seg_nodes[seg_idx, 1]
            rows.append(node_from)
            cols.append(seg_idx)
            data.append(1.0)
            rows.append(node_to)
            cols.append(seg_idx)
            data.append(-1.0)

        A_sparse = sparse.csr_matrix(
            (data, (rows, cols)), shape=(n_nodes, n_seg))

        return {
            'incidence_data': np.array(A_sparse.data),
            'incidence_indices': np.array(A_sparse.indices),
            'incidence_indptr': np.array(A_sparse.indptr),
            'n_junction': n_nodes,
        }

    def to_circuit_solver(self):
        """Create PEECCircuitSolver from ngbem matrices.

        Convenience method that calls to_topology_dict() and creates
        a PEECCircuitSolver instance.

        Returns:
            PEECCircuitSolver instance
        """
        from peec_topology import PEECCircuitSolver
        topo = self.to_topology_dict()
        return PEECCircuitSolver(topo)

    def print_summary(self):
        """Print summary of the bridge mapping."""
        if self._virtual_topo is None:
            self._build_virtual_topology()
        if self._ports is None:
            self._assign_ports()

        n_nodes = self._virtual_topo['n_nodes']
        n_loop = self.solver.n_loop
        n_star = self.solver.n_star

        print("NGBEMBridge Summary:")
        print(f"  Backend: ngbem (order={self.solver.order})")
        print(f"  Mesh vertices (virtual nodes): {n_nodes}")
        print(f"  Mesh edges (virtual segments / Loop DOFs): {n_loop}")
        print(f"  Mesh triangles (Star DOFs): {n_star}")
        print(f"  Edge lengths: [{np.min(self._edge_geom['lengths']):.4e}, "
              f"{np.max(self._edge_geom['lengths']):.4e}] m")

        for port_pos, port_neg, port_id in self._ports:
            pos_coord = self._virtual_topo['node_positions'][port_pos]
            neg_coord = self._virtual_topo['node_positions'][port_neg]
            print(f"  Port {port_id}: node {port_pos} "
                  f"({pos_coord[0]:.4f}, {pos_coord[1]:.4f}, "
                  f"{pos_coord[2]:.4f})")
            print(f"           -> node {port_neg} "
                  f"({neg_coord[0]:.4f}, {neg_coord[1]:.4f}, "
                  f"{neg_coord[2]:.4f})")
