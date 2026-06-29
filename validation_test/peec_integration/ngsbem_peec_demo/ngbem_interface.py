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

Loop-Star Decomposition:
    Decomposes HDivSurface (edge-based) DOFs into:
    - Loop (face-based): solenoidal currents, div(J_loop) = 0  -> eddy currents
    - Star (vertex-based): irrotational currents, curl(J_star) = 0  -> charge

    T_loop: (n_edges, n_loop) face circulation basis
    T_star: (n_edges, n_star) vertex star basis
    T = [T_loop | T_star]: full transformation (square, invertible)

    For BEM matrix A_edge:
        A_LS = T^T * A_edge * T  (congruence transform, preserves symmetry)

    Properties:
        T_star^T * T_loop = 0  (orthogonality)
        M_LS * T_loop = 0      (loops are divergence-free)

Usage:
    from ngbem_peec import NGBEMPEECSolver, create_plate_mesh
    from ngbem_interface import NGBEMBridge, LoopStarTransform

    mesh = create_plate_mesh(0.01, 0.001, 0.003)
    solver = NGBEMPEECSolver(mesh, sigma=5.8e7, thickness=1e-3)
    solver.assemble()

    # Option 1: Direct bridge to PEEC
    bridge = NGBEMBridge(solver, port_spec=([0,0,0], [0.1,0,0]))
    topo = bridge.to_topology_dict()

    # Option 2: Loop-Star decomposition
    ls = LoopStarTransform(mesh)
    L_LS = ls.transform_matrix(solver.L)
    L_LL, L_LS_block, L_SL, L_SS = ls.extract_blocks(L_LS)

    # Decompose BEM solution into eddy current + charge components
    J_loop, J_star, alpha_loop, alpha_star = ls.decompose_current(J_edge)

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


class LoopStarTransform:
    """Loop-Star basis transformation for surface BEM mesh.

    Decomposes HDivSurface (edge-based) DOFs into:
    - Loop (face-based): solenoidal currents, div(J_loop) = 0
    - Star (vertex-based): irrotational currents, curl(J_star) = 0

    Physical interpretation:
    - Loop = eddy currents (no charge accumulation at vertices)
    - Star = charge-driven currents (charge accumulates at vertices)

    Mathematical basis:
    - T_loop columns = face circulation vectors (incidence matrix B^T)
    - T_star columns = vertex star vectors (node-edge incidence A^T)
    - Orthogonality: T_star^T * T_loop = 0 (exact, from graph theory)
    - Euler formula: n_loop + n_star = n_edges

    Surface type handling:
    - Open surface (chi=1): remove 1 reference vertex from T_star
    - Closed surface (chi=2): remove 1 vertex + 1 face

    Supported element orders:
    - order=0: Full Loop-Star (1 DOF per edge, exact decomposition)
    - order>0: Not yet supported (requires extended Helmholtz decomposition)

    Usage:
        ls = LoopStarTransform(mesh)

        # Transform BEM matrix to Loop-Star basis
        L_LS = ls.transform_matrix(L_edge)
        L_LL, L_LS_block, L_SL, L_SS = ls.extract_blocks(L_LS)

        # Decompose BEM solution
        J_loop, J_star, alpha_loop, alpha_star = ls.decompose_current(J_edge)
    """

    def __init__(self, mesh):
        """Build Loop-Star transformation from NGSolve surface mesh.

        Args:
            mesh: NGSolve Mesh (surface mesh with triangular elements)

        Raises:
            ValueError: If mesh has no surface elements
        """
        self.mesh = mesh
        self._build()

    def _build(self):
        """Build T_loop and T_star matrices from mesh topology."""
        try:
            from ngsolve import BND
        except ImportError:
            raise ImportError(
                "LoopStarTransform requires NGSolve. "
                "Install with: pip install ngsolve>=6.2.2601")

        mesh = self.mesh

        # --- Collect vertices ---
        vertex_nrs = []
        vertex_positions = {}
        for v in mesh.vertices:
            vertex_nrs.append(v.nr)
            pt = v.point
            vertex_positions[v.nr] = np.array([pt[0], pt[1], pt[2]])

        vertex_nrs = sorted(vertex_nrs)
        vert_to_idx = {v: i for i, v in enumerate(vertex_nrs)}
        n_vertices = len(vertex_nrs)

        # --- Collect edges with global direction ---
        edge_list = []        # [(v0_nr, v1_nr)] in mesh.edges iteration order
        edge_to_idx = {}     # (min_v, max_v) -> sequential index

        for idx, edge in enumerate(mesh.edges):
            verts = list(edge.vertices)
            v0, v1 = verts[0].nr, verts[1].nr
            edge_list.append((v0, v1))
            key = (min(v0, v1), max(v0, v1))
            edge_to_idx[key] = idx

        n_edges = len(edge_list)

        # --- Collect faces (surface elements) ---
        face_vertex_list = []
        for el in mesh.Elements(BND):
            face_verts = [v.nr for v in el.vertices]
            face_vertex_list.append(face_verts)

        n_faces = len(face_vertex_list)

        if n_faces == 0:
            raise ValueError(
                "No surface elements found. "
                "LoopStarTransform requires a surface mesh.")

        # --- Euler characteristic ---
        chi = n_vertices - n_edges + n_faces

        # --- Build T_loop: face-edge incidence (n_edges x n_faces) ---
        # T_loop[e, f] = +1 if edge e's global direction matches face f's
        #                     local circulation (va -> vb), -1 otherwise
        T_loop_full = np.zeros((n_edges, n_faces))

        for f_idx, face_verts in enumerate(face_vertex_list):
            n_v = len(face_verts)
            for k in range(n_v):
                va = face_verts[k]
                vb = face_verts[(k + 1) % n_v]
                key = (min(va, vb), max(va, vb))

                if key not in edge_to_idx:
                    continue

                e_idx = edge_to_idx[key]
                global_v0, global_v1 = edge_list[e_idx]

                # +1 if global direction matches face local (va -> vb)
                if global_v0 == va:
                    sign = +1.0
                else:
                    sign = -1.0

                T_loop_full[e_idx, f_idx] = sign

        # --- Build T_star: node-edge incidence (n_edges x n_vertices) ---
        # T_star[e, v] = +1 if edge e leaves vertex v (v is v0)
        #               -1 if edge e enters vertex v (v is v1)
        T_star_full = np.zeros((n_edges, n_vertices))

        for e_idx, (v0, v1) in enumerate(edge_list):
            v0_idx = vert_to_idx[v0]
            v1_idx = vert_to_idx[v1]
            T_star_full[e_idx, v0_idx] = +1.0
            T_star_full[e_idx, v1_idx] = -1.0

        # --- Remove dependent columns to make T square ---
        # Euler: V - E + F = chi
        # For open surface (chi=1): remove 1 star column (reference vertex)
        # For closed surface (chi=2): remove 1 star + 1 loop column
        T_star = T_star_full
        T_loop = T_loop_full

        if chi >= 1:
            T_star = T_star_full[:, :-1]
        if chi >= 2:
            T_loop = T_loop_full[:, :-1]

        # --- Store results ---
        self.T_loop = T_loop          # (n_edges, n_loop)
        self.T_star = T_star          # (n_edges, n_star)
        self.T = np.hstack([T_loop, T_star])  # (n_edges, n_edges)

        self._n_loop = T_loop.shape[1]
        self._n_star = T_star.shape[1]
        self._n_edges = n_edges
        self._n_faces = n_faces
        self._n_vertices = n_vertices
        self._chi = chi

        self._vertex_positions = vertex_positions
        self._vert_to_idx = vert_to_idx
        self._edge_list = edge_list
        self._face_vertex_list = face_vertex_list

        # Verify T is square
        expected_cols = n_edges
        actual_cols = self._n_loop + self._n_star
        if actual_cols != expected_cols:
            raise RuntimeError(
                f"Loop-Star dimension mismatch: "
                f"n_loop({self._n_loop}) + n_star({self._n_star}) = "
                f"{actual_cols} != n_edges({expected_cols}). "
                f"Euler chi={chi}")

        # Compute inverse
        try:
            cond = np.linalg.cond(self.T)
            if cond > 1e14:
                self.T_inv = np.linalg.pinv(self.T)
                self._invertible = False
            else:
                self.T_inv = np.linalg.inv(self.T)
                self._invertible = True
        except np.linalg.LinAlgError:
            self.T_inv = np.linalg.pinv(self.T)
            self._invertible = False

    @property
    def n_loop(self):
        """Number of Loop (face-based, solenoidal) DOFs."""
        return self._n_loop

    @property
    def n_star(self):
        """Number of Star (vertex-based, irrotational) DOFs."""
        return self._n_star

    @property
    def n_edge(self):
        """Number of edge DOFs (HDivSurface order=0)."""
        return self._n_edges

    def transform_matrix(self, A_edge):
        """Transform BEM matrix from edge basis to Loop-Star basis.

        Congruence transformation: A_LS = T^T * A_edge * T
        Preserves symmetry: if A_edge is symmetric, A_LS is symmetric.

        Args:
            A_edge: BEM matrix in edge basis (n_edges x n_edges)

        Returns:
            A_LS: Matrix in Loop-Star basis (n_edges x n_edges)
                  Block structure: [[A_LL, A_LS], [A_SL, A_SS]]
        """
        A_edge = np.asarray(A_edge)
        if A_edge.shape != (self._n_edges, self._n_edges):
            raise ValueError(
                f"Matrix shape {A_edge.shape} doesn't match "
                f"n_edges={self._n_edges}")
        return self.T.T @ A_edge @ self.T

    def extract_blocks(self, A_LS):
        """Extract Loop-Loop, Loop-Star, Star-Loop, Star-Star blocks.

        Args:
            A_LS: Matrix in Loop-Star basis (n_edges x n_edges)

        Returns:
            A_LL: (n_loop, n_loop) Loop-Loop block (inductance-like)
            A_LS_block: (n_loop, n_star) Loop-Star coupling
            A_SL: (n_star, n_loop) Star-Loop coupling
            A_SS: (n_star, n_star) Star-Star block (capacitance-like)
        """
        nL = self._n_loop
        A_LL = A_LS[:nL, :nL]
        A_LS_block = A_LS[:nL, nL:]
        A_SL = A_LS[nL:, :nL]
        A_SS = A_LS[nL:, nL:]
        return A_LL, A_LS_block, A_SL, A_SS

    def transform_vector(self, b_edge):
        """Transform RHS vector from edge basis to Loop-Star basis.

        b_LS = T^T * b_edge

        Args:
            b_edge: Vector in edge basis (n_edges,) or (n_edges, k)

        Returns:
            b_LS: Vector in Loop-Star basis
        """
        return self.T.T @ np.asarray(b_edge)

    def decompose_current(self, J_edge):
        """Decompose edge-basis current into Loop and Star components.

        Solves: T * [alpha_loop; alpha_star] = J_edge

        Physical meaning:
        - J_loop = T_loop * alpha_loop: eddy currents (div-free)
        - J_star = T_star * alpha_star: charge currents (curl-free)

        Verification: J_loop + J_star = J_edge (exact if T is invertible)

        Args:
            J_edge: Current coefficients in edge basis (n_edges,)

        Returns:
            J_loop: Loop (eddy) component in edge basis (n_edges,)
            J_star: Star (charge) component in edge basis (n_edges,)
            alpha_loop: Loop coefficients (n_loop,)
            alpha_star: Star coefficients (n_star,)
        """
        J_edge = np.asarray(J_edge, dtype=float)

        if self._invertible:
            alpha = self.T_inv @ J_edge
        else:
            alpha = np.linalg.lstsq(self.T, J_edge, rcond=None)[0]

        nL = self._n_loop
        alpha_loop = alpha[:nL]
        alpha_star = alpha[nL:]

        J_loop = self.T_loop @ alpha_loop
        J_star = self.T_star @ alpha_star

        return J_loop, J_star, alpha_loop, alpha_star

    def verify_orthogonality(self):
        """Verify T_star^T * T_loop = 0 (graph-theoretic orthogonality).

        Returns:
            max_abs: Maximum absolute value (should be ~0 to machine precision)
        """
        product = self.T_star.T @ self.T_loop
        return np.max(np.abs(product))

    def verify_divergence_free(self, M_LS):
        """Verify that Loop basis functions are divergence-free.

        M_LS is the divergence coupling matrix from ngbem:
            M_LS[i, j] = integral div(basis_j) * phi_i dS

        For loop basis: M_LS * T_loop should be ~0

        Args:
            M_LS: Divergence coupling matrix (n_star_ngbem x n_edges)

        Returns:
            max_abs: Maximum absolute value of M_LS * T_loop
        """
        M_LS = np.asarray(M_LS)
        product = M_LS @ self.T_loop
        return np.max(np.abs(product))

    def print_summary(self):
        """Print Loop-Star transform summary."""
        print("LoopStarTransform Summary:")
        print(f"  Mesh: {self._n_vertices} vertices, "
              f"{self._n_edges} edges, {self._n_faces} faces")
        print(f"  Euler characteristic chi = {self._chi} "
              f"({'closed' if self._chi == 2 else 'open'} surface)")
        print(f"  Loop DOFs (face-based, solenoidal): {self._n_loop}")
        print(f"  Star DOFs (vertex-based, irrotational): {self._n_star}")
        print(f"  Total: {self._n_loop} + {self._n_star} = "
              f"{self._n_loop + self._n_star} (edges: {self._n_edges})")
        print(f"  T matrix: {self.T.shape}, "
              f"invertible: {self._invertible}")

        ortho = self.verify_orthogonality()
        print(f"  Orthogonality |T_star^T * T_loop|_max = {ortho:.2e}")

        if not self._invertible:
            print("  WARNING: T matrix is near-singular, using pseudoinverse")


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

    def to_topology_dict(self, decompose_loop_star=False):
        """Convert ngbem matrices to topology_dict format.

        This is the main entry point. Produces a dict compatible
        with PEECCircuitSolver.__init__().

        Args:
            decompose_loop_star: If True, apply Loop-Star decomposition
                to transform edge-basis matrices into Loop-Star blocks.
                L_LL becomes the PEEC inductance, L_SS the capacitance-like.
                Only supported for order=0.

        Returns:
            topology_dict: dict with L, R, segment_nodes, n_nodes,
                          n_loop, ports, and optional P, M_LS, n_star.
                          If decompose_loop_star=True, also includes
                          'loop_star_transform' with the LoopStarTransform
                          object and Loop-Star block matrices.
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

            # Optional: full Loop-Star (from ngbem product space)
            'P': mats['P'],
            'M_LS': mats['M_LS'],
            'n_star': mats['n_star'],

            # Metadata
            'backend': 'ngbem',
            'order': self.solver.order,
        }

        # Build incidence data for MNA solver
        topo.update(self._build_incidence_data())

        # --- Loop-Star decomposition ---
        if decompose_loop_star:
            if self.solver.order > 0:
                raise NotImplementedError(
                    "Loop-Star decomposition is only supported for order=0. "
                    "High-order (order>0) requires extended Helmholtz "
                    "decomposition (Andriulli et al., 2008).")

            ls = LoopStarTransform(self.mesh)

            # Transform L matrix to Loop-Star basis
            L_LS = ls.transform_matrix(mats['L'])
            L_LL, L_LS_block, L_SL, L_SS = ls.extract_blocks(L_LS)

            # Transform R to Loop-Star basis
            R_diag = np.diag(mats['R_loop'])
            R_LS = ls.transform_matrix(R_diag)
            R_LL, R_LS_block, R_SL, R_SS = ls.extract_blocks(R_LS)

            # Verify divergence-free property
            div_free_error = 0.0
            if mats['M_LS'] is not None:
                div_free_error = ls.verify_divergence_free(mats['M_LS'])

            topo['loop_star'] = {
                'transform': ls,
                'L_LL': L_LL,           # Loop inductance (n_loop x n_loop)
                'L_LS': L_LS_block,     # Loop-Star coupling
                'L_SL': L_SL,           # Star-Loop coupling
                'L_SS': L_SS,           # Star "inductance" (small)
                'R_LL': R_LL,           # Loop resistance
                'R_SS': R_SS,           # Star resistance
                'n_loop_ls': ls.n_loop, # Loop DOFs (face-based)
                'n_star_ls': ls.n_star, # Star DOFs (vertex-based)
                'T_loop': ls.T_loop,    # Transformation matrices
                'T_star': ls.T_star,
                'T': ls.T,
                'T_inv': ls.T_inv,
                'orthogonality_error': ls.verify_orthogonality(),
                'div_free_error': div_free_error,
            }

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
