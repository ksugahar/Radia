"""
peec_topology.py

PEEC Circuit Solver with Node-Segment Topology

Computes port impedance Z(f) using Modified Nodal Analysis (MNA)
with full node incidence matrix from PEECBuilder.build_topology().

Supports:
- Series connections (wire segments)
- Parallel connections (Litz wire strands, multi-filament)
- Multi-port extraction
- Frequency-dependent surface impedance (Bessel SIBC)

Usage:
    from peec_matrices import PEECBuilder
    from peec_topology import PEECCircuitSolver

    builder = PEECBuilder()
    n1 = builder.add_node_at(0, 0, 0)
    n2 = builder.add_node_at(0.05, 0, 0)
    n3 = builder.add_node_at(0.1, 0, 0)
    builder.add_connected_segment(n1, n2, 1e-3, 1e-3)
    builder.add_connected_segment(n2, n3, 1e-3, 1e-3)
    builder.add_port(n1, n3)

    topo = builder.build_topology()
    solver = PEECCircuitSolver(topo)
    Z = solver.compute_port_impedance(1e6)

Part of Radia project
"""

import numpy as np


class PEECCircuitSolver:
    """
    PEEC port impedance solver using nodal admittance (MNA).

    Formulation:
      Full node incidence: A_full (n_nodes x n_filaments)
      A_full[node, fil] = +1 if filament leaves node (node_from)
      A_full[node, fil] = -1 if filament enters node (node_to)

      Branch impedance: Z_branch = diag(R_dc + Zs) + jw*L

      Nodal admittance: Y_node = A_full * Z_branch^{-1} * A_full^T

      Port injection: I_ext[pos] = +1, I_ext[neg] = -1

      Solve: Y_node * V_node = I_ext
      Z_port = V_node[pos] - V_node[neg]

    For grounding, fix one node voltage to 0 (delete that row/col).
    """

    def __init__(self, topology_dict):
        """
        Initialize from build_topology() result dict.

        Args:
            topology_dict: Dict from PEECBuilder.build_topology() containing:
                L, R, segment_nodes, n_nodes, ports, n_loop
        """
        self.L = np.array(topology_dict['L'])
        self.R_dc = np.array(topology_dict['R'])
        self.n_loop = topology_dict['n_loop']
        self.n_nodes = topology_dict['n_nodes']

        # Segment connectivity: (n_filaments, 2) array of [node_from, node_to]
        self.segment_nodes = np.array(topology_dict['segment_nodes'])

        # Port definitions: list of (node_positive, node_negative, port_id)
        self.ports = topology_dict['ports']

        # Build full incidence matrix A_full (n_nodes x n_filaments)
        self._build_full_incidence()

    def _build_full_incidence(self):
        """Build full node incidence matrix from segment connectivity."""
        n_nodes = self.n_nodes
        n_fil = self.n_loop

        self.A_full = np.zeros((n_nodes, n_fil))

        for f in range(n_fil):
            nf = self.segment_nodes[f, 0]  # node_from
            nt = self.segment_nodes[f, 1]  # node_to

            if nf >= 0 and nf < n_nodes:
                self.A_full[nf, f] = +1.0  # filament leaves node_from
            if nt >= 0 and nt < n_nodes:
                self.A_full[nt, f] = -1.0  # filament enters node_to

    def compute_port_impedance(self, freq, Zs=None):
        """
        Compute port impedance at a given frequency.

        Method:
        1. Build Z_branch = diag(R_dc + Zs) + jw*L
        2. Compute Y_branch = Z_branch^{-1}
        3. Build Y_node = A * Y_branch * A^T
        4. Remove ground node (port negative), solve Y * V = I_ext
        5. Z_port = V[pos] - V[neg]

        Args:
            freq: Frequency in Hz
            Zs: Surface impedance array (n_loop,) complex, or None for DC

        Returns:
            Complex port impedance [Ohm]
        """
        omega = 2.0 * np.pi * freq

        # Branch impedance: Z_branch = diag(R_dc) + jw*L + diag(Zs)
        Z_branch = np.diag(self.R_dc.astype(complex)) + 1j * omega * self.L

        if Zs is not None:
            Zs = np.asarray(Zs)
            if Zs.shape[0] == self.n_loop:
                Z_branch += np.diag(Zs)

        if len(self.ports) == 0:
            return 0.0 + 0.0j

        port = self.ports[0]
        node_pos, node_neg = port[0], port[1]

        # Branch admittance
        Y_branch = np.linalg.inv(Z_branch)

        # Nodal admittance: Y_node = A * Y_branch * A^T
        A = self.A_full
        Y_node = A @ Y_branch @ A.T  # (n_nodes x n_nodes)

        # External current injection: I_ext[pos] = +1, I_ext[neg] = -1
        # Ground the negative terminal: V[neg] = 0
        # Remove the ground node row/col from the system

        # Reorder nodes: ground node last
        n = self.n_nodes
        node_order = [i for i in range(n) if i != node_neg] + [node_neg]
        inv_order = [0] * n
        for i, j in enumerate(node_order):
            inv_order[j] = i

        # Permute Y_node
        Y_perm = Y_node[np.ix_(node_order, node_order)]

        # Remove ground node (last row/col)
        Y_reduced = Y_perm[:n-1, :n-1]

        # Current injection (ground node removed)
        I_ext = np.zeros(n-1, dtype=complex)
        pos_idx = inv_order[node_pos]
        if pos_idx < n-1:
            I_ext[pos_idx] = 1.0  # 1A into positive terminal

        # Solve Y_reduced * V = I_ext
        V_reduced = np.linalg.solve(Y_reduced, I_ext)

        # Port voltage = V[pos] - V[neg] = V[pos] - 0 = V[pos]
        Z_port = V_reduced[pos_idx]

        return Z_port

    def frequency_sweep(self, freqs, Zs_func=None):
        """
        Compute port impedance over a frequency range.

        Args:
            freqs: Array of frequencies [Hz]
            Zs_func: Optional callable(freq) returning Zs array (n_loop,)
                     for surface impedance at each frequency.
                     If None, only DC resistance is used.

        Returns:
            Z_port: Complex array of port impedances, shape (len(freqs),)
        """
        freqs = np.asarray(freqs)
        Z_port = np.zeros(len(freqs), dtype=complex)

        for i, f in enumerate(freqs):
            Zs = Zs_func(f) if Zs_func is not None else None
            Z_port[i] = self.compute_port_impedance(f, Zs)

        return Z_port
