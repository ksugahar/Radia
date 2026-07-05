"""
peec_shielded.py

PEEC solver with conducting shield coupling via ShieldBEMSIBC.

Computes port impedance Z(f) for PEEC coils near conducting shields.
The shield's eddy current effect is captured via Delta_Z_shield(freq),
a complex frequency-dependent impedance coupling matrix computed from
the BEM solution of surface currents.

At each frequency:
    Z_branch(f) = diag(R_dc + Zs_fil) + jw*L_air + Delta_Z_shield(f)

Key difference from retired magnetic-coupled PEEC solver:
    - Delta_Z_shield is COMPLEX and FREQUENCY-DEPENDENT (BEM solve per freq)
    - The retired magnetic-coupled PEEC Delta_L is REAL and FREQUENCY-INDEPENDENT (static)

Usage:
    from radia.ngsbem_eddy import ShieldBEMSIBC
    from radia.peec_shielded import ShieldedPEECSolver
    from radia.peec_matrices import PEECBuilder

    # Build PEEC topology
    builder = PEECBuilder()
    n1 = builder.add_node_at(0, 0, 0)
    n2 = builder.add_node_at(0.1, 0, 0)
    builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
    builder.add_port(n1, n2)
    topo = builder.build_topology()

    # Build shield BEM solver (NGSolve mesh)
    shield = ShieldBEMSIBC(mesh, sigma=3.7e7)
    shield.assemble(intorder=4)

    # Coupled solve
    solver = ShieldedPEECSolver(topo, shield)
    Z = solver.compute_port_impedance(1e6)

Part of Radia project
"""

import numpy as np
from radia.peec_topology import PEECCircuitSolver


class ShieldedPEECSolver(PEECCircuitSolver):
    """PEEC solver with conducting shield coupling via ShieldBEMSIBC.

    At each frequency, computes Delta_Z_shield from BEM solution and adds
    it to the PEEC branch impedance before MNA port extraction.

    Delta_Z_shield[i,j] captures the back-EMF on PEEC filament i due to
    shield eddy currents induced by unit current in filament j.
    """

    def __init__(self, topology_dict, shield_solver):
        """Initialize shielded PEEC solver.

        Args:
            topology_dict: Dict from PEECBuilder.build_topology() containing:
                L, R, segment_nodes, n_nodes, ports, n_loop,
                segment_centers, segment_directions, segment_lengths
            shield_solver: ShieldBEMSIBC instance (must be assembled)
        """
        super().__init__(topology_dict)

        self.shield = shield_solver
        self._topology_dict = topology_dict

        # Build full node incidence matrix A_full (n_nodes x n_loop)
        # A[node, fil] = +1 if current leaves node (node_from)
        # A[node, fil] = -1 if current enters node (node_to)
        self.A_full = np.zeros((self.n_nodes, self.n_loop))
        for j in range(self.n_loop):
            nf = self.segment_nodes[j, 0]  # node_from
            nt = self.segment_nodes[j, 1]  # node_to
            self.A_full[nf, j] = +1.0
            self.A_full[nt, j] = -1.0

        # Cache for Delta_Z at each frequency (avoid recomputing if same freq)
        self._cached_freq = None
        self._cached_Delta_Z = None

    def _get_delta_z_shield(self, freq):
        """Get Delta_Z_shield at given frequency, using cache if available.

        Args:
            freq: Frequency [Hz]

        Returns:
            Delta_Z_shield: (n_loop x n_loop) complex impedance matrix [Ohm]
        """
        if self._cached_freq is not None and abs(freq - self._cached_freq) < 1e-6:
            return self._cached_Delta_Z

        Delta_Z = self.shield.compute_impedance_matrix(freq, self._topology_dict)
        self._cached_freq = freq
        self._cached_Delta_Z = Delta_Z
        return Delta_Z

    def compute_port_impedance(self, freq, Zs=None):
        """Compute port impedance with shield coupling at given frequency.

        Overrides base class to add Delta_Z_shield to the branch impedance
        before the MNA solve.

        Z_branch = diag(R_dc + Zs) + jw*L_air + Delta_Z_shield(freq)

        Args:
            freq: Frequency in Hz
            Zs: Surface impedance array (n_loop,) complex, or None for DC

        Returns:
            Complex port impedance [Ohm]
        """
        omega = 2.0 * np.pi * freq

        # Build base branch impedance
        Z_branch = np.diag(self.R_dc.astype(complex)) + 1j * omega * self.L

        if Zs is not None:
            Zs = np.asarray(Zs)
            if Zs.shape[0] == self.n_loop:
                Z_branch += np.diag(Zs)

        # Add shield coupling (full matrix, complex, frequency-dependent)
        if freq > 0:
            Delta_Z = self._get_delta_z_shield(freq)
            Z_branch += Delta_Z

        if len(self.ports) == 0:
            return 0.0 + 0.0j

        port = self.ports[0]
        node_pos, node_neg = port[0], port[1]

        # Branch admittance
        Y_branch = np.linalg.inv(Z_branch)

        # Nodal admittance: Y_node = A * Y_branch * A^T
        A = self.A_full
        Y_node = A @ Y_branch @ A.T

        # Reorder nodes: ground node last
        n = self.n_nodes
        node_order = [i for i in range(n) if i != node_neg] + [node_neg]
        inv_order = [0] * n
        for i, j in enumerate(node_order):
            inv_order[j] = i

        Y_perm = Y_node[np.ix_(node_order, node_order)]
        Y_reduced = Y_perm[:n-1, :n-1]

        I_ext = np.zeros(n-1, dtype=complex)
        pos_idx = inv_order[node_pos]
        if pos_idx < n-1:
            I_ext[pos_idx] = 1.0

        V_reduced = np.linalg.solve(Y_reduced, I_ext)
        Z_port = V_reduced[pos_idx]

        return Z_port

    def frequency_sweep(self, freqs, Zs_func=None):
        """Compute port impedance over a frequency range with shield coupling.

        At each frequency, computes Delta_Z_shield from BEM solution.

        Args:
            freqs: Array of frequencies [Hz]
            Zs_func: Optional callable(freq) returning Zs array (n_loop,)

        Returns:
            Z_port: Complex array of port impedances, shape (len(freqs),)
        """
        freqs = np.asarray(freqs)
        Z_port = np.zeros(len(freqs), dtype=complex)

        for i, f in enumerate(freqs):
            # Invalidate cache for each new frequency
            self._cached_freq = None
            Zs = Zs_func(f) if Zs_func is not None else None
            Z_port[i] = self.compute_port_impedance(f, Zs)

        return Z_port
