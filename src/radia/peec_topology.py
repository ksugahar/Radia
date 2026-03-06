"""
peec_topology.py

PEEC Circuit Solver with Node-Segment Topology

Computes port impedance Z(f) using Modified Nodal Analysis (MNA)
with full node incidence matrix from PEECBuilder.build_topology().

All linear algebra is delegated to C++ LAPACK via MNASolver
(zgesv_/zgetrf_/zgetrs_ + cblas_zgemm, shared MKL infrastructure with MSC).

Supports:
- Series connections (wire segments)
- Parallel connections (Litz wire strands, multi-filament)
- Multi-port extraction
- Frequency-dependent surface impedance (Bessel SIBC)
- Full RLCM system with P matrix (capacitive panels) for resonance

Full PEEC system (when panels present):

  [Z_LL   Z_LS] [I_L]   [V_L]       Z_LL = diag(R + Zs) + jw*L
  [Z_SL   Z_SS] [I_S] = [V_S]       Z_SS = P / (jw)
                                      Z_LS = jw * M_LS
                                      Z_SL = Z_LS^T

  Resonance: wL = 1/(wC) at f_res = 1/(2*pi*sqrt(LC))

Loop-only system (no panels, backward compatible):

  Z_branch = diag(R + Zs) + jw*L
  Y_node = A * Z_branch^{-1} * A^T
  Y_node * V = I_ext

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
from peec_matrices import MNASolver


class PEECCircuitSolver:
    """
    PEEC port impedance solver using nodal admittance (MNA).

    All linear algebra delegated to C++ LAPACK via MNASolver.

    Supports two modes:

    1. Loop-only (no panels): Classic PEEC with L, R, Zs
       Z_branch = diag(R_dc + Zs) + jw*L
       Y_node = A * Z_branch^{-1} * A^T

    2. Full RLCM (with panels): Includes capacitive coupling
       [Z_LL   Z_LS] [I_L]   [V_L]
       [Z_SL   Z_SS] [I_S] = [V_S]

       Z_LL = diag(R + Zs) + jw*L
       Z_SS = P / (jw)
       Z_LS = jw * M_LS
       Z_SL = Z_LS^T

    Mode is auto-detected from topology_dict contents.
    """

    def __init__(self, topology_dict):
        """
        Initialize from build_topology() result dict.

        Args:
            topology_dict: Dict from PEECBuilder.build_topology() containing:
                L, R, segment_nodes, n_nodes, ports, n_loop
                Optional: P, M_LS, n_star (for full RLCM mode)
                Optional: C_gathered (for proper MNA with resonance)
        """
        self.L = np.array(topology_dict['L'])
        self.R_dc = np.array(topology_dict['R'])
        self.n_loop = topology_dict['n_loop']
        self.n_nodes = topology_dict['n_nodes']

        # Segment connectivity: (n_filaments, 2) array of [node_from, node_to]
        self.segment_nodes = np.array(topology_dict['segment_nodes'])

        # Port definitions: list of (node_positive, node_negative, port_id)
        self.ports = topology_dict['ports']

        # Panel/Star data (optional - enables full RLCM mode)
        P_raw = topology_dict.get('P', None)
        M_LS_raw = topology_dict.get('M_LS', None)
        self.n_star = topology_dict.get('n_star', 0)

        if P_raw is not None and not isinstance(P_raw, type(None)):
            self.P = np.array(P_raw)
            self.has_panels = self.P.size > 0 and self.n_star > 0
        else:
            self.P = None
            self.has_panels = False

        if M_LS_raw is not None and not isinstance(M_LS_raw, type(None)):
            self.M_LS = np.array(M_LS_raw)
        else:
            self.M_LS = None

        # Gathered capacitance (optional - enables proper MNA with resonance)
        C_gathered_raw = topology_dict.get('C_gathered', None)
        if C_gathered_raw is not None and not isinstance(C_gathered_raw, type(None)):
            self.C_gathered = np.array(C_gathered_raw)
        else:
            self.C_gathered = None

        # Create C++ MNA solver from raw arrays
        seg_nodes_int = np.ascontiguousarray(self.segment_nodes, dtype=np.int32)
        port_tuples = [(int(p[0]), int(p[1]), int(p[2])) for p in self.ports]

        P_arg = self.P if self.has_panels else None
        M_LS_arg = self.M_LS if (self.has_panels and self.M_LS is not None) else None
        C_gath_arg = (np.ascontiguousarray(self.C_gathered, dtype=np.float64)
                      if self.C_gathered is not None else None)

        self._solver = MNASolver(
            np.ascontiguousarray(self.L, dtype=np.float64),
            np.ascontiguousarray(self.R_dc, dtype=np.float64),
            seg_nodes_int,
            self.n_nodes,
            port_tuples,
            P_arg,
            M_LS_arg,
            C_gath_arg,
        )

    def set_solver_method(self, method):
        """Set solver method for Z_branch inversion.

        Args:
            method: 0 = LU (default, LAPACK zgesv_), 1 = BiCGSTAB
        """
        self._solver.set_solver_method(method)

    def set_bicgstab_params(self, tol=1e-10, max_iter=1000):
        """Set BiCGSTAB solver parameters.

        Args:
            tol: Convergence tolerance (default 1e-10)
            max_iter: Maximum iterations (default 1000)
        """
        self._solver.set_bicgstab_params(tol, max_iter)

    def compute_port_impedance(self, freq, Zs=None):
        """
        Compute port impedance at a given frequency.

        Args:
            freq: Frequency in Hz
            Zs: Surface impedance array (n_loop,) complex, or None for DC

        Returns:
            Complex port impedance [Ohm]
        """
        if len(self.ports) == 0:
            return 0.0 + 0.0j

        Zs_arg = np.asarray(Zs, dtype=complex) if Zs is not None else None
        Z_mat = np.array(self._solver.solve_z_matrix(freq, Zs_arg))
        return Z_mat[0, 0]

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
        freqs = np.asarray(freqs, dtype=np.float64)

        if Zs_func is None:
            # Batch sweep in C++ (no per-frequency callback)
            result = self._solver.frequency_sweep(freqs)
            Z_all = np.array(result['Z_matrix'])  # (n_freq, n_ports, n_ports)
            return Z_all[:, 0, 0]

        # Per-frequency Zs: compute in Python, batch in C++
        n_freq = len(freqs)
        Z_port = np.zeros(n_freq, dtype=complex)
        for i, f in enumerate(freqs):
            Zs = Zs_func(f)
            Zs_arg = np.asarray(Zs, dtype=complex) if Zs is not None else None
            Z_mat = np.array(self._solver.solve_z_matrix(f, Zs_arg))
            Z_port[i] = Z_mat[0, 0]

        return Z_port

    def compute_Z_matrix(self, freq, Zs=None):
        """Compute full n_port x n_port Z-parameter matrix.

        For a 2-port system:
            Z = [[Z11, Z12],
                 [Z21, Z22]]

        where Z_ij = V_i / I_j with all other ports open-circuited.

        Args:
            freq: Frequency in Hz
            Zs: Surface impedance array (n_loop,) complex, or None

        Returns:
            Z_mat: (n_ports x n_ports) complex Z-parameter matrix
        """
        Zs_arg = np.asarray(Zs, dtype=complex) if Zs is not None else None
        return np.array(self._solver.solve_z_matrix(freq, Zs_arg))

    def compute_coupling_coefficient(self, freq, Zs=None):
        """Compute coupling coefficient k between ports from Z-parameters.

        Uses the relation:
            L_ij = Im(Z_ij) / omega
            k_ij = M_ij / sqrt(L_ii * L_jj)

        For a 2-port transformer:
            k = M12 / sqrt(L1 * L2)

        Args:
            freq: Frequency in Hz
            Zs: Surface impedance array (n_loop,) complex, or None

        Returns:
            dict with keys:
                Z_matrix: (n_ports x n_ports) complex Z-parameter matrix
                L_matrix: (n_ports x n_ports) inductance matrix [H]
                k_matrix: (n_ports x n_ports) coupling coefficient matrix
                k: Scalar coupling coefficient k_12 (for 2-port)
                L1, L2: Self-inductances of port 1 and 2 [H]
                M: Mutual inductance between port 1 and 2 [H]
        """
        Z_mat = self.compute_Z_matrix(freq, Zs)
        n_ports = Z_mat.shape[0]
        omega = 2.0 * np.pi * freq

        if n_ports < 2 or omega < 1e-10:
            return {
                'Z_matrix': Z_mat,
                'L_matrix': np.zeros_like(Z_mat, dtype=float),
                'k_matrix': np.zeros_like(Z_mat, dtype=float),
                'k': 0.0, 'L1': 0.0, 'L2': 0.0, 'M': 0.0,
            }

        # Extract inductance matrix from imaginary part of Z
        L = np.imag(Z_mat) / omega

        # Coupling coefficient matrix
        k_mat = np.zeros((n_ports, n_ports))
        for i in range(n_ports):
            for j in range(n_ports):
                if L[i, i] > 0 and L[j, j] > 0:
                    k_mat[i, j] = L[i, j] / np.sqrt(L[i, i] * L[j, j])

        return {
            'Z_matrix': Z_mat,
            'L_matrix': L,
            'k_matrix': k_mat,
            'k': k_mat[0, 1] if n_ports >= 2 else 0.0,
            'L1': L[0, 0],
            'L2': L[1, 1] if n_ports >= 2 else 0.0,
            'M': L[0, 1] if n_ports >= 2 else 0.0,
        }

    def frequency_sweep_multiport(self, freqs, Zs_func=None):
        """Compute Z-parameter matrix and coupling over a frequency range.

        Args:
            freqs: Array of frequencies [Hz]
            Zs_func: Optional callable(freq) returning Zs array (n_loop,)

        Returns:
            dict with keys:
                freqs: Frequency array [Hz]
                Z_matrix: (n_freq, n_ports, n_ports) complex Z-parameters
                L_matrix: (n_freq, n_ports, n_ports) inductance [H]
                k: (n_freq,) coupling coefficient k_12
                L1, L2: (n_freq,) self-inductances [H]
                M: (n_freq,) mutual inductance [H]
        """
        freqs = np.asarray(freqs, dtype=np.float64)
        n_freq = len(freqs)
        n_ports = len(self.ports)

        if Zs_func is None:
            # Batch sweep in C++
            result = self._solver.frequency_sweep(freqs)
            Z_all = np.array(result['Z_matrix'])
        else:
            # Pre-compute Zs array, then batch in C++
            n_loop = self.n_loop
            Zs_all = np.zeros((n_freq, n_loop), dtype=complex)
            for i, f in enumerate(freqs):
                Zs_val = Zs_func(f)
                if Zs_val is not None:
                    Zs_all[i] = Zs_val
            result = self._solver.frequency_sweep(freqs, Zs_all)
            Z_all = np.array(result['Z_matrix'])

        # Post-process: extract L, k from Z (numpy only, no scipy)
        L_all = np.zeros((n_freq, n_ports, n_ports))
        k_arr = np.zeros(n_freq)
        L1_arr = np.zeros(n_freq)
        L2_arr = np.zeros(n_freq)
        M_arr = np.zeros(n_freq)

        for i in range(n_freq):
            omega = 2.0 * np.pi * freqs[i]
            if omega > 1e-10:
                L_all[i] = np.imag(Z_all[i]) / omega

            L1_val = L_all[i, 0, 0]
            L2_val = L_all[i, 1, 1] if n_ports >= 2 else 0.0
            M_val = L_all[i, 0, 1] if n_ports >= 2 else 0.0

            L1_arr[i] = L1_val
            L2_arr[i] = L2_val
            M_arr[i] = M_val
            k_arr[i] = M_val / np.sqrt(L1_val * L2_val) if (L1_val > 0 and L2_val > 0) else 0.0

        return {
            'freqs': freqs,
            'Z_matrix': Z_all,
            'L_matrix': L_all,
            'k': k_arr,
            'L1': L1_arr,
            'L2': L2_arr,
            'M': M_arr,
        }

    def get_resonant_frequencies(self, freq_min=1e3, freq_max=1e9, n_points=1000):
        """
        Find self-resonant frequencies by scanning Im(Z) zero crossings.

        Only meaningful when panels are present (has_panels=True).
        Resonance occurs where Im(Z) changes sign (inductive -> capacitive).

        Args:
            freq_min: Scan start frequency [Hz]
            freq_max: Scan end frequency [Hz]
            n_points: Number of scan points

        Returns:
            List of resonant frequencies [Hz], or empty list if no resonance found
        """
        if not self.has_panels:
            return []

        freqs = np.logspace(np.log10(freq_min), np.log10(freq_max), n_points)
        Z_sweep = self.frequency_sweep(freqs)
        im_Z = Z_sweep.imag

        # Find zero crossings of Im(Z)
        resonances = []
        for i in range(len(freqs) - 1):
            if im_Z[i] * im_Z[i + 1] < 0:
                # Linear interpolation for zero crossing
                f1, f2 = freqs[i], freqs[i + 1]
                z1, z2 = im_Z[i], im_Z[i + 1]
                f_res = f1 - z1 * (f2 - f1) / (z2 - z1)
                resonances.append(float(f_res))

        return resonances
