"""
peec_topology.py

PEEC Circuit Solver with Node-Segment Topology

Computes port impedance Z(f) using Modified Nodal Analysis (MNA)
with full node incidence matrix from PEECBuilder.build_topology().

All linear algebra is delegated to C++ LAPACK via MNASolver
(zgesv_/zgetrf_/zgetrs_ + cblas_zgemm, shared MKL infrastructure with face integration).

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
    from radia.peec_matrices import PEECBuilder
    from radia.peec_topology import PEECCircuitSolver

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
from radia.peec_matrices import MNASolver


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

    def __init__(self, topology_dict, *,
                 use_hacapk=False,
                 hacapk_params=None,
                 bicgstab_outer_tol=1e-8,
                 bicgstab_inner_tol=1e-10,
                 bicgstab_maxiter=1000,
                 outer_method="saddle",
                 outer_m=50):
        """
        Initialize from build_topology() result dict.

        Args:
            topology_dict: Dict from PEECBuilder.build_topology() containing:
                L, R, segment_nodes, n_nodes, ports, n_loop
                Optional: P, M_LS, n_star (for full RLCM mode)
                Optional: C_gathered (for proper MNA with resonance)
                Optional (required when use_hacapk=True):
                    segment_centers, segment_directions, segment_lengths,
                    segment_widths, segment_heights, segment_sigmas
            use_hacapk: If True, replace the C++ dense MNA solver with a
                Python MNA driver that uses a HACApK H-matrix for L and
                a flexible outer Krylov for the branch Z_branch^{-1}
                matvec.  Only compute_Z_matrix / compute_port_impedance /
                frequency_sweep / compute_coupling_coefficient are
                supported in this mode (others raise NotImplementedError).
            hacapk_params: Optional dict forwarded to PEECHACApKSolver
                (keys: aca_eps, leaf_size, eta, max_rank, print_level).
            bicgstab_outer_tol: Relative tol for the nodal outer Krylov.
            bicgstab_inner_tol: Relative tol for the branch (HACApK) solve.
                Should be tighter than bicgstab_outer_tol.
            bicgstab_maxiter: Max iterations for both outer and inner
                Krylov calls.
            outer_method: Nodal solver for the HACApK-MNA path.
                ``"saddle"`` (default, stage-6b/7) solves the combined
                2x2 block system
                ``[Z_branch A^T; A 0] [I; V] = [0; -i_ext]`` with a
                single lgmres call and a block-diagonal preconditioner
                (``diag(Z)^{-1}`` + sparse LU of the Schur complement
                ``A diag(Z)^{-1} A^T``).  No nested Krylov; machine-
                precision accuracy and crossover with dense LU at
                N ~ 3000 (see validation_test/solver_benchmarks/
                findings_peec_mna_crossover.md).
                ``"lgmres"`` and ``"gcrotmk"`` use the older nested
                (outer Krylov + inner HACApK BiCGSTAB) path; ``"lgmres"``
                is robust but 100x slower at small N, ``"gcrotmk"``
                NaN-divides at high omega.  ``"bicgstab"`` is the
                stage-4 legacy and stalls above N ~ 200.
            outer_m: Krylov subspace dimension for lgmres / gcrotmk.
                Also used as the ``inner_m`` of the saddle-point lgmres.
        """
        self.n_loop = topology_dict['n_loop']
        self.n_nodes = topology_dict['n_nodes']

        # Segment connectivity: (n_filaments, 2) array of [node_from, node_to]
        self.segment_nodes = np.array(topology_dict['segment_nodes'])

        # Port definitions: list of (node_positive, node_negative, port_id)
        self.ports = topology_dict['ports']

        # HACApK switch + parameters
        self._use_hacapk = bool(use_hacapk)
        hp = dict(hacapk_params) if hacapk_params else {}
        # Override the C++ PEEC default aca_eps=1e-4 with 1e-8 unless
        # the caller explicitly set it.  1e-8 is the sweet spot for the
        # Ruehli mutual-inductance kernel: ACA rank jumps from 1 (16%
        # L matvec error) to ~100 (< 0.4% error) between 1e-7 and 1e-8,
        # while memory stays within 2x of the 1e-4 setting.  See
        # findings_peec_mna_crossover.md "Stage-8" section.
        hp.setdefault("aca_eps", 1.0e-8)
        self._hacapk_params = hp
        self._bicgstab_outer_tol = float(bicgstab_outer_tol)
        self._bicgstab_inner_tol = float(bicgstab_inner_tol)
        self._bicgstab_maxiter = int(bicgstab_maxiter)
        method = str(outer_method).lower()
        if method not in ("bicgstab", "gcrotmk", "lgmres", "saddle"):
            raise ValueError(
                "outer_method must be one of 'bicgstab', 'gcrotmk', "
                "'lgmres', 'saddle', got %r" % outer_method
            )
        self._outer_method = method
        self._outer_m = int(outer_m)
        self._hacapk_solver = None   # lazy-built on first call

        if self._use_hacapk:
            required = ('segment_centers', 'segment_directions',
                        'segment_lengths', 'segment_widths',
                        'segment_heights', 'segment_sigmas')
            missing = [k for k in required if topology_dict.get(k) is None]
            if missing:
                raise KeyError(
                    "use_hacapk=True requires per-filament geometry in "
                    "topology_dict (missing: %s). Rebuild via "
                    "PEECBuilder.build_topology() from a recent Radia."
                    % missing
                )
            self._seg_centers = np.ascontiguousarray(
                topology_dict['segment_centers'], dtype=np.float64)
            self._seg_directions = np.ascontiguousarray(
                topology_dict['segment_directions'], dtype=np.float64)
            self._seg_lengths = np.ascontiguousarray(
                topology_dict['segment_lengths'], dtype=np.float64)
            self._seg_widths = np.ascontiguousarray(
                topology_dict['segment_widths'], dtype=np.float64)
            self._seg_heights = np.ascontiguousarray(
                topology_dict['segment_heights'], dtype=np.float64)
            self._seg_sigmas = np.ascontiguousarray(
                topology_dict['segment_sigmas'], dtype=np.float64)
            # No dense L needed; R is recomputed inside PEECHACApKSolver.
            self.L = None
            self.R_dc = None
        else:
            self.L = np.array(topology_dict['L'])
            self.R_dc = np.array(topology_dict['R'])

        # Per-segment geometry kept for post-processing exports.  Only
        # O(n_seg) arrays are retained (never the dense L), so this is
        # free even in the HACApK path where L is deliberately not
        # materialized.  Consumed by export_gmsh().
        R_geom = topology_dict.get('R', None)
        if R_geom is not None:
            R_geom = np.asarray(R_geom, dtype=np.float64)
            if R_geom.ndim == 2:
                R_geom = np.diag(R_geom)
        self._export_geometry = {
            k: topology_dict.get(k, None)
            for k in ('segment_centers', 'segment_directions',
                      'segment_lengths', 'segment_widths',
                      'segment_heights')
        }
        self._export_geometry['R'] = R_geom

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

        # Create C++ MNA solver from raw arrays (dense path only; HACApK
        # path goes through _compute_Z_matrix_hacapk and skips this).
        if self._use_hacapk:
            if self.has_panels:
                raise NotImplementedError(
                    "use_hacapk=True does not yet support P/M_LS (panel "
                    "capacitance / full RLCM). Drop include_star=True when "
                    "calling build_topology()."
                )
            self._solver = None
        else:
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
        if self._use_hacapk:
            raise NotImplementedError(
                "set_solver_method is not applicable when use_hacapk=True"
            )
        self._solver.set_solver_method(method)

    def set_bicgstab_params(self, tol=1e-10, max_iter=1000):
        """Set BiCGSTAB solver parameters.

        Args:
            tol: Convergence tolerance (default 1e-10)
            max_iter: Maximum iterations (default 1000)
        """
        if self._use_hacapk:
            self._bicgstab_inner_tol = float(tol)
            self._bicgstab_maxiter = int(max_iter)
            return
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

        if self._use_hacapk:
            Z_mat = self._compute_Z_matrix_hacapk(freq, Zs)
            return Z_mat[0, 0]

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

        if self._use_hacapk:
            # Python loop; HACApK H-matrix is reused across freqs.
            n_freq = len(freqs)
            Z_port = np.zeros(n_freq, dtype=complex)
            for i, f in enumerate(freqs):
                Zs = Zs_func(f) if Zs_func is not None else None
                Z_mat = self._compute_Z_matrix_hacapk(float(f), Zs)
                Z_port[i] = Z_mat[0, 0]
            return Z_port

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

    def compute_branch_currents(self, freq, port_currents, Zs=None):
        """Solve for per-filament complex branch currents under user-specified
        port current injection.

        Use case: all N filaments wired in parallel between a single port's
        positive and negative nodes. Pass port_currents=[I_total] to get the
        AC current distribution across the N filaments (skin + proximity).

        Args:
            freq: Frequency in Hz
            port_currents: array-like of complex port currents, length n_ports
            Zs: Optional surface impedance array (n_loop,) complex

        Returns:
            I_branch: (n_loop,) complex array. I_branch[f] > 0 means current
                flows from segment_nodes[f, 0] to segment_nodes[f, 1].
                For a parallel filament bundle, sum(I_branch) == port_currents[0].
        """
        if self._use_hacapk:
            pc = np.asarray(port_currents, dtype=complex).reshape(-1)
            return self._compute_branch_currents_hacapk(freq, pc, Zs)
        pc = np.asarray(port_currents, dtype=complex).reshape(-1)
        Zs_arg = np.asarray(Zs, dtype=complex) if Zs is not None else None
        return np.array(self._solver.solve_branch_currents(freq, pc, Zs_arg))

    def export_gmsh(self, output_msh, freq, port_currents, *, Zs=None,
                    **kwargs):
        """Solve at one frequency and write the result as a GMSH post file.

        One hop from a solved PEEC circuit to a viewable .msh: solves the
        branch currents at ``freq`` and hands them, together with the
        per-segment geometry this solver was built from, to
        ``radia.gmsh_post_export.export_peec_topology_msh``.  Views:
        ``|I| [A]``, ``|J| [A/m^2]``, ``P [W]`` (plus the opt-in arrow
        and complex two-step views).

        Overlay the conductor CAD in the same picture with the radia-mcp
        gmsh server::

            solver.export_gmsh("coil_I.msh", 50e3, [1.0],
                               direction_view=True, complex_steps=True)
            gmsh_render("coil_I.msh", "coil_I.png",
                        merge_files=["coil.step"])

        Args:
            output_msh: output .msh path.
            freq: frequency in Hz.
            port_currents: array-like of complex port currents
                (length n_ports), as for compute_branch_currents.
            Zs: optional per-segment surface impedance; used BOTH in the
                solve and in the dissipation (its real part).
            **kwargs: forwarded to export_peec_topology_msh
                (current_convention, groups, group_names, label,
                direction_view, complex_steps).

        Returns:
            The export summary dict, with the frequency added.
        """
        from radia.gmsh_post_export import export_peec_topology_msh

        I_branch = self.compute_branch_currents(freq, port_currents, Zs=Zs)
        summary = export_peec_topology_msh(
            self._export_geometry, output_msh, currents=I_branch, Zs=Zs,
            **kwargs)
        summary["frequency_Hz"] = float(freq)
        return summary

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
        if self._use_hacapk:
            return self._compute_Z_matrix_hacapk(freq, Zs)
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

    # ------------------------------------------------------------------
    # HACApK MNA path (use_hacapk=True)
    # ------------------------------------------------------------------

    def _get_hacapk_solver(self):
        """Lazily instantiate the PEECHACApKSolver (H-matrix of L)."""
        if self._hacapk_solver is None:
            from radia.peec_hacapk_solver import PEECHACApKSolver
            hp = self._hacapk_params
            self._hacapk_solver = PEECHACApKSolver(
                self._seg_centers,
                self._seg_directions,
                self._seg_lengths,
                self._seg_widths,
                self._seg_heights,
                self._seg_sigmas,
                aca_eps=hp.get('aca_eps', -1.0),
                leaf_size=hp.get('leaf_size', -1),
                eta=hp.get('eta', -1.0),
                max_rank=hp.get('max_rank', -1),
                print_level=hp.get('print_level', 0),
            )
        return self._hacapk_solver

    def _build_reduced_incidence(self):
        """Return (A_hat, keep, gnd) for the nodal MNA system.

        A_hat is a scipy.sparse CSR matrix of shape (n_nodes-1, n_loop)
        representing the reduced incidence after grounding one node.
        ``keep`` is the list of original node indices retained (row order).
        ``gnd`` is the index of the grounded node.

        Convention: A[i, j] = -1 if node i is branch j's source
        (``segment_nodes[j, 0]``), +1 if node i is the sink
        (``segment_nodes[j, 1]``). Branch KVL is ``Z_branch I = -A^T V``
        so that the nodal admittance ``Y = A Z^{-1} A^T`` satisfies
        ``Y V = i_ext`` with ``i_ext[pos] = +I, i_ext[neg] = -I``.
        """
        import scipy.sparse as sp

        N = self.n_loop
        n_nodes = self.n_nodes
        src = self.segment_nodes[:, 0].astype(np.int64)
        dst = self.segment_nodes[:, 1].astype(np.int64)

        rows = np.concatenate([src, dst])
        cols = np.concatenate([np.arange(N, dtype=np.int64),
                               np.arange(N, dtype=np.int64)])
        data = np.concatenate([-np.ones(N, dtype=np.float64),
                               np.ones(N, dtype=np.float64)])
        A_inc = sp.csr_matrix((data, (rows, cols)),
                              shape=(n_nodes, N))

        # Ground port 0's negative terminal (must exist since len(ports)>0
        # is checked by the caller).
        gnd = int(self.ports[0][1])
        keep = np.array([i for i in range(n_nodes) if i != gnd],
                        dtype=np.int64)
        A_hat = A_inc[keep, :]
        return A_hat, keep, gnd

    def _build_outer_jacobi(self, freq, Zs, A_hat):
        """Jacobi preconditioner for the nodal admittance Y = A Z^{-1} A^T.

        Since A_hat has entries in {-1, 0, +1}, (A_hat[i,j])**2 = |A_hat[i,j]|
        and the exact diagonal of Y is

            diag_Y[i] = sum_{branch j incident to i} 1 / Z_branch[j,j]

        which is one sparse matvec of |A_hat| against the reciprocal of
        ``diag(Z_branch) = R + Zs + j*omega*diag(L)``.  ``diag(L)`` is
        extracted once via ``PEECHACApKSolver.get_L_diag()`` (cached).
        """
        from scipy.sparse.linalg import LinearOperator

        omega = 2.0 * np.pi * float(freq)
        solver = self._get_hacapk_solver()
        R_diag = np.asarray(solver._R_dc_default, dtype=np.float64)
        Zs_diag = (np.zeros(self.n_loop, dtype=np.complex128)
                   if Zs is None else np.asarray(Zs, dtype=np.complex128))
        diag_L = solver.get_L_diag()
        diag_Z = (R_diag + Zs_diag.real) + 1j * (omega * diag_L + Zs_diag.imag)
        A_abs = abs(A_hat)
        diag_Y = A_abs @ (1.0 / diag_Z)
        tiny = np.max(np.abs(diag_Y)) * 1.0e-30
        diag_Y_safe = np.where(np.abs(diag_Y) > tiny, diag_Y, 1.0)

        n_red = A_hat.shape[0]
        return LinearOperator(
            shape=(n_red, n_red),
            matvec=lambda r: np.asarray(r) / diag_Y_safe,
            dtype=np.complex128,
        )

    def _solve_nodal_system(self, freq, Zs, A_hat, i_ext_red):
        """Solve ``A_hat Z_branch^{-1} A_hat^T V = i_ext_red`` via an
        outer Krylov (gcrotmk by default) + inner HACApK BiCGSTAB, with
        a diagonal Jacobi preconditioner on diag(Y).

        Raises RuntimeError on non-convergence (fail-loud policy).

        The default outer method is gcrotmk because plain BiCGSTAB
        stalls above N ~ 200 on strongly coupled helical coils: the
        inner HACApK BiCGSTAB breaks down at residual ~1e-7 and the
        resulting varying-operator perturbation destroys BiCGSTAB's
        three-term recurrence (see stage-5 findings).  gcrotmk's
        augmented subspace + FGMRES-style outer loop absorbs this noise.
        """
        from scipy.sparse.linalg import (
            LinearOperator, bicgstab, gcrotmk, lgmres,
        )

        solver = self._get_hacapk_solver()
        n_red = A_hat.shape[0]

        counters = {"outer": 0, "inner_iters": 0, "inner_matvec_real": 0}
        M_inner = solver.build_jacobi_preconditioner(float(freq), Zs=Zs)
        M_outer = self._build_outer_jacobi(float(freq), Zs, A_hat)

        def matvec_Y(v):
            counters["outer"] += 1
            b = A_hat.T @ np.asarray(v, dtype=np.complex128)
            b = np.ascontiguousarray(b.reshape(-1))
            res = solver.solve(
                float(freq), b, Zs=Zs,
                rtol=self._bicgstab_inner_tol,
                maxiter=self._bicgstab_maxiter,
                M=M_inner,
                track_residual=False,
            )
            if not res.converged:
                raise RuntimeError(
                    "HACApK inner BiCGSTAB did not converge at f=%g Hz "
                    "(info=%d, final_residual=%.2e)"
                    % (freq, res.info, res.final_residual)
                )
            counters["inner_iters"] += res.iterations
            counters["inner_matvec_real"] += res.n_matvec_real
            return A_hat @ res.x

        Y_op = LinearOperator(shape=(n_red, n_red),
                              matvec=matvec_Y,
                              dtype=np.complex128)

        m = min(self._outer_m, max(2, n_red - 1))
        if self._outer_method == "gcrotmk":
            V_red, info = gcrotmk(
                Y_op, i_ext_red,
                rtol=self._bicgstab_outer_tol,
                maxiter=self._bicgstab_maxiter,
                m=m,
                M=M_outer,
            )
        elif self._outer_method == "lgmres":
            V_red, info = lgmres(
                Y_op, i_ext_red,
                rtol=self._bicgstab_outer_tol,
                maxiter=self._bicgstab_maxiter,
                inner_m=m,
                M=M_outer,
            )
        else:  # "bicgstab"
            V_red, info = bicgstab(
                Y_op, i_ext_red,
                rtol=self._bicgstab_outer_tol,
                maxiter=self._bicgstab_maxiter,
                M=M_outer,
            )

        # Residual-based convergence: scipy info alone is not reliable
        # under a varying-operator preconditioner; recompute once.
        b_norm = np.linalg.norm(i_ext_red)
        final_res = (np.linalg.norm(i_ext_red - Y_op @ V_red) / b_norm
                     if b_norm > 0 else 0.0)
        counters["outer"] = max(0, counters["outer"] - 1)
        if final_res > self._bicgstab_outer_tol * 10.0:
            raise RuntimeError(
                "Nodal %s did not converge at f=%g Hz "
                "(scipy info=%d, outer matvecs=%d, final_res=%.2e > %.2e)"
                % (self._outer_method, freq, info, counters["outer"],
                   final_res, self._bicgstab_outer_tol)
            )
        counters["final_outer_res"] = float(final_res)

        self._last_nodal_stats = counters
        return V_red

    def _solve_saddle_point(self, freq, Zs, A_hat, i_ext_red):
        """Solve the 2x2 block saddle-point system directly::

            [Z_branch   A_hat^T] [I]   [0        ]
            [A_hat       0     ] [V] = [i_ext_red]

        No nested Krylov: ONE outer lgmres iteration on the combined
        (N + n_red)-dimensional complex system, block-diagonal
        preconditioned with ``diag(Z_branch)^{-1}`` on the upper block
        and a sparse LU of the Schur complement
        ``S = A_hat diag(Z)^{-1} A_hat^T`` on the lower block.

        ``Z_branch = diag(R + Zs) + j*omega*L``; the L matvec uses the
        HACApK H-matrix, so each block matvec costs two real H-matrix
        calls plus O(|A_hat|) sparse-vector ops.

        Returns V_red (length n_red).
        """
        from scipy.sparse.linalg import LinearOperator, lgmres
        from scipy.sparse.linalg import splu
        import scipy.sparse as sp

        omega = 2.0 * np.pi * float(freq)
        solver = self._get_hacapk_solver()
        N = self.n_loop
        n_red = A_hat.shape[0]
        dim_total = N + n_red

        R_diag = np.asarray(solver._R_dc_default, dtype=np.float64)
        Zs_diag = (np.zeros(N, dtype=np.complex128)
                   if Zs is None else np.asarray(Zs, dtype=np.complex128))
        diag_L = solver.get_L_diag()
        Z_diag = (R_diag + Zs_diag.real) + 1j * (omega * diag_L + Zs_diag.imag)
        diag_scal = R_diag.astype(np.complex128) + Zs_diag

        A_hat_T = A_hat.T
        counters = {"matvec": 0}

        def block_matvec(xy):
            counters["matvec"] += 1
            I = xy[:N]
            V = xy[N:]
            I_re = np.ascontiguousarray(I.real)
            I_im = np.ascontiguousarray(I.imag)
            LI_re = solver.apply_L(I_re)
            LI_im = solver.apply_L(I_im)
            # Z_branch I = (R + Re Zs) I + j (Im Zs) I + j*omega*L I
            # j*omega*(LI_re + j*LI_im) = -omega*LI_im + j*omega*LI_re
            ZI = diag_scal * I + (-omega * LI_im + 1j * omega * LI_re)
            top = ZI + A_hat_T @ V
            bot = A_hat @ I
            out = np.empty(dim_total, dtype=np.complex128)
            out[:N] = top
            out[N:] = bot
            return out

        A_block = LinearOperator(shape=(dim_total, dim_total),
                                 matvec=block_matvec,
                                 dtype=np.complex128)

        # Block-diagonal preconditioner (Murphy-Golub-Wathen, 2000).
        inv_Z = 1.0 / Z_diag
        Diag_invZ = sp.diags(inv_Z, format='csc')
        # Schur complement S = A_hat diag(Z)^{-1} A_hat^T has the same
        # sparsity as A_hat A_hat^T (graph Laplacian for a connected
        # topology); for a chain of n_red=N nodes this is tridiagonal,
        # so splu is O(N).
        S_mat = (A_hat @ Diag_invZ @ A_hat_T).tocsc()
        S_lu = splu(S_mat)

        def precond(r):
            r_top = r[:N]
            r_bot = r[N:]
            out = np.empty(dim_total, dtype=np.complex128)
            out[:N] = inv_Z * r_top
            out[N:] = S_lu.solve(r_bot)
            return out

        M_block = LinearOperator(shape=(dim_total, dim_total),
                                 matvec=precond,
                                 dtype=np.complex128)

        # Sign convention matches _solve_nodal_system: the incidence
        # matrix is built so that ``Z_branch I = -A_hat^T V`` (see
        # _build_reduced_incidence docstring), hence
        # ``I = -Z^{-1} A_hat^T V`` and ``A_hat I = -A_hat Z^{-1} A_hat^T V
        # = -Y V``.  Requiring ``Y V = i_ext_red`` (the stage-4
        # convention) therefore gives ``A_hat I = -i_ext_red``; in block
        # form the RHS is ``[0; -i_ext_red]``.
        rhs = np.zeros(dim_total, dtype=np.complex128)
        rhs[N:] = -i_ext_red

        m = min(self._outer_m, max(2, dim_total - 1))
        xy, info = lgmres(
            A_block, rhs,
            rtol=self._bicgstab_outer_tol,
            maxiter=self._bicgstab_maxiter,
            inner_m=m,
            M=M_block,
        )

        b_norm = np.linalg.norm(rhs)
        final_res = (np.linalg.norm(rhs - A_block @ xy) / b_norm
                     if b_norm > 0 else 0.0)
        # A_block @ xy above adds one matvec for the residual check.
        n_solve_matvec = max(0, counters["matvec"] - 1)
        if final_res > self._bicgstab_outer_tol * 10.0:
            raise RuntimeError(
                "Saddle-point lgmres did not converge at f=%g Hz "
                "(info=%d, block matvecs=%d, final_res=%.2e > %.2e)"
                % (freq, info, n_solve_matvec,
                   final_res, self._bicgstab_outer_tol)
            )

        V_red = xy[N:]
        self._last_nodal_stats = {
            "outer": n_solve_matvec,
            "inner_iters": 0,
            "inner_matvec_real": 2 * n_solve_matvec,
            "final_outer_res": float(final_res),
        }
        return V_red

    def _dispatch_nodal_solve(self, freq, Zs, A_hat, i_ext_red):
        """Route to the configured nodal solver."""
        if self._outer_method == "saddle":
            return self._solve_saddle_point(freq, Zs, A_hat, i_ext_red)
        return self._solve_nodal_system(freq, Zs, A_hat, i_ext_red)

    def _compute_Z_matrix_hacapk(self, freq, Zs=None):
        """Nodal-MNA Z-parameter matrix via HACApK + nested BiCGSTAB."""
        n_ports = len(self.ports)
        if n_ports == 0:
            return np.zeros((0, 0), dtype=np.complex128)

        A_hat, keep, gnd = self._build_reduced_incidence()
        keep_map = {int(k): int(i) for i, k in enumerate(keep)}

        Zs_arr = None
        if Zs is not None:
            Zs_arr = np.asarray(Zs, dtype=np.complex128).reshape(-1)
            if Zs_arr.shape != (self.n_loop,):
                raise ValueError(
                    "Zs must have shape (%d,), got %s"
                    % (self.n_loop, Zs_arr.shape)
                )

        Z_mat = np.zeros((n_ports, n_ports), dtype=np.complex128)
        for p_idx, port in enumerate(self.ports):
            pos, neg, _ = int(port[0]), int(port[1]), int(port[2])
            i_ext = np.zeros(self.n_nodes, dtype=np.complex128)
            i_ext[pos] += 1.0
            i_ext[neg] -= 1.0
            i_ext_red = i_ext[keep]

            V_red = self._dispatch_nodal_solve(freq, Zs_arr, A_hat, i_ext_red)

            V_full = np.zeros(self.n_nodes, dtype=np.complex128)
            V_full[keep] = V_red  # V[gnd] = 0

            for q_idx, q_port in enumerate(self.ports):
                qpos = int(q_port[0])
                qneg = int(q_port[1])
                Z_mat[q_idx, p_idx] = V_full[qpos] - V_full[qneg]

        return Z_mat

    def _compute_branch_currents_hacapk(self, freq, port_currents, Zs=None):
        """Per-filament currents under port_current injection (HACApK MNA)."""
        n_ports = len(self.ports)
        if n_ports == 0:
            return np.zeros(self.n_loop, dtype=np.complex128)
        if port_currents.shape != (n_ports,):
            raise ValueError(
                "port_currents must have shape (%d,), got %s"
                % (n_ports, port_currents.shape)
            )

        A_hat, keep, gnd = self._build_reduced_incidence()

        Zs_arr = None
        if Zs is not None:
            Zs_arr = np.asarray(Zs, dtype=np.complex128).reshape(-1)
            if Zs_arr.shape != (self.n_loop,):
                raise ValueError(
                    "Zs must have shape (%d,), got %s"
                    % (self.n_loop, Zs_arr.shape)
                )

        # Build external current injection from the user's port_currents
        i_ext = np.zeros(self.n_nodes, dtype=np.complex128)
        for (pos, neg, _), I_p in zip(self.ports, port_currents):
            i_ext[int(pos)] += I_p
            i_ext[int(neg)] -= I_p
        i_ext_red = i_ext[keep]

        V_red = self._solve_nodal_system(freq, Zs_arr, A_hat, i_ext_red)
        V_full = np.zeros(self.n_nodes, dtype=np.complex128)
        V_full[keep] = V_red

        # I = -Z_branch^{-1} A^T V  (branch KVL: Z I = -A^T V)
        # Reuse the inner HACApK solve for the final branch-current back-sub.
        solver = self._get_hacapk_solver()
        # A_full^T V computed directly from segment endpoints (A_hat drops gnd
        # row, but V[gnd]=0 so the missing row contributes nothing).
        A_full_T_V = A_hat.T @ V_red
        res = solver.solve(
            float(freq), np.ascontiguousarray(A_full_T_V),
            Zs=Zs_arr,
            rtol=self._bicgstab_inner_tol,
            maxiter=self._bicgstab_maxiter,
            M=solver.build_jacobi_preconditioner(float(freq), Zs=Zs_arr),
            track_residual=False,
        )
        if not res.converged:
            raise RuntimeError(
                "HACApK back-sub BiCGSTAB did not converge at f=%g Hz "
                "(info=%d, final_residual=%.2e)"
                % (freq, res.info, res.final_residual)
            )
        return -res.x
