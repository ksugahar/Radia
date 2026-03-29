"""
peec_coupled.py

Coupled PEEC + MMM Solver

Computes port impedance Z(f) for conductors near magnetic materials.
Uses Biot-Savart for H-field from conductors, Radia Solve + A-field
for magnetic material response, yielding a coupling matrix Delta_L.

Effective impedance with complex permeability mu = mu'_r - j*mu"_r:
    Z_branch(f) = diag(R + Zs(f)) + jw * (L_air + Delta_L)
                  + omega * tan_delta_m * Delta_L

where tan_delta_m = mu"_r / mu'_r is the magnetic loss tangent.

The Delta_L matrix is computed from mu'_r (real part) via Radia Solve.
The mu"_r loss term adds frequency-proportional resistance (magnetic loss).

For linear materials, Delta_L is computed ONCE and reused for all frequencies.

Usage:
    import radia as rad
    from peec_matrices import PEECBuilder
    from peec_coupled import CoupledPEECSolver

    # Build conductor topology
    builder = PEECBuilder()
    n1 = builder.add_node_at(0, 0, 0)
    n2 = builder.add_node_at(0.1, 0, 0)
    builder.add_connected_segment(n1, n2, 1e-3, 1e-3)
    builder.add_port(n1, n2)
    topo = builder.build_topology()

    # Create magnetic object (coordinates in meters)
    verts = [[-0.02,-0.02,-0.02],[0.02,-0.02,-0.02],[0.02,0.02,-0.02],[-0.02,0.02,-0.02],
             [-0.02,-0.02,0.02],[0.02,-0.02,0.02],[0.02,0.02,0.02],[-0.02,0.02,0.02]]
    core = rad.ObjHexahedron(verts, [0, 0, 0])
    mat = rad.MatLin(999)
    rad.MatApl(core, mat)

    # Coupled solver with magnetic loss
    solver = CoupledPEECSolver(topo, [core], mu_r_imag=10.0)
    solver.compute_coupling_matrix()
    Z = solver.compute_port_impedance(1e6)

Part of Radia project
"""

import numpy as np
from peec_topology import PEECCircuitSolver

# Physical constants
MU_0 = 4.0 * np.pi * 1e-7  # H/m


class CoupledPEECSolver(PEECCircuitSolver):
    """
    Coupled PEEC + MMM solver for conductors near magnetic materials.

    Extends PEECCircuitSolver with magnetic material coupling via Delta_L matrix.

    The coupling matrix Delta_L[i][j] represents the additional inductance
    on segment i due to unit current in segment j, caused by the magnetized
    magnetic material.

    Physics:
        1. Unit current in segment j produces H-field via Biot-Savart
        2. H-field magnetizes the material (Radia Solve with ObjBckg)
        3. Magnetized material produces vector potential A at other segments
        4. Delta_L[i][j] = dot(A(center_i), dir_i) * length_i

    Complex permeability mu = mu'_r - j*mu"_r:
        - Delta_L is computed using mu'_r (real part) via Radia MatLin
        - mu"_r adds magnetic loss as frequency-proportional resistance:
          R_mag(f) = omega * tan_delta_m * Delta_L
          where tan_delta_m = mu"_r / mu'_r

    For linear materials, Delta_L is frequency-independent and computed once.
    """

    def __init__(self, topology_dict, magnetic_objects=None, mu_r_imag=0.0):
        """
        Initialize coupled solver.

        Args:
            topology_dict: Dict from PEECBuilder.build_topology() containing:
                L, R, segment_nodes, n_nodes, ports, n_loop,
                segment_centers, segment_directions, segment_lengths
            magnetic_objects: List of Radia object handles (hex/tet with material)
                             or None for uncoupled mode
            mu_r_imag: Imaginary part of relative permeability (positive value).
                       Represents magnetic loss: mu = mu'_r - j*mu"_r.
                       Set to 0 for lossless magnetic material.
        """
        super().__init__(topology_dict)

        # Store segment geometry for Biot-Savart and coupling
        self.seg_centers = np.array(topology_dict['segment_centers'])
        self.seg_directions = np.array(topology_dict['segment_directions'])
        self.seg_lengths = np.array(topology_dict['segment_lengths'])

        # Segment endpoints for Biot-Savart
        # Always reconstruct from center/direction/length (reliable)
        self._reconstruct_endpoints()

        # Magnetic objects
        self.magnetic_objects = magnetic_objects or []
        self.Delta_L = None  # Coupling matrix (computed on demand)

        # Complex permeability: mu = mu'_r - j*mu"_r
        # mu_r_imag is the POSITIVE imaginary part (loss factor)
        self.mu_r_imag = float(mu_r_imag)

        # mu'_r will be detected from the material applied to magnetic objects
        # (used to compute tan_delta = mu"_r / mu'_r)
        self._mu_r_real = None  # Set during compute_coupling_matrix or externally

    def _reconstruct_endpoints(self):
        """Reconstruct segment endpoints from center/direction/length."""
        n_seg = self.n_loop
        self.seg_p1 = np.zeros((n_seg, 3))
        self.seg_p2 = np.zeros((n_seg, 3))

        for i in range(n_seg):
            half = self.seg_directions[i] * self.seg_lengths[i] * 0.5
            self.seg_p1[i] = self.seg_centers[i] - half
            self.seg_p2[i] = self.seg_centers[i] + half

    @property
    def tan_delta_m(self):
        """Magnetic loss tangent: mu"_r / mu'_r."""
        if self.mu_r_imag == 0.0:
            return 0.0
        if self._mu_r_real is not None and self._mu_r_real > 0:
            return self.mu_r_imag / self._mu_r_real
        return 0.0

    def compute_coupling_matrix(self, solver_method=0, solver_prec=0.0001,
                                solver_maxiter=1000, mu_r_real=None):
        """
        Compute Delta_L matrix (n_seg x n_seg) from magnetic material response.

        For each column j:
        1. Compute H-field from unit current in segment j (Biot-Savart)
        2. Set as background field via rad.ObjBckg()
        3. Solve for induced magnetization M via rad.Solve()
        4. Compute A-field at each segment center
        5. Delta_L[i][j] = dot(A(center_i), dir_i) * length_i

        Args:
            solver_method: Radia solver method (0=LU, 1=BiCGSTAB, 2=HACApK)
            solver_prec: Solver precision
            solver_maxiter: Maximum iterations
            mu_r_real: Real part of mu_r (for tan_delta computation).
                       If None, auto-detected from Radia material.

        Returns:
            Delta_L matrix (n_seg x n_seg) in Henries
        """
        try:
            import sys
            import os
            sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
            from _radia_pybind import (ObjBckg, ObjCnt, Solve, Fld, UtiDel)
        except ImportError:
            import radia as rad
            ObjBckg = rad.ObjBckg
            ObjCnt = rad.ObjCnt
            Solve = rad.Solve
            Fld = rad.Fld
            UtiDel = rad.UtiDel

        # Store mu_r_real for tan_delta computation
        if mu_r_real is not None:
            self._mu_r_real = float(mu_r_real)

        if len(self.magnetic_objects) == 0:
            self.Delta_L = np.zeros((self.n_loop, self.n_loop))
            return self.Delta_L

        n_seg = self.n_loop
        Delta_L = np.zeros((n_seg, n_seg))

        # Create container of magnetic objects
        mag_container = ObjCnt(self.magnetic_objects)

        for j in range(n_seg):
            p1_j = self.seg_p1[j]
            p2_j = self.seg_p2[j]

            # H-field from unit current in segment j
            # TODO: Replace with ngsolve.bem MaxwellDL(J*dC) when available
            def h_field_from_seg_j(point, _p1=p1_j, _p2=p2_j):
                raise NotImplementedError(
                    "biot_savart_finite_filament removed. "
                    "Use ngsolve.bem MaxwellDL(J*dC) instead.")

            # Create background field
            bkg = ObjBckg(h_field_from_seg_j)

            # Create container with magnetic objects + background
            container = ObjCnt([mag_container, bkg])

            # Solve for induced magnetization
            Solve(container, solver_prec, solver_maxiter, solver_method)

            # Compute A-field at each segment center from magnetized material
            for i in range(n_seg):
                center_i = self.seg_centers[i].tolist()
                A_vec = Fld(mag_container, 'a', center_i)
                A_vec = np.array(A_vec)

                # Delta_L[i][j] = dot(A, dir_i) * length_i
                Delta_L[i, j] = np.dot(A_vec, self.seg_directions[i]) * self.seg_lengths[i]

            # Clean up background field
            UtiDel(bkg)

        self.Delta_L = Delta_L

        # Store original L_air and update L to include coupling
        self._L_air = self.L.copy()
        self.L = self._L_air + Delta_L

        # Rebuild C++ MNASolver with coupled L
        # (The C++ solver holds its own copy of L, so we must rebuild it)
        self._rebuild_mna_solver()

        return Delta_L

    def _rebuild_mna_solver(self):
        """Rebuild C++ MNASolver with current self.L matrix."""
        from peec_matrices import MNASolver

        seg_nodes_int = np.ascontiguousarray(self.segment_nodes, dtype=np.int32)
        port_tuples = [(int(p[0]), int(p[1]), int(p[2])) for p in self.ports]

        P_arg = self.P if self.has_panels else None
        M_LS_arg = self.M_LS if (self.has_panels and self.M_LS is not None) else None

        self._solver = MNASolver(
            np.ascontiguousarray(self.L, dtype=np.float64),
            np.ascontiguousarray(self.R_dc, dtype=np.float64),
            seg_nodes_int,
            self.n_nodes,
            port_tuples,
            P_arg,
            M_LS_arg,
        )

    def compute_port_impedance(self, freq, Zs=None):
        """
        Compute port impedance at a given frequency with magnetic coupling.

        Uses L_total = L_air + Delta_L (set during compute_coupling_matrix).

        For lossy materials (mu"_r > 0), adds frequency-proportional
        magnetic loss resistance: R_mag(f) = omega * tan_delta * diag(Delta_L).

        Args:
            freq: Frequency in Hz
            Zs: Surface impedance array (n_loop,) complex, or None for DC

        Returns:
            Complex port impedance [Ohm]
        """
        if self.Delta_L is not None and self.tan_delta_m > 0:
            omega = 2.0 * np.pi * freq
            if omega > 0:
                # R_mag = omega * tan_delta * diag(Delta_L)
                R_mag_diag = omega * self.tan_delta_m * np.diag(self.Delta_L)
                Zs_total = R_mag_diag.astype(complex)
                if Zs is not None:
                    Zs_total += np.asarray(Zs, dtype=complex)
                return super().compute_port_impedance(freq, Zs_total)

        return super().compute_port_impedance(freq, Zs)

    def frequency_sweep(self, freqs, Zs_func=None):
        """
        Compute port impedance over a frequency range with magnetic coupling.

        For linear materials, Delta_L is frequency-independent. Call
        compute_coupling_matrix() once before frequency_sweep().

        Args:
            freqs: Array of frequencies [Hz]
            Zs_func: Optional callable(freq) returning Zs array (n_loop,)

        Returns:
            Z_port: Complex array of port impedances, shape (len(freqs),)
        """
        if self.Delta_L is not None and self.tan_delta_m > 0:
            # Lossy: per-frequency Zs with R_mag
            freqs = np.asarray(freqs)
            Z_port = np.zeros(len(freqs), dtype=complex)
            for i, f in enumerate(freqs):
                Zs = Zs_func(f) if Zs_func is not None else None
                Z_port[i] = self.compute_port_impedance(f, Zs)
            return Z_port

        # Lossless: use C++ batch sweep (fast path)
        return super().frequency_sweep(freqs, Zs_func)

    def get_effective_inductance(self):
        """
        Get effective inductance matrix L_total = L_air + Delta_L.

        Returns:
            L_total: (n_seg x n_seg) inductance matrix [H]
        """
        # self.L already contains L_air + Delta_L after compute_coupling_matrix()
        return self.L.copy()

    def get_magnetic_loss_resistance(self, freq):
        """
        Get frequency-dependent magnetic loss resistance from mu"_r.

        R_mag(f) = omega * tan_delta_m * Delta_L

        Args:
            freq: Frequency in Hz

        Returns:
            R_mag: (n_seg x n_seg) magnetic loss resistance matrix [Ohm]
        """
        omega = 2.0 * np.pi * freq
        if self.Delta_L is not None and self.tan_delta_m > 0:
            return omega * self.tan_delta_m * self.Delta_L
        return np.zeros((self.n_loop, self.n_loop))
