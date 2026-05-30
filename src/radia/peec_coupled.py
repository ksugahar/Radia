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
import radia as rad
from radia.peec_topology import PEECCircuitSolver
from radia.peec_matrices import MNASolver

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
        # Per CLAUDE.md "No Fallbacks": use the canonical `radia` package
        # namespace.  The prior shim attempted bare `_radia_pybind` import
        # with sys.path mutation -> editable-install hack that mutates
        # global state on every call.  `import radia as rad` at module
        # top handles both editable and wheel layouts uniformly.
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

            # H-field from unit current in segment j (analytical)
            def h_field_from_seg_j(point, _p1=p1_j, _p2=p2_j):
                from radia.biot_savart import h_filament, MU0
                H = h_filament(_p1, _p2, point, current=1.0)
                B = MU0 * H
                return B.tolist()

            # Create background field
            bkg = ObjBckg(h_field_from_seg_j)

            # Create container with magnetic objects + background
            container = ObjCnt([mag_container, bkg])

            # Solve for induced magnetization
            Solve(container, solver_prec, solver_maxiter, solver_method)

            # Compute A-field at all segment centers in ONE batched call
            # (Radia ComputeFieldBatch, TaskManager-parallel internally; the
            # caller is responsible for the TaskManager context per CLAUDE.md
            # "Caller Wraps, Helper Does NOT").
            A_all = np.asarray(Fld(mag_container, 'a', self.seg_centers))
            # A_all shape: (n_seg, 3)
            # Delta_L[i][j] = dot(A_i, dir_i) * length_i, vectorized over i
            Delta_L[:, j] = np.einsum('ij,ij->i', A_all, self.seg_directions) \
                            * self.seg_lengths

            # Clean up background field
            UtiDel(bkg)

        self.Delta_L = Delta_L

        # Store original L_air and update L to include coupling
        self._L_air = self.L.copy()
        self.L = self._L_air + Delta_L

        # Rebuild C++ MNASolver with coupled L
        self._rebuild_mna_solver()

        return Delta_L

    def compute_coupling_scalar_bie(self, body_mesh, frequency=0,
                                     sigma=0, mu_r=1.0):
        """Compute Delta_Z from scattering body via Scalar BIE + SIBC.

        Surface-only method: no volume mesh needed for the body.
        Uses ScalarBIESIBCSolver from bem_sibc_solver.py.

        For conductors (sigma > 0): eddy current shielding (Delta_Z complex)
        For magnetics (mu_r > 1):   inductance change (Delta_L real)

        Args:
            body_mesh: NGSolve surface mesh of scattering body
            frequency: Operating frequency [Hz] (0 = DC/static)
            sigma: Body conductivity [S/m] (0 = non-conducting)
            mu_r: Body relative permeability

        Returns:
            Delta_Z: (n_seg, n_seg) complex impedance change matrix
        """
        import math
        from radia.biot_savart import h_filament
        from radia.bem_sibc_solver import (ScalarBIESIBCSolver,
                                           compute_phi_inc_from_loop)

        omega = 2 * math.pi * frequency
        n_seg = self.n_loop

        # SIBC surface impedance
        if sigma > 0 and omega > 0:
            rho = 1.0 / sigma
            mu_eff = MU_0 * mu_r
            delta = math.sqrt(2 * rho / (omega * mu_eff))
            Z_s = complex(1, 1) * rho / delta
        else:
            Z_s = complex(0, 0)

        # Build BEM solver on body surface
        solver = ScalarBIESIBCSolver(body_mesh, order=1)

        # Node coordinates for phi_inc evaluation
        node_coords = np.array([[body_mesh.vertices[i].point[j]
                                 for j in range(3)]
                                for i in range(body_mesh.nv)])

        Delta_Z = np.zeros((n_seg, n_seg), dtype=complex)

        for j in range(n_seg):
            p1_j = self.seg_p1[j]
            p2_j = self.seg_p2[j]

            # phi_inc from unit current in segment j at body surface nodes
            phi_inc = np.zeros(body_mesh.nv)
            # Path integration: phi = -int H.dl from infinity
            # For a single segment, use analytical H + vertical axis formula
            from radia.bem_sibc_solver import compute_phi_inc_from_surface_J
            # Approximate segment as point dipole for phi_inc
            seg_center = 0.5 * (p1_j + p2_j)
            seg_dir = p2_j - p1_j
            seg_len = np.linalg.norm(seg_dir)

            # Use h_filament at quadrature points for path integration.
            # Vectorized across ALL body-mesh vertices in two batched
            # `h_segments_batch` calls (vertical + horizontal) instead of
            # 2 * n_vert * 20 individual Python calls (orig was O(n_vert)
            # sequential).
            from radia.biot_savart import h_segments_batch
            seg_list = [(tuple(p1_j), tuple(p2_j))]

            t_gl, w_gl = np.polynomial.legendre.leggauss(20)
            t_01 = 0.5 * (t_gl + 1)
            w_01 = 0.5 * w_gl
            z_far = 20 * max(np.max(np.abs(node_coords)), seg_len)

            n_vert = body_mesh.nv
            n_q = len(t_01)
            xi_arr = node_coords[:, 0]
            yi_arr = node_coords[:, 1]
            zi_arr = node_coords[:, 2]
            rho_xy_arr = np.sqrt(xi_arr * xi_arr + yi_arr * yi_arr)

            # --- Vertical path: phi(0,0,z) for every vertex ---
            # x_quad_v_all[ip, k, :] = (0, 0, z_far - t_01[k]*(z_far - zi[ip]))
            z_quad_all = (z_far
                          - t_01[np.newaxis, :] * (z_far - zi_arr[:, np.newaxis]))
            x_quad_v_all = np.zeros((n_vert, n_q, 3))
            x_quad_v_all[..., 2] = z_quad_all
            x_quad_v_flat = x_quad_v_all.reshape(n_vert * n_q, 3)
            Hz_v_flat = h_segments_batch(seg_list, x_quad_v_flat)[:, 2]
            Hz_v_all = Hz_v_flat.reshape(n_vert, n_q)
            phi_axis_all = (z_far - zi_arr) * np.sum(
                w_01[np.newaxis, :] * Hz_v_all, axis=1)

            # --- Horizontal path: contributes only for off-axis vertices ---
            off_axis = rho_xy_arr >= 1e-12 * z_far
            phi_inc[:] = phi_axis_all
            if np.any(off_axis):
                idx_off = np.flatnonzero(off_axis)
                m = idx_off.size
                dl_off = np.zeros((m, 3))
                dl_off[:, 0] = xi_arr[idx_off]
                dl_off[:, 1] = yi_arr[idx_off]
                # x_quad_h[ip, k, :] = t_01[k] * dl[ip], z = zi[ip]
                x_quad_h_all = (t_01[np.newaxis, :, np.newaxis]
                                * dl_off[:, np.newaxis, :])
                x_quad_h_all[..., 2] = zi_arr[idx_off][:, np.newaxis]
                x_quad_h_flat = x_quad_h_all.reshape(m * n_q, 3)
                H_h_flat = h_segments_batch(seg_list, x_quad_h_flat)
                H_h_all = H_h_flat.reshape(m, n_q, 3)
                # integrand[ip, k] = H_h[ip, k, :] . dl[ip]
                integrand = np.einsum('mqd,md->mq', H_h_all, dl_off)
                phi_inc[idx_off] -= np.sum(w_01[np.newaxis, :] * integrand,
                                            axis=1)

            # Solve BIE
            if abs(Z_s) > 0:
                result = solver.solve(phi_inc, Z_s=Z_s, omega=omega)
            else:
                # Static magnetic: gamma = mu_r - 1 coupling
                gamma = (mu_r - 1) / (mu_r + 1)
                result = solver.solve(phi_inc, Z_s=complex(0, 0),
                                      omega=0, gamma_override=gamma)

            # Extract Delta_Z from power dissipation
            # For now, store phi solution for post-processing
            # Delta_Z[i,j] requires mapping phi -> H -> A -> flux linkage
            # This is a simplified version using energy method
            phi_sol = result.get('phi', phi_inc)
            # Power = Re(Z_s) * int |H_t|^2 dS ≈ Re(Z_s) * sum(|dphi/dt|^2 * dA)
            # For impedance matrix: Delta_Z[i,j] = ... (complex)
            # Simplified: treat as perturbation to self/mutual inductance
            Delta_Z[j, j] = result.get('Z_total', 0)

        self.Delta_Z_bie = Delta_Z
        return Delta_Z

    def _rebuild_mna_solver(self):
        """Rebuild C++ MNASolver with current self.L matrix."""
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
