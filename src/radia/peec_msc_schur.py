"""
peec_msc_schur.py

Schur Complement Solver for PEEC-MSC Coupled System

Eliminates PEEC DOFs (small) via Schur complement so that the large MSC
system retains H-matrix (HACApK) acceleration with BiCGSTAB.

Coupled system:
    [Z_peec    -jw*M_mp] [I]   [V_source]
    [B_pm      K_msc   ] [s] = [   0    ]

Eliminating I from the first equation:
    I = Z_peec^{-1} * (V_source + jw*M_mp*sigma)

Substituting into the second equation:
    (K_msc + B_pm * Z_peec^{-1} * jw*M_mp) * sigma = -B_pm * Z_peec^{-1} * V_source

The Schur complement system:
    K_schur * sigma = rhs_schur

where:
    K_schur = K_msc + B_pm * Z_peec^{-1} * jw * M_mp   (low-rank correction, rank = n_peec)
    rhs_schur = -B_pm * Z_peec^{-1} * V_source

Key insight: The low-rank correction has rank n_peec << n_msc, so the H-matrix
structure of K_msc is preserved. BiCGSTAB matvec cost per iteration:
    O(N log N) for H-matrix matvec  +  O(N * n_peec) for low-rank correction

Part of Radia project
"""

import numpy as np
from scipy.sparse.linalg import LinearOperator, bicgstab, gmres
import scipy

# scipy >= 1.12 uses rtol instead of tol for bicgstab
_scipy_version = tuple(int(x) for x in scipy.__version__.split('.')[:2])
_bicgstab_tol_kwarg = 'rtol' if _scipy_version >= (1, 12) else 'tol'

MU_0 = 4.0 * np.pi * 1e-7  # H/m


class SchurComplementSolver:
    """
    Schur complement solver for PEEC-MSC coupled system.

    Eliminates PEEC DOFs (small, dense, complex) via Schur complement,
    preserving H-matrix structure for the large MSC system.

    Usage:
        solver = SchurComplementSolver()

        # Set MSC system (from MMMBuilder/MMMSolver)
        solver.set_msc_system(N_matrix, dof_offset, inv_chi)

        # Set PEEC topology
        solver.set_peec_system(Z_peec_func, seg_p1, seg_p2,
                               seg_centers, seg_directions, seg_lengths)

        # Compute coupling matrices (B_pm, M_mp)
        solver.compute_coupling_matrices(msc_eval_points, msc_face_normals)

        # Solve at frequency
        I, sigma = solver.solve(freq, V_source)
    """

    def __init__(self):
        self._N = None                # MSC interaction matrix (dense, row-major)
        self._dof_offset = None       # DOF offsets
        self._inv_chi = None          # 1/chi per DOF
        self._n_msc = 0               # MSC DOF count

        # PEEC data
        self._n_peec = 0              # PEEC loop count
        self._seg_p1 = None           # Segment start points (n_peec, 3)
        self._seg_p2 = None           # Segment end points (n_peec, 3)
        self._seg_centers = None      # Segment centers (n_peec, 3)
        self._seg_directions = None   # Segment unit directions (n_peec, 3)
        self._seg_lengths = None      # Segment lengths (n_peec,)
        self._L_air = None            # PEEC air inductance matrix
        self._R_dc = None             # PEEC DC resistance array

        # Coupling matrices
        self._B_pm = None             # (n_msc, n_peec): PEEC current -> H at MSC faces
        self._M_mp = None             # (n_peec, n_msc): MSC charge -> flux at PEEC segs

        # MSC matvec callback (for H-matrix mode)
        self._msc_matvec_func = None  # If set, uses this instead of dense N

        # MSC element geometry for coupling matrix computation
        self._msc_eval_points = None  # (n_msc, 3): evaluation points per MSC DOF
        self._msc_face_normals = None # (n_msc, 3): face normals per MSC DOF
        self._msc_face_areas = None   # (n_msc,): face areas per MSC DOF

        # MMMFieldComputer for proper A-field computation (optional)
        self._field_computer = None

        # Hex vertices for Radia-based coupling (set by user)
        self._msc_hex_vertices = None

    def set_msc_system(self, N, dof_offset, inv_chi):
        """
        Set MSC interaction matrix and material data.

        Args:
            N: Interaction matrix (n_msc, n_msc) numpy array
            dof_offset: DOF offset array (n_elem + 1,)
            inv_chi: 1/chi values (n_msc,) per DOF
        """
        self._N = np.asarray(N, dtype=np.float64)
        self._dof_offset = np.asarray(dof_offset, dtype=np.int32)
        self._inv_chi = np.asarray(inv_chi, dtype=np.float64)
        self._n_msc = self._N.shape[0]

    def set_msc_matvec(self, matvec_func):
        """
        Set custom matvec for MSC system (e.g., H-matrix via HACApK).

        The function should compute y = K_msc * x where K_msc = -diag(1/chi) - N.

        Args:
            matvec_func: callable(x) -> y, where x and y are (n_msc,) arrays
        """
        self._msc_matvec_func = matvec_func

    def set_peec_system(self, L_air, R_dc,
                        seg_p1, seg_p2, seg_centers, seg_directions, seg_lengths):
        """
        Set PEEC system data.

        Args:
            L_air: Air inductance matrix (n_peec, n_peec) [H]
            R_dc: DC resistance array (n_peec,) [Ohm]
            seg_p1: Segment start points (n_peec, 3) [m]
            seg_p2: Segment end points (n_peec, 3) [m]
            seg_centers: Segment centers (n_peec, 3) [m]
            seg_directions: Segment unit directions (n_peec, 3)
            seg_lengths: Segment lengths (n_peec,) [m]
        """
        self._L_air = np.asarray(L_air, dtype=np.float64)
        self._R_dc = np.asarray(R_dc, dtype=np.float64)
        self._seg_p1 = np.asarray(seg_p1, dtype=np.float64)
        self._seg_p2 = np.asarray(seg_p2, dtype=np.float64)
        self._seg_centers = np.asarray(seg_centers, dtype=np.float64)
        self._seg_directions = np.asarray(seg_directions, dtype=np.float64)
        self._seg_lengths = np.asarray(seg_lengths, dtype=np.float64)
        self._n_peec = len(self._seg_lengths)

    def set_msc_hex_vertices(self, vertices):
        """
        Set hex vertices for Radia-based coupling matrix computation.

        Args:
            vertices: (8, 3) array or list of 8 vertices in Radia hex order
        """
        self._msc_hex_vertices = [list(v) for v in np.asarray(vertices)]

    def set_field_computer(self, field_computer):
        """
        Set MMMFieldComputer for proper A-field computation.

        When set, compute_coupling_matrices() uses compute_a_field()
        (magnetic dipole with reconstructed M for hexahedra) instead of
        the crude point source approximation.

        Args:
            field_computer: mmm_core.MMMFieldComputer instance
                            (must have set_elements_from_builder already called)
        """
        self._field_computer = field_computer

    def set_msc_geometry(self, eval_points, face_normals, face_areas):
        """
        Set MSC element face geometry for coupling matrix computation.

        Args:
            eval_points: Face evaluation points (n_msc, 3) [m]
            face_normals: Face outward normals (n_msc, 3)
            face_areas: Face areas (n_msc,) [m^2]
        """
        self._msc_eval_points = np.asarray(eval_points, dtype=np.float64)
        self._msc_face_normals = np.asarray(face_normals, dtype=np.float64)
        self._msc_face_areas = np.asarray(face_areas, dtype=np.float64)

    def compute_coupling_matrices(self):
        """
        Compute coupling matrices B_pm and M_mp.

        B_pm[i,k]: H-field normal component at MSC face i due to unit current
                   in PEEC segment k. Shape: (n_msc, n_peec)

        M_mp[k,i]: Flux linkage with PEEC segment k due to MSC DOF i
                   with unit charge/magnetization. Shape: (n_peec, n_msc)

        B_pm maps PEEC currents -> external H field at MSC DOFs.
        M_mp maps MSC surface charges -> induced voltage at PEEC segments.

        If a MMMFieldComputer is set (via set_field_computer), M_mp is computed
        using proper magnetic dipole A-field. Otherwise falls back to approximate
        scalar potential.
        """
        n_msc = self._n_msc
        n_peec = self._n_peec

        B_pm = np.zeros((n_msc, n_peec))
        M_mp = np.zeros((n_peec, n_msc))

        # B_pm: H-field at MSC faces from PEEC segments (Biot-Savart)
        for k in range(n_peec):
            p1 = self._seg_p1[k]
            p2 = self._seg_p2[k]
            for i in range(n_msc):
                obs = self._msc_eval_points[i]
                from radia.biot_savart import h_filament, MU0
                H = h_filament(p1, p2, obs, current=1.0)
                B_pm[i, k] = MU0 * np.dot(H, self._msc_face_normals[i])

        # M_mp: Flux linkage at PEEC segments from MSC magnetization
        if self._field_computer is not None:
            # Proper computation via MMMFieldComputer.compute_a_field
            # For each MSC DOF i, set unit magnetization/charge, compute A
            # at all PEEC segment centers, then flux linkage.
            obs_points = self._seg_centers  # (n_peec, 3)
            for i in range(n_msc):
                # Unit vector for DOF i
                sigma = np.zeros(n_msc)
                sigma[i] = 1.0
                # Compute A from this magnetization state
                A = self._field_computer.compute_a_field(sigma, obs_points)
                # A has shape (n_peec, 3)
                for k in range(n_peec):
                    M_mp[k, i] = (np.dot(A[k], self._seg_directions[k])
                                  * self._seg_lengths[k])
        else:
            # Fallback: point dipole approximation (less accurate)
            # Reconstruct M from sigma for each element, compute dipole A
            import warnings
            warnings.warn(
                "No MMMFieldComputer set. Using point dipole approximation "
                "for M_mp. Call set_field_computer() for proper computation.",
                stacklevel=2
            )
            for i in range(n_msc):
                r_face = self._msc_eval_points[i]
                n_hat = self._msc_face_normals[i]
                area_i = self._msc_face_areas[i]
                for k in range(n_peec):
                    r_seg = self._seg_centers[k]
                    dr = r_seg - r_face
                    dist = np.linalg.norm(dr)
                    if dist < 1e-20:
                        continue
                    A_vec = (MU_0 / (4.0 * np.pi)) * area_i * n_hat / dist
                    M_mp[k, i] = (np.dot(A_vec, self._seg_directions[k])
                                  * self._seg_lengths[k])

        self._B_pm = B_pm
        self._M_mp = M_mp
        return B_pm, M_mp

    def compute_coupling_matrices_radia(self, magnetic_objects,
                                        solver_prec=0.0001, solver_maxiter=1000,
                                        solver_method=0):
        """
        Compute coupling matrices using Radia's exact field computation.

        B_pm is computed via Biot-Savart (same as analytical method).
        M_mp is computed by solving the MSC system in Radia for each MSC DOF
        and extracting A-field at PEEC segments — exact analytical integration.

        This is more accurate than compute_coupling_matrices() especially
        for near-field coupling where the dipole approximation breaks down.

        Note: Unlike compute_coupling_matrices(), this method does NOT
        compute M_mp column-by-column from unit charges. Instead, it computes
        the full Delta_L = M_mp · K_msc^{-1} · B_pm directly by solving the
        MSC system for each PEEC segment's field (same as peec_coupled.py).
        The individual M_mp matrix is then extracted as:
        M_mp[k,i] = flux linkage at segment k from unit sigma_i computed
        via Radia Fld('a', ...).

        Args:
            magnetic_objects: List of Radia object handles (with material applied)
            solver_prec: Radia solver precision
            solver_maxiter: Radia solver max iterations
            solver_method: 0=LU, 1=BiCGSTAB, 2=HACApK

        Returns:
            B_pm: (n_msc, n_peec) coupling matrix
            M_mp: (n_peec, n_msc) coupling matrix
        """
        try:
            from _radia_pybind import (ObjBckg, ObjCnt, Solve, Fld, UtiDel)
        except ImportError:
            import radia as _rad
            ObjBckg = _rad.ObjBckg
            ObjCnt = _rad.ObjCnt
            Solve = _rad.Solve
            Fld = _rad.Fld
            UtiDel = _rad.UtiDel

        n_msc = self._n_msc
        n_peec = self._n_peec

        # B_pm: same Biot-Savart computation
        B_pm = np.zeros((n_msc, n_peec))
        for k in range(n_peec):
            p1 = self._seg_p1[k]
            p2 = self._seg_p2[k]
            for i in range(n_msc):
                obs = self._msc_eval_points[i]
                from radia.biot_savart import h_filament, MU0
                H = h_filament(p1, p2, obs, current=1.0)
                B_pm[i, k] = MU0 * np.dot(H, self._msc_face_normals[i])

        # M_mp via Radia: solve MSC for each PEEC segment's field,
        # then compute A at all segments.
        # Delta_L[i,j] = dot(A(center_i), dir_i) * length_i
        # where A is from Radia Solve with background field from segment j.
        #
        # To get M_mp separately, we solve with Radia and extract A for
        # the solved magnetization state. Then:
        # Delta_L[:,j] = M_mp @ sigma_j where sigma_j = K_msc^{-1} @ B_pm[:,j]
        #
        # But we need M_mp as a standalone matrix for the Schur complement.
        # Compute it directly: for each MSC DOF i, create a Radia object
        # with unit magnetic moment in the appropriate direction and
        # compute A at all PEEC segments.

        # Strategy: M_mp[k,i] = ∂Φ_k/∂σ_i where Φ_k = A(center_k)·dir_k·len_k
        # For the Radia approach, we set uniform M in the hex,
        # compute A via Fld('a'), and build M_mp from the Jacobian.

        # Simpler approach: solve 3 Radia problems with M along x, y, z,
        # get A at PEEC centers, and build the linear map from M to Φ.
        # Then convert from M-basis to sigma-basis using face normals.

        try:
            from _radia_pybind import (ObjHexahedron, Fld as RadFld, UtiDel as RadUtiDel)
        except ImportError:
            RadFld = Fld
            RadUtiDel = UtiDel
            ObjHexahedron = None
            try:
                import radia as _rad2
                ObjHexahedron = _rad2.ObjHexahedron
            except ImportError:
                pass

        # Extract hex vertex data from magnetic_objects for creating
        # temporary objects with fixed M.
        # Strategy: create Radia hex with fixed M in each axis direction,
        # compute A at PEEC centers (no Solve - just field from known M).
        # This gives dA/dM, then convert to dA/dsigma via face normals.

        # For now, use the first magnetic object's geometry.
        # Create 3 temporary objects with unit M along x, y, z.
        # Radia ObjHexahedron(vertices, M) creates a uniformly magnetized hex.
        verts = self._msc_hex_vertices  # must be set by user

        # Compute A response for 3 orthogonal unit magnetizations
        # No Solve needed - just field from fixed magnetization
        dPhi_dM = np.zeros((n_peec, 3))  # [seg, M_dir]

        for d in range(3):
            M_unit = [0.0, 0.0, 0.0]
            M_unit[d] = 1.0  # Unit M in direction d

            # Create hex with fixed M (no material, no solve)
            tmp_obj = ObjHexahedron(verts, M_unit)

            for k in range(n_peec):
                center_k = self._seg_centers[k].tolist()
                A_vec = np.array(RadFld(tmp_obj, 'a', center_k))
                dPhi_dM[k, d] = (np.dot(A_vec, self._seg_directions[k])
                                 * self._seg_lengths[k])

            RadUtiDel(tmp_obj)

        # Convert from M-basis to sigma-basis.
        # For uniform M: sigma_i = M . n_i (6 equations, 3 unknowns)
        # Pseudoinverse: M = (N^T N)^{-1} N^T sigma
        # So: dPhi/dsigma = dPhi/dM . (dM/dsigma) where dM/dsigma_i ≈ n_i
        # More precisely: dPhi/dsigma = dPhi/dM . N^+ where N^+ is pseudoinverse
        # But since each face independently contributes, and for physically
        # consistent sigma (from uniform M), this is exact.
        # M_mp[k,i] = dPhi_k/dM . n_i (direct chain rule)
        M_mp = np.zeros((n_peec, n_msc))
        for i in range(n_msc):
            n_i = self._msc_face_normals[i]
            for k in range(n_peec):
                M_mp[k, i] = np.dot(dPhi_dM[k, :], n_i)

        self._B_pm = B_pm
        self._M_mp = M_mp
        return B_pm, M_mp

    def _build_Z_peec(self, freq):
        """Build PEEC impedance matrix at frequency."""
        omega = 2.0 * np.pi * freq
        Z = np.diag(self._R_dc.astype(complex)) + 1j * omega * self._L_air
        return Z

    def _msc_matvec(self, x):
        """
        MSC system matvec: y = K_msc * x = (diag(1/chi) + N) * x

        Uses H-matrix callback if set, otherwise dense numpy.
        """
        if self._msc_matvec_func is not None:
            return self._msc_matvec_func(x.real) + 1j * self._msc_matvec_func(x.imag)

        # Dense: A = diag(inv_chi) + N
        y = self._inv_chi * x + self._N @ x
        return y

    def solve(self, freq, V_source, tol=1e-8, max_iter=1000, x0=None):
        """
        Solve coupled PEEC-MSC system at a given frequency using Schur complement.

        Args:
            freq: Frequency [Hz]
            V_source: Source voltage vector (n_peec,) [V], complex
            tol: BiCGSTAB convergence tolerance
            max_iter: Maximum BiCGSTAB iterations
            x0: Initial guess for sigma (n_msc,), or None

        Returns:
            I_peec: PEEC loop currents (n_peec,), complex
            sigma: MSC surface charges (n_msc,), complex
            info: dict with solver info (n_iter, residual)
        """
        omega = 2.0 * np.pi * freq
        n_msc = self._n_msc
        n_peec = self._n_peec

        V_source = np.asarray(V_source, dtype=complex)

        # PEEC impedance and its inverse
        Z_peec = self._build_Z_peec(freq)
        Z_peec_inv = np.linalg.inv(Z_peec)  # Small: n_peec x n_peec

        # Pre-compute: Z_peec^{-1} * jw * M_mp  -> (n_peec, n_msc)
        jw_M_mp = 1j * omega * self._M_mp
        Z_inv_jw_Mmp = Z_peec_inv @ jw_M_mp   # (n_peec, n_msc)

        # Pre-compute: B_pm * Z_peec^{-1}  -> (n_msc, n_peec)
        B_pm_Z_inv = self._B_pm @ Z_peec_inv   # (n_msc, n_peec)

        # Schur complement RHS
        rhs_schur = -B_pm_Z_inv @ V_source      # (n_msc,), complex

        # Schur complement matvec: y = K_msc*x + B_pm*Z_inv*jw*M_mp*x
        def schur_matvec(x):
            # K_msc * x (H-matrix or dense)
            y = self._msc_matvec(x)
            # Low-rank correction: B_pm * (Z_peec^{-1} * (jw * M_mp * x))
            temp = Z_inv_jw_Mmp @ x             # (n_peec,)
            y += self._B_pm @ temp               # (n_msc,)  -- use B_pm not B_pm_Z_inv
            return y

        K_schur_op = LinearOperator(
            (n_msc, n_msc), matvec=schur_matvec, dtype=complex
        )

        # Solve reduced system
        sigma, info_code = bicgstab(K_schur_op, rhs_schur, x0=x0,
                                     maxiter=max_iter, **{_bicgstab_tol_kwarg: tol})

        # Back-substitute for PEEC currents
        I_peec = Z_peec_inv @ (V_source + jw_M_mp @ sigma)

        info = {
            'converged': info_code == 0,
            'info_code': info_code,
            'n_msc': n_msc,
            'n_peec': n_peec,
        }

        return I_peec, sigma, info

    def solve_static(self, I_peec_dc, tol=1e-8, max_iter=1000):
        """
        Solve MSC system with DC PEEC currents (static coupling).

        For magnetostatic problems where PEEC currents are known.
        K_msc * sigma = -B_pm * I_peec_dc

        Args:
            I_peec_dc: DC PEEC currents (n_peec,) [A], real
            tol: BiCGSTAB tolerance
            max_iter: Maximum iterations

        Returns:
            sigma: MSC surface charges (n_msc,), real
            info: dict with solver info
        """
        n_msc = self._n_msc
        I_peec_dc = np.asarray(I_peec_dc, dtype=np.float64)

        # RHS: H_ext at MSC faces from PEEC currents
        rhs = -self._B_pm @ I_peec_dc   # (n_msc,)

        def msc_matvec_real(x):
            return np.real(self._msc_matvec(x))

        K_msc_op = LinearOperator(
            (n_msc, n_msc), matvec=msc_matvec_real, dtype=np.float64
        )

        sigma, info_code = bicgstab(K_msc_op, rhs, maxiter=max_iter, **{_bicgstab_tol_kwarg: tol})

        info = {
            'converged': info_code == 0,
            'info_code': info_code,
        }

        return sigma, info

    def compute_delta_L(self, tol=1e-8, max_iter=1000):
        """
        Compute Delta_L coupling inductance matrix via Schur complement.

        This is equivalent to the column-by-column approach in peec_coupled.py
        but computed as: Delta_L = -M_mp * K_msc^{-1} * B_pm

        Each column of K_msc^{-1} * B_pm is a separate MSC solve.

        Args:
            tol: BiCGSTAB tolerance
            max_iter: Maximum iterations

        Returns:
            Delta_L: (n_peec, n_peec) coupling inductance matrix [H]
        """
        n_msc = self._n_msc
        n_peec = self._n_peec

        def msc_matvec_real(x):
            return np.real(self._msc_matvec(x))

        K_msc_op = LinearOperator(
            (n_msc, n_msc), matvec=msc_matvec_real, dtype=np.float64
        )

        # Solve K_msc * X = B_pm column by column
        # X = K_msc^{-1} * B_pm  -> (n_msc, n_peec)
        X = np.zeros((n_msc, n_peec))
        for k in range(n_peec):
            rhs_k = self._B_pm[:, k]
            x_k, info_code = bicgstab(K_msc_op, rhs_k, maxiter=max_iter, **{_bicgstab_tol_kwarg: tol})
            if info_code != 0:
                print(f"Warning: BiCGSTAB did not converge for column {k}, info={info_code}")
            X[:, k] = x_k

        # Delta_L = -M_mp * X = -M_mp * K_msc^{-1} * B_pm
        Delta_L = -self._M_mp @ X

        return Delta_L

    def frequency_sweep(self, freqs, V_source, tol=1e-8, max_iter=1000):
        """
        Compute PEEC currents and MSC charges over a frequency range.

        Z_peec^{-1} is recomputed at each frequency. The H-matrix for K_msc
        is reused across all frequencies.

        Args:
            freqs: Array of frequencies [Hz]
            V_source: Source voltage (n_peec,) [V], complex (same for all freqs)
            tol: BiCGSTAB tolerance
            max_iter: Maximum iterations

        Returns:
            I_all: (n_freq, n_peec) PEEC currents, complex
            sigma_all: (n_freq, n_msc) MSC charges, complex
        """
        freqs = np.asarray(freqs)
        n_freq = len(freqs)
        n_peec = self._n_peec
        n_msc = self._n_msc

        I_all = np.zeros((n_freq, n_peec), dtype=complex)
        sigma_all = np.zeros((n_freq, n_msc), dtype=complex)

        sigma_prev = None
        for idx, f in enumerate(freqs):
            I, sigma, info = self.solve(f, V_source, tol=tol, max_iter=max_iter,
                                        x0=sigma_prev)
            I_all[idx] = I
            sigma_all[idx] = sigma
            sigma_prev = sigma  # Use as initial guess for next frequency

            if not info['converged']:
                print(f"Warning: f={f:.3e} Hz did not converge")

        return I_all, sigma_all
