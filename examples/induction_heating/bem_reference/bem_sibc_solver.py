"""
Scalar potential BIE with SIBC for MQS eddy current on conducting surfaces.

Solves for the scalar magnetic potential phi on a conductor surface:
    (1/2*M - DL + gamma * SL * M^{-1} * K) phi = rhs
    gamma = Z_s / (jw * mu_0)

where M = H1 surface mass, K = H1 surface stiffness (Laplace-Beltrami),
DL = Laplace double layer, SL = Laplace single layer (all from ngsolve.bem).

Usage:
    from radia.bem_sibc_solver import ScalarBIESIBCSolver

    solver = ScalarBIESIBCSolver(mesh, order=1)
    result = solver.solve(phi_inc_cf, Z_s=Z_s, omega=omega)

    J_rms = result['H_t_rms']       # Surface current RMS [A/m]
    P = result['P_density']          # Power loss density [W/m^2]
    phi_gf = result['phi']           # Solved potential GridFunction
"""

import math
import time
import numpy as np
from scipy.linalg import solve as scipy_solve

MU_0 = 4e-7 * np.pi


class ScalarBIESIBCSolver:
    """Scalar potential BIE + SIBC solver for conducting surfaces.

    The solver assembles BEM operators once, then can be called repeatedly
    with different Z_s (frequency sweep) or phi_inc (different sources).
    """

    def __init__(self, mesh, order=1, assemble_dense=True):
        """Initialize solver and assemble BEM operators.

        Args:
            mesh: NGSolve Mesh (surface mesh, dim=2, or volume mesh with BND).
            order: H1 polynomial order on surface (default 1).
            assemble_dense: if True (default, backward-compat), extract
                ``DL`` and ``SL`` to dense ``ndof x ndof`` numpy arrays
                via N column matvecs.  This is O(N^3) total (~67 min for
                N=2474) but lets ``solve()`` use a single dense scipy
                solve.
                If False, KEEP ``DL_bf`` and ``SL_bf`` as NGSolve
                bilinear forms only and use ``solve_iterative()`` (GMRES
                with LinearOperator wrappers).  ~50x faster for typical
                IH wp meshes -- recommended for new code.
        """
        from ngsolve import (H1, BilinearForm, GridFunction, ds, grad,
                             TaskManager, InnerProduct)
        from ngsolve.bem import LaplaceDL, LaplaceSL

        self.mesh = mesh
        self.order = order

        # H1 on surface
        self.fes = H1(mesh, order=order)
        u, v = self.fes.TnT()
        self.ndof = self.fes.ndof

        t0 = time.perf_counter()

        # BEM operators.  use_fmm=True selects FMM (Fast Multipole)
        # matvec, dropping per-matvec cost from O(N^2) to O(N log N).
        # For typical IH wp meshes (N=2000-5000) this is 100-1000x
        # faster than dense matvec, both for the legacy column-by-
        # column extraction (assemble_dense=True) and the matrix-free
        # GMRES path (solve_iterative).  Verified 2026-04-29: 67 min
        # assembly -> seconds with FMM on the same 2474-dof wp.
        with TaskManager():
            DL_bf = LaplaceDL(u.Trace() * ds, use_fmm=True) * v.Trace() * ds
            SL_bf = LaplaceSL(u.Trace() * ds, use_fmm=True) * v.Trace() * ds
        # Always retain the bilinear forms so solve_iterative() can
        # call DL_bf.mat * vec without re-assembling.
        self._DL_bf = DL_bf
        self._SL_bf = SL_bf

        ndof = self.ndof
        if assemble_dense:
            self.DL = np.zeros((ndof, ndof))
            self.SL = np.zeros((ndof, ndof))
            for j in range(ndof):
                ej = GridFunction(self.fes)
                ej.vec[:] = 0
                ej.vec[j] = 1.0
                r1 = ej.vec.CreateVector()
                r1.data = DL_bf.mat * ej.vec
                self.DL[:, j] = r1.FV().NumPy().copy()
                r2 = ej.vec.CreateVector()
                r2.data = SL_bf.mat * ej.vec
                self.SL[:, j] = r2.FV().NumPy().copy()
        else:
            # Dense matrices are NOT extracted -- caller must use
            # solve_iterative().  These attributes are set None so that
            # any accidental call to solve() (which expects dense .DL
            # and .SL) fails loudly instead of producing garbage.
            self.DL = None
            self.SL = None

        # Surface mass M
        mass_bf = BilinearForm(self.fes)
        mass_bf += u.Trace() * v.Trace() * ds
        mass_bf.Assemble()
        self.M = np.zeros((ndof, ndof))
        rows, cols, vals = mass_bf.mat.COO()
        for r_, c_, val in zip(rows, cols, vals):
            self.M[int(r_), int(c_)] = val

        # Surface stiffness K (Laplace-Beltrami)
        stiff_bf = BilinearForm(self.fes)
        stiff_bf += InnerProduct(grad(u).Trace(), grad(v).Trace()) * ds
        stiff_bf.Assemble()
        self.K = np.zeros((ndof, ndof))
        rows, cols, vals = stiff_bf.mat.COO()
        for r_, c_, val in zip(rows, cols, vals):
            self.K[int(r_), int(c_)] = val

        self.M_inv = np.linalg.inv(self.M)

        # Gauge vector: <1, v>_S for Lagrange multiplier
        self._c_gauge = self.M @ np.ones(ndof)

        self.t_assembly = time.perf_counter() - t0

    def solve(self, phi_inc_cf, Z_s, omega):
        """Solve scalar BIE + SIBC for given incident potential and impedance.

        Args:
            phi_inc_cf: NGSolve CoefficientFunction for incident scalar potential
                        on the surface, OR ndarray of length ndof (nodal values).
            Z_s: Surface impedance.
                 - **complex scalar**: legacy global SIBC (uniform Z_s
                   over the workpiece)
                 - **ndarray of length ndof (complex)**: per-node Z_s
                   (per-panel curvature SIBC; build the array from
                   panel-level Z_s values via vertex averaging or any
                   other H1 projection).
            omega: Angular frequency [rad/s].

        Returns:
            dict with keys:
                phi: GridFunction (complex) - solved surface potential
                phi_vec: ndarray (complex) - coefficient vector
                H_t_rms: float - RMS tangential H [A/m] (= surface current density)
                P_density: float - time-averaged power loss density [W/m^2]
                gamma: complex - Z_s / (jw * mu_0) (or its mean if per-node)
                t_solve: float - solve time [s]
        """
        from ngsolve import (LinearForm, GridFunction, Integrate, CF, ds,
                             grad, BND, InnerProduct)

        t0 = time.perf_counter()
        ndof = self.ndof

        # RHS: <phi_inc, v>_S
        if isinstance(phi_inc_cf, np.ndarray):
            rhs_vec = self.M @ phi_inc_cf
        else:
            v_h1 = self.fes.TestFunction()
            lf = LinearForm(self.fes)
            lf += phi_inc_cf * v_h1.Trace() * ds
            lf.Assemble()
            rhs_vec = lf.vec.FV().NumPy().copy()

        # System matrix.
        # Original global-Z_s formulation:
        #   gamma = Z_s / (jw * mu_0)         (scalar)
        #   A = 1/2 M - DL + gamma * (SL M^{-1} K)
        # Per-node Z_s formulation:
        #   Gamma_ii = Z_s_i / (jw * mu_0)    (diagonal)
        #   The Robin term gamma * Delta_s phi has the per-node coefficient
        #   on the test side: <gamma_i grad_s phi, grad_s v> -> the SL M^-1 K
        #   block is left-multiplied by diag(gamma) (per-row scaling).
        # That is, the discrete form (1/2 M - DL) phi + diag(gamma) (SL M^-1 K)
        # phi = M phi_inc.
        if isinstance(Z_s, np.ndarray):
            if Z_s.shape != (ndof,):
                raise ValueError(
                    f"Per-node Z_s must have shape ({ndof},), got {Z_s.shape}")
            if omega <= 0:
                gamma_vec = np.zeros(ndof, dtype=complex)
            else:
                gamma_vec = Z_s.astype(complex) / (1j * omega * MU_0)
            gamma_for_log = complex(np.mean(gamma_vec))
            # diag(gamma_vec) @ M = gamma_vec[:, None] * M (row-scaling)
            robin_block = gamma_vec[:, None] * (self.SL @ self.M_inv @ self.K)
            A_sys = (0.5 * self.M - self.DL + robin_block).astype(complex)
        else:
            gamma = Z_s / (1j * omega * MU_0) if omega > 0 and Z_s != 0 else 0
            gamma_for_log = complex(gamma)
            A_sys = (0.5 * self.M - self.DL
                     + gamma * self.SL @ self.M_inv @ self.K).astype(complex)

        # Solve with gauge (Lagrange multiplier for int phi dS = 0)
        phi_vec = self._solve_with_gauge(A_sys, rhs_vec.astype(complex))
        t_solve = time.perf_counter() - t0

        # Extract H_t_rms = sqrt(<|grad_s phi|^2> / area)
        gf = GridFunction(self.fes)
        gf.vec.FV().NumPy()[:] = phi_vec.real
        Hsq_re = Integrate(InnerProduct(grad(gf), grad(gf)),
                           self.mesh, BND)
        gf.vec.FV().NumPy()[:] = phi_vec.imag
        Hsq_im = Integrate(InnerProduct(grad(gf), grad(gf)),
                           self.mesh, BND)
        area = Integrate(CF(1), self.mesh, BND)
        H_t_rms = math.sqrt((abs(Hsq_re) + abs(Hsq_im)) / abs(area))

        # Power density: P' = (1/2) Re(Z_s) |J_s|^2 = (1/2) Re(Z_s) H_t_rms^2
        # (time-averaged). For per-node Z_s, the area-averaged Re(Z_s) is
        # the right scalar to report; the caller already integrates the
        # local power per panel via the panel-level R if needed.
        if isinstance(Z_s, np.ndarray):
            Z_s_avg_re = float(np.mean(Z_s.real))
            P_density = 0.5 * Z_s_avg_re * H_t_rms ** 2
        else:
            P_density = 0.5 * Z_s.real * H_t_rms ** 2 if Z_s != 0 else 0

        # GridFunction output
        gf_phi = GridFunction(self.fes)
        gf_phi.vec.FV().NumPy()[:] = phi_vec.real  # real part

        return {
            'phi': gf_phi,
            'phi_vec': phi_vec,
            'H_t_rms': float(H_t_rms),
            'P_density': float(P_density),
            'area': float(abs(area)),
            'gamma': gamma_for_log,
            't_solve': round(t_solve, 3),
        }

    def _solve_with_gauge(self, A_mat, rhs):
        """Solve A*phi = rhs with gauge constraint int(phi)dS = 0."""
        n = len(rhs)
        A_aug = np.zeros((n + 1, n + 1), dtype=complex)
        A_aug[:n, :n] = A_mat
        A_aug[:n, n] = self._c_gauge
        A_aug[n, :n] = self._c_gauge
        rhs_aug = np.zeros(n + 1, dtype=complex)
        rhs_aug[:n] = rhs
        return scipy_solve(A_aug, rhs_aug)[:n]

    def solve_iterative(self, phi_inc_cf, Z_s, omega, *,
                         tol=1e-8, maxiter=200, restart=50):
        """GMRES-based solve that AVOIDS dense DL/SL extraction.

        Drop-in replacement for ``solve()`` with ~50x speedup on the
        typical IH wp mesh (N=2000-5000 dofs).  The expensive
        column-by-column matvec extraction in ``__init__`` is skipped
        (you pass ``assemble_dense=False``) and the BEM operators
        ``DL_bf`` and ``SL_bf`` are wrapped as scipy ``LinearOperator``
        objects for matrix-free GMRES.

        Math is identical to ``solve()``:
            A_sys = (1/2) M - DL + gamma * SL @ M^{-1} @ K
            A_sys @ phi = M @ phi_inc
            gauge: int(phi) dS = 0  via Lagrange multiplier
        Same convergence properties; the iterative path just doesn't
        store ``A_sys`` as dense numpy.
        """
        from scipy.sparse.linalg import LinearOperator, gmres
        from ngsolve import (LinearForm, GridFunction, BND, ds,
                              Integrate, CF, InnerProduct, grad)

        t0 = time.perf_counter()
        ndof = self.ndof

        # RHS: <phi_inc, v>_S
        if isinstance(phi_inc_cf, np.ndarray):
            rhs_vec = self.M @ phi_inc_cf
        else:
            v_h1 = self.fes.TestFunction()
            lf = LinearForm(self.fes)
            lf += phi_inc_cf * v_h1.Trace() * ds
            lf.Assemble()
            rhs_vec = lf.vec.FV().NumPy().copy()

        # Per-node Z_s NOT supported in iterative path yet (the dense
        # solve has its own row-scaling logic).
        if isinstance(Z_s, np.ndarray):
            raise NotImplementedError(
                "solve_iterative does not support per-node Z_s yet -- "
                "use solve() with assemble_dense=True for that case.")

        gamma = (Z_s / (1j * omega * MU_0)
                 if (omega > 0 and Z_s != 0) else complex(0))

        # Helper: complex matvec on a real bem bilinear form.  NGSolve
        # bem operators are real-valued, so we dispatch real and imag
        # parts separately and recombine.
        fes = self.fes
        DL_bf = self._DL_bf
        SL_bf = self._SL_bf

        def _bem_matvec(bf, x_complex):
            gf_re = GridFunction(fes)
            gf_im = GridFunction(fes)
            gf_re.vec.FV().NumPy()[:] = x_complex.real
            gf_im.vec.FV().NumPy()[:] = x_complex.imag
            out_re = gf_re.vec.CreateVector()
            out_im = gf_im.vec.CreateVector()
            out_re.data = bf.mat * gf_re.vec
            out_im.data = bf.mat * gf_im.vec
            return (out_re.FV().NumPy().astype(complex)
                    + 1j * out_im.FV().NumPy().astype(complex))

        # A_sys @ x  (matrix-free)
        def A_matvec(x):
            x = np.asarray(x, dtype=complex).ravel()
            t1 = 0.5 * (self.M @ x)
            t2 = _bem_matvec(DL_bf, x)
            kx = self.K @ x
            mkx = self.M_inv @ kx
            slmkx = _bem_matvec(SL_bf, mkx)
            return t1 - t2 + gamma * slmkx

        # Augmented system to enforce int(phi) dS = 0:
        #   [ A    c ] [phi]   [b]
        #   [ c^T  0 ] [lam] = [0]
        c_gauge = self._c_gauge.astype(complex)

        def A_aug_matvec(x_aug):
            x = x_aug[:ndof]
            lam = x_aug[ndof]
            top = A_matvec(x) + lam * c_gauge
            bot = c_gauge @ x
            out = np.empty(ndof + 1, dtype=complex)
            out[:ndof] = top
            out[ndof] = bot
            return out

        A_op = LinearOperator((ndof + 1, ndof + 1),
                               matvec=A_aug_matvec, dtype=complex)

        rhs_aug = np.zeros(ndof + 1, dtype=complex)
        rhs_aug[:ndof] = rhs_vec

        sol_aug, info = gmres(A_op, rhs_aug, rtol=tol,
                               maxiter=maxiter, restart=restart)
        phi_vec = sol_aug[:ndof]

        t_solve = time.perf_counter() - t0

        if info > 0:
            print(f"[bem_sibc] WARN: GMRES did not converge in {info} iters")

        # Extract H_t_rms / P_density (same as solve())
        gf = GridFunction(self.fes)
        gf.vec.FV().NumPy()[:] = phi_vec.real
        Hsq_re = Integrate(InnerProduct(grad(gf), grad(gf)),
                            self.mesh, BND)
        gf.vec.FV().NumPy()[:] = phi_vec.imag
        Hsq_im = Integrate(InnerProduct(grad(gf), grad(gf)),
                            self.mesh, BND)
        area = Integrate(CF(1), self.mesh, BND)
        H_t_rms = math.sqrt((abs(Hsq_re) + abs(Hsq_im)) / abs(area))
        P_density = (0.5 * Z_s.real * H_t_rms ** 2
                     if Z_s != 0 else 0.0)

        gf_phi = GridFunction(self.fes)
        gf_phi.vec.FV().NumPy()[:] = phi_vec.real
        return {
            'phi': gf_phi,
            'phi_vec': phi_vec,
            'H_t_rms': float(H_t_rms),
            'P_density': float(P_density),
            'area': float(abs(area)),
            'gamma': complex(gamma),
            't_solve': round(t_solve, 3),
            'gmres_info': int(info),
        }



def compute_phi_inc_from_loop(obs_points, loop_center, loop_radius, current,
                              n_quad=30, gap_deg=0):
    """Compute incident scalar magnetic potential from a circular current loop.

    Uses path integration: phi(P) = phi_axis(z) - int_axis^P H.dl
    H field from radia.biot_savart.h_segments (analytical formula).
    """
    from radia.biot_savart import h_segments_batch

    obs = np.asarray(obs_points, dtype=float)
    center = np.asarray(loop_center, dtype=float)
    a = float(loop_radius)
    I = float(current)
    if obs.ndim == 1:
        obs = obs.reshape(1, 3)

    obs_local = obs - center[np.newaxis, :]

    # Build coil wire segments as array (N_seg, 2, 3)
    arc_deg = 360 - gap_deg
    n_seg = max(200, int(arc_deg))
    theta = np.linspace(0, np.radians(arc_deg), n_seg + 1)
    coil_segs = np.zeros((n_seg, 2, 3))
    coil_segs[:, 0, 0] = a * np.cos(theta[:-1])
    coil_segs[:, 0, 1] = a * np.sin(theta[:-1])
    coil_segs[:, 1, 0] = a * np.cos(theta[1:])
    coil_segs[:, 1, 1] = a * np.sin(theta[1:])

    frac = arc_deg / 360.0

    # Gauss-Legendre quadrature for horizontal integration
    t_gl, w_gl = np.polynomial.legendre.leggauss(n_quad)
    t_01 = 0.5 * (t_gl + 1)
    w_01 = 0.5 * w_gl

    n_pts = len(obs_local)
    phi = np.zeros(n_pts)

    xy = obs_local[:, :2]
    z_arr = obs_local[:, 2]
    rho = np.sqrt(xy[:, 0]**2 + xy[:, 1]**2)

    # Analytical phi on z-axis for circular loop
    r_za = np.sqrt(z_arr**2 + a * a)
    phi_axis = (I / 2.0) * (1.0 - z_arr / r_za) * frac

    on_axis = rho < 1e-12 * a
    phi[on_axis] = phi_axis[on_axis]

    off_axis = ~on_axis
    if np.any(off_axis):
        idx_off = np.where(off_axis)[0]
        n_off = len(idx_off)

        # Build all quadrature points at once: (n_off * n_quad, 3)
        dl_vecs = np.zeros((n_off, 3))
        dl_vecs[:, 0] = obs_local[idx_off, 0]
        dl_vecs[:, 1] = obs_local[idx_off, 1]

        all_quad = np.zeros((n_off * n_quad, 3))
        for iq in range(n_quad):
            sl = slice(iq * n_off, (iq + 1) * n_off)
            all_quad[sl, 0] = t_01[iq] * dl_vecs[:, 0]
            all_quad[sl, 1] = t_01[iq] * dl_vecs[:, 1]
            all_quad[sl, 2] = z_arr[idx_off]

        # Vectorized H-field at all quadrature points
        H_all = h_segments_batch(coil_segs, all_quad, current=I)

        # Dot with dl_vec and integrate
        for iq in range(n_quad):
            sl = slice(iq * n_off, (iq + 1) * n_off)
            H_dot_dl = (H_all[sl, 0] * dl_vecs[:, 0]
                        + H_all[sl, 1] * dl_vecs[:, 1])
            phi[idx_off] -= w_01[iq] * H_dot_dl

        phi[idx_off] += phi_axis[idx_off]

    return phi



def compute_phi_inc_from_surface_J(obs_points, src_centroids, src_areas,
                                    src_J_vecs, n_quad=20, chunk_size=512):
    """Compute phi_inc from solved surface current via path integration.

    H field computed as sum of dipole contributions from surface elements.
    Fully vectorised over both the axis sweep (Stage 1) and the
    horizontal path sweep (Stage 2). For ih_sample.vol (170 obs nodes,
    4006 source panels, n_quad=20) this brings the runtime from ~160 s
    (per-node loop) down to ~1 s.

    Args:
        obs_points: (N, 3) observation node coordinates.
        src_centroids: (M, 3) source panel centroids.
        src_areas: (M,) source panel areas.
        src_J_vecs: (M, 3) source panel surface current densities (A/m).
        n_quad: Gauss-Legendre order for the path integrals.
        chunk_size: maximum number of obs nodes processed per chunk
            (memory cap; the inner H field array is
            chunk_size * n_quad * M * 3 floats).

    Returns:
        (N,) phi_inc values.
    """
    obs = np.asarray(obs_points, dtype=float)
    if obs.ndim == 1:
        obs = obs.reshape(1, 3)

    centers = np.asarray(src_centroids, dtype=float)
    areas = np.asarray(src_areas, dtype=float)
    J = np.asarray(src_J_vecs, dtype=float)

    INV_4PI = 1.0 / (4.0 * np.pi)

    def _H_at_points(pts):
        """Vectorized H field at *pts* (shape (N, 3)) from all sources.

        Returns (N, 3) H values.  This is the per-chunk hot loop:
        memory ~ N * M * 3 floats for the dx and cross arrays.
        """
        dx = pts[:, None, :] - centers[None, :, :]  # (N, M, 3)
        r2 = np.sum(dx * dx, axis=2)
        r2 = np.maximum(r2, 1e-60)
        r3_inv = areas[None, :] / (r2 * np.sqrt(r2))  # (N, M)
        cross = np.cross(J[None, :, :], dx)  # (N, M, 3)
        return INV_4PI * np.sum(cross * r3_inv[:, :, None], axis=1)

    t_gl, w_gl = np.polynomial.legendre.leggauss(n_quad)
    t_01 = 0.5 * (t_gl + 1)
    w_01 = 0.5 * w_gl

    src_extent = float(np.max(np.abs(centers)))
    z_far = 20.0 * src_extent

    # ------------------------------------------------------------------
    # Stage 1: phi(0, 0, z) for every unique z used by an observation.
    # ------------------------------------------------------------------
    z_unique = np.unique(obs[:, 2])
    n_z = len(z_unique)

    # Quadrature points along the (0,0,z) -> (0,0,z_far) ray:
    #   z(t) = z_far - t * (z_far - z_unique[k])
    # x_quad shape: (n_z, n_quad, 3); we then flatten to feed _H_at_points.
    x_axis_quad = np.zeros((n_z, n_quad, 3))
    x_axis_quad[:, :, 2] = (z_far
                            - t_01[None, :] * (z_far - z_unique[:, None]))

    # Memory-bounded chunking on n_z * n_quad observation rows.
    rows_per_chunk = max(1, chunk_size // max(1, n_quad))
    H_axis = np.empty((n_z, n_quad, 3))
    flat = x_axis_quad.reshape(-1, 3)
    H_flat = np.empty_like(flat)
    n_rows = flat.shape[0]
    for s in range(0, n_rows, chunk_size):
        e = min(s + chunk_size, n_rows)
        H_flat[s:e] = _H_at_points(flat[s:e])
    H_axis = H_flat.reshape(n_z, n_quad, 3)

    # phi_axis[k] = (z_far - z_unique[k]) * int_0^1 H_z(z(t)) dt
    phi_axis_arr = (z_far - z_unique) * np.sum(
        w_01[None, :] * H_axis[:, :, 2], axis=1)
    phi_axis_lookup = dict(zip(z_unique, phi_axis_arr))

    # ------------------------------------------------------------------
    # Stage 2: horizontal path (0,0,z) -> (x,y,z) per observation node,
    # vectorised across all off-axis nodes.
    # ------------------------------------------------------------------
    n_obs = len(obs)
    phi = np.empty(n_obs)

    rho2 = obs[:, 0] ** 2 + obs[:, 1] ** 2
    eps_axis = (1e-12 * src_extent) ** 2
    on_axis = rho2 < eps_axis

    # On-axis nodes: phi(P) = phi_axis(z_P) directly.
    if np.any(on_axis):
        for k in np.where(on_axis)[0]:
            phi[k] = phi_axis_lookup[obs[k, 2]]

    off_idx = np.where(~on_axis)[0]
    if off_idx.size:
        x_off = obs[off_idx, 0]
        y_off = obs[off_idx, 1]
        z_off = obs[off_idx, 2]

        # dl vector per off-axis node (z component is zero).
        dl_x = x_off  # (n_off,)
        dl_y = y_off

        # Phi at the axis end of each path (z = z_off[k]).
        phi_axis_off = np.array([phi_axis_lookup[z] for z in z_off])

        # Process off-axis nodes in chunks to keep peak memory bounded.
        max_nodes_per_chunk = max(1, chunk_size // max(1, n_quad))
        for s in range(0, off_idx.size, max_nodes_per_chunk):
            e = min(s + max_nodes_per_chunk, off_idx.size)
            n_chunk = e - s

            # Build (n_chunk, n_quad, 3) quadrature points along
            # (0,0,z) -> (x,y,z): position(t) = (t*x, t*y, z)
            quad_pts = np.empty((n_chunk, n_quad, 3))
            quad_pts[:, :, 0] = t_01[None, :] * dl_x[s:e, None]
            quad_pts[:, :, 1] = t_01[None, :] * dl_y[s:e, None]
            quad_pts[:, :, 2] = z_off[s:e, None]
            flat_pts = quad_pts.reshape(-1, 3)

            H_flat = _H_at_points(flat_pts)
            H_quad = H_flat.reshape(n_chunk, n_quad, 3)

            # integrand[k, q] = H_x[k, q] * dl_x[k] + H_y[k, q] * dl_y[k]
            integrand = (H_quad[:, :, 0] * dl_x[s:e, None]
                         + H_quad[:, :, 1] * dl_y[s:e, None])
            path_int = np.sum(w_01[None, :] * integrand, axis=1)  # (n_chunk,)
            phi[off_idx[s:e]] = phi_axis_off[s:e] - path_int

    return phi


def _h_segments_complex(segments, obs_points, currents,
                         seg_chunk=4096):
    """Finite-segment Biot-Savart H at obs_points with complex per-segment currents.

    segments:  (N_seg, 2, 3) real endpoints
    obs_points: (N_obs, 3) real
    currents:  (N_seg,) complex
    seg_chunk: max segments processed at once (memory cap
        ~ n_obs * seg_chunk * 3 complex floats)
    returns:   (N_obs, 3) complex

    Uses the exact analytical line-segment formula:
        H = I/(4 pi d) * (cos(alpha1) - cos(alpha2)) * e_perp
    mirroring biot_savart.h_segments_batch but with per-segment
    complex current. Fully vectorized across segments (previous
    Python loop over N_seg was the phi_inc bottleneck).
    """
    INV_4PI = 1.0 / (4.0 * np.pi)
    obs = np.asarray(obs_points, dtype=float)
    if obs.ndim == 1:
        obs = obs.reshape(1, 3)
    n_obs = len(obs)

    seg = np.asarray(segments, dtype=float)
    p1s_all = seg[:, 0, :]
    p2s_all = seg[:, 1, :]
    I_all = np.asarray(currents, dtype=complex)

    H_total = np.zeros((n_obs, 3), dtype=complex)

    # Chunk on segments to cap peak memory at
    # n_obs * seg_chunk * 3 (* 16 B complex).
    for s0 in range(0, len(p1s_all), seg_chunk):
        s1 = min(s0 + seg_chunk, len(p1s_all))
        p1s = p1s_all[s0:s1]  # (M, 3)
        p2s = p2s_all[s0:s1]
        I = I_all[s0:s1]      # (M,)
        M = s1 - s0

        dl = p2s - p1s                             # (M, 3)
        L = np.linalg.norm(dl, axis=1)             # (M,)
        valid_seg = L > 1e-30
        if not np.any(valid_seg):
            continue
        e_l = np.zeros_like(dl)
        e_l[valid_seg] = dl[valid_seg] / L[valid_seg, np.newaxis]

        # r1, r2 shape: (n_obs, M, 3)
        r1 = obs[:, None, :] - p1s[None, :, :]
        r2 = obs[:, None, :] - p2s[None, :, :]

        # cross(e_l, r1) per (obs, seg) -> (n_obs, M, 3)
        cross = np.cross(e_l[None, :, :], r1)
        d = np.linalg.norm(cross, axis=2)          # (n_obs, M)
        r1_mag = np.linalg.norm(r1, axis=2)
        r2_mag = np.linalg.norm(r2, axis=2)

        ok = (d > 1e-30) & (r1_mag > 1e-30) & (r2_mag > 1e-30)
        ok &= valid_seg[None, :]

        # Safe denominators to avoid /0 warnings outside 'ok'
        d_safe = np.where(ok, d, 1.0)
        r1_safe = np.where(ok, r1_mag, 1.0)
        r2_safe = np.where(ok, r2_mag, 1.0)

        cos_a1 = np.sum(r1 * e_l[None, :, :], axis=2) / r1_safe
        cos_a2 = np.sum(r2 * e_l[None, :, :], axis=2) / r2_safe
        geom = np.where(ok, INV_4PI / d_safe * (cos_a1 - cos_a2), 0.0)

        # e_perp = cross / d (zero where not ok)
        e_perp = cross / d_safe[..., None]

        # H contribution per (obs, seg): I[seg] * geom[obs, seg] * e_perp[obs, seg, :]
        weight = I[None, :] * geom                  # (n_obs, M) complex
        H_total += np.sum(weight[..., None] * e_perp, axis=1)

    return H_total


def compute_phi_inc_from_filaments(obs_points, filament_paths, currents,
                                    n_quad=20, chunk_size=512, z_far=None):
    """Complex phi_inc from a bundle of filaments with per-filament currents.

    Path-integration follows the same two-stage pattern as
    ``compute_phi_inc_from_surface_J`` (axis ray from z_far + horizontal
    sweep) but evaluates H via the exact finite-segment Biot-Savart
    formula and supports COMPLEX per-filament currents, so the output
    is a complex phi_inc suitable for ScalarBIESIBCSolver driven by a
    coil with non-trivial AC current distribution (e.g. Tier A skin
    effect).

    The method assumes the z-axis ray (0, 0, z_far) -> (0, 0, z_obs)
    and the horizontal ray (0, 0, z_obs) -> (x_obs, y_obs, z_obs) do
    not pierce any filament wire — true for a ring coil centered on
    the z-axis with the workpiece inside the coil bore.

    Args:
        obs_points: (N, 3) observation points [m].
        filament_paths: list of K filaments; each filament is a list of
            ``(p1, p2)`` endpoint tuples, as returned by
            ``CoilBuilder.to_filaments``.
        currents: length-K array of per-filament currents, real or
            complex [A].
        n_quad: Gauss-Legendre order for the 1-D path integrals.
        chunk_size: memory cap on ``n_obs * n_quad`` rows per H batch.
        z_far: reference height where phi ~ 0. Default: 20x max extent.

    Returns:
        (N,) complex phi values.
    """
    obs = np.asarray(obs_points, dtype=float)
    if obs.ndim == 1:
        obs = obs.reshape(1, 3)
    n_obs = len(obs)

    flat_segs = []
    flat_I = []
    for fil_segs, Ik in zip(filament_paths, currents):
        Ik_c = complex(Ik)
        for (p1, p2) in fil_segs:
            flat_segs.append((p1, p2))
            flat_I.append(Ik_c)

    if not flat_segs:
        return np.zeros(n_obs, dtype=complex)

    segs = np.asarray(flat_segs, dtype=float)
    Iseg = np.asarray(flat_I, dtype=complex)

    extent = float(np.max(np.abs(segs)))
    if z_far is None:
        z_far = 20.0 * max(extent, 1e-3)

    # --- Dipole tail correction: phi(0,0,z_far) ~ m_z / (4 pi z_far^2) ---
    # For the gauge phi(infinity) = 0 to be numerically honored, we seed the
    # axis integration with the dipole value at z_far instead of zero. This
    # cancels the O(1/z_far^2) offset from truncating the axis integral at
    # a finite z_far; residual is quadrupole-order (1/z_far^4).
    r_mid = 0.5 * (segs[:, 0, :] + segs[:, 1, :])
    dl_vec = segs[:, 1, :] - segs[:, 0, :]
    m_z = 0.5 * np.sum(
        Iseg * (r_mid[:, 0] * dl_vec[:, 1] - r_mid[:, 1] * dl_vec[:, 0]))
    phi_tail = m_z / (4.0 * np.pi * z_far ** 2)

    t_gl, w_gl = np.polynomial.legendre.leggauss(n_quad)
    t_01 = 0.5 * (t_gl + 1)
    w_01 = 0.5 * w_gl

    # ---- Stage 1: axis phi at unique z ----
    z_unique = np.unique(obs[:, 2])
    n_z = len(z_unique)
    axis_quad = np.zeros((n_z, n_quad, 3))
    axis_quad[:, :, 2] = z_far - t_01[None, :] * (z_far - z_unique[:, None])
    flat = axis_quad.reshape(-1, 3)

    max_rows = max(1, chunk_size)
    H_axis_flat = np.empty((flat.shape[0], 3), dtype=complex)
    for s in range(0, flat.shape[0], max_rows):
        e = min(s + max_rows, flat.shape[0])
        H_axis_flat[s:e] = _h_segments_complex(segs, flat[s:e], Iseg)
    H_axis = H_axis_flat.reshape(n_z, n_quad, 3)
    phi_axis_arr = phi_tail + (z_far - z_unique) * np.sum(
        w_01[None, :] * H_axis[:, :, 2], axis=1)
    phi_axis_lookup = dict(zip(z_unique, phi_axis_arr))

    # ---- Stage 2: horizontal path (0,0,z) -> (x,y,z) ----
    phi = np.empty(n_obs, dtype=complex)
    rho2 = obs[:, 0] ** 2 + obs[:, 1] ** 2
    eps_axis = (1e-12 * max(extent, 1e-3)) ** 2
    on_axis = rho2 < eps_axis
    if np.any(on_axis):
        for k in np.where(on_axis)[0]:
            phi[k] = phi_axis_lookup[obs[k, 2]]
    off_idx = np.where(~on_axis)[0]
    if off_idx.size:
        x_off = obs[off_idx, 0]
        y_off = obs[off_idx, 1]
        z_off = obs[off_idx, 2]
        phi_axis_off = np.array([phi_axis_lookup[z] for z in z_off],
                                dtype=complex)
        max_nodes_per_chunk = max(1, chunk_size // max(1, n_quad))
        for s in range(0, off_idx.size, max_nodes_per_chunk):
            e = min(s + max_nodes_per_chunk, off_idx.size)
            n_chunk = e - s
            quad_pts = np.empty((n_chunk, n_quad, 3))
            quad_pts[:, :, 0] = t_01[None, :] * x_off[s:e, None]
            quad_pts[:, :, 1] = t_01[None, :] * y_off[s:e, None]
            quad_pts[:, :, 2] = z_off[s:e, None]
            flat_pts = quad_pts.reshape(-1, 3)
            H_flat = _h_segments_complex(segs, flat_pts, Iseg)
            H_quad = H_flat.reshape(n_chunk, n_quad, 3)
            integrand = (H_quad[:, :, 0] * x_off[s:e, None]
                         + H_quad[:, :, 1] * y_off[s:e, None])
            path_int = np.sum(w_01[None, :] * integrand, axis=1)
            phi[off_idx[s:e]] = phi_axis_off[s:e] - path_int

    return phi


def compute_phi_inc_from_filaments_arbitrary_axis(
        obs_points, filament_paths, currents, *,
        n_quad=20, chunk_size=512,
        far_ref=None, axis_dir=None):
    """Complex phi_inc with arbitrary integration axis (generalises
    compute_phi_inc_from_filaments to non-z-symmetric coils).

    The original ``compute_phi_inc_from_filaments`` assumes the coil is
    a ring centered on the +z axis with the workpiece inside the bore --
    its 2-stage path (axis ray (0,0,z_far)->(0,0,z_obs) + horizontal
    sweep) only stays clear of filament wires under that geometry.

    This variant lets the caller specify:
      * ``axis_dir``  : unit vector for the "axis" leg of the path
                         (default (0,0,1) = z-axis, matching the
                         original behaviour).
      * ``far_ref``    : point at which phi is taken to be ~0
                         (default: ``axis_dir * 20*max_extent``,
                         matching the original ``z_far``).

    Path topology (matches the original 2-stage form, just rotated):
      Stage 1: ``far_ref`` --(along axis_dir)--> ``(obs projected onto
                 the axis line through far_ref)``
      Stage 2: that intermediate --(perpendicular to axis_dir)-->
                 obs

    The dipole tail correction at far_ref is the projection of the
    coil's magnetic dipole moment onto axis_dir (so ``axis_dir = +z``
    recovers the m_z formula).

    For a 3-turn coil whose helix axis is z but is offset to +x (e.g.
    ``3turnCoil_work``), the original function still works because
    z-axis is INSIDE the coil bore.  This generalised form is needed
    when the coil's helix axis is not aligned with z, OR when the
    workpiece sticks out of the coil bore so far that the original
    horizontal sweep would graze a filament wire.

    Args:
        obs_points: (N, 3) observation points [m].
        filament_paths: list of K filaments; each filament is a list of
            ``(p1, p2)`` endpoint tuples.
        currents: length-K array of per-filament currents [A] (real or
            complex).
        n_quad: Gauss-Legendre order for the 1-D path integrals.
        chunk_size: memory cap on ``n_obs * n_quad`` rows per H batch.
        far_ref: (3,) point [m] where phi ~ 0.  Default: 20x max extent
            along axis_dir.
        axis_dir: (3,) unit vector for stage-1 integration direction.
            Default (0, 0, 1).

    Returns:
        (N,) complex phi values.
    """
    obs = np.asarray(obs_points, dtype=float)
    if obs.ndim == 1:
        obs = obs.reshape(1, 3)
    n_obs = len(obs)

    flat_segs = []
    flat_I = []
    for fil_segs, Ik in zip(filament_paths, currents):
        Ik_c = complex(Ik)
        for (p1, p2) in fil_segs:
            flat_segs.append((p1, p2))
            flat_I.append(Ik_c)
    if not flat_segs:
        return np.zeros(n_obs, dtype=complex)
    segs = np.asarray(flat_segs, dtype=float)
    Iseg = np.asarray(flat_I, dtype=complex)

    if axis_dir is None:
        axis_dir = np.array([0.0, 0.0, 1.0])
    else:
        axis_dir = np.asarray(axis_dir, dtype=float)
        n = np.linalg.norm(axis_dir)
        if n < 1e-30:
            raise ValueError("axis_dir must be non-zero")
        axis_dir = axis_dir / n

    extent = float(np.max(np.abs(segs)))
    if far_ref is None:
        far_ref = axis_dir * (20.0 * max(extent, 1e-3))
    far_ref = np.asarray(far_ref, dtype=float)

    # Dipole tail correction at far_ref:
    #   phi_dipole(r) = (m . r_hat) / (4 pi |r|^2)
    # m = (1/2) sum_seg I_seg * (r_mid x dl)
    r_mid = 0.5 * (segs[:, 0, :] + segs[:, 1, :])
    dl_vec = segs[:, 1, :] - segs[:, 0, :]
    cross_md = np.cross(r_mid, dl_vec)  # (N_seg, 3) real
    # m_complex per cartesian component: shape (3,) complex
    m_vec = 0.5 * np.sum(Iseg[:, None] * cross_md, axis=0)
    r_far_mag = np.linalg.norm(far_ref)
    if r_far_mag > 1e-30:
        r_hat = far_ref / r_far_mag
        phi_far = float(np.dot(m_vec.real, r_hat)) / (4.0 * np.pi * r_far_mag ** 2) \
                  + 1j * float(np.dot(m_vec.imag, r_hat)) / (4.0 * np.pi * r_far_mag ** 2)
    else:
        phi_far = 0.0 + 0.0j

    t_gl, w_gl = np.polynomial.legendre.leggauss(n_quad)
    t_01 = 0.5 * (t_gl + 1)  # (n_quad,)
    w_01 = 0.5 * w_gl

    # ---- Stage 1: along axis_dir, from far_ref to (obs projected onto
    # axis line through far_ref).
    # Param: r(t) = far_ref + t * (mid_obs - far_ref)
    # where mid_obs = far_ref + ((obs - far_ref) . axis_dir) * axis_dir
    proj = np.dot(obs - far_ref[None, :], axis_dir)  # (n_obs,) along-axis distance
    mid_obs = far_ref[None, :] + proj[:, None] * axis_dir[None, :]  # (n_obs, 3)
    # Segment vector from far_ref to mid_obs (length = proj along axis_dir)
    s1_vec = mid_obs - far_ref[None, :]  # (n_obs, 3)

    # Quadrature points for stage 1
    quad_pts_s1 = far_ref[None, None, :] + t_01[None, :, None] * s1_vec[:, None, :]
    flat_s1 = quad_pts_s1.reshape(-1, 3)
    H_s1_flat = np.empty((flat_s1.shape[0], 3), dtype=complex)
    max_rows = max(1, chunk_size)
    for s in range(0, flat_s1.shape[0], max_rows):
        e = min(s + max_rows, flat_s1.shape[0])
        H_s1_flat[s:e] = _h_segments_complex(segs, flat_s1[s:e], Iseg)
    H_s1 = H_s1_flat.reshape(n_obs, n_quad, 3)
    integrand_s1 = np.sum(H_s1 * s1_vec[:, None, :], axis=2)  # H . dr/dt
    # phi at mid_obs = phi_far - integral of H.dl from far_ref to mid_obs
    phi_mid = phi_far - np.sum(w_01[None, :] * integrand_s1, axis=1)

    # ---- Stage 2: from mid_obs to obs (perpendicular to axis_dir).
    s2_vec = obs - mid_obs  # (n_obs, 3)
    quad_pts_s2 = mid_obs[:, None, :] + t_01[None, :, None] * s2_vec[:, None, :]
    flat_s2 = quad_pts_s2.reshape(-1, 3)
    H_s2_flat = np.empty((flat_s2.shape[0], 3), dtype=complex)
    for s in range(0, flat_s2.shape[0], max_rows):
        e = min(s + max_rows, flat_s2.shape[0])
        H_s2_flat[s:e] = _h_segments_complex(segs, flat_s2[s:e], Iseg)
    H_s2 = H_s2_flat.reshape(n_obs, n_quad, 3)
    integrand_s2 = np.sum(H_s2 * s2_vec[:, None, :], axis=2)
    phi = phi_mid - np.sum(w_01[None, :] * integrand_s2, axis=1)

    return phi


def extract_surface_J_from_phi(mesh, phi_vec_complex, order=1):
    """Extract per-element complex surface current density J_s.

    SIBC scalar BIE: phi on surface is magnetic scalar potential. The
    tangential magnetic field is H_t = -grad_s(phi), and the surface
    current density is J_s = n x H_t = -n x grad_s(phi).

    Args:
        mesh: NGSolve surface mesh (2D or 3D with BND).
        phi_vec_complex: (ndof,) complex coefficient vector for phi.
        order: H1 polynomial order matching the BIE solve (default 1).

    Returns:
        centroids: (M, 3) float, element area-centroids.
        areas: (M,) float, element areas.
        J_s: (M, 3) complex, element-averaged surface current density.
    """
    from ngsolve import (H1, GridFunction, Integrate, CF, BND,
                         grad, specialcf, Cross)

    fes = H1(mesh, order=order)
    gf_re = GridFunction(fes)
    gf_im = GridFunction(fes)
    gf_re.vec.FV().NumPy()[:] = phi_vec_complex.real
    gf_im.vec.FV().NumPy()[:] = phi_vec_complex.imag

    n_cf = specialcf.normal(3)
    Jre_cf = -Cross(n_cf, grad(gf_re).Trace())
    Jim_cf = -Cross(n_cf, grad(gf_im).Trace())

    elem_A = Integrate(CF(1), mesh, VOL_or_BND=BND, element_wise=True)
    Jre = [Integrate(Jre_cf[i], mesh, VOL_or_BND=BND, element_wise=True)
           for i in range(3)]
    Jim = [Integrate(Jim_cf[i], mesh, VOL_or_BND=BND, element_wise=True)
           for i in range(3)]

    centroids, areas, J_s = [], [], []
    for el in mesh.Elements(BND):
        a = abs(elem_A[el.nr])
        if a < 1e-30:
            continue
        verts = [mesh.vertices[v.nr].point for v in el.vertices]
        c = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)
        Jvec = np.array([(Jre[i][el.nr] + 1j * Jim[i][el.nr]) / a
                         for i in range(3)], dtype=complex)
        centroids.append(c)
        areas.append(a)
        J_s.append(Jvec)
    return np.array(centroids), np.array(areas), np.array(J_s)


def A_from_surface_J(obs_points, src_centroids, src_areas, src_J,
                     chunk_size=4096):
    """Magnetic vector potential at obs points from surface current density.

        A(r) = (mu_0 / 4 pi) * sum_j  J_s[j] * area[j] / |r - c_j|

    Args:
        obs_points: (N, 3) float observation coordinates.
        src_centroids: (M, 3) panel centroids.
        src_areas: (M,) panel areas.
        src_J: (M, 3) complex (or real) surface current density.
        chunk_size: rows of the pair-wise distance matrix per chunk
            (memory cap; each chunk is N_chunk x M).

    Returns:
        (N, 3) complex array (real if src_J is real).
    """
    r = np.asarray(obs_points, dtype=float)
    c = np.asarray(src_centroids, dtype=float)
    a = np.asarray(src_areas, dtype=float)
    J = np.asarray(src_J)
    if r.ndim == 1:
        r = r.reshape(1, 3)
    N = r.shape[0]
    dtype = complex if np.iscomplexobj(J) else float
    MU_0 = 4e-7 * np.pi
    pref = MU_0 / (4.0 * np.pi)
    A = np.zeros((N, 3), dtype=dtype)
    for s in range(0, N, chunk_size):
        e = min(s + chunk_size, N)
        dx = r[s:e, None, 0] - c[None, :, 0]
        dy = r[s:e, None, 1] - c[None, :, 1]
        dz = r[s:e, None, 2] - c[None, :, 2]
        dist = np.sqrt(dx * dx + dy * dy + dz * dz)
        dist[dist < 1e-12] = 1e-12
        weight = pref * (a / dist)                    # (n, M)
        A[s:e] = np.einsum('nm,mc->nc', weight, J)
    return A


def flux_linkage_in_filaments(filament_paths, src_centroids, src_areas, src_J):
    """Phi_k = integral of A_src . dl along each filament polyline.

    Uses midpoint rule on each polyline segment (A evaluated at midpoint,
    dotted with the segment's dl vector).

    Args:
        filament_paths: list of K filament polylines, each a list of
            ``(p1, p2)`` endpoint tuples (3-vectors in meters).
        src_centroids, src_areas, src_J: panel data for A_from_surface_J.

    Returns:
        (K,) complex (or real) array of flux linkages.
    """
    K = len(filament_paths)
    dtype = complex if np.iscomplexobj(src_J) else float
    Phi = np.zeros(K, dtype=dtype)
    for k, path in enumerate(filament_paths):
        mids = np.array([0.5 * (np.asarray(p1, float) + np.asarray(p2, float))
                         for p1, p2 in path])
        dls = np.array([np.asarray(p2, float) - np.asarray(p1, float)
                        for p1, p2 in path])
        A_mid = A_from_surface_J(mids, src_centroids, src_areas, src_J)
        Phi[k] = np.sum(A_mid * dls)
    return Phi


# Legacy alias
def compute_phi_inc_from_coil(obs_points, coil_path, current, n_quad=50):
    """Compute phi_inc using filament approximation (legacy interface)."""
    path = np.asarray(coil_path, dtype=float)
    center = np.mean(path, axis=0)
    center[2] = path[0, 2]
    radius = np.mean(np.linalg.norm(path[:, :2] - center[np.newaxis, :2], axis=1))
    theta_start = np.arctan2(path[0, 1], path[0, 0])
    theta_end = np.arctan2(path[-1, 1], path[-1, 0])
    arc = (theta_end - theta_start) % (2 * np.pi)
    gap_deg = max(0, 360 - np.degrees(arc))
    return compute_phi_inc_from_loop(obs_points, center, radius, current,
                                     n_quad=n_quad, gap_deg=gap_deg)
