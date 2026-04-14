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

    def __init__(self, mesh, order=1):
        """Initialize solver and assemble BEM operators.

        Args:
            mesh: NGSolve Mesh (surface mesh, dim=2, or volume mesh with BND).
            order: H1 polynomial order on surface (default 1).
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

        # BEM operators (dense extraction via COO)
        with TaskManager():
            DL_bf = LaplaceDL(u.Trace() * ds) * v.Trace() * ds
            SL_bf = LaplaceSL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds

        ndof = self.ndof
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


def _h_segments_complex(segments, obs_points, currents):
    """Finite-segment Biot-Savart H at obs_points with complex per-segment currents.

    segments:  (N_seg, 2, 3) real endpoints
    obs_points: (N_obs, 3) real
    currents:  (N_seg,) complex
    returns:   (N_obs, 3) complex

    Uses the exact analytical line-segment formula:
        H = I/(4 pi d) * (cos(alpha1) - cos(alpha2)) * e_perp
    mirroring biot_savart.h_segments_batch but with per-segment current.
    """
    INV_4PI = 1.0 / (4.0 * np.pi)
    obs = np.asarray(obs_points, dtype=float)
    if obs.ndim == 1:
        obs = obs.reshape(1, 3)
    n_obs = len(obs)

    seg = np.asarray(segments, dtype=float)
    p1s = seg[:, 0, :]
    p2s = seg[:, 1, :]
    I = np.asarray(currents, dtype=complex)

    dl = p2s - p1s
    L = np.linalg.norm(dl, axis=1)
    valid = L > 1e-30
    e_l = np.zeros_like(dl)
    e_l[valid] = dl[valid] / L[valid, np.newaxis]

    H_total = np.zeros((n_obs, 3), dtype=complex)

    for si in range(len(p1s)):
        if not valid[si]:
            continue
        r1 = obs - p1s[si]
        r2 = obs - p2s[si]
        el = e_l[si]
        cross = np.cross(el, r1)
        d = np.linalg.norm(cross, axis=1)
        ok = d > 1e-30
        if not np.any(ok):
            continue
        e_perp = np.zeros_like(cross)
        e_perp[ok] = cross[ok] / d[ok, np.newaxis]
        r1_mag = np.linalg.norm(r1, axis=1)
        r2_mag = np.linalg.norm(r2, axis=1)
        ok &= (r1_mag > 1e-30) & (r2_mag > 1e-30)
        if not np.any(ok):
            continue
        cos_a1 = np.zeros(n_obs)
        cos_a2 = np.zeros(n_obs)
        cos_a1[ok] = np.sum(r1[ok] * el, axis=1) / r1_mag[ok]
        cos_a2[ok] = np.sum(r2[ok] * el, axis=1) / r2_mag[ok]
        geom = np.zeros(n_obs)
        geom[ok] = INV_4PI / d[ok] * (cos_a1[ok] - cos_a2[ok])
        H_total[ok] += (I[si] * geom[ok])[:, np.newaxis] * e_perp[ok]

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
