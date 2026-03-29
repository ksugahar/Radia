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
            Z_s: Surface impedance [Ohm] (complex scalar).
            omega: Angular frequency [rad/s].

        Returns:
            dict with keys:
                phi: GridFunction (complex) - solved surface potential
                phi_vec: ndarray (complex) - coefficient vector
                H_t_rms: float - RMS tangential H [A/m] (= surface current density)
                P_density: float - time-averaged power loss density [W/m^2]
                gamma: complex - Z_s / (jw * mu_0)
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

        # System matrix
        gamma = Z_s / (1j * omega * MU_0) if omega > 0 and Z_s != 0 else 0
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
        # (time-averaged)
        P_density = 0.5 * Z_s.real * H_t_rms**2 if Z_s != 0 else 0

        # GridFunction output
        gf_phi = GridFunction(self.fes)
        gf_phi.vec.FV().NumPy()[:] = phi_vec.real  # real part

        return {
            'phi': gf_phi,
            'phi_vec': phi_vec,
            'H_t_rms': float(H_t_rms),
            'P_density': float(P_density),
            'area': float(abs(area)),
            'gamma': complex(gamma),
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
    from radia.biot_savart import h_segments

    obs = np.asarray(obs_points, dtype=float)
    center = np.asarray(loop_center, dtype=float)
    a = float(loop_radius)
    I = float(current)
    if obs.ndim == 1:
        obs = obs.reshape(1, 3)

    obs_local = obs - center[np.newaxis, :]

    # Build coil wire segments
    arc_deg = 360 - gap_deg
    n_seg = max(200, int(arc_deg))
    theta = np.linspace(0, np.radians(arc_deg), n_seg + 1)
    coil_segs = [(
        (a * np.cos(theta[i]), a * np.sin(theta[i]), 0),
        (a * np.cos(theta[i + 1]), a * np.sin(theta[i + 1]), 0),
    ) for i in range(n_seg)]

    frac = arc_deg / 360.0

    # Gauss-Legendre quadrature for horizontal integration
    t_gl, w_gl = np.polynomial.legendre.leggauss(n_quad)
    t_01 = 0.5 * (t_gl + 1)
    w_01 = 0.5 * w_gl

    n_pts = len(obs_local)
    phi = np.zeros(n_pts)

    for ip in range(n_pts):
        x, y, z = obs_local[ip]
        rho = math.sqrt(x * x + y * y)

        # Analytical phi on z-axis for circular loop
        r_za = math.sqrt(z * z + a * a)
        phi_axis = (I / 2.0) * (1.0 - z / r_za) * frac

        if rho < 1e-12 * a:
            phi[ip] = phi_axis
        else:
            # Horizontal path integration from (0,0,z) to (x,y,z)
            dl_vec = np.array([x, y, 0.0])
            x_quad = np.outer(t_01, dl_vec)
            x_quad[:, 2] = z

            H_z_integrand = np.zeros(n_quad)
            for iq in range(n_quad):
                H = h_segments(coil_segs, x_quad[iq], current=I)
                H_z_integrand[iq] = np.dot(H, dl_vec)

            phi[ip] = phi_axis - np.sum(w_01 * H_z_integrand)

    return phi



def compute_phi_inc_from_surface_J(obs_points, src_centroids, src_areas,
                                    src_J_vecs, n_quad=20):
    """Compute phi_inc from solved surface current via path integration.

    H field computed as sum of dipole contributions from surface elements.
    Uses vectorized NumPy for efficiency.
    """
    obs = np.asarray(obs_points, dtype=float)
    if obs.ndim == 1:
        obs = obs.reshape(1, 3)

    centers = np.asarray(src_centroids, dtype=float)
    areas = np.asarray(src_areas, dtype=float)
    J = np.asarray(src_J_vecs, dtype=float)

    INV_4PI = 1.0 / (4.0 * np.pi)

    def _H_at_points(pts):
        """Vectorized H field from surface current elements."""
        dx = pts[:, None, :] - centers[None, :, :]  # (N, M, 3)
        r = np.sqrt(np.maximum(np.sum(dx**2, axis=2), 1e-60))  # (N, M)
        r3_inv = areas[None, :] / (r ** 3)  # (N, M)
        cross = np.cross(J[None, :, :], dx)  # (N, M, 3)
        return INV_4PI * np.sum(cross * r3_inv[:, :, None], axis=1)

    t_gl, w_gl = np.polynomial.legendre.leggauss(n_quad)
    t_01 = 0.5 * (t_gl + 1)
    w_01 = 0.5 * w_gl

    src_extent = np.max(np.abs(centers))
    z_far = 20 * src_extent

    # Stage 1: phi(0,0,z) for each unique z
    z_vals = np.unique(obs[:, 2])
    phi_axis = {}
    for z_i in z_vals:
        z_quad = z_far - t_01 * (z_far - z_i)
        x_quad = np.zeros((n_quad, 3))
        x_quad[:, 2] = z_quad
        H_quad = _H_at_points(x_quad)
        phi_axis[z_i] = (z_far - z_i) * np.sum(w_01 * H_quad[:, 2])

    # Stage 2: horizontal path for each node
    phi = np.zeros(len(obs))
    for ip in range(len(obs)):
        x_i, y_i, z_i = obs[ip]
        rho = math.sqrt(x_i * x_i + y_i * y_i)

        if rho < 1e-12 * src_extent:
            phi[ip] = phi_axis[z_i]
        else:
            dl_vec = np.array([x_i, y_i, 0.0])
            x_quad = np.outer(t_01, dl_vec)
            x_quad[:, 2] = z_i
            H_quad = _H_at_points(x_quad)
            integrand = np.sum(H_quad * dl_vec[np.newaxis, :], axis=1)
            phi[ip] = phi_axis[z_i] - np.sum(w_01 * integrand)

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
