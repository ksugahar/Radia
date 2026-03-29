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



# biot_savart_filament removed. Use ngsolve.bem MaxwellDL(J*dC) instead.


def compute_phi_inc_from_loop(obs_points, loop_center, loop_radius, current,
                              n_quad=30, gap_deg=0):
    """Compute incident scalar magnetic potential from a circular current loop.

    TODO: Reimplement with ngsolve.bem MaxwellDL(J*dC).
    """
    raise NotImplementedError(
        "Biot-Savart removed. Use ngsolve.bem MaxwellDL(J*dC) for H field, "
        "then path-integrate for phi_inc.")



# biot_savart_surface_current removed. Use ngsolve.bem MaxwellDL(J*ds) instead.


def compute_phi_inc_from_surface_J(obs_points, src_centroids, src_areas,
                                    src_J_vecs, n_quad=20):
    """Compute phi_inc from solved surface current.

    TODO: Reimplement with ngsolve.bem MaxwellDL(J*ds).
    """
    raise NotImplementedError(
        "Biot-Savart removed. Use ngsolve.bem MaxwellDL(J*ds) for H field, "
        "then path-integrate for phi_inc.")


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
