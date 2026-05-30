"""
Scalar Potential BIE with SIBC for MQS eddy current computation.

Solves for eddy currents on a conducting workpiece surface using a scalar
potential boundary integral equation with Surface Impedance BC.

Formulation (verified on sphere, <0.1% error for all Z_s):

    (1/2*M - DL + gamma * SL * M_inv * K) phi = rhs

    M  = H1 surface mass matrix
    K  = H1 surface stiffness (Laplace-Beltrami)
    DL = scalar Laplace double layer (ngsolve.bem)
    SL = scalar Laplace single layer (ngsolve.bem)
    gamma = Z_s / (jw * mu0)

    SIBC: dphi/dn = -(Z_s/(jw*mu0)) * Delta_s(phi)
    Derived from Faraday: curl_s(E_t).n = -jw*n.B, E_t = Z_s*J_s

Requires: NGSolve with ngsolve.bem, scipy, numpy
"""

import time
import math
import numpy as np
from scipy.linalg import solve as scipy_solve
from scipy.sparse import coo_matrix

MU_0 = 4e-7 * np.pi


def _to_dense_coo(mat):
    """Extract dense NumPy array from NGSolve BaseMatrix via COO."""
    rows, cols, vals = mat.COO()
    return coo_matrix((vals, (rows, cols)),
                      shape=(mat.height, mat.width)).toarray()


def _to_dense_matvec(mat, ndof):
    """Extract dense matrix via column-by-column MatVec (for BEM operators)."""
    from ngsolve import GridFunction
    dense = np.zeros((ndof, ndof))
    for j in range(ndof):
        ej = mat.CreateColVector()
        ej[:] = 0
        ej[j] = 1.0
        res = mat.CreateColVector()
        res.data = mat * ej
        dense[:, j] = res.FV().NumPy().copy()
    return dense


class ScalarBIE_SIBC:
    """Scalar potential BIE solver with SIBC for workpiece eddy currents.

    Usage:
        solver = ScalarBIE_SIBC(mesh, order=1)
        result = solver.solve(H_inc_cf, Zs, omega)
        # result['phi'], result['H_rms'], result['P_loss'], ...
    """

    def __init__(self, mesh, order=1):
        """Initialize solver and assemble frequency-independent operators.

        Args:
            mesh: NGSolve surface mesh (dim=2, from OCCGeometry or Cubit export)
            order: H1 polynomial order (default 1)
        """
        from ngsolve import (H1, BilinearForm, GridFunction, Integrate,
                             CF, ds, grad, BND, TaskManager, InnerProduct)
        from ngsolve.bem import LaplaceDL, LaplaceSL

        self.mesh = mesh
        self.order = order

        self.fes = H1(mesh, order=order)
        u, v = self.fes.TnT()
        self.ndof = self.fes.ndof
        self.ne = mesh.GetNE(BND)

        t0 = time.perf_counter()

        # BEM operators (frequency-independent).  Caller MUST be
        # inside `with TaskManager():` per CLAUDE.md
        # "Caller Wraps, Helper Does NOT" (2026-05-27).
        DL_bf = LaplaceDL(u.Trace() * ds) * v.Trace() * ds
        SL_bf = LaplaceSL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds

        # Try COO extraction first (fast), fall back to matvec
        try:
            self.DL = _to_dense_coo(DL_bf.mat)
            self.SL = _to_dense_coo(SL_bf.mat)
        except Exception:
            self.DL = _to_dense_matvec(DL_bf.mat, self.ndof)
            self.SL = _to_dense_matvec(SL_bf.mat, self.ndof)

        # Surface mass M
        mass_bf = BilinearForm(self.fes)
        mass_bf += u.Trace() * v.Trace() * ds
        mass_bf.Assemble()
        self.M = _to_dense_coo(mass_bf.mat)

        # Surface stiffness K (Laplace-Beltrami)
        stiff_bf = BilinearForm(self.fes)
        stiff_bf += InnerProduct(grad(u).Trace(), grad(v).Trace()) * ds
        stiff_bf.Assemble()
        self.K = _to_dense_coo(stiff_bf.mat)

        self.M_inv = np.linalg.inv(self.M)

        # Gauge vector for int(phi)dS = 0 constraint
        self._c_gauge = self.M @ np.ones(self.ndof)

        # Surface area
        self.area = float(Integrate(CF(1), mesh, BND))

        self.t_assembly = time.perf_counter() - t0

    def _solve_with_gauge(self, A_mat, rhs):
        """Solve A*phi = rhs with constraint int(phi)dS = 0."""
        n = self.ndof
        dtype = complex if np.iscomplexobj(A_mat) or np.iscomplexobj(rhs) else float
        A_aug = np.zeros((n + 1, n + 1), dtype=dtype)
        A_aug[:n, :n] = A_mat
        A_aug[:n, n] = self._c_gauge
        A_aug[n, :n] = self._c_gauge
        rhs_aug = np.zeros(n + 1, dtype=dtype)
        rhs_aug[:n] = rhs
        return scipy_solve(A_aug, rhs_aug)[:n]

    def compute_phi_inc(self, H_inc_cf):
        """Compute incident scalar potential on surface from H_inc field.

        Given H_inc (CoefficientFunction, 3-vector), recovers phi_inc
        on the surface by solving the surface Poisson equation:
            K * phi_inc = -integral(H_inc . grad_s(v)) dS

        Args:
            H_inc_cf: NGSolve CoefficientFunction (dim=3) for incident H field

        Returns:
            phi_inc: ndarray (ndof,) nodal values of incident potential
        """
        from ngsolve import (LinearForm, ds, grad, InnerProduct)

        v = self.fes.TestFunction()
        lf = LinearForm(self.fes)
        lf += InnerProduct(H_inc_cf, grad(v).Trace()) * ds
        lf.Assemble()
        h_rhs = lf.vec.FV().NumPy().copy()

        # Solve K * phi_inc = -h_rhs with gauge
        # Note: K has null space (constants), gauge fixes it
        return self._solve_with_gauge(self.K, -h_rhs)

    def compute_rhs_direct(self, phi_inc_cf):
        """Compute RHS from known incident potential CoefficientFunction.

        Args:
            phi_inc_cf: NGSolve CoefficientFunction (scalar) for phi_inc

        Returns:
            rhs: ndarray (ndof,) mass-weighted nodal values
        """
        from ngsolve import LinearForm, ds

        v = self.fes.TestFunction()
        lf = LinearForm(self.fes)
        lf += phi_inc_cf * v.Trace() * ds
        lf.Assemble()
        return lf.vec.FV().NumPy().copy()

    def solve(self, H_inc_cf, Zs, omega, phi_inc_cf=None):
        """Solve scalar BIE with SIBC for given incident field and impedance.

        Args:
            H_inc_cf: NGSolve CoefficientFunction (dim=3) for incident H field.
                      Used to compute phi_inc via surface Poisson if phi_inc_cf
                      is not provided.
            Zs: Complex surface impedance [Ohm]. For linear conductor:
                Zs = (1+1j)/sqrt(2) * sqrt(rho*omega*mu/2) where rho=resistivity.
            omega: Angular frequency [rad/s]
            phi_inc_cf: Optional CoefficientFunction for phi_inc directly
                        (e.g., -H0*z for uniform field). If given, H_inc_cf
                        is used only for post-processing.

        Returns:
            dict with keys:
                phi: complex ndarray (ndof,) scalar potential on surface
                H_rms: float, RMS tangential H field [A/m]
                P_loss: float, total power loss [W] (= Re(Z_s) * integral |J|^2 dS)
                Q_reactive: float, reactive power [var]
                Zs_used: complex, surface impedance used
                omega: float, angular frequency
                gamma: complex, Z_s / (jw*mu0)
                t_solve: float, solve time [s]
                ndof: int, number of DOFs
                gf_phi_re: GridFunction for real part of phi
                gf_phi_im: GridFunction for imaginary part of phi
        """
        from ngsolve import (GridFunction, Integrate, CF, ds, grad, BND,
                             InnerProduct)

        t0 = time.perf_counter()

        # RHS
        if phi_inc_cf is not None:
            rhs = self.compute_rhs_direct(phi_inc_cf)
        else:
            phi_inc_nodal = self.compute_phi_inc(H_inc_cf)
            rhs = self.M @ phi_inc_nodal

        # System matrix
        gamma = Zs / (1j * omega * MU_0)
        A_sys = (0.5 * self.M - self.DL
                 + gamma * self.SL @ self.M_inv @ self.K).astype(complex)

        # Solve
        phi = self._solve_with_gauge(A_sys, rhs.astype(complex))

        t_solve = time.perf_counter() - t0

        # Post-processing: H_t = -grad_s(phi), J = n x H_t
        # |J|^2 = |H_t|^2 = |grad_s(phi)|^2
        gf_re = GridFunction(self.fes)
        gf_im = GridFunction(self.fes)
        gf_re.vec.FV().NumPy()[:] = phi.real
        gf_im.vec.FV().NumPy()[:] = phi.imag

        Hsq_re = float(Integrate(
            InnerProduct(grad(gf_re), grad(gf_re)), self.mesh, BND))
        Hsq_im = float(Integrate(
            InnerProduct(grad(gf_im), grad(gf_im)), self.mesh, BND))
        Jsq_total = abs(Hsq_re) + abs(Hsq_im)

        H_rms = math.sqrt(Jsq_total / self.area)

        # Power: P + jQ = (1/2) * Z_s * integral |J_s|^2 dS
        # (factor 1/2 for time average of sinusoidal fields)
        S_complex = 0.5 * Zs * Jsq_total
        P_loss = S_complex.real
        Q_reactive = S_complex.imag

        return {
            'phi': phi,
            'H_rms': float(H_rms),
            'P_loss': float(P_loss),
            'Q_reactive': float(Q_reactive),
            'Zs_used': complex(Zs),
            'omega': float(omega),
            'gamma': complex(gamma),
            't_solve': round(t_solve, 3),
            't_assembly': round(self.t_assembly, 2),
            'ndof': self.ndof,
            'gf_phi_re': gf_re,
            'gf_phi_im': gf_im,
        }

    def frequency_sweep(self, H_inc_cf, freqs, sigma, mu_r=1.0,
                        phi_inc_cf=None):
        """Sweep frequency with automatic Z_s computation.

        The system matrix changes with frequency only through gamma.
        DL and SL are frequency-independent (Laplace kernel).

        Args:
            H_inc_cf: Incident H field CoefficientFunction
            freqs: array of frequencies [Hz]
            sigma: Conductivity [S/m]
            mu_r: Relative permeability (default 1.0)
            phi_inc_cf: Optional direct phi_inc CoefficientFunction

        Returns:
            list of result dicts (one per frequency)
        """
        results = []
        rho = 1.0 / sigma
        mu = MU_0 * mu_r

        for f in freqs:
            omega = 2 * math.pi * f
            delta = math.sqrt(2 * rho / (omega * mu))
            Zs = (1 + 1j) * rho / delta  # = (1+1j)/sqrt(2) * sqrt(rho*omega*mu)

            result = self.solve(H_inc_cf, Zs, omega, phi_inc_cf=phi_inc_cf)
            result['frequency'] = float(f)
            result['skin_depth'] = float(delta)
            results.append(result)

        return results
