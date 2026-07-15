"""Planar (2D) HDiv-VIM soft-iron demag layer -- the motor-cross-section twin of the tet RT1 solve.

The C++ 2D log-kernel charge Gram (``build_charge_gram`` auto-routes ``HDiv(mesh2d, order=1)``:
charges = -div M on tri/quad cells + M.n on boundary edges, kernel -ln(r)/(2 pi)) supplies the
matrix-free demag operator N = B^T G B.

* ``PlanarDemagBody`` -- C++ matrix-free N, a sparse per-element averaging operator, the
  per-element secant-chi weighted mass, linear and nonlinear (scalar-chi Picard + safeguarded Anderson(1),
  the 2D twin of the C++ tet ``SolveNonlinearPicard``) solves, analytic exterior field evaluation
  from the charge quadrature clouds, and volume-average magnetization.
* ``maxwell_torque_circle`` -- Maxwell-stress torque on a circle in air (real fields, or complex
  phasors -> the TIME-AVERAGED torque).
* ``vim.PlanarSolve`` -- the ``vim.Solve`` 2D dispatch target (one-call solve).

Validation lineage (2026-07-03, distilled from the planar HDiv-VIM research layer; details in
memory ``hdiv-vim-tri-quad-motor``): disk deep-saturation vs the analytic fixed point
M = Mof(H0 - M/2) to 3e-4 (mesh floor) .. 2e-6; ellipse 2:1 nonlinear at 0/45/90 deg vs the
secant-bisection reference 1e-4..3e-3; salient-bar + 6-wire-stator torque-angle sweep vs an
exact-Newton nonlinear FEM mean 0.58%; rotating-field induction gates (conducting cylinder vs the
Bessel closed form 0.19%; mini-cage torque-slip vs an all-in-one FEM 0.5%).

TaskManager: NGSolve assembly uses the caller's region; the C++ solve/apply kernels also self-wrap so
direct ``PlanarDemagBody`` use remains parallel and deterministic.
"""
from __future__ import annotations

import numpy as np
import ngsolve as ng
import scipy.sparse as sp

from ._vim import build_charge_gram
from radia.planar_materials import law_from_table as _law_from_table

MU0 = 4e-7 * np.pi

class PlanarDemagBody:
    """One planar soft-iron body on the matrix-free C++ 2D charge Gram."""

    def __init__(self, mesh, eta=2.0, cg_tol=1e-10, cg_maxit=5000,
                 image_masks=None, image_signs=None):
        if mesh.dim != 2:
            raise ValueError("PlanarDemagBody: mesh.dim must be 2 (got %d)" % mesh.dim)
        self.mesh = mesh
        self.fes = ng.HDiv(mesh, order=1)
        self.B, self.G, self.Mm = build_charge_gram(
            self.fes, eta=eta, image_masks=image_masks, image_signs=image_signs)
        chk = self.G.hex_state_check()
        if chk["ctor"] != chk["now"]:
            raise RuntimeError("PlanarDemagBody: 2D Gram state canary mismatch %r" % (chk,))
        self.B = self.B.tocsr()
        self.Mm = self.Mm.tocsr()
        self.ndof = self.B.shape[1]
        self.n_charge = int(self.G.ndof())
        self.cg_tol = float(cg_tol)
        self.cg_maxit = int(cg_maxit)
        self.last_linear_iterations = 0
        self.last_solve_timings = {}
        # ---- sparse per-element averaging operator E [2*nel, ndof] + areas ----
        fesL = ng.VectorL2(mesh, order=0)
        u = self.fes.TrialFunction()
        vL = fesL.TestFunction()
        mixed = ng.BilinearForm(trialspace=self.fes, testspace=fesL)
        mixed += ng.InnerProduct(u, vL) * ng.dx
        mixed.Assemble()
        areas = ng.Integrate(ng.CoefficientFunction(1.0), mesh, element_wise=True)
        rows, cols, vals = mixed.mat.COO()
        self.nel = mesh.ne
        self.areas = np.array([areas[k] for k in range(self.nel)])
        row_to_out = np.full(fesL.ndof, -1, dtype=np.int64)
        for el in mesh.Elements(ng.VOL):
            dn = fesL.GetDofNrs(el)
            if len(dn) != 2:
                raise RuntimeError("PlanarDemagBody: VectorL2(order=0) element dof count != 2")
            row_to_out[dn[0]] = 2 * el.nr
            row_to_out[dn[1]] = 2 * el.nr + 1
        mixed_rows = row_to_out[np.asarray(rows, dtype=np.int64)]
        if np.any(mixed_rows < 0):
            raise RuntimeError("PlanarDemagBody: VectorL2 row was not assigned to an element")
        mixed_vals = np.asarray(vals, dtype=float) / self.areas[mixed_rows // 2]
        self._E = sp.coo_matrix(
            (mixed_vals, (mixed_rows, np.asarray(cols, dtype=np.int64))),
            shape=(2 * self.nel, self.ndof),
        ).tocsr()
        # averaging-operator gate: the constant (1,0) must average to (1,0) on every element
        gf = ng.GridFunction(self.fes)
        gf.Set(ng.CoefficientFunction((1.0, 0.0)))
        chk1 = self._element_average(gf.vec.FV().NumPy().copy())
        if not (np.allclose(chk1[:, 0], 1.0, atol=1e-8) and np.allclose(chk1[:, 1], 0.0, atol=1e-8)):
            raise RuntimeError("PlanarDemagBody: per-element averaging operator failed the "
                               "constant gate")
        # element -> L2(order=0) scalar dof map (for the secant-chi weighted mass)
        self._fes0 = ng.L2(mesh, order=0)
        self._el2dof0 = np.zeros(self.nel, dtype=int)
        for el in mesh.Elements(ng.VOL):
            dn = self._fes0.GetDofNrs(el)
            self._el2dof0[el.nr] = dn[0]
        self._gfchi = ng.GridFunction(self._fes0)
        self._field_coefficients = None
        self._field_evaluator = None

    # ---------------- projections / solves ----------------
    @staticmethod
    def _coo_arrays(matrix):
        matrix = matrix.tocoo()
        return (
            np.ascontiguousarray(matrix.row, dtype=np.int32),
            np.ascontiguousarray(matrix.col, dtype=np.int32),
            np.ascontiguousarray(matrix.data, dtype=np.float64),
        )

    def _element_average(self, coefficients):
        return np.asarray(self._E @ np.asarray(coefficients, dtype=float)).reshape(self.nel, 2)

    def apply_demag(self, m):
        """Apply N m = B^T G B m in C++ without materializing N."""
        return self.G.apply_configured_demag(
            np.ascontiguousarray(m, dtype=np.float64), True)

    def _mass_riesz(self, rhs):
        return self.G.apply_configured_mass_riesz(
            np.ascontiguousarray(rhs, dtype=np.float64))

    def _solve_mass_system(self, mass, inv_chi, rhs, x0=None):
        if sp.issparse(mass):
            mI, mJ, mV = self._coo_arrays(mass)
            self.G.configure_mass_matrix(mI, mJ, mV, self.ndof)
        else:
            self.G.configure_mass_matrix_ngsolve(mass)
        result = self.G.solve_configured_linear_material_mass_riesz(
            float(inv_chi),
            np.ascontiguousarray(rhs, dtype=np.float64), self.cg_tol,
            self.cg_maxit, True,
            None if x0 is None else np.ascontiguousarray(x0, dtype=np.float64))
        self.last_linear_iterations = int(result["iters"])
        self.last_solve_timings = dict(result.get("timings", {}))
        if self.last_linear_iterations >= self.cg_maxit:
            raise RuntimeError(
                "PlanarDemagBody: C++ mass-Riesz CG did not converge in %d iterations"
                % self.cg_maxit)
        return np.asarray(result["m"], dtype=float)

    def project(self, H_cf):
        """RT-interpolate a (2-component) CoefficientFunction -> coefficient vector."""
        gf = ng.GridFunction(self.fes)
        gf.Set(H_cf)
        return gf.vec.FV().NumPy().copy()

    def weighted_mass(self, invchi_e):
        """W = INT (1/chi(x)) u.v dx with per-element 1/chi (L2 order-0 coefficient)."""
        arr = self._gfchi.vec.FV().NumPy()
        arr[self._el2dof0] = invchi_e
        u, v = self.fes.TnT()
        W = ng.BilinearForm(self.fes)
        W += self._gfchi * ng.InnerProduct(u, v) * ng.dx
        W.Assemble()
        return W.mat

    def solve_linear(self, chi, mu_ext):
        """(M/chi + B^T G B) m = M mu_ext via C++ mass-Riesz CG."""
        return self._solve_mass_system(self.Mm, 1.0 / chi, self.Mm @ mu_ext)

    def solve_weighted(self, invchi_e, mu_ext, x0=None):
        """Solve a per-element secant-susceptibility system without a dense matrix."""
        return self._solve_mass_system(
            self.weighted_mass(invchi_e), 1.0, self.Mm @ mu_ext, x0=x0)

    def elem_H(self, m, mu_ext):
        """Per-element average of the PROJECTED local field H = H_ext + H_dem  [nel, 2]."""
        hdem = -self._mass_riesz(self.apply_demag(m))
        return self._element_average(mu_ext + hdem)

    def solve_nonlinear(self, M_of_h, chi_sec, mu_ext, tol=1e-6, maxit=300, damp=0.6,
                        chi_floor=1e-12):
        """Scalar-chi Picard with damping + safeguarded Anderson(1) on the chi vector.

        ``M_of_h(h)`` and ``chi_sec(h)`` take |H| >= 0 (vectorized).  RAISES on non-convergence
        (fail-loud; no silent partial result).  The engineering default tol is 1e-6 to match the
        tet path; ~1e-3 is the lab's engineering standard when speed matters.
        Returns (m, chi_e, iters, res)."""
        He0 = self._element_average(mu_ext)
        chi_e = np.maximum(chi_sec(np.linalg.norm(He0, axis=1)), chi_floor)
        prev = None
        res = np.inf
        m = None
        for it in range(maxit):
            m = self.solve_weighted(1.0 / chi_e, mu_ext, x0=m)
            He = self.elem_H(m, mu_ext)
            nH = np.maximum(np.linalg.norm(He, axis=1), 1e-300)
            chi_star = np.maximum(M_of_h(nH) / nH, chi_floor)
            r = chi_star - chi_e
            res = np.linalg.norm(r) / max(np.linalg.norm(chi_star), 1e-300)
            if res < tol:
                return m, chi_e, it + 1, res
            chi_next = chi_e + damp * r
            if prev is not None:
                chi_p, r_p = prev
                dr = r - r_p
                den = float(dr @ dr)
                if den > 1e-300:
                    th = float(r @ dr) / den
                    cand = (1 - th) * (chi_e + r) + th * (chi_p + r_p)
                    if np.all(cand > 0):
                        chi_next = np.maximum(cand, chi_floor)
            prev = (chi_e.copy(), r.copy())
            chi_e = np.maximum(chi_next, chi_floor)
        raise RuntimeError(
            "PlanarDemagBody.solve_nonlinear: Picard NOT converged (res=%.2e after %d iters). "
            "Refine the mesh, soften the drive, or lower damp." % (res, maxit))

    # ---------------- postprocessing ----------------
    def M_elem(self, m):
        """Per-element average magnetization [nel, 2]."""
        return self._element_average(m)

    def M_avg(self, m):
        """Area-averaged magnetization (Mx, My)."""
        Me = self._element_average(m)
        w = self.areas / self.areas.sum()
        return float(w @ Me[:, 0]), float(w @ Me[:, 1])

    def demag_factors(self):
        """Rayleigh-quotient demag factors (Dx, Dy) of the uniform modes."""
        out = []
        for comp in ((1.0, 0.0), (0.0, 1.0)):
            gf = ng.GridFunction(self.fes)
            gf.Set(ng.CoefficientFunction(comp))
            mu = gf.vec.FV().NumPy().copy()
            out.append(float(mu @ self.apply_demag(mu)) / float(mu @ (self.Mm @ mu)))
        return tuple(out)

    def H_at(self, P, m):
        """Exterior H of the body's charges at P [n,2] (branch-free; valid outside the body).

        The solved source is materialized and cached in the persistent C++ evaluator."""
        evaluator = self.field_evaluator(m)
        return np.asarray(evaluator.field(np.ascontiguousarray(P, dtype=np.float64).reshape(-1, 2)), float)

    def field_evaluator(self, m):
        """Return the immutable C++ evaluator for ``m``, rebuilding only when coefficients change."""
        coefficients = np.ascontiguousarray(m, dtype=np.float64).reshape(-1)
        if (self._field_evaluator is None or self._field_coefficients is None
                or not np.array_equal(coefficients, self._field_coefficients)):
            self._field_evaluator = self.G.create_planar_field_evaluator(coefficients)
            self._field_coefficients = coefficients.copy()
        return self._field_evaluator

    def Az_at(self, P, m):
        """Exterior A_z of the body's charges: A = +mu0 q/(2 pi) atan2(dy, dx) summed over the
        charge clouds (dA/dy = mu0 H_x, -dA/dx = mu0 H_y; single-valued as a FIELD because the
        total charge is zero).

        BRANCH-CUT CAVEAT: the atan2 FORMULA has a cut along the -x ray of EVERY charge.  It is
        safe when the whole evaluation set sees the charges from one side (e.g. strictly below /
        above / to the +x side of the body).  For points that SURROUND the body (e.g. a bar ring
        around a rotor core) use the polar-integrated single-valued construction
        (dA/dphi = mu0 r H_r anchored on the cut-free +x axis; closure over 2 pi is exact by
        Gauss/zero-total-charge) -- see docs/electric_machine's helper module.

        The solved source and all IMA copies are retained by the persistent C++ evaluator."""
        evaluator = self.field_evaluator(m)
        return np.asarray(evaluator.az(np.ascontiguousarray(P, dtype=np.float64).reshape(-1, 2)), float)


def maxwell_torque_circle(H_total_at, Rc, n=1440, center=(0.0, 0.0)):
    """Maxwell-stress torque per unit length about ``center`` on a circle of radius Rc in AIR
    (B = mu0 H).  ``H_total_at(P)`` returns the TOTAL H at P [n,2] -- real for a static field,
    or complex phasors, in which case the TIME-AVERAGED torque (1/2) Re(H_r H_phi^*) is returned.
    T_static = mu0 Rc^2 oint H_r H_phi dphi."""
    phi = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    P = np.stack([center[0] + Rc * np.cos(phi), center[1] + Rc * np.sin(phi)], axis=1)
    H = H_total_at(P)
    er = np.stack([np.cos(phi), np.sin(phi)], axis=1)
    et = np.stack([-np.sin(phi), np.cos(phi)], axis=1)
    Hr = (H * er).sum(axis=1)
    Ht = (H * et).sum(axis=1)
    if np.iscomplexobj(H):
        acc = 0.5 * float((Hr * np.conj(Ht)).real.sum())
    else:
        acc = float((Hr * Ht).sum())
    return MU0 * Rc * Rc * (2 * np.pi / n) * acc


def solve_planar_demag(mesh, mu_r=None, H_ext=None, bh_table=None, *, magnets=None, eta=2.0,
                       nl_tol=1e-6, nl_maxit=300, image_masks=None, image_signs=None):
    """The ``vim.PlanarSolve`` / ``vim.Solve`` 2D dispatch target: single-region planar soft-iron demag solve.

    ``magnets`` is an optional list of SEPARATE-body PERMANENT MAGNETS [(pm_mesh, M_fixed), ...]
    whose RIGID field (the shared planar_charges.magnet_field_cf) is added to the applied field before
    the H(div) projection -- a hard PM does not demagnetize, so it is a one-way source (no iteration).
    Embedded-PM regions (design B, ``pm=``) are not yet wired here (they need a soft/hard partition of
    the PlanarDemagBody).

    Returns dict: M (n_el,2) per-element magnetization, M_avg (2,), demag_factors (Dx, Dy),
    iters, residual, ndof, n_el, n_charge, nonlinear (bool), linear_solver='mass-riesz-cg-2d', and
    body (the PlanarDemagBody -- reuse it for field evaluation / Maxwell torque / sweeps: G is
    built ONCE, a rigid rotation of the body is only a new H_ext).  The caller wraps TaskManager.
    """
    if H_ext is None:
        raise ValueError("vim.PlanarSolve: H_ext (2-component CoefficientFunction) is required")
    if (mu_r is None) == (bh_table is None):
        raise ValueError("vim.PlanarSolve: provide EXACTLY ONE of mu_r (linear) or bh_table "
                         "(nonlinear)")
    if isinstance(mu_r, dict) or isinstance(bh_table, dict):
        raise NotImplementedError(
            "vim.PlanarSolve: per-region (dict) materials are not wired for the 2D layer yet; "
            "the first increment is a single soft-iron region")
    if magnets:
        from radia.planar_charges import magnet_field_cf
        H_ext = H_ext + magnet_field_cf(magnets)             # rigid PM source (design A), shared CF
    body = PlanarDemagBody(
        mesh, eta=eta, image_masks=image_masks, image_signs=image_signs)
    mu_ext = body.project(H_ext)
    if mu_r is not None:
        if not mu_r > 1.0:
            raise ValueError("vim.PlanarSolve: mu_r must be > 1 (got %r)" % (mu_r,))
        m = body.solve_linear(mu_r - 1.0, mu_ext)
        iters, res, nonlinear = 1, 0.0, False
    else:
        M_of_h, chi_sec, _ = _law_from_table(bh_table)
        m, _, iters, res = body.solve_nonlinear(M_of_h, chi_sec, mu_ext,
                                                tol=nl_tol, maxit=nl_maxit)
        nonlinear = True
    Mx, My = body.M_avg(m)
    field_evaluator = body.field_evaluator(m)
    return {
        "M": body.M_elem(m),
        "m": m,
        "M_avg": np.array([Mx, My]),
        "demag_factors": body.demag_factors(),
        "iters": iters,
        "residual": res,
        "ndof": body.ndof,
        "n_el": body.nel,
        "n_charge": body.n_charge,
        "nonlinear": nonlinear,
        "linear_solver": "mass-riesz-cg-2d",
        "linear_iterations": body.last_linear_iterations,
        "solve_timings": body.last_solve_timings,
        "body": body,
        "image_masks": [] if image_masks is None else list(image_masks),
        "image_signs": [] if image_signs is None else list(image_signs),
        "_field_evaluator": field_evaluator,
        "field_evaluator_stats": dict(field_evaluator.stats()),
    }
