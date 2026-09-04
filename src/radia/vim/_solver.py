"""Persistent NGSolve-style owner for HDiv-VIM geometry operators."""

from __future__ import annotations

from math import pi

from ._hysteresis import PlayHysteresisMaterial, SolveHysteresis
from ._solve import hdiv_demag_solve
from ._vim import _curve_mesh


class HDivSolver:
    """Own one 3D mesh and reuse its HDiv space and ChargeGram.

    ``vim.Solve`` remains the convenient one-shot API.  Use this object for
    load sweeps, nonlinear studies, history stepping, and coupled-body solves;
    geometry-only operators stay private to the object instead of travelling
    through result dictionaries.

    The caller owns the surrounding ``ngsolve.TaskManager`` region.
    """

    def __init__(self, mesh, *, order=1, image=None, image_cyclic=None,
                 image_cyclic_alternating=False,
                 cyclic_periodic_boundaries=None,
                 gram_eps=None, leaf=64,
                 eta=2.0, far_quad=None, curve_order=None, curve_gauss=8,
                 ho_far_factor=None, gram_backend="hmat",
                 exact_dense_memory_mb=None):
        if int(getattr(mesh, "dim", -1)) != 3:
            raise ValueError(
                "vim.HDivSolver requires a 3D mesh; use PlanarDemagBody for 2D")
        order = int(order)
        if order not in (1, 2):
            raise ValueError("vim.HDivSolver order must be 1 or 2")
        if curve_order is not None and int(curve_order) != 2:
            raise ValueError(
                "vim.HDivSolver curve_order must be None or 2 for production 3D solves")
        self.mesh = mesh
        self.order = order
        self.image = image
        self.image_cyclic = image_cyclic
        self.image_cyclic_alternating = bool(image_cyclic_alternating)
        if isinstance(cyclic_periodic_boundaries, str):
            raise TypeError(
                "cyclic_periodic_boundaries must be a pair of labels")
        self.cyclic_periodic_boundaries = (
            None if cyclic_periodic_boundaries is None
            else tuple(str(name) for name in cyclic_periodic_boundaries))
        self.gram_eps = gram_eps
        self.leaf = int(leaf)
        self.eta = float(eta)
        self.far_quad = far_quad
        self.curve_order = curve_order
        self.curve_gauss = int(curve_gauss)
        self.ho_far_factor = ho_far_factor
        self.gram_backend = str(gram_backend)
        self.exact_dense_memory_mb = exact_dense_memory_mb
        self._material_caches = {}
        self._history_caches = {}
        self._last_result = None
        self._operator_build_count = 0

    def _geometry_options(self):
        return dict(
            image=self.image, image_cyclic=self.image_cyclic,
            image_cyclic_alternating=self.image_cyclic_alternating,
            cyclic_periodic_boundaries=self.cyclic_periodic_boundaries,
            gram_eps=self.gram_eps, leaf=self.leaf,
            eta=self.eta, far_quad=self.far_quad,
            curve_order=self.curve_order, curve_gauss=self.curve_gauss,
            ho_far_factor=self.ho_far_factor,
            gram_backend=self.gram_backend,
            exact_dense_memory_mb=self.exact_dense_memory_mb,
        )

    @property
    def operator_build_count(self):
        """Number of ChargeGram builds owned by this solver."""
        return int(self._operator_build_count)

    @property
    def last_result(self):
        return self._last_result

    @property
    def fes(self):
        if self._last_result is None:
            raise RuntimeError("vim.HDivSolver has not solved a problem yet")
        return self._last_result["gfM"].space

    @property
    def charge_gram(self):
        if self._last_result is None:
            raise RuntimeError("vim.HDivSolver has not solved a problem yet")
        return self._last_result["_charge_gram"]

    def Solve(self, mu_r=None, H_ext=None, *, B_r=None, bh_table=None,
              magnetization_sources=None, tol=1e-8, maxit=4000,
              nl_maxit=300, nl_tol=1e-6,
              nonlinear_solver="energy-newton", preconditioner="auto",
              linear_solver="auto", newton_inner_tol="auto",
              newton_warmstart="linear", newton_continuation=1,
              newton_reuse_tangent_steps=1, newton_cg_x0=False):
        """Solve one material/load case while reusing this geometry."""
        cache_key = "nonlinear" if bh_table is not None else "linear"
        cache = self._material_caches.get(cache_key)
        options = self._geometry_options()
        options.update(
            B_r=B_r, bh_table=bh_table,
            magnetization_sources=magnetization_sources,
            tol=tol, maxit=maxit, nl_maxit=nl_maxit, nl_tol=nl_tol,
            nonlinear_solver=nonlinear_solver, preconditioner=preconditioner,
            linear_solver=linear_solver, order=self.order,
            newton_inner_tol=newton_inner_tol,
            newton_warmstart=newton_warmstart,
            newton_continuation=newton_continuation,
            newton_reuse_tangent_steps=newton_reuse_tangent_steps,
            newton_cg_x0=newton_cg_x0, _operator_cache=cache,
        )
        result = hdiv_demag_solve(
            self.mesh, mu_r=mu_r, H_ext=H_ext, **options)
        new_cache = result._operator_cache
        if cache is None:
            self._operator_build_count += 1
        self._material_caches[cache_key] = new_cache
        self._last_result = result
        return result

    @staticmethod
    def _material_and_nu0(play, material, nu0):
        if (play is None) == (material is None):
            raise ValueError(
                "vim.HDivSolver.SolveHysteresis requires exactly one of play or material")
        if play is not None:
            material = PlayHysteresisMaterial(*play)
        if nu0 is None:
            mu0 = 4.0e-7*pi
            scaled = mu0*float(material.nu_bound())
            if not (0.0 < scaled < 1.0):
                raise ValueError(
                    "material.nu_bound() must satisfy 0 < mu0*dH/dB < 1")
            nu0 = scaled/(1.0-scaled)
        return material, float(nu0)

    def SolveHysteresis(self, h_steps, play=None, material=None, *, nu0=None,
                        tol=1e-8, maxit=4000, nl_maxit=200,
                        nl_tol=1e-3, initial_b_path=None,
                        initial_state=None, state_quadrature_order=None):
        """Advance a B-input material history on this persistent geometry."""
        if self.gram_backend != "hmat":
            raise NotImplementedError(
                "vim.HDivSolver.SolveHysteresis currently requires gram_backend='hmat'")
        if self.image is not None or self.image_cyclic is not None:
            raise NotImplementedError(
                "vim.HDivSolver.SolveHysteresis does not yet support image symmetry")
        if self.curve_gauss != 8:
            raise NotImplementedError(
                "vim.HDivSolver.SolveHysteresis currently uses curve_gauss=8")
        if self.curve_order == 2 and int(self.mesh.GetCurveOrder()) < 2:
            _curve_mesh(self.mesh, 2)
        material, nu0 = self._material_and_nu0(play, material, nu0)
        cache = self._history_caches.get(nu0)
        result = SolveHysteresis(
            self.mesh, h_steps, material=material, nu0=nu0,
            gram_eps=self.gram_eps, leaf=self.leaf, eta=self.eta,
            far_quad=self.far_quad, ho_far_factor=self.ho_far_factor,
            tol=tol, maxit=maxit, nl_maxit=nl_maxit, nl_tol=nl_tol,
            initial_b_path=initial_b_path, initial_state=initial_state,
            order=self.order, state_quadrature_order=state_quadrature_order,
            _operator_cache=cache,
        )
        if cache is None:
            self._operator_build_count += 1
        self._history_caches[nu0] = result._operator_cache
        self._last_result = result
        return result


__all__ = ["HDivSolver"]
