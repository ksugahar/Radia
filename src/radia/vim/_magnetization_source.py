"""Prescribed magnetization projected into an independent HDiv source space."""

import time

import ngsolve as ng
import numpy as np
import radia._radia_pybind as _rp

from . import _image
from ._field_batch import _create_field_evaluator
from ._vim import _curve_mesh, build_charge_gram


def _field_coefficient_algorithm(evaluator_stats, requested):
    """Choose the native field path used while NGSolve assembles a solve.

    NGSolve evaluates a coefficient at short SIMD batches.  The evaluator's
    generic ``auto`` threshold is therefore not reached for a large prescribed
    source even when the surrounding form has millions of quadrature calls.
    Use the C++ tree explicitly once the charge cloud is large; scalar callers
    retain the exact direct path through :meth:`MagnetizationSource.Field`.
    """
    if requested is not None:
        requested = str(requested)
        if requested not in {"direct", "tree"}:
            raise ValueError("field_cf_algorithm must be 'direct' or 'tree'")
        return requested
    return "tree" if int(evaluator_stats.get("source_count", 0)) >= 512 else "direct"


class MagnetizationSource:
    """Fixed 3D magnetization source for HDiv-VIM coupling.

    ``magnetization`` is L2-projected into a source-owned HDiv space.  The
    resulting normal-continuous field is converted once into an immutable C++
    charge evaluator; no ChargeGram H-matrix is built.  A separate iron mesh
    and HDiv space must be passed to :func:`vim.Solve`, so a physical jump of
    normal magnetization at a PM/iron interface remains representable.

    The caller follows the ordinary NGSolve convention and constructs the
    source inside ``with ng.TaskManager():``.
    """

    permanent_magnet_model = "fixed-given"
    permanent_magnet_level = 1

    def __init__(self, mesh, magnetization, *, order=1, curve_order=None,
                 image=None, curve_gauss=8, field_cf_algorithm=None,
                 field_tree_options=None):
        started = time.perf_counter()
        if mesh.dim != 3:
            raise ValueError("vim.MagnetizationSource: the prescribed-source API is 3D; "
                             "use the planar magnets= path for 2D.")
        order = int(order)
        if order not in (1, 2):
            raise ValueError("vim.MagnetizationSource: order must be 1 or 2")
        if getattr(magnetization, "dim", None) is None:
            uniform = np.asarray(magnetization, dtype=float)
            if uniform.shape != (3,):
                raise ValueError("vim.MagnetizationSource: a numeric magnetization must have length 3")
            magnetization = ng.CoefficientFunction(tuple(float(value) for value in uniform))
        if getattr(magnetization, "dim", None) != 3:
            raise ValueError("vim.MagnetizationSource: magnetization must be a 3-vector "
                             "NGSolve CoefficientFunction/GridFunction")

        if curve_order is None:
            current_curve_order = int(mesh.GetCurveOrder())
            curve_order = current_curve_order if current_curve_order >= 2 else None
        elif int(curve_order) == 0:
            curve_order = None
        elif int(curve_order) != 2:
            raise NotImplementedError(
                "vim.MagnetizationSource: only curve_order=2 is supported")
        else:
            curve_order = 2
            _curve_mesh(mesh, 2)

        image_masks, image_signs = [], []
        if image is not None:
            for axes, sign in _image.image_group(_image.parse_image_string(image)):
                image_masks.append(int(sum(1 << axis for axis in axes)))
                image_signs.append(float(sign))

        space = ng.HDiv(mesh, order=order)
        trial, test = space.TnT()
        mass = ng.BilinearForm(space)
        mass += trial * test * ng.dx
        mass.Assemble()
        rhs = ng.LinearForm(space)
        rhs += magnetization * test * ng.dx
        rhs.Assemble()
        gf_m = ng.GridFunction(space)
        inverse = mass.mat.Inverse(space.FreeDofs(), inverse="sparsecholesky")
        gf_m.vec.data = inverse * rhs.vec
        residual = mass.mat * gf_m.vec - rhs.vec
        relative_residual = float(ng.Norm(residual) / max(ng.Norm(rhs.vec), 1.0e-300))
        projection_done = time.perf_counter()

        _, gram, _ = build_charge_gram(
            space, curve_order=curve_order, curve_gauss=curve_gauss,
            image_masks=image_masks, image_signs=image_signs,
            _materialize_mass=False, _build_hmatrix=False)
        coefficients = np.ascontiguousarray(gf_m.vec.FV().NumPy(), dtype=np.float64).copy()
        evaluator, evaluator_stats = _create_field_evaluator(
            gram, coefficients, order, tree_options=field_tree_options
        )
        coefficient_algorithm = _field_coefficient_algorithm(
            evaluator_stats, field_cf_algorithm
        )

        self.mesh = mesh
        self.space = space
        self.magnetization = gf_m
        self.field_cf = _rp._HDivFieldCoefficient(evaluator, coefficient_algorithm)
        self.image = image
        self.order = order
        self.curve_order = curve_order
        self._field_evaluator = evaluator
        self._coefficients = coefficients
        self.field_cf_algorithm = coefficient_algorithm
        self.stats = dict(
            permanent_magnet_model=self.permanent_magnet_model,
            permanent_magnet_level=self.permanent_magnet_level,
            ndof=int(space.ndof),
            order=order,
            curve_order=curve_order,
            image=image,
            projection_relative_residual=relative_residual,
            projection_wall_s=projection_done-started,
            geometry_wall_s=time.perf_counter()-projection_done,
            total_wall_s=time.perf_counter()-started,
            hmatrix_built=False,
            field_cf_algorithm=coefficient_algorithm,
            field_tree_options={
                name: evaluator_stats[name]
                for name in (
                    "leaf_size", "theta", "tree_min_sources", "auto_min_work",
                    "tree_relative_tolerance", "probe_count",
                )
            },
            field_evaluator=evaluator_stats,
        )

    def Field(self, points, algorithm="auto"):
        """Evaluate the prescribed source field H (A/m) at ``points`` (N,3)."""
        points = np.ascontiguousarray(np.asarray(points, dtype=float).reshape(-1, 3))
        return np.asarray(self._field_evaluator.field(points, str(algorithm)), dtype=float) / (4.0*np.pi)
