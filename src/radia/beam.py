"""Canonical six-dimensional transfer-map analysis.

The numerical implementation lives in the Radia C++ core.  This module only
normalizes NumPy inputs and exposes the checked pybind11 boundary.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from . import _radia_pybind as _native


def _real_finite_array(value, name: str) -> np.ndarray:
    raw = np.asarray(value)
    if np.iscomplexobj(raw):
        raise TypeError(f"{name} must be real")
    array = np.ascontiguousarray(raw, dtype=np.float64)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _positive_integer(value, name: str) -> int:
    if isinstance(value, bool) or int(value) != value or int(value) < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def propagate_variational_map(
    lengths_m,
    A_per_m,
    F2_per_m=None,
    F3_per_m=None,
    names: Sequence[str] | None = None,
    *,
    maximum_order: int = 3,
    maximum_step_m: float = 1.0e-3,
    maximum_steps: int = 1_000_000,
    maximum_region_pairs: int = 100_000,
    input_symmetry_tolerance: float = 1.0e-12,
) -> dict:
    """Propagate a canonical Taylor map through piecewise-constant jets.

    ``A_per_m``, ``F2_per_m``, and ``F3_per_m`` have shapes ``(n,6,6)``,
    ``(n,6,6,6)``, and ``(n,6,6,6,6)``.  The returned map uses
    ``u_out = R*u + T[u,u]/2 + U[u,u,u]/6`` and includes region-resolved
    quadratic, direct cubic, local-cascade, and ordered region-pair terms.
    Coordinates are ``(x, px/p0, y, py/p0, sigma, delta)`` in SI units.
    """
    lengths = _real_finite_array(lengths_m, "lengths_m")
    if lengths.ndim != 1 or lengths.size == 0:
        raise ValueError("lengths_m must have shape (n_segment,)")
    if np.any(lengths <= 0.0):
        raise ValueError("lengths_m must be positive")

    a = _real_finite_array(A_per_m, "A_per_m")
    f2 = None if F2_per_m is None else _real_finite_array(F2_per_m, "F2_per_m")
    f3 = None if F3_per_m is None else _real_finite_array(F3_per_m, "F3_per_m")
    region_names = None if names is None else [str(name) for name in names]
    if region_names is not None and len(region_names) != lengths.size:
        raise ValueError("names must contain one entry per segment")

    order = _positive_integer(maximum_order, "maximum_order")
    if order > 3:
        raise ValueError("maximum_order must be 1, 2, or 3")
    step = float(maximum_step_m)
    tolerance = float(input_symmetry_tolerance)
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("maximum_step_m must be finite and positive")
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError(
            "input_symmetry_tolerance must be finite and nonnegative"
        )

    return _native._beam_variational_map(
        lengths,
        a,
        f2,
        f3,
        region_names,
        maximum_order=order,
        maximum_step_m=step,
        maximum_steps=_positive_integer(maximum_steps, "maximum_steps"),
        maximum_region_pairs=_positive_integer(
            maximum_region_pairs, "maximum_region_pairs"
        ),
        input_symmetry_tolerance=tolerance,
    )


def propagate_grid_function_linear_map(
    field,
    lengths_m,
    reference_positions_m,
    reference_tangents,
    magnetic_rigidity_t_m: float,
    *,
    sample_radius_m: float = 1.0e-3,
    initial_horizontal=(1.0, 0.0, 0.0),
    names: Sequence[str] | None = None,
    curvature_sign: float = 1.0,
    gradient_sign: float = 1.0,
    maximum_step_m: float = 1.0e-3,
    maximum_steps: int = 1_000_000,
) -> dict:
    """Build a linear transfer map directly from an NGSolve GridFunction.

    ``reference_positions_m`` and ``reference_tangents`` contain one local
    reference station per positive entry in ``lengths_m``.  At every station
    the C++ adapter asks NGSolve to evaluate the real three-component field at
    nine transverse points.  A transparent affine fit yields curvature,
    normal/skew quadrupole gradients, Maxwell-residual diagnostics, and the
    accumulated six-dimensional ``R`` map.  No regular-grid field map is
    created; NGSolve retains ownership of element lookup, transformations, and
    GridFunction evaluation.

    This entry point is deliberately first order.  Its returned ``T`` and
    ``U`` tensors are zero; nonlinear transfer attribution remains available
    through :func:`propagate_variational_map` once a physical ``F2/F3`` jet is
    supplied.
    """
    lengths = _real_finite_array(lengths_m, "lengths_m")
    if lengths.ndim != 1 or lengths.size == 0:
        raise ValueError("lengths_m must have shape (n_segment,)")
    if np.any(lengths <= 0.0):
        raise ValueError("lengths_m must be positive")
    positions = _real_finite_array(
        reference_positions_m, "reference_positions_m"
    )
    tangents = _real_finite_array(reference_tangents, "reference_tangents")
    expected_shape = (lengths.size, 3)
    if positions.shape != expected_shape:
        raise ValueError(
            "reference_positions_m must have shape (n_segment, 3)"
        )
    if tangents.shape != expected_shape:
        raise ValueError(
            "reference_tangents must have shape (n_segment, 3)"
        )
    horizontal = _real_finite_array(initial_horizontal, "initial_horizontal")
    if horizontal.shape != (3,):
        raise ValueError("initial_horizontal must have shape (3,)")

    rigidity = float(magnetic_rigidity_t_m)
    radius = float(sample_radius_m)
    curvature = float(curvature_sign)
    gradient = float(gradient_sign)
    step = float(maximum_step_m)
    if not np.isfinite(rigidity) or rigidity == 0.0:
        raise ValueError("magnetic_rigidity_t_m must be finite and nonzero")
    if not np.isfinite(radius) or radius <= 0.0:
        raise ValueError("sample_radius_m must be finite and positive")
    if not np.all(np.isfinite([curvature, gradient])):
        raise ValueError("curvature_sign and gradient_sign must be finite")
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("maximum_step_m must be finite and positive")
    region_names = None if names is None else [str(name) for name in names]
    if region_names is not None and len(region_names) != lengths.size:
        raise ValueError("names must contain one entry per segment")

    return _native._beam_grid_function_linear_map(
        field,
        lengths,
        positions,
        tangents,
        rigidity,
        horizontal,
        sample_radius_m=radius,
        names=region_names,
        curvature_sign=curvature,
        gradient_sign=gradient,
        maximum_step_m=step,
        maximum_steps=_positive_integer(maximum_steps, "maximum_steps"),
    )


__all__ = [
    "propagate_grid_function_linear_map",
    "propagate_variational_map",
]
