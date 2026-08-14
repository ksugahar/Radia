"""Native charged-particle tracking and six-dimensional transfer maps.

The numerical implementation lives in the Radia C++ core.  This module only
normalizes NumPy inputs and exposes the checked pybind11 boundary.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from . import _radia_pybind as _native

# Thin owners of the native C++ tracking objects.  These aliases intentionally
# keep all equations, field evaluation, stepping, and trajectory construction
# below the pybind11 boundary.
ParticleSpecies = _native.BeamParticleSpecies
ReferenceParticle = _native.BeamReferenceParticle
CartesianState = _native.BeamCartesianState
Field = _native.BeamField
FieldSample = _native.BeamFieldSample
ZeroField = _native.BeamZeroField
UniformField = _native.BeamUniformField
StateDerivative = _native.BeamStateDerivative
InvariantReport = _native.BeamInvariantReport
Equation = _native.BeamEquation
LorentzEquation = _native.BeamLorentzEquation
Stepper = _native.BeamStepper
ClassicalRK4 = _native.BeamClassicalRK4
Boris2 = _native.BeamBoris2
StepResult = _native.BeamStepResult
TrackPlan = _native.BeamTrackPlan
StepRecord = _native.BeamStepRecord
TrajectorySummary = _native.BeamTrajectorySummary
Trajectory = _native.BeamTrajectory
Tracker = _native.BeamTracker
GridFunctionField = _native.BeamNGSolveGridFunctionField


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


def canonical_body_hamiltonian_jet(
    coefficients,
    magnetic_rigidity_t_m: float,
    *,
    curvature_sign: float = 1.0,
    gradient_sign: float = 1.0,
    reference_beta: float = 1.0,
) -> dict:
    """Build the native fourth-degree canonical body-multipole jet.

    ``coefficients`` contains dipole, normal/skew quadrupole,
    normal/skew sextupole, and normal/skew octupole coefficients.  Returned
    ``H2/H3/H4`` and ``A/F2/F3`` tensors use coordinates
    ``(x, px/p0, y, py/p0, ell, delta)`` with longitudinal Poisson sign -1.
    """
    values = _real_finite_array(coefficients, "coefficients")
    if values.shape != (7,):
        raise ValueError("coefficients must have shape (7,)")
    rigidity = float(magnetic_rigidity_t_m)
    curvature = float(curvature_sign)
    gradient = float(gradient_sign)
    beta = float(reference_beta)
    if not np.isfinite(rigidity) or rigidity == 0.0:
        raise ValueError("magnetic_rigidity_t_m must be finite and nonzero")
    if not np.all(np.isfinite([curvature, gradient, beta])):
        raise ValueError("curvature_sign, gradient_sign, and beta must be finite")
    if beta <= 0.0 or beta > 1.0:
        raise ValueError("reference_beta must be in (0, 1]")
    return _native._beam_canonical_hamiltonian_jet(
        values,
        rigidity,
        curvature_sign=curvature,
        gradient_sign=gradient,
        reference_beta=beta,
    )


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
    nine transverse points.  A first-order harmonic multipole fit yields
    curvature, normal/skew quadrupole gradients, Maxwell-residual diagnostics,
    and the accumulated six-dimensional ``R`` map.  No regular-grid field map
    is created; NGSolve retains ownership of element lookup, transformations,
    and GridFunction evaluation.

    This entry point is deliberately first order.  Its returned ``T`` and
    ``U`` tensors are zero.  Use
    :func:`propagate_grid_function_multipole_map` for the declared paraxial
    field expansion through ``A/F2/F3`` and region-attributed ``R/T/U``.
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


def propagate_grid_function_multipole_map(
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
    multipole_order: int = 3,
    maximum_map_order: int = 3,
    maximum_step_m: float = 1.0e-3,
    maximum_steps: int = 1_000_000,
) -> dict:
    """Expand a solved field in a moving frame and propagate ``R/T/U``.

    At each supplied reference station the native adapter evaluates the live
    NGSolve ``GridFunction`` at the center and eight angles on a transverse
    ring.  It fits the source-free expansion
    ``By + 1j*Bx = sum(C[n] * (x + 1j*y)**n)`` through cubic order, converts
    it to the declared paraxial canonical ``A/F2/F3`` jet (including the
    chromatic expansion of ``1/(1+delta)``), and returns region-attributed
    ``R/T/U`` maps.  The raw samples, frames, coefficients, jets, and fit
    residuals remain inspectable.

    This is a local paraxial body-field model.  The transverse expansion
    assumes a source-free, current-free cross-section and treats each segment
    jet as constant over its supplied length.  It does not supply a complete
    curved-coordinate Hamiltonian, longitudinal/fringe or edge derivatives,
    closed-orbit finding, or a general high-order symplectic Lie map.  Use
    :class:`GridFunctionField` with the native tracker as the independent
    point-evaluation check when those omitted effects may matter.
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
    if horizontal.shape != (3,) or np.linalg.norm(horizontal) == 0.0:
        raise ValueError("initial_horizontal must be a nonzero three-vector")
    rigidity = float(magnetic_rigidity_t_m)
    radius = float(sample_radius_m)
    step = float(maximum_step_m)
    if not np.isfinite(rigidity) or rigidity == 0.0:
        raise ValueError("magnetic_rigidity_t_m must be finite and nonzero")
    if not np.isfinite(radius) or radius <= 0.0:
        raise ValueError("sample_radius_m must be finite and positive")
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("maximum_step_m must be finite and positive")
    field_order = _positive_integer(multipole_order, "multipole_order")
    map_order = _positive_integer(maximum_map_order, "maximum_map_order")
    if field_order > 3 or map_order > 3:
        raise ValueError("multipole_order and maximum_map_order must be <= 3")
    region_names = None if names is None else [str(name) for name in names]
    if region_names is not None and len(region_names) != lengths.size:
        raise ValueError("names must contain one entry per segment")

    return _native._beam_grid_function_multipole_map(
        field,
        lengths,
        positions,
        tangents,
        rigidity,
        horizontal,
        sample_radius_m=radius,
        names=region_names,
        curvature_sign=float(curvature_sign),
        gradient_sign=float(gradient_sign),
        multipole_order=field_order,
        maximum_map_order=map_order,
        maximum_step_m=step,
        maximum_steps=_positive_integer(maximum_steps, "maximum_steps"),
    )


__all__ = [
    "Boris2",
    "CartesianState",
    "ClassicalRK4",
    "Equation",
    "Field",
    "FieldSample",
    "GridFunctionField",
    "InvariantReport",
    "LorentzEquation",
    "ParticleSpecies",
    "ReferenceParticle",
    "StateDerivative",
    "StepRecord",
    "StepResult",
    "Stepper",
    "TrackPlan",
    "Tracker",
    "Trajectory",
    "TrajectorySummary",
    "UniformField",
    "ZeroField",
    "canonical_body_hamiltonian_jet",
    "propagate_grid_function_linear_map",
    "propagate_grid_function_multipole_map",
    "propagate_variational_map",
]
