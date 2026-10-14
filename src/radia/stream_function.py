"""(ACA+)+TSVD accelerated least-norm solver -- generic and kernel-agnostic.

Solves the underdetermined least-norm system  A phi = B  (M field points x
N basis sources, M < N) via a TSVD-regularized pseudo-inverse.  A ~= C D^T is
factored with ACA+ (HACApK's cHACApK_acaplus), then only the small factors are
TSVD'd -> ~ (M/k)^2 faster than naive O(N M^2) TSVD.

KERNEL-AGNOSTIC.  The matrix entry A(i,j) is supplied by the caller as a
callable ``entry(i, j) -> float`` (0-based).  The SAME machinery therefore
serves any Radia source family using Radia's *already-implemented* field
computation -- there is no field kernel baked into this module:

    - coils                : Biot-Savart H/A from filaments
    - permanent magnets,   : fixed-magnet field from magnetization
      soft iron

For convenience, ``radia_field_kernel`` builds such a callback directly from
Radia object handles via ``radia.Fld`` (works for coils AND magnetic materials).

C++ core: src/core/rad_stream_function.cpp (ACA+ delegated to HACApK; only the
TSVD recompression, manuscript Method 2/3 of IEEJ SA-25-020, lives there).

Example
-------
    import numpy as np, radia as rad
    from radia.stream_function import aca_tsvd, pseudo_inverse_solve, radia_field_kernel

    obs = np.random.rand(M, 3)          # observation points
    sources = [rad.ObjRecMag(...), ...] # N Radia objects (coils or magnets)

    entry = radia_field_kernel(obs, sources, component=2)   # A(i,j) = Bz(obs_i, src_j)
    res = aca_tsvd(M, N, entry, modes=40)
    phi = pseudo_inverse_solve(res, B_target, k_mode=30)
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from radia._radia_pybind import _stream_aca_tsvd as _cpp_aca_tsvd

__all__ = [
    "AbeBoundedCurrentPotentialSolution",
    "AbeCurrentPotentialModeDiagnostics",
    "AbeCurrentPotentialSolution",
    "RegularizedTSVD",
    "StreamTSVD",
    "abe_nearest_field_distance_scales",
    "abe_reduce_node_potential_scales",
    "aca_tsvd",
    "pseudo_inverse_solve",
    "pseudo_inverse_solve_regularized",
    "radia_field_kernel",
    "solve",
    "solve_abe_bounded_current_potential",
    "solve_abe_current_potential",
]


@dataclass
class StreamTSVD:
    """Recompressed truncated SVD of A:  A ~= U diag(S) V^T (truncated to modes).

    Attributes
    ----------
    U : ndarray, shape (M, modes), row-major (C-contiguous)
    S : ndarray, shape (modes,)
    V : ndarray, shape (N, modes), row-major (C-contiguous)
    k_aca : int
        ACA+ rank found before TSVD truncation.
    method : str
        Which method produced this SVD: ``"aca_qr_tsvd"`` (ACA + QR + TSVD, the
        fast default) or ``"dense"`` (direct dense TSVD, the exact reference).
        (The legacy manuscript Method 2/3 were removed 2026-07-05, JIAM-2026-36.)
    """

    U: np.ndarray
    S: np.ndarray
    V: np.ndarray
    k_aca: int
    method: str

    @property
    def M(self) -> int:
        return int(self.U.shape[0])

    @property
    def N(self) -> int:
        return int(self.V.shape[0])

    @property
    def modes(self) -> int:
        return int(self.S.shape[0])


def _validated_stream_tsvd_arrays(factor, expected_m, expected_n):
    """Validate a caller-supplied factor before using it as an inverse."""
    if not isinstance(factor, StreamTSVD):
        raise TypeError("precomputed_factor must be a StreamTSVD")
    U = np.asarray(factor.U, dtype=float)
    S = np.asarray(factor.S, dtype=float).reshape(-1)
    V = np.asarray(factor.V, dtype=float)
    modes = int(S.size)
    if modes < 1:
        raise ValueError("precomputed_factor must retain at least one mode")
    if U.shape != (int(expected_m), modes) or V.shape != (
            int(expected_n), modes):
        raise ValueError(
            "precomputed_factor dimensions do not match the weighted response")
    if not np.all(np.isfinite(U)) or not np.all(np.isfinite(S)) or not np.all(
            np.isfinite(V)):
        raise ValueError("precomputed_factor arrays must be finite")
    if S[0] <= 0.0 or np.any(S < 0.0):
        raise ValueError(
            "precomputed_factor singular values must start positive and be "
            "nonnegative")
    ordering_tolerance = 16.0 * np.finfo(float).eps * max(1.0, float(S[0]))
    if np.any(np.diff(S) > ordering_tolerance):
        raise ValueError(
            "precomputed_factor singular values must be nonincreasing")
    return U, S, V


def _positive_integer(value, name):
    numeric = float(value)
    if (not np.isfinite(numeric) or numeric < 1.0
            or numeric != np.floor(numeric)):
        raise ValueError(f"{name} must be a positive integer")
    return int(numeric)


@dataclass
class AbeCurrentPotentialModeDiagnostics:
    """Per-mode evidence used by the DUCAS/Abe selection rule.

    The arrays use Python's zero-based mode numbering.  ``mode_strength`` is
    ``u_i.T @ (W_B B_target)`` and ``normalized_mode_strength`` is its
    magnitude divided by ``sqrt(M)`` (Abe's magnetic-field strength).  Modes
    are *not* required to form one contiguous prefix: symmetry, target-field
    significance, or engineering judgement can be supplied through
    ``allowed_modes`` in :func:`solve_abe_current_potential`.
    """

    singular_values: np.ndarray
    mode_strength: np.ndarray
    normalized_mode_strength: np.ndarray
    target_correlation: np.ndarray
    peak_potential_correction: np.ndarray
    selected: np.ndarray
    allowed: np.ndarray
    rejection_reason: tuple[str, ...]
    residual_peak_to_peak_history: np.ndarray
    residual_rms_history: np.ndarray
    residual_max_abs_history: np.ndarray


@dataclass
class AbeCurrentPotentialSolution:
    """Weighted node-current-potential inverse-design result.

    ``potential`` contains all physical node potentials ``T`` and
    ``independent_potential`` contains ``T_IN`` in ``T = R T_IN``.  The
    initial potential is retained exactly in unselected high-order directions;
    selected modes add only the correction required by the target field.
    """

    potential: np.ndarray
    independent_potential_correction: np.ndarray
    reconstructed_field: np.ndarray
    residual_field: np.ndarray
    factor: StreamTSVD
    diagnostics: AbeCurrentPotentialModeDiagnostics
    selected_modes: np.ndarray
    converged: bool
    stop_reason: str
    residual_peak_to_peak: float
    residual_rms: float
    residual_max_abs: float
    peak_abs_potential: float


@dataclass
class AbeBoundedCurrentPotentialSolution:
    """Positive/box-bounded repeated DUCAS correction result.

    ``solution`` is the final bounded physical potential.  Histories contain
    the checked field residual *after clipping* at every iteration.  The SVD
    factor in ``solution`` is reused throughout; no refactorisation occurs.
    """

    solution: AbeCurrentPotentialSolution
    iterations: int
    clipped_dof_history: np.ndarray
    potential_change_history: np.ndarray
    residual_peak_to_peak_history: np.ndarray
    residual_rms_history: np.ndarray
    residual_max_abs_history: np.ndarray
    converged: bool
    stop_reason: str


def _as_positive_vector(value, size, name, *, default=1.0):
    if value is None:
        result = np.full(int(size), float(default), dtype=float)
    else:
        result = np.asarray(value, dtype=float).reshape(-1)
    if result.shape != (int(size),):
        raise ValueError(f"{name} must have length {size}, got {result.shape}")
    if not np.all(np.isfinite(result)) or np.any(result <= 0.0):
        raise ValueError(f"{name} must contain finite positive values")
    return result


def _reduction_matrix_product(response, reduction):
    """Return ``response @ reduction`` without densifying ``reduction``."""
    if hasattr(reduction, "tocsr"):
        return np.asarray((reduction.T @ response.T).T, dtype=float)
    return np.asarray(response @ np.asarray(reduction, dtype=float), dtype=float)


def _reduction_apply(reduction, value):
    return np.asarray(reduction @ np.asarray(value, dtype=float),
                      dtype=float).reshape(-1)


def abe_nearest_field_distance_scales(node_points, field_points, *,
                                      distance_floor=None, normalize=True):
    """Return Abe's squared nearest-field-point node weights.

    Improved DUCAS compensates the Biot--Savart distance dependence with
    ``delta_i proportional to d_i**2``, where ``d_i`` is the distance from
    current-potential node ``i`` to the nearest magnetic-field evaluation
    point (Abe 2013, eq. 20).  These are *not* finite-element area weights:
    the response matrix already contains the element-area contribution.

    ``distance_floor`` is an optional physical lower bound used only to avoid
    a zero scale when the source and evaluation meshes touch.  The default
    rejects that geometry instead of silently inventing a length scale.
    """
    nodes = np.asarray(node_points, dtype=float)
    fields = np.asarray(field_points, dtype=float)
    if (nodes.ndim != 2 or fields.ndim != 2 or nodes.shape[1] != fields.shape[1]
            or min(nodes.shape[0], fields.shape[0], nodes.shape[1]) <= 0
            or not np.all(np.isfinite(nodes))
            or not np.all(np.isfinite(fields))):
        raise ValueError(
            "node_points and field_points must be finite non-empty (N, dim) "
            "and (M, dim) arrays")
    if distance_floor is not None:
        distance_floor = float(distance_floor)
        if not np.isfinite(distance_floor) or distance_floor <= 0.0:
            raise ValueError("distance_floor must be finite and positive")

    # Chunk the evaluation set: the public API remains dense-array simple,
    # while its temporary storage is bounded for realistic surface meshes.
    nearest_squared = np.full(nodes.shape[0], np.inf, dtype=float)
    chunk = max(1, min(fields.shape[0], 4096))
    for start in range(0, fields.shape[0], chunk):
        difference = nodes[:, None, :] - fields[None, start:start + chunk, :]
        nearest_squared = np.minimum(
            nearest_squared, np.min(np.sum(difference * difference, axis=2),
                                    axis=1))
    if distance_floor is not None:
        nearest_squared = np.maximum(nearest_squared, distance_floor ** 2)
    elif np.any(nearest_squared <= 0.0):
        raise ValueError(
            "a current-potential node coincides with a field evaluation "
            "point; provide a positive distance_floor or separate the meshes")
    if normalize:
        nearest_squared /= float(np.max(nearest_squared))
    return np.ascontiguousarray(nearest_squared)


def abe_reduce_node_potential_scales(reduction, node_potential_scales,
                                     *, normalize=True):
    """Average full-node DUCAS scales into the independent-potential space.

    This is Abe's node-weight reduction (2013, eqs. 24--25) for
    ``T = R T_IN``.  Classical DUCAS ``R`` contains only zero/one equality
    constraints.  Absolute coefficients are used for the same physically
    meaningful average when a caller supplies signed connection constraints.

    Parameters
    ----------
    reduction : array_like or scipy sparse matrix, shape (N, K)
        Independent-to-full potential map ``R``.
    node_potential_scales : array_like, shape (N,)
        Positive full-node scales ``delta_i``.  Larger values make potential
        on that part of the current-carrying surface less expensive.
    normalize : bool, optional
        Divide by the maximum independent scale, as in Abe eq. 25.

    Returns
    -------
    ndarray, shape (K,)
        Positive independent scales ``delta'_j``.
    """
    shape = getattr(reduction, "shape", None)
    if shape is None or len(shape) != 2 or min(shape) <= 0:
        raise ValueError("reduction must be a non-empty two-dimensional matrix")
    node = _as_positive_vector(
        node_potential_scales, int(shape[0]), "node_potential_scales")
    magnitude = abs(reduction) if hasattr(reduction, "tocsr") else np.abs(
        np.asarray(reduction, dtype=float))
    numerator = np.asarray(magnitude.T @ node, dtype=float).reshape(-1)
    denominator = np.asarray(
        magnitude.T @ np.ones(int(shape[0])), dtype=float).reshape(-1)
    if np.any(denominator <= 0.0):
        raise ValueError("every independent potential must map to at least one node")
    result = numerator / denominator
    if normalize:
        result /= float(np.max(result))
    return np.ascontiguousarray(result)


def solve_abe_current_potential(response, target_field, *, reduction=None,
        field_weights=None, node_potential_scales=None,
        independent_potential_scales=None, initial_potential=None,
        external_field=None, allowed_modes=None, minimum_mode_strength=0.0,
        minimum_target_correlation=0.0, relative_singular_threshold=1.0e-12,
        residual_peak_to_peak=None, residual_rms=None,
        maximum_abs_potential=None, modes=None, kmax=None, aca_eps=1.0e-8,
        method="aca_qr_tsvd",
        precomputed_factor=None) -> AbeCurrentPotentialSolution:
    """Solve the weighted DUCAS/Abe node-current-potential inverse problem.

    The numerical contract is the improved DUCAS formulation

    ``W_B B_TG = (W_B A R diag(delta')) q`` and
    ``T = T0 + R diag(delta') q``,

    where ``B_TG = B0 - A T0 - B_external``.  The factorisation is Radia's
    row-major ACA--QR--TSVD kernel.  Selection follows Abe's design rule: a
    mode must have a safe singular value, significant target-field strength,
    an allowed field distribution, and it is accumulated in decreasing
    singular-value order only until the requested physical residual is met.
    This deliberately differs from blindly retaining the first ``k`` modes.

    ``initial_potential`` is the improved-DUCAS manufacturing lever: its
    unselected high-order components remain in the result and the selected
    modes correct only the field error.  ``node_potential_scales`` is the
    geometric/distance lever (larger values permit a larger potential there).

    Parameters
    ----------
    response : array_like, shape (M, N)
        Field response ``A`` from full node potentials to field samples.
    target_field : array_like, shape (M,)
        Desired physical field ``B0``.
    reduction : array_like or sparse matrix, shape (N, K), optional
        Abe constraint map ``T = R T_IN``.  Identity by default.
    field_weights : array_like, shape (M,), optional
        Positive least-squares multipliers ``W_B`` (typically inverse field
        uncertainty or inverse local tolerance).
    node_potential_scales, independent_potential_scales : array_like, optional
        Supply either positive full-node ``delta_i`` (reduced by Abe eq. 24)
        or positive independent scales ``delta'_j`` directly.
    allowed_modes : array_like, optional
        Boolean mask of length ``modes`` or zero-based mode indices.  This is
        how symmetry and inspection of each magnetic-field eigen-distribution
        are imposed; non-contiguous modes are supported.
    minimum_mode_strength : float, optional
        Minimum ``abs(u_i.T W_B B_TG) / sqrt(M)``.
    minimum_target_correlation : float, optional
        Minimum ``abs(u_i.T W_B B_TG) / ||W_B B_TG||``.
    residual_peak_to_peak, residual_rms : float, optional
        Physical-field stopping tolerances.  If both are omitted, every
        significant allowed mode is accumulated.
    maximum_abs_potential : float, optional
        Engineering feasibility limit checked on the final full potential.
    precomputed_factor : StreamTSVD, optional
        Factor of the identically weighted/reduced response.  Supplying it
        skips ACA--QR--TSVD and is the intended repeated-shimming path.

    Returns
    -------
    AbeCurrentPotentialSolution
        Potentials, reconstructed field, residuals, complete mode evidence,
        and a fail-loud convergence reason.
    """
    A = np.asarray(response, dtype=float)
    if A.ndim != 2 or min(A.shape) <= 0 or not np.all(np.isfinite(A)):
        raise ValueError("response must be a finite non-empty two-dimensional array")
    M, N = A.shape
    target = np.asarray(target_field, dtype=float).reshape(-1)
    if target.shape != (M,) or not np.all(np.isfinite(target)):
        raise ValueError(f"target_field must be finite with length {M}")

    if reduction is None:
        reduction = np.eye(N, dtype=float)
    shape = getattr(reduction, "shape", None)
    if shape is None or len(shape) != 2 or shape[0] != N:
        raise ValueError(f"reduction must have shape ({N}, K)")
    K = int(shape[1])
    if K <= 0:
        raise ValueError("reduction must contain at least one independent potential")

    wb = _as_positive_vector(field_weights, M, "field_weights")
    if (node_potential_scales is not None and
            independent_potential_scales is not None):
        raise ValueError("supply node_potential_scales or independent_potential_scales, not both")
    if node_potential_scales is not None:
        delta = abe_reduce_node_potential_scales(
            reduction, node_potential_scales, normalize=True)
    else:
        delta = _as_positive_vector(
            independent_potential_scales, K,
            "independent_potential_scales")

    initial = (np.zeros(N, dtype=float) if initial_potential is None else
               np.asarray(initial_potential, dtype=float).reshape(-1))
    external = (np.zeros(M, dtype=float) if external_field is None else
                np.asarray(external_field, dtype=float).reshape(-1))
    if initial.shape != (N,) or not np.all(np.isfinite(initial)):
        raise ValueError(f"initial_potential must be finite with length {N}")
    if external.shape != (M,) or not np.all(np.isfinite(external)):
        raise ValueError(f"external_field must be finite with length {M}")

    for value, name in ((minimum_mode_strength, "minimum_mode_strength"),
                        (minimum_target_correlation,
                         "minimum_target_correlation"),
                        (relative_singular_threshold,
                         "relative_singular_threshold")):
        if not np.isfinite(value) or float(value) < 0.0:
            raise ValueError(f"{name} must be finite and nonnegative")
    if float(minimum_target_correlation) > 1.0 + 1.0e-15:
        raise ValueError("minimum_target_correlation cannot exceed one")
    for value, name in ((residual_peak_to_peak, "residual_peak_to_peak"),
                        (residual_rms, "residual_rms"),
                        (maximum_abs_potential, "maximum_abs_potential")):
        if value is not None and (not np.isfinite(value) or float(value) < 0.0):
            raise ValueError(f"{name} must be finite and nonnegative")

    reduced = _reduction_matrix_product(A, reduction)
    weighted = wb[:, None] * reduced * delta[None, :]
    field_to_correct = target - external - A @ initial
    weighted_target = wb * field_to_correct
    max_rank = (min(M, K) if kmax is None else
                min(_positive_integer(kmax, "kmax"), M, K))
    mode_count = (max_rank if modes is None else
                  min(_positive_integer(modes, "modes"), max_rank))
    if max_rank <= 0 or mode_count <= 0:
        raise ValueError("modes and kmax must retain at least one mode")
    aca_eps = float(aca_eps)
    if not np.isfinite(aca_eps) or aca_eps <= 0.0:
        raise ValueError("aca_eps must be finite and positive")
    if precomputed_factor is None:
        factor = aca_tsvd(
            M, K, lambda i, j: float(weighted[i, j]), modes=mode_count,
            kmax=max_rank, aca_eps=aca_eps, method=method)
    else:
        factor = precomputed_factor

    factor_U, factor_S, factor_V = _validated_stream_tsvd_arrays(factor, M, K)
    nm = int(factor_S.size)
    strengths = np.asarray(factor_U.T @ weighted_target, dtype=float)
    normalized = np.abs(strengths) / np.sqrt(float(M))
    target_norm = float(np.linalg.norm(weighted_target))
    correlation = np.abs(strengths) / max(target_norm, np.finfo(float).tiny)
    singular_safe = factor_S > (float(relative_singular_threshold) *
                                max(float(factor_S[0]), np.finfo(float).tiny))
    allowed = np.ones(nm, dtype=bool)
    if allowed_modes is not None:
        supplied = np.asarray(allowed_modes)
        if supplied.dtype == bool:
            supplied = supplied.reshape(-1)
            if supplied.shape != (nm,):
                raise ValueError(f"boolean allowed_modes must have length {nm}")
            allowed = supplied.copy()
        else:
            try:
                numeric_indices = np.asarray(
                    allowed_modes, dtype=float).reshape(-1)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    "allowed_modes must contain integer mode indices") from exc
            if (not np.all(np.isfinite(numeric_indices))
                    or np.any(numeric_indices != np.floor(numeric_indices))):
                raise ValueError(
                    "allowed_modes must contain finite integer mode indices")
            indices = numeric_indices.astype(np.int64)
            if np.any(indices < 0) or np.any(indices >= nm):
                raise ValueError("allowed_modes contains an out-of-range index")
            allowed[:] = False
            allowed[indices] = True
    def metrics(residual):
        residual = np.asarray(residual, dtype=float)
        return (float(np.ptp(residual)),
                float(np.sqrt(np.mean(residual * residual))),
                float(np.max(np.abs(residual))))

    require_target = (residual_peak_to_peak is not None or
                      residual_rms is not None)

    def target_met(values):
        pp, rms, _ = values
        return ((residual_peak_to_peak is None or
                 pp <= float(residual_peak_to_peak)) and
                (residual_rms is None or rms <= float(residual_rms)))

    correction_independent = np.zeros(K, dtype=float)
    reconstructed = external + A @ initial
    residual = target - reconstructed
    history = [metrics(residual)]
    selected = np.zeros(nm, dtype=bool)
    reasons = ["" for _ in range(nm)]
    peak_correction = np.zeros(nm, dtype=float)
    already_met = require_target and target_met(history[-1])
    stopped = already_met
    for index in range(nm):
        if not allowed[index]:
            reasons[index] = "field_distribution_not_allowed"
            continue
        if not singular_safe[index]:
            reasons[index] = "singular_value_below_threshold"
            continue
        if normalized[index] <= float(minimum_mode_strength):
            reasons[index] = "mode_strength_below_threshold"
            continue
        if correlation[index] < float(minimum_target_correlation):
            reasons[index] = "target_correlation_below_threshold"
            continue
        mode_independent = (delta * factor_V[:, index] *
                            (strengths[index] / factor_S[index]))
        mode_full = _reduction_apply(reduction, mode_independent)
        peak_correction[index] = float(np.max(np.abs(mode_full)))
        if stopped:
            reasons[index] = "not_needed_after_residual_target"
            continue
        correction_independent += mode_independent
        reconstructed += A @ mode_full
        residual = target - reconstructed
        selected[index] = True
        reasons[index] = "selected"
        history.append(metrics(residual))
        if require_target and target_met(history[-1]):
            stopped = True

    full_correction = _reduction_apply(reduction, correction_independent)
    potential = initial + full_correction
    # Re-evaluate once from the complete result rather than relying on the
    # incremental sum, so the returned residual is an auditable A @ T value.
    reconstructed = external + A @ potential
    residual = target - reconstructed
    final_metrics = metrics(residual)
    peak_potential = float(np.max(np.abs(potential)))
    potential_ok = (maximum_abs_potential is None or
                    peak_potential <= float(maximum_abs_potential))
    residual_ok = (target_met(final_metrics) if require_target else True)
    if not potential_ok:
        stop_reason = "maximum_abs_potential_exceeded"
    elif require_target and residual_ok:
        stop_reason = ("initial_potential_already_met_target" if already_met
                       else "residual_target_met")
    elif require_target:
        stop_reason = "significant_allowed_modes_exhausted"
    else:
        stop_reason = "significant_allowed_modes_accumulated"
    converged = bool(potential_ok and residual_ok)
    histories = np.asarray(history, dtype=float)
    diagnostics = AbeCurrentPotentialModeDiagnostics(
        singular_values=np.ascontiguousarray(factor_S),
        mode_strength=np.ascontiguousarray(strengths),
        normalized_mode_strength=np.ascontiguousarray(normalized),
        target_correlation=np.ascontiguousarray(correlation),
        peak_potential_correction=np.ascontiguousarray(peak_correction),
        selected=np.ascontiguousarray(selected),
        allowed=np.ascontiguousarray(allowed),
        rejection_reason=tuple(reasons),
        residual_peak_to_peak_history=np.ascontiguousarray(histories[:, 0]),
        residual_rms_history=np.ascontiguousarray(histories[:, 1]),
        residual_max_abs_history=np.ascontiguousarray(histories[:, 2]))
    return AbeCurrentPotentialSolution(
        potential=np.ascontiguousarray(potential),
        independent_potential_correction=np.ascontiguousarray(
            correction_independent),
        reconstructed_field=np.ascontiguousarray(reconstructed),
        residual_field=np.ascontiguousarray(residual), factor=factor,
        diagnostics=diagnostics,
        selected_modes=np.flatnonzero(selected).astype(np.int64),
        converged=converged, stop_reason=stop_reason,
        residual_peak_to_peak=final_metrics[0],
        residual_rms=final_metrics[1],
        residual_max_abs=final_metrics[2],
        peak_abs_potential=peak_potential)


def solve_abe_bounded_current_potential(response, target_field, *,
        lower_potential=None, upper_potential=None, max_iterations=64,
        relaxation=1.0, stagnation_tolerance=1.0e-12, **solve_options
        ) -> AbeBoundedCurrentPotentialSolution:
    """Repeat DUCAS correction with physical potential/iron bounds.

    This implements the positive-only passive-shimming loop described by Abe:

    1. solve the error field in the precomputed SVD eigenmodes;
    2. force negative or over-capacity material to its physical bound;
    3. recompute the error field caused by that bounded placement; and
    4. solve the new error again with the *same* SVD factorisation.

    Use ``lower_potential=0`` for positive-only saturated iron.  For an
    add/remove topology layer, give signed lower/upper capacity arrays.  This
    remains a continuous current-potential planning solve; the caller must
    still map signed material demand to full-strength topology changes and
    re-solve the physical HDiv-MMM system.

    All keyword arguments accepted by :func:`solve_abe_current_potential` are
    accepted through ``solve_options`` except ``precomputed_factor`` (managed
    here), ``maximum_abs_potential`` (use explicit bounds instead), and a
    nontrivial ``reduction``.  Apply ``R`` first and pass the reduced response
    when bounds belong to independent potentials/material cells; clipping full
    nodes independently would otherwise violate the equality constraints.
    """
    if "precomputed_factor" in solve_options:
        raise ValueError("bounded solve owns precomputed_factor reuse")
    if "maximum_abs_potential" in solve_options:
        raise ValueError("bounded solve uses lower_potential/upper_potential")
    if solve_options.get("reduction", None) is not None:
        raise ValueError(
            "bounded solve requires an already reduced response; clipping "
            "full potentials can violate reduction constraints")
    A = np.asarray(response, dtype=float)
    if A.ndim != 2 or min(A.shape) <= 0 or not np.all(np.isfinite(A)):
        raise ValueError("response must be a finite non-empty two-dimensional array")
    M, N = A.shape
    target = np.asarray(target_field, dtype=float).reshape(-1)
    if target.shape != (M,) or not np.all(np.isfinite(target)):
        raise ValueError(f"target_field must be finite with length {M}")
    max_iterations = _positive_integer(max_iterations, "max_iterations")
    if max_iterations < 1:
        raise ValueError("max_iterations must be positive")
    relaxation = float(relaxation)
    stagnation_tolerance = float(stagnation_tolerance)
    if (not np.isfinite(relaxation) or not 0.0 < relaxation <= 1.0 or
            not np.isfinite(stagnation_tolerance) or
            stagnation_tolerance < 0.0):
        raise ValueError("relaxation and stagnation_tolerance are invalid")

    lower = (np.full(N, -np.inf) if lower_potential is None else
             np.broadcast_to(np.asarray(lower_potential, dtype=float),
                             (N,)).copy())
    upper = (np.full(N, np.inf) if upper_potential is None else
             np.broadcast_to(np.asarray(upper_potential, dtype=float),
                             (N,)).copy())
    if (np.any(np.isnan(lower)) or np.any(np.isnan(upper)) or
            np.any(lower > upper)):
        raise ValueError("potential bounds must be ordered and not NaN")
    if np.all(np.isneginf(lower)) and np.all(np.isposinf(upper)):
        raise ValueError("bounded solve requires at least one finite bound")

    initial_option = solve_options.pop("initial_potential", None)
    current = (np.zeros(N, dtype=float) if initial_option is None else
               np.asarray(initial_option, dtype=float).reshape(-1))
    if current.shape != (N,) or not np.all(np.isfinite(current)):
        raise ValueError(f"initial_potential must be finite with length {N}")
    current = np.clip(current, lower, upper)
    external = solve_options.get("external_field", None)
    external = (np.zeros(M, dtype=float) if external is None else
                np.asarray(external, dtype=float).reshape(-1))
    if external.shape != (M,) or not np.all(np.isfinite(external)):
        raise ValueError(f"external_field must be finite with length {M}")

    residual_pp_target = solve_options.get("residual_peak_to_peak", None)
    residual_rms_target = solve_options.get("residual_rms", None)
    require_target = (residual_pp_target is not None or
                      residual_rms_target is not None)
    if not require_target:
        raise ValueError(
            "bounded solve requires residual_peak_to_peak and/or residual_rms "
            "as an acceptance criterion")

    def residual_metrics(potential):
        residual = target - external - A @ potential
        return residual, (float(np.ptp(residual)),
                          float(np.sqrt(np.mean(residual * residual))),
                          float(np.max(np.abs(residual))))

    def target_met(values):
        return ((residual_pp_target is None or
                 values[0] <= float(residual_pp_target)) and
                (residual_rms_target is None or
                 values[1] <= float(residual_rms_target)))

    factor = None
    clipped_history = []
    change_history = []
    residual_history = []
    final = None
    stop_reason = "bounded_max_iterations"
    converged = False
    for _ in range(max_iterations):
        unconstrained = solve_abe_current_potential(
            A, target, initial_potential=current,
            precomputed_factor=factor, **solve_options)
        if factor is None:
            factor = unconstrained.factor
        trial = current + relaxation * (unconstrained.potential - current)
        bounded = np.clip(trial, lower, upper)
        clipped = int(np.count_nonzero(bounded != trial))
        change = float(np.max(np.abs(bounded - current)))
        residual, values = residual_metrics(bounded)
        clipped_history.append(clipped)
        change_history.append(change)
        residual_history.append(values)
        residual_ok = target_met(values) if require_target else False
        final = replace(
            unconstrained, potential=np.ascontiguousarray(bounded),
            reconstructed_field=np.ascontiguousarray(target - residual),
            residual_field=np.ascontiguousarray(residual),
            converged=bool(residual_ok),
            stop_reason=("bounded_residual_target_met" if residual_ok else
                         "bounded_iteration"),
            residual_peak_to_peak=values[0], residual_rms=values[1],
            residual_max_abs=values[2],
            peak_abs_potential=float(np.max(np.abs(bounded))))
        current = bounded
        if residual_ok:
            converged = True
            stop_reason = "bounded_residual_target_met"
            break
        if change <= stagnation_tolerance:
            stop_reason = "bounded_stagnation"
            break
    if final is None:  # defensive: max_iterations validation makes this unreachable
        raise RuntimeError("bounded current-potential iteration did not start")
    if not converged:
        final = replace(final,converged=False,stop_reason=stop_reason)
    histories = np.asarray(residual_history, dtype=float)
    return AbeBoundedCurrentPotentialSolution(
        solution=final, iterations=len(clipped_history),
        clipped_dof_history=np.asarray(clipped_history, dtype=np.int64),
        potential_change_history=np.asarray(change_history, dtype=float),
        residual_peak_to_peak_history=np.ascontiguousarray(histories[:, 0]),
        residual_rms_history=np.ascontiguousarray(histories[:, 1]),
        residual_max_abs_history=np.ascontiguousarray(histories[:, 2]),
        converged=converged, stop_reason=stop_reason)


def aca_tsvd(M, N, entry, modes=None, kmax=None,
             aca_eps=1.0e-4, method="aca_qr_tsvd") -> StreamTSVD:
    """Truncated SVD of an M x N matrix A supplied by the callback ``entry(i,j)``.

    Two methods -- and ONLY two (peer review JIAM-2026-36):

    - ``method="aca_qr_tsvd"`` (default): ACA+ factors ``A ~= C D^T`` (only
      O(k_aca*(M+N)) entry evaluations), then the standard "SVD of a low-rank
      product" recompression -- QR each tall-skinny factor + ONE small
      ``k_aca x k_aca`` TSVD.  Fast; the production path.  (Name = the full
      pipeline: ACA + QR + TSVD.)
    - ``method="dense"`` (alias ``"tsvd"``): the plain/direct TSVD -- materialise
      the full A via the callback and take its dense SVD (``numpy.linalg.svd``),
      then truncate.  EXACT (no ACA approximation), the trusted reference /
      validation baseline, but O(M*N) entry calls + an O(N*M^2) SVD -- use for
      SMALL problems only.

    Parameters
    ----------
    M, N : int
        Rows (field / observation points) and columns (basis sources); usually N>M.
    entry : callable
        ``entry(i, j) -> float`` returning A(i,j) for 0-based i in [0,M), j in [0,N).
    modes : int, optional
        Singular triplets to return (clamped to the rank).  Defaults to ``kmax``.
    kmax : int, optional
        Maximum ACA+ rank (``aca_qr_tsvd`` only).  Defaults to ``min(M, N)``.
    aca_eps : float, optional
        ACA+ stopping tolerance (``aca_qr_tsvd`` only).  Default 1e-4.
    method : {"aca_qr_tsvd", "dense"}, optional
        "aca_qr_tsvd" (default) = ACA + QR + TSVD; "dense" (alias "tsvd") =
        direct dense SVD (exact reference).  NO backward compatibility: the
        legacy integers 2/3 and the terse "qr"/"aca" now RAISE (the manuscript
        Method 2/3 were removed 2026-07-05; see
        ``memory/aca_tsvd_qr_recompression.md``).

    Returns
    -------
    StreamTSVD
    """
    M = int(M)
    N = int(N)
    if M <= 0 or N <= 0:
        raise ValueError(f"M and N must be positive, got M={M}, N={N}")
    if not callable(entry):
        raise TypeError("entry must be callable: entry(i, j) -> float")

    # Exactly TWO canonical methods -- NO backward-compat shims.  MCP-delivered
    # software carries no legacy aliases: unknown / legacy values (incl. the old
    # ints 2/3 and the terse "qr"/"aca") RAISE.
    _m = "aca_qr_tsvd" if method is None else str(method).strip().lower()
    if _m in ("aca_qr_tsvd", "aca+qr+tsvd"):
        _m = "aca_qr_tsvd"
    elif _m in ("dense", "tsvd"):
        _m = "dense"
    else:
        raise ValueError(
            f"method must be 'aca_qr_tsvd' (default, ACA+QR+TSVD) or 'dense' "
            f"(direct TSVD); got {method!r}")

    if kmax is None:
        kmax = min(M, N)
    kmax = int(min(kmax, M, N))
    if modes is None:
        modes = kmax
    modes = int(min(modes, kmax))

    if _m == "dense":
        # 通常のTSVD: materialise A via the callback, take its dense SVD, truncate.
        A = np.fromiter((float(entry(i, j)) for i in range(M) for j in range(N)),
                        dtype=float, count=M * N).reshape(M, N)
        Uf, Sf, Vtf = np.linalg.svd(A, full_matrices=False)
        m = int(min(modes, Sf.shape[0]))
        return StreamTSVD(U=np.ascontiguousarray(Uf[:, :m]),
                          S=np.ascontiguousarray(Sf[:m]),
                          V=np.ascontiguousarray(Vtf[:m, :].T),
                          k_aca=int(min(M, N)), method="dense")

    # method="aca_qr_tsvd": ACA+ then the QR-of-a-low-rank-product recompression (C++).
    def _entry(i, j):                       # C++ side always gets plain floats
        return float(entry(i, j))
    U, S, V, k_aca = _cpp_aca_tsvd(
        M, N, _entry, int(modes), int(kmax), float(aca_eps))
    return StreamTSVD(U=U, S=S, V=V, k_aca=int(k_aca), method="aca_qr_tsvd")


def pseudo_inverse_solve(result: StreamTSVD, B, k_mode=None) -> np.ndarray:
    """Least-norm pseudo-inverse solve:  phi = V diag(1/S) U^T B.

    Uses the first ``k_mode`` singular triplets (TSVD regularization).

    Parameters
    ----------
    result : StreamTSVD
    B : array_like, shape (M,)
        Target field values.
    k_mode : int, optional
        Number of modes to use (<= result.modes).  Default = result.modes.

    Returns
    -------
    ndarray, shape (N,)
        Basis coefficients phi.
    """
    B = np.asarray(B, dtype=float).ravel()
    if B.shape[0] != result.M:
        raise ValueError(f"B length {B.shape[0]} != M {result.M}")
    k = result.modes if k_mode is None else int(min(k_mode, result.modes))
    k = max(1, k)
    U = result.U[:, :k]   # (M, k)
    S = result.S[:k]      # (k,)
    V = result.V[:, :k]   # (N, k)
    with np.errstate(divide="ignore", invalid="ignore"):
        c = (U.T @ B) / S
    c = np.where(S > 0.0, c, 0.0)
    return V @ c


# --------------------------------------------------------------------------
# Regularisation-folded ACA+TSVD pseudo-inverse
# --------------------------------------------------------------------------
@dataclass
class RegularizedTSVD:
    """Folds an SPD stiffness ``S`` into the ACA+TSVD pseudo-inverse.

    Solves the constrained min-norm problem

        min  phi^T S phi    s.t.    A phi = B

    with the closed-form Lagrangian solution

        phi = S^-1 A^T (A S^-1 A^T)^-1 B.

    Substituting the truncated SVD ``A = U Sigma V^T`` (orthonormal V
    columns from ACA+TSVD) gives the FACTORED form

        phi = (S^-1 V) . W^-1 . Sigma^-1 . U^T B,    W = V^T S^-1 V (k x k).

    The factorisation ``(S^-1 V, W^-1)`` is precomputed once; subsequent
    ``solve(B)`` calls cost only ``O(k * (M + N))`` floating-point ops
    (no further sparse / dense solves).  This makes the Path-A
    compensated iteration's inner solve cheap regardless of the FE
    DOF count.

    ``solve(B, alpha)`` with ``alpha > 0`` is the **generalised Tikhonov**
    (soft data fit) ``min ||A phi - B||^2 + alpha . phi^T S phi`` via the
    SAME cached factorisation -- only an ``alpha I`` is added inside the
    ``k x k`` core ``(alpha I + Sigma^2 W)``.  So an entire L-curve /
    Pareto alpha-sweep re-uses one factorisation (the ``+ alpha I`` is the
    whole cost of moving along the front); ``alpha = 0`` recovers the
    exact-fit form above.

    Special case ``S = I``: ``W = V^T V = I_k`` since ACA+TSVD V columns
    are orthonormal, so ``phi = V Sigma^-1 U^T B`` -- reduces exactly to
    the standard L2 pseudo-inverse (see ``pseudo_inverse_solve``).

    Parameters
    ----------
    base : StreamTSVD
        ACA+TSVD factorisation of A (M x N, M <= k).
    Sinv_V : ndarray, shape (N, k)
        Precomputed ``S^-1 V`` (one k-RHS solve with S).
    W_inv : ndarray, shape (k, k)
        Inverse of ``W = V^T Sinv_V``.

    Notes
    -----
    For an N x N stiffness ``S`` restricted to free DOFs, the caller
    must pass ``S_free`` and the matching ``V`` (= base.V on free
    DOFs).  Padding back to full DOF (e.g. Dirichlet zero on boundary)
    is the caller's responsibility.

    See ``docs/stream_function/theory.ipynb`` and the docs-local
    ``demo_planar_uniform_fem_psi.py`` helper for the canonical Path-A usage
    pattern.
    """

    base: StreamTSVD
    Sinv_V: np.ndarray         # (N, k)
    W_inv: np.ndarray          # (k, k)
    W: np.ndarray = None       # (k, k) = V^T S^-1 V; needed for alpha > 0

    @classmethod
    def from_stiffness(cls, base: "StreamTSVD", S) -> "RegularizedTSVD":
        """Precompute ``Sinv_V`` and ``W_inv`` from an SPD stiffness ``S``.

        Parameters
        ----------
        base : StreamTSVD
            ACA+TSVD factorisation of the system matrix A.
        S : ndarray or scipy.sparse matrix, shape (N, N)
            SPD regularisation matrix (e.g. H1 stiffness on free DOFs).
            Dense numpy arrays go through ``np.linalg.solve``; scipy sparse
            matrices are factored ONCE with ``scipy.sparse.linalg.splu`` and
            all columns of ``V`` are back-solved together.

        Raises
        ------
        ValueError
            If ``S`` is singular / not SPD (No-Fallback: the dense
            ``LinAlgError`` and the sparse ``RuntimeError`` are normalised into
            one clear error), or if the inverted core ``W = V^T S^-1 V`` is
            catastrophically ill-conditioned (``cond > 1e12``, where the
            exact-fit solve would silently lose most of its digits).  Add a
            ridge / mass term to ``S`` or use a better-conditioned seminorm.

        Returns
        -------
        RegularizedTSVD
            Cached factorisation; pass ``B`` to ``.solve(B)``.
        """
        V = base.V                                 # (N, k)
        try:
            if hasattr(S, "tocsc"):
                # scipy sparse path -- factor S ONCE (sparse LU) then back-solve
                # all k columns of V together.  (Was k separate spsolve calls,
                # each re-factorising S; this is O(1) factorisation + k solves --
                # the win that lets N=15k-DOF designs build in ~1 s, not ~13 s.)
                from scipy.sparse.linalg import splu
                Sinv_V = splu(S.tocsc()).solve(np.asarray(V, dtype=float))
            else:
                S_dense = np.asarray(S, dtype=float)
                if S_dense.shape != (V.shape[0], V.shape[0]):
                    raise ValueError(
                        f"S shape {S_dense.shape} must be ({V.shape[0]}, "
                        f"{V.shape[0]}) to match base.V")
                Sinv_V = np.linalg.solve(S_dense, V)
        except (np.linalg.LinAlgError, RuntimeError) as e:
            # No-Fallback: the dense path raises numpy.linalg.LinAlgError and the
            # sparse splu raises RuntimeError on a singular S -- normalise BOTH
            # into ONE clear, actionable error rather than leaking two types.
            raise ValueError(
                f"regularisation matrix S is singular / not SPD "
                f"({type(e).__name__}: {e}). Add a ridge / mass term to S (the "
                f"_seminorm helper adds +1e-10 mass for exactly this), or use "
                f"--confine abe/on to remove the constant-psi null space.") from e
        W = V.T @ Sinv_V                            # (k, k) -- the inverted core
        # No-Fallback: a catastrophically ill-conditioned core makes the exact-
        # fit solve SILENTLY lose most of its digits.  Production paths sit at
        # cond(W) ~ 1e2 (l2/abe) .. 1e7 (h1/off); only a genuinely broken S
        # exceeds 1e12 (< ~4 good digits) -- raise rather than return garbage.
        cond_W = float(np.linalg.cond(W))
        if not np.isfinite(cond_W) or cond_W > 1e12:
            raise ValueError(
                f"regularised core W = V^T S^-1 V is too ill-conditioned "
                f"(cond ~ {cond_W:.1e}); the exact-fit solve would lose most "
                f"digits. Add a ridge / mass term to S, use --confine abe, or a "
                f"better-conditioned --regularize.")
        W_inv = np.linalg.inv(W)               # back-compat; W is sound past here
        return cls(base=base, Sinv_V=Sinv_V, W_inv=W_inv, W=W)

    def solve(self, B, k_mode=None, alpha=0.0) -> np.ndarray:
        """Apply the cached regularised pseudo-inverse to ``B``.

        ``alpha = 0`` (default) -- exact data fit, minimum ``S`` seminorm::

            phi = Sinv_V . W_inv . diag(1/Sigma) . U^T B

        ``alpha > 0`` -- generalised **Tikhonov** (soft data fit), trading
        misfit against the seminorm,
        ``min ||A phi - B||^2 + alpha . phi^T S phi``::

            phi(alpha) = Sinv_V . (alpha I + Sigma^2 W)^-1 . Sigma . U^T B
                       == (A^T A + alpha S)^-1 A^T B

        This is the SAME cached factorisation with a single ``alpha I``
        added inside the ``k x k`` core, so an L-curve / Pareto alpha-sweep
        re-solves only the small core (no re-factorisation).  Special case
        ``S = I``: the core reduces to the classic Tikhonov filter factors
        ``sigma / (sigma^2 + alpha)``.  See
        ``docs/stream_function/regularization.md`` and its docs-local
        ``demo_pareto_tikhonov_aca.py`` helper.

        Parameters
        ----------
        B : array_like, shape (M,)
            Target field values.
        k_mode : int, optional
            Truncate to first ``k_mode`` SVD modes (re-inverts the k_mode
            x k_mode top-left block of W since W_inv depends on k).
            Default = base.modes (use full cached factorisation).
        alpha : float, optional
            Tikhonov weight.  ``0`` (default) = exact-fit min-seminorm;
            ``> 0`` = soft-fit Tikhonov in the ``S`` metric.  Composes
            with ``k_mode`` (hard spectral truncation + smooth damping).

        Returns
        -------
        ndarray, shape (N,)
            ``alpha = 0``: ``phi`` satisfying ``A phi = B`` with minimum
            ``phi^T S phi``.  ``alpha > 0``: the Tikhonov solution
            ``(A^T A + alpha S)^-1 A^T B``.
        """
        B = np.asarray(B, dtype=float).ravel()
        if B.shape[0] != self.base.M:
            raise ValueError(f"B length {B.shape[0]} != M {self.base.M}")
        k_full = self.base.modes
        k = k_full if k_mode is None else int(min(k_mode, k_full))
        k = max(1, k)

        U = self.base.U[:, :k]                      # (M, k)
        Sigma = self.base.S[:k]                     # (k,)

        if k == k_full:
            Sinv_V = self.Sinv_V                    # (N, k_full)
            W_inv = self.W_inv                      # (k_full, k_full)
            W = self.W if self.W is not None else np.linalg.inv(W_inv)
        else:
            Sinv_V = self.Sinv_V[:, :k]             # (N, k)
            V_k = self.base.V[:, :k]                # (N, k)
            W = V_k.T @ Sinv_V                      # (k, k)

        UtB = U.T @ B                               # (k,)

        if alpha <= 0.0:
            # exact-fit min-seminorm: phi = Sinv_V . W^-1 . Sigma^-1 . U^T B.
            # Solve W y = c (factored) instead of caching/applying inv(W): a
            # factored solve is strictly more stable on an ill-conditioned W
            # (No-Fallback: do not silently lose digits to an explicit inverse).
            with np.errstate(divide="ignore", invalid="ignore"):
                c = UtB / Sigma                     # (k,)
            c = np.where(Sigma > 0.0, c, 0.0)
            return Sinv_V @ np.linalg.solve(W, c)

        # generalised Tikhonov: (alpha I + Sigma^2 W) y = Sigma U^T B
        #   phi = Sinv_V . y == (A^T A + alpha S)^-1 A^T B
        core = alpha * np.eye(k) + (Sigma ** 2)[:, None] * W
        y = np.linalg.solve(core, Sigma * UtB)
        return Sinv_V @ y


def pseudo_inverse_solve_regularized(result: StreamTSVD, B, S, k_mode=None):
    """One-shot regularised pseudo-inverse solve.

    Builds ``RegularizedTSVD.from_stiffness(result, S)`` and applies
    it once.  Use ``RegularizedTSVD`` directly when the same ``S`` is
    reused across multiple ``B`` (= Path-A compensated iteration).

    Parameters
    ----------
    result : StreamTSVD
        ACA+TSVD factorisation of A.
    B : array_like, shape (M,)
        Target field values.
    S : ndarray or scipy.sparse matrix, shape (N, N)
        SPD regularisation matrix (e.g. H1 stiffness on free DOFs).
    k_mode : int, optional
        Truncate to first ``k_mode`` SVD modes (default = result.modes).

    Returns
    -------
    ndarray, shape (N,)
        ``phi`` satisfying ``A phi = B`` with minimum ``phi^T S phi``.
    """
    reg = RegularizedTSVD.from_stiffness(result, S)
    return reg.solve(B, k_mode=k_mode)


def solve(M, N, entry, B, modes=None, k_mode=None,
          kmax=None, aca_eps=1.0e-4, method="aca_qr_tsvd"):
    """Convenience: (ACA+)+TSVD decompose then pseudo-inverse solve.
    ``method`` selects "aca_qr_tsvd" (ACA+QR+TSVD, default) or "dense" (direct
    TSVD); see :func:`aca_tsvd`.

    Returns
    -------
    (phi, result) : (ndarray shape (N,), StreamTSVD)
    """
    result = aca_tsvd(M, N, entry, modes=modes, kmax=kmax,
                      aca_eps=aca_eps, method=method)
    phi = pseudo_inverse_solve(result, B, k_mode=k_mode)
    return phi, result


def radia_field_kernel(obs_points, sources, component=2, field="b"):
    """Build a matrix-entry callback from Radia's existing field computation.

    ``A(i, j)`` = (``component`` of the ``field``) at observation point i
    produced by Radia source object ``sources[j]``, evaluated via ``radia.Fld``.
    This is the generic reuse of Radia's already-implemented kernels and works
    for ANY source family -- coils (Biot-Savart), fixed-magnet objects, or
    solved magnetic materials -- without embedding a new field formula.

    Parameters
    ----------
    obs_points : array_like, shape (M, 3)
        Observation points.
    sources : sequence of int
        Length-N sequence of Radia object handles.
    component : int, optional
        Field component index 0/1/2 for x/y/z (default 2 = z).
    field : str, optional
        Radia field id passed to ``radia.Fld`` (default "b"; "h", "a", ...).

    Returns
    -------
    callable
        ``entry(i, j) -> float`` for use with ``aca_tsvd`` / ``solve``
        (M = len(obs_points), N = len(sources)).
    """
    import radia as rad

    obs = np.ascontiguousarray(obs_points, dtype=float).reshape(-1, 3)
    src = [int(s) for s in sources]
    comp = int(component)

    def entry(i, j):
        B = rad.Fld(src[j], field, obs[i].tolist())
        return float(np.asarray(B, dtype=float).ravel()[comp])

    return entry
