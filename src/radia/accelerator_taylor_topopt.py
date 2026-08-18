"""High-order Taylor-map topology optimization for accelerator magnets.

This module extends the planar first-order transfer-matrix objective to the
factorial Taylor convention

``u_out = R @ u + T[u,u]/2 + U[u,u,u]/6 + O(u**4)``.

The second-order contract contains five source-free transverse multipoles per
fixed-orbit segment.  The third-order contract adds normal/skew octupoles and
propagates direct cubic, chromatic, and lower-order cascade terms.  Forward-mode
algorithmic differentiation follows the same piecewise-constant RK4 Taylor-map
composition used by the native C++ beam kernel.  The native result remains the
value source of truth; finite differences are used only by regression tests.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations
from math import factorial

import numpy as np

from .accelerator_magnet_topopt import (
    PlanarDesignOrbit,
    TransferMatrixAutomaticDifferentiation,
    TransferMatrixFieldCorrection,
    differentiate_transfer_matrix_field_response,
    solve_transfer_matrix_field_correction,
)
from .isochronous_topopt import MU0
from .topology_optimization import (
    GrowthTopologyReport,
    HDivMMMGenerationResult,
    TSVDElementCandidateSelection,
    grow_hdiv_mmm_by_superposition,
    ngsolve_growth_topology,
    select_tsvd_element_candidates,
)

_EINSUM_PATHS: dict = {}


def _cached_einsum(subscripts, *operands, optimize=True, **kwargs):
    """``np.einsum`` with the contraction path cached per (subscripts, shapes).

    The Taylor/Lie map machinery issues tens of thousands of einsum calls over
    a handful of small fixed-shape tensor contractions; numpy re-plans the
    greedy path on every call, which costs more than the contractions
    themselves.  The first call per signature computes the exact path
    ``optimize=True`` would have used and every later call replays it, so
    results stay bit-identical while the planning cost vanishes.  The
    ``optimize`` argument is accepted for signature compatibility and ignored.
    """
    del optimize
    key = (subscripts, tuple(np.shape(operand) for operand in operands))
    path = _EINSUM_PATHS.get(key)
    if path is None:
        path = np.einsum_path(subscripts, *operands, optimize=True)[0]
        _EINSUM_PATHS[key] = path
    return np.einsum(subscripts, *operands, optimize=path, **kwargs)

SECOND_ORDER_MULTIPOLE_COMPONENTS = (
    "normal_dipole",
    "normal_quadrupole",
    "skew_quadrupole",
    "normal_sextupole",
    "skew_sextupole",
)
THIRD_ORDER_MULTIPOLE_COMPONENTS = SECOND_ORDER_MULTIPOLE_COMPONENTS + (
    "normal_octupole",
    "skew_octupole",
)
_ALL_R_ENTRIES = tuple((row, column) for row in range(6) for column in range(6))
_ALL_SYMMETRIC_T_ENTRIES = tuple(
    (output, first, second)
    for output in range(6)
    for first in range(6)
    for second in range(first, 6)
)
_ALL_SYMMETRIC_U_ENTRIES = tuple(
    (output, first, second, third)
    for output in range(6)
    for first in range(6)
    for second in range(first, 6)
    for third in range(second, 6)
)


def _finite_array(value, *, shape=None, name):
    result = np.asarray(value, dtype=float)
    if shape is not None and result.shape != shape:
        raise ValueError(f"{name} must have shape {shape}")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be finite")
    return result


def _checked_entries(entries, *, rank, name):
    result = tuple(tuple(int(value) for value in item) for item in entries)
    if (
        not result
        or len(set(result)) != len(result)
        or any(len(item) != rank for item in result)
        or any(index < 0 or index >= 6 for item in result for index in item)
    ):
        raise ValueError(f"{name} must contain unique zero-based rank-{rank} indices")
    if rank == 3 and any(first > second for _, first, second in result):
        raise ValueError(
            "T entries must use canonical symmetric input order first<=second"
        )
    if rank == 4 and any(
        first > second or second > third for _, first, second, third in result
    ):
        raise ValueError(
            "U entries must use canonical symmetric input order first<=second<=third"
        )
    return result


@dataclass(frozen=True)
class Symplectic2x2KAN:
    """Global KAN coordinates for one real two-dimensional symplectic map.

    The convention is ``M = K(rotation) A(log_scale) N(shear)`` with

    ``K = [[cos, sin],[-sin, cos]]``, ``A = diag(exp(e), exp(-e))``, and
    ``N = [[1, shear],[0, 1]]``.  Iwasawa/KAN coordinates cover every real
    ``2 x 2`` matrix of determinant one; the rotation is periodic modulo
    ``2*pi``.
    """

    rotation_rad: float = 0.0
    log_scale: float = 0.0
    shear: float = 0.0

    def __post_init__(self):
        values = np.asarray(
            [self.rotation_rad, self.log_scale, self.shear], dtype=float
        )
        if not np.all(np.isfinite(values)):
            raise ValueError("KAN parameters must be finite")
        object.__setattr__(self, "rotation_rad", float(values[0]))
        object.__setattr__(self, "log_scale", float(values[1]))
        object.__setattr__(self, "shear", float(values[2]))

    @property
    def parameters(self) -> np.ndarray:
        return np.asarray(
            [self.rotation_rad, self.log_scale, self.shear], dtype=float
        )

    @property
    def matrix(self) -> np.ndarray:
        cosine = np.cos(self.rotation_rad)
        sine = np.sin(self.rotation_rad)
        expansion = np.exp(self.log_scale)
        contraction = np.exp(-self.log_scale)
        return np.asarray(
            [
                [
                    cosine * expansion,
                    cosine * expansion * self.shear + sine * contraction,
                ],
                [
                    -sine * expansion,
                    -sine * expansion * self.shear + cosine * contraction,
                ],
            ],
            dtype=float,
        )

    @classmethod
    def from_matrix(cls, matrix, *, tolerance=1.0e-10):
        """Recover KAN coordinates after checking ``det(matrix) == 1``."""
        value = _finite_array(matrix, shape=(2, 2), name="symplectic block")
        tolerance = float(tolerance)
        if not np.isfinite(tolerance) or tolerance < 0.0:
            raise ValueError("symplectic tolerance must be finite and nonnegative")
        determinant = float(np.linalg.det(value))
        scale = max(1.0, float(np.max(np.abs(value))) ** 2)
        if abs(determinant - 1.0) > tolerance * scale:
            raise ValueError(
                "symplectic 2x2 block must have determinant one; "
                f"residual={abs(determinant - 1.0):.3e}"
            )
        first_column_squared = float(value[0, 0] ** 2 + value[1, 0] ** 2)
        if first_column_squared <= np.finfo(float).tiny:
            raise ValueError("symplectic block has a degenerate first column")
        rotation = float(np.arctan2(-value[1, 0], value[0, 0]))
        log_scale = 0.5 * float(np.log(first_column_squared))
        shear = float(
            (value[0, 0] * value[0, 1] + value[1, 0] * value[1, 1])
            / first_column_squared
        )
        return cls(rotation, log_scale, shear)


@dataclass(frozen=True)
class DecoupledFirstOrderTarget:
    """Physical block-diagonal first-order target for ``(x,px,y,py,ell,delta)``.

    Static magnetic optics preserve ``delta``.  The longitudinal block is
    therefore the one-parameter shear ``[[1, R56], [0, 1]]``; an arbitrary
    longitudinal ``Sp(2,R)`` block needs an RF/electric energy kick and is not
    represented by this target type.
    """

    horizontal: Symplectic2x2KAN
    vertical: Symplectic2x2KAN
    longitudinal_shear: float = 0.0

    def __post_init__(self):
        if not isinstance(self.horizontal, Symplectic2x2KAN) or not isinstance(
            self.vertical, Symplectic2x2KAN
        ):
            raise TypeError("horizontal and vertical must be Symplectic2x2KAN")
        shear = float(self.longitudinal_shear)
        if not np.isfinite(shear):
            raise ValueError("longitudinal_shear must be finite")
        object.__setattr__(self, "longitudinal_shear", shear)

    @property
    def parameters(self) -> np.ndarray:
        return np.r_[
            self.horizontal.parameters,
            self.vertical.parameters,
            self.longitudinal_shear,
        ]

    @property
    def matrix(self) -> np.ndarray:
        result = np.zeros((6, 6), dtype=float)
        result[0:2, 0:2] = self.horizontal.matrix
        result[2:4, 2:4] = self.vertical.matrix
        result[4:6, 4:6] = np.asarray(
            [[1.0, self.longitudinal_shear], [0.0, 1.0]]
        )
        return result

    @classmethod
    def from_matrix(cls, matrix, *, tolerance=1.0e-10):
        """Validate and parameterize a decoupled magnetostatic ``6 x 6`` map."""
        value = _finite_array(
            matrix, shape=(6, 6), name="decoupled first-order target"
        )
        tolerance = float(tolerance)
        if not np.isfinite(tolerance) or tolerance < 0.0:
            raise ValueError("target tolerance must be finite and nonnegative")
        scale = max(1.0, float(np.max(np.abs(value))))
        off_block = value.copy()
        for begin in (0, 2, 4):
            off_block[begin : begin + 2, begin : begin + 2] = 0.0
        off_block_residual = float(np.max(np.abs(off_block), initial=0.0))
        if off_block_residual > tolerance * scale:
            raise ValueError(
                "first-order target must decouple x, y, and longitudinal blocks; "
                f"residual={off_block_residual:.3e}"
            )
        longitudinal = value[4:6, 4:6]
        expected_longitudinal = np.asarray(
            [[1.0, longitudinal[0, 1]], [0.0, 1.0]]
        )
        longitudinal_residual = float(
            np.max(np.abs(longitudinal - expected_longitudinal), initial=0.0)
        )
        if longitudinal_residual > tolerance * scale:
            raise ValueError(
                "magnetostatic longitudinal target must be [[1,R56],[0,1]]; "
                f"residual={longitudinal_residual:.3e}"
            )
        return cls(
            Symplectic2x2KAN.from_matrix(value[0:2, 0:2], tolerance=tolerance),
            Symplectic2x2KAN.from_matrix(value[2:4, 2:4], tolerance=tolerance),
            float(longitudinal[0, 1]),
        )


def _unpack_multipoles(response, segment_count):
    values = _finite_array(response, name="multipole response").reshape(-1)
    expected = len(SECOND_ORDER_MULTIPOLE_COMPONENTS) * int(segment_count)
    if values.shape != (expected,):
        raise ValueError(
            "multipole response must contain five component blocks per orbit segment"
        )
    return values.reshape(len(SECOND_ORDER_MULTIPOLE_COMPONENTS), int(segment_count))


def _unpack_third_order_multipoles(response, segment_count):
    values = _finite_array(response, name="multipole response").reshape(-1)
    expected = len(THIRD_ORDER_MULTIPOLE_COMPONENTS) * int(segment_count)
    if values.shape != (expected,):
        raise ValueError(
            "multipole response must contain seven component blocks per orbit segment"
        )
    return values.reshape(len(THIRD_ORDER_MULTIPOLE_COMPONENTS), int(segment_count))


def _planar_orbit_multipole_observations(
    orbit: PlanarDesignOrbit, *, sample_radius, maximum_degree
) -> tuple[np.ndarray, np.ndarray]:
    """Return the shared nine-point harmonic stencil through ``maximum_degree``."""
    if maximum_degree not in (2, 3, 4):
        raise ValueError("maximum_degree must be two, three, or four")
    if not isinstance(orbit, PlanarDesignOrbit):
        raise TypeError("orbit must be a PlanarDesignOrbit")
    radius = float(sample_radius)
    scale = max(1.0, float(np.max(orbit.segment_lengths)))
    if not np.isfinite(radius) or radius <= 1.0e-12 * scale:
        raise ValueError("sample_radius must be a positive physical length")

    centers = orbit.sample_positions
    tangents = orbit.tangents[:-1] + orbit.tangents[1:]
    tangents /= np.linalg.norm(tangents, axis=1)[:, None]
    horizontal = np.cross(orbit.bend_axis[None, :], tangents)
    horizontal /= np.linalg.norm(horizontal, axis=1)[:, None]
    vertical = np.cross(tangents, horizontal)
    vertical /= np.linalg.norm(vertical, axis=1)[:, None]
    diagonal = radius / np.sqrt(2.0)
    offsets = np.asarray(
        [
            [0.0, 0.0],
            [radius, 0.0],
            [-radius, 0.0],
            [0.0, radius],
            [0.0, -radius],
            [diagonal, diagonal],
            [diagonal, -diagonal],
            [-diagonal, diagonal],
            [-diagonal, -diagonal],
        ]
    )
    count = len(centers)
    component_count = 1 + 2 * maximum_degree
    points = np.empty((9 * count, 3), dtype=float)
    weights = np.zeros((component_count * count, 9 * count, 3), dtype=float)
    for segment in range(count):
        begin = 9 * segment
        points[begin : begin + 9] = (
            centers[segment]
            + offsets[:, :1] * horizontal[segment]
            + offsets[:, 1:] * vertical[segment]
        )
        weights[segment, begin] = vertical[segment]
        for degree in range(1, maximum_degree + 1):
            normal_block = 2 * degree - 1
            skew_block = 2 * degree
            denominator = 8.0 * radius**degree
            for sample in range(1, 9):
                normalized = complex(
                    offsets[sample, 0] / radius, offsets[sample, 1] / radius
                )
                projection = np.conj(normalized**degree)
                real = float(projection.real) / denominator
                imag = float(projection.imag) / denominator
                weights[normal_block * count + segment, begin + sample] = (
                    real * vertical[segment] - imag * horizontal[segment]
                )
                weights[skew_block * count + segment, begin + sample] = (
                    imag * vertical[segment] + real * horizontal[segment]
                )
    return (np.ascontiguousarray(points), np.ascontiguousarray(weights))


def planar_orbit_multipole_observations(
    orbit: PlanarDesignOrbit, *, sample_radius
) -> tuple[np.ndarray, np.ndarray]:
    """Return native field points/weights for five multipole blocks.

    The nine-point stencil and complex harmonic projection match the native
    GridFunction multipole map.  Rows are ordered by component block according
    to :data:`SECOND_ORDER_MULTIPOLE_COMPONENTS`, then by segment.  Since the
    projection is linear in the sampled field, the returned weights are valid
    HDiv response functionals and require no finite-difference material step.
    """
    return _planar_orbit_multipole_observations(
        orbit, sample_radius=sample_radius, maximum_degree=2
    )


def planar_orbit_cubic_multipole_observations(
    orbit: PlanarDesignOrbit, *, sample_radius
) -> tuple[np.ndarray, np.ndarray]:
    """Return native field points/weights through normal/skew octupole."""
    return _planar_orbit_multipole_observations(
        orbit, sample_radius=sample_radius, maximum_degree=3
    )


def build_planar_orbit_multipole_response_matrix(
    charge_gram, orbit: PlanarDesignOrbit, *, sample_radius, field_scale=MU0
) -> np.ndarray:
    """Build native HDiv rows for dipole through skew-sextupole response."""
    native = getattr(charge_gram, "configured_field_functional_rows", None)
    if native is None:
        raise TypeError("charge_gram must expose configured_field_functional_rows")
    scale = float(field_scale)
    if not np.isfinite(scale) or scale == 0.0:
        raise ValueError("field_scale must be finite and nonzero")
    points, weights = planar_orbit_multipole_observations(
        orbit, sample_radius=sample_radius
    )
    rows = scale * np.asarray(native(points, weights), dtype=float)
    if (
        rows.ndim != 2
        or rows.shape[0] != 5 * len(orbit.segment_lengths)
        or not np.all(np.isfinite(rows))
    ):
        raise RuntimeError(
            "native configured-field API returned invalid multipole rows"
        )
    return np.ascontiguousarray(rows)


def build_planar_orbit_cubic_multipole_response_matrix(
    charge_gram, orbit: PlanarDesignOrbit, *, sample_radius, field_scale=MU0
) -> np.ndarray:
    """Build native HDiv rows for dipole through skew-octupole response."""
    native = getattr(charge_gram, "configured_field_functional_rows", None)
    if native is None:
        raise TypeError("charge_gram must expose configured_field_functional_rows")
    scale = float(field_scale)
    if not np.isfinite(scale) or scale == 0.0:
        raise ValueError("field_scale must be finite and nonzero")
    points, weights = planar_orbit_cubic_multipole_observations(
        orbit, sample_radius=sample_radius
    )
    rows = scale * np.asarray(native(points, weights), dtype=float)
    expected_rows = len(THIRD_ORDER_MULTIPOLE_COMPONENTS) * len(orbit.segment_lengths)
    if (
        rows.ndim != 2
        or rows.shape[0] != expected_rows
        or not np.all(np.isfinite(rows))
    ):
        raise RuntimeError(
            "native configured-field API returned invalid cubic multipole rows"
        )
    return np.ascontiguousarray(rows)


def _add_monomial(value, tangent, output, powers, coefficient, coefficient_tangent):
    degree = int(sum(powers))
    derivative = float(coefficient)
    for power in powers:
        derivative *= factorial(int(power))
    derivative_tangent = np.asarray(coefficient_tangent, dtype=float).copy()
    for power in powers:
        derivative_tangent *= factorial(int(power))
    inputs = []
    for coordinate, power in enumerate(powers):
        inputs.extend([coordinate] * int(power))
    if degree == 1:
        value[output, inputs[0]] += derivative
        tangent[:, output, inputs[0]] += derivative_tangent
    elif degree in (2, 3):
        for ordered in set(permutations(inputs)):
            value[(output,) + ordered] += derivative
            tangent[(slice(None), output) + ordered] += derivative_tangent


def _paraxial_taylor_jet(
    coefficients,
    magnetic_rigidity,
    curvature_sign,
    gradient_sign,
    *,
    maximum_order,
):
    maximum_order = int(maximum_order)
    if maximum_order not in (2, 3):
        raise ValueError("maximum_order must be two or three")
    parameter_count = 1 + 2 * maximum_order
    coefficients = _finite_array(
        coefficients,
        shape=(parameter_count,),
        name="segment multipoles",
    )
    rigidity = float(magnetic_rigidity)
    curvature_sign = float(curvature_sign)
    gradient_sign = float(gradient_sign)
    if (
        not np.isfinite(rigidity)
        or rigidity == 0.0
        or not np.all(np.isfinite([curvature_sign, gradient_sign]))
    ):
        raise ValueError("rigidity and multipole signs must be finite")

    A = np.zeros((6, 6), dtype=float)
    F2 = np.zeros((6, 6, 6), dtype=float)
    F3 = np.zeros((6, 6, 6, 6), dtype=float)
    dA = np.zeros((parameter_count, 6, 6), dtype=float)
    dF2 = np.zeros((parameter_count, 6, 6, 6), dtype=float)
    dF3 = np.zeros((parameter_count, 6, 6, 6, 6), dtype=float)

    normal_dipole = coefficients[0]
    curvature = curvature_sign * normal_dipole / rigidity
    dcurvature = np.zeros(parameter_count)
    dcurvature[0] = curvature_sign / rigidity
    A[0, 1] = 1.0
    A[1, 0] = -curvature * curvature
    dA[:, 1, 0] = -2.0 * curvature * dcurvature
    A[1, 5] = curvature
    dA[:, 1, 5] = dcurvature
    A[2, 3] = 1.0
    A[4, 0] = curvature
    dA[:, 4, 0] = dcurvature

    zero_tangent = np.zeros(parameter_count)
    _add_monomial(F2, dF2, 0, (0, 1, 0, 0, 0, 1), -1.0, zero_tangent)
    _add_monomial(F2, dF2, 2, (0, 0, 0, 1, 0, 1), -1.0, zero_tangent)
    if maximum_order >= 3:
        _add_monomial(F3, dF3, 0, (0, 1, 0, 0, 0, 2), 1.0, zero_tangent)
        _add_monomial(F3, dF3, 2, (0, 0, 0, 1, 0, 2), 1.0, zero_tangent)

    for order in range(1, maximum_order + 1):
        normal_parameter = 2 * order - 1
        skew_parameter = 2 * order
        multipole = complex(
            coefficients[normal_parameter], coefficients[skew_parameter]
        )
        multipole_tangent = np.zeros(parameter_count, dtype=complex)
        multipole_tangent[normal_parameter] = 1.0
        multipole_tangent[skew_parameter] = 1.0j
        for y_power in range(order + 1):
            polynomial_scale = (
                factorial(order)
                / (factorial(y_power) * factorial(order - y_power))
                * (1.0j**y_power)
            )
            polynomial = multipole * polynomial_scale
            polynomial_tangent = multipole_tangent * polynomial_scale
            for delta_power in range(maximum_order - order + 1):
                powers = [0] * 6
                powers[0] = order - y_power
                powers[2] = y_power
                powers[5] = delta_power
                chromatic = 1.0 if delta_power % 2 == 0 else -1.0
                horizontal_scale = -gradient_sign * chromatic / rigidity
                vertical_scale = gradient_sign * chromatic / rigidity
                degree = sum(powers)
                if degree == 1:
                    target_value, target_tangent = A, dA
                elif degree == 2:
                    target_value, target_tangent = F2, dF2
                else:
                    target_value, target_tangent = F3, dF3
                _add_monomial(
                    target_value,
                    target_tangent,
                    1,
                    tuple(powers),
                    horizontal_scale * polynomial.real,
                    horizontal_scale * polynomial_tangent.real,
                )
                _add_monomial(
                    target_value,
                    target_tangent,
                    3,
                    tuple(powers),
                    vertical_scale * polynomial.imag,
                    vertical_scale * polynomial_tangent.imag,
                )
    return A, F2, F3, dA, dF2, dF3


def _paraxial_second_order_jet(
    coefficients, magnetic_rigidity, curvature_sign, gradient_sign
):
    A, F2, _, dA, dF2, _ = _paraxial_taylor_jet(
        coefficients,
        magnetic_rigidity,
        curvature_sign,
        gradient_sign,
        maximum_order=2,
    )
    return A, F2, dA, dF2


def _paraxial_third_order_jet(
    coefficients, magnetic_rigidity, curvature_sign, gradient_sign
):
    return _paraxial_taylor_jet(
        coefficients,
        magnetic_rigidity,
        curvature_sign,
        gradient_sign,
        maximum_order=3,
    )


def _second_order_rhs(A, F2, dA, dF2, R, T, dR, dT):
    R_dot = _cached_einsum("ia,aj->ij", A, R, optimize=True)
    T_dot = _cached_einsum("ia,ajk->ijk", A, T, optimize=True) + _cached_einsum(
        "iab,aj,bk->ijk", F2, R, R, optimize=True
    )
    dR_dot = _cached_einsum("ia,paj->pij", A, dR, optimize=True) + _cached_einsum(
        "pia,aj->pij", dA, R, optimize=True
    )
    dT_dot = (
        _cached_einsum("ia,pajk->pijk", A, dT, optimize=True)
        + _cached_einsum("pia,ajk->pijk", dA, T, optimize=True)
        + _cached_einsum("piab,aj,bk->pijk", dF2, R, R, optimize=True)
        + _cached_einsum("iab,paj,bk->pijk", F2, dR, R, optimize=True)
        + _cached_einsum("iab,aj,pbk->pijk", F2, R, dR, optimize=True)
    )
    return R_dot, T_dot, dR_dot, dT_dot


def _second_order_rk4_step(A, F2, dA, dF2, length):
    parameter_count = dA.shape[0]
    state = (
        np.eye(6),
        np.zeros((6, 6, 6)),
        np.zeros((parameter_count, 6, 6)),
        np.zeros((parameter_count, 6, 6, 6)),
    )

    def rhs(value):
        return _second_order_rhs(A, F2, dA, dF2, *value)

    def shifted(value, derivative, scale):
        return tuple(item + scale * rate for item, rate in zip(value, derivative))

    k1 = rhs(state)
    k2 = rhs(shifted(state, k1, 0.5 * length))
    k3 = rhs(shifted(state, k2, 0.5 * length))
    k4 = rhs(shifted(state, k3, length))
    return tuple(
        item + length * (first + 2.0 * second + 2.0 * third + fourth) / 6.0
        for item, first, second, third, fourth in zip(state, k1, k2, k3, k4)
    )


def _compose_second_order(outer, inner):
    Ro, To, dRo, dTo = outer
    Ri, Ti, dRi, dTi = inner
    R = Ro @ Ri
    T = _cached_einsum("ia,ajk->ijk", Ro, Ti, optimize=True) + _cached_einsum(
        "iab,aj,bk->ijk", To, Ri, Ri, optimize=True
    )
    dR = _cached_einsum("pia,aj->pij", dRo, Ri, optimize=True) + _cached_einsum(
        "ia,paj->pij", Ro, dRi, optimize=True
    )
    dT = (
        _cached_einsum("pia,ajk->pijk", dRo, Ti, optimize=True)
        + _cached_einsum("ia,pajk->pijk", Ro, dTi, optimize=True)
        + _cached_einsum("piab,aj,bk->pijk", dTo, Ri, Ri, optimize=True)
        + _cached_einsum("iab,paj,bk->pijk", To, dRi, Ri, optimize=True)
        + _cached_einsum("iab,aj,pbk->pijk", To, Ri, dRi, optimize=True)
    )
    return R, T, dR, dT


def _identity_second_order(parameter_count):
    return (
        np.eye(6),
        np.zeros((6, 6, 6)),
        np.zeros((int(parameter_count), 6, 6)),
        np.zeros((int(parameter_count), 6, 6, 6)),
    )


def _cross_second_order(F2, R, T):
    return (
        _cached_einsum("iab,aj,bkl->ijkl", F2, R, T, optimize=True)
        + _cached_einsum("iab,ak,bjl->ijkl", F2, R, T, optimize=True)
        + _cached_einsum("iab,al,bjk->ijkl", F2, R, T, optimize=True)
    )


def _cross_second_order_tangent(F2, R, T, dF2, dR, dT):
    return (
        _cached_einsum("piab,aj,bkl->pijkl", dF2, R, T, optimize=True)
        + _cached_einsum("piab,ak,bjl->pijkl", dF2, R, T, optimize=True)
        + _cached_einsum("piab,al,bjk->pijkl", dF2, R, T, optimize=True)
        + _cached_einsum("iab,paj,bkl->pijkl", F2, dR, T, optimize=True)
        + _cached_einsum("iab,pak,bjl->pijkl", F2, dR, T, optimize=True)
        + _cached_einsum("iab,pal,bjk->pijkl", F2, dR, T, optimize=True)
        + _cached_einsum("iab,aj,pbkl->pijkl", F2, R, dT, optimize=True)
        + _cached_einsum("iab,ak,pbjl->pijkl", F2, R, dT, optimize=True)
        + _cached_einsum("iab,al,pbjk->pijkl", F2, R, dT, optimize=True)
    )


def _transform_cubic(F3, R):
    return _cached_einsum("iabc,aj,bk,cl->ijkl", F3, R, R, R, optimize=True)


def _transform_cubic_tangent(F3, R, dF3, dR):
    return (
        _cached_einsum("piabc,aj,bk,cl->pijkl", dF3, R, R, R, optimize=True)
        + _cached_einsum("iabc,paj,bk,cl->pijkl", F3, dR, R, R, optimize=True)
        + _cached_einsum("iabc,aj,pbk,cl->pijkl", F3, R, dR, R, optimize=True)
        + _cached_einsum("iabc,aj,bk,pcl->pijkl", F3, R, R, dR, optimize=True)
    )


def _third_order_rhs(A, F2, F3, dA, dF2, dF3, R, T, U, dR, dT, dU):
    R_dot, T_dot, dR_dot, dT_dot = _second_order_rhs(A, F2, dA, dF2, R, T, dR, dT)
    U_dot = (
        _cached_einsum("ia,ajkl->ijkl", A, U, optimize=True)
        + _cross_second_order(F2, R, T)
        + _transform_cubic(F3, R)
    )
    dU_dot = (
        _cached_einsum("ia,pajkl->pijkl", A, dU, optimize=True)
        + _cached_einsum("pia,ajkl->pijkl", dA, U, optimize=True)
        + _cross_second_order_tangent(F2, R, T, dF2, dR, dT)
        + _transform_cubic_tangent(F3, R, dF3, dR)
    )
    return R_dot, T_dot, U_dot, dR_dot, dT_dot, dU_dot


def _third_order_rk4_step(A, F2, F3, dA, dF2, dF3, length):
    parameter_count = dA.shape[0]
    state = _identity_third_order(parameter_count)

    def rhs(value):
        return _third_order_rhs(A, F2, F3, dA, dF2, dF3, *value)

    def shifted(value, derivative, scale):
        return tuple(item + scale * rate for item, rate in zip(value, derivative))

    k1 = rhs(state)
    k2 = rhs(shifted(state, k1, 0.5 * length))
    k3 = rhs(shifted(state, k2, 0.5 * length))
    k4 = rhs(shifted(state, k3, length))
    return tuple(
        item + length * (first + 2.0 * second + 2.0 * third + fourth) / 6.0
        for item, first, second, third, fourth in zip(state, k1, k2, k3, k4)
    )


def _compose_third_order(outer, inner):
    Ro, To, Uo, dRo, dTo, dUo = outer
    Ri, Ti, Ui, dRi, dTi, dUi = inner
    R = Ro @ Ri
    T = _cached_einsum("ia,ajk->ijk", Ro, Ti, optimize=True) + _cached_einsum(
        "iab,aj,bk->ijk", To, Ri, Ri, optimize=True
    )
    U = (
        _cached_einsum("ia,ajkl->ijkl", Ro, Ui, optimize=True)
        + _cross_second_order(To, Ri, Ti)
        + _transform_cubic(Uo, Ri)
    )
    dR = _cached_einsum("pia,aj->pij", dRo, Ri, optimize=True) + _cached_einsum(
        "ia,paj->pij", Ro, dRi, optimize=True
    )
    dT = (
        _cached_einsum("pia,ajk->pijk", dRo, Ti, optimize=True)
        + _cached_einsum("ia,pajk->pijk", Ro, dTi, optimize=True)
        + _cached_einsum("piab,aj,bk->pijk", dTo, Ri, Ri, optimize=True)
        + _cached_einsum("iab,paj,bk->pijk", To, dRi, Ri, optimize=True)
        + _cached_einsum("iab,aj,pbk->pijk", To, Ri, dRi, optimize=True)
    )
    dU = (
        _cached_einsum("pia,ajkl->pijkl", dRo, Ui, optimize=True)
        + _cached_einsum("ia,pajkl->pijkl", Ro, dUi, optimize=True)
        + _cross_second_order_tangent(To, Ri, Ti, dTo, dRi, dTi)
        + _transform_cubic_tangent(Uo, Ri, dUo, dRi)
    )
    return R, T, U, dR, dT, dU


def _identity_third_order(parameter_count):
    return (
        np.eye(6),
        np.zeros((6, 6, 6)),
        np.zeros((6, 6, 6, 6)),
        np.zeros((int(parameter_count), 6, 6)),
        np.zeros((int(parameter_count), 6, 6, 6)),
        np.zeros((int(parameter_count), 6, 6, 6, 6)),
    )


@dataclass(frozen=True)
class SecondOrderTaylorMap:
    """Native-valued ``R/T`` map and forward-mode multipole Jacobians."""

    R: np.ndarray
    T: np.ndarray
    R_jacobian: np.ndarray
    T_jacobian: np.ndarray
    multipole_components: tuple[str, ...]
    derivative_backend: str = "forward-mode-rk4-taylor-ad"
    value_backend: str = "native-cpp-variational-map"


def second_order_taylor_map_from_multipoles(
    multipole_response,
    segment_lengths,
    magnetic_rigidity,
    *,
    curvature_sign=1.0,
    gradient_sign=1.0,
    maximum_step_m=1.0e-3,
    maximum_steps=1_000_000,
) -> SecondOrderTaylorMap:
    """Propagate an AD-differentiated second-order canonical Taylor map."""
    lengths = _finite_array(segment_lengths, name="segment_lengths").reshape(-1)
    if lengths.size == 0 or np.any(lengths <= 0.0):
        raise ValueError("segment_lengths must be non-empty and positive")
    coefficients = _unpack_multipoles(multipole_response, lengths.size)
    step_limit = float(maximum_step_m)
    step_cap = int(maximum_steps)
    if not np.isfinite(step_limit) or step_limit <= 0.0 or step_cap < 1:
        raise ValueError("Taylor integration limits must be positive")
    parameter_count = coefficients.size
    accumulated = _identity_second_order(parameter_count)
    A_values = np.empty((lengths.size, 6, 6))
    F2_values = np.empty((lengths.size, 6, 6, 6))
    total_steps = 0
    for segment, length in enumerate(lengths):
        local_coefficients = coefficients[:, segment]
        A, F2, dA, dF2 = _paraxial_second_order_jet(
            local_coefficients, magnetic_rigidity, curvature_sign, gradient_sign
        )
        A_values[segment] = A
        F2_values[segment] = F2
        step_count = int(np.ceil(length / step_limit))
        total_steps += step_count
        if total_steps > step_cap:
            raise ValueError("Taylor integration exceeds maximum_steps")
        step = _second_order_rk4_step(A, F2, dA, dF2, length / step_count)
        local = _identity_second_order(5)
        for _ in range(step_count):
            local = _compose_second_order(step, local)
        embedded = _identity_second_order(parameter_count)
        embedded = (
            local[0],
            local[1],
            np.zeros_like(embedded[2]),
            np.zeros_like(embedded[3]),
        )
        indexes = np.asarray(
            [block * lengths.size + segment for block in range(5)], dtype=np.int64
        )
        embedded[2][indexes] = local[2]
        embedded[3][indexes] = local[3]
        accumulated = _compose_second_order(embedded, accumulated)

    # The native C++ implementation remains the value source of truth.  The
    # Python path above exists to propagate exact tangents through the same
    # RK4/composition algebra and must agree before a result is exposed.
    from .beam import propagate_variational_map

    native = propagate_variational_map(
        lengths,
        A_values,
        F2_values,
        maximum_order=2,
        maximum_step_m=step_limit,
        maximum_steps=step_cap,
    )
    native_R = np.asarray(native["R"], dtype=float)
    native_T = np.asarray(native["T"], dtype=float)
    if not np.allclose(
        accumulated[0], native_R, rtol=2.0e-12, atol=2.0e-13
    ) or not np.allclose(accumulated[1], native_T, rtol=3.0e-11, atol=3.0e-12):
        raise RuntimeError(
            "forward-AD Taylor propagation disagrees with the native value"
        )
    return SecondOrderTaylorMap(
        R=native_R,
        T=native_T,
        R_jacobian=np.asarray(accumulated[2], dtype=float),
        T_jacobian=np.asarray(accumulated[3], dtype=float),
        multipole_components=SECOND_ORDER_MULTIPOLE_COMPONENTS,
    )


@dataclass(frozen=True)
class ThirdOrderTaylorMap:
    """Native-valued ``R/T/U`` map and forward-mode multipole Jacobians."""

    R: np.ndarray
    T: np.ndarray
    U: np.ndarray
    R_jacobian: np.ndarray
    T_jacobian: np.ndarray
    U_jacobian: np.ndarray
    multipole_components: tuple[str, ...]
    derivative_backend: str = "forward-mode-rk4-taylor-ad"
    value_backend: str = "native-cpp-variational-map"


def third_order_taylor_map_from_multipoles(
    multipole_response,
    segment_lengths,
    magnetic_rigidity,
    *,
    curvature_sign=1.0,
    gradient_sign=1.0,
    maximum_step_m=1.0e-3,
    maximum_steps=1_000_000,
) -> ThirdOrderTaylorMap:
    """Propagate an AD-differentiated third-order canonical Taylor map."""
    lengths = _finite_array(segment_lengths, name="segment_lengths").reshape(-1)
    if lengths.size == 0 or np.any(lengths <= 0.0):
        raise ValueError("segment_lengths must be non-empty and positive")
    coefficients = _unpack_third_order_multipoles(multipole_response, lengths.size)
    step_limit = float(maximum_step_m)
    step_cap = int(maximum_steps)
    if not np.isfinite(step_limit) or step_limit <= 0.0 or step_cap < 1:
        raise ValueError("Taylor integration limits must be positive")
    parameter_count = coefficients.size
    accumulated = _identity_third_order(parameter_count)
    A_values = np.empty((lengths.size, 6, 6))
    F2_values = np.empty((lengths.size, 6, 6, 6))
    F3_values = np.empty((lengths.size, 6, 6, 6, 6))
    total_steps = 0
    for segment, length in enumerate(lengths):
        local_coefficients = coefficients[:, segment]
        A, F2, F3, dA, dF2, dF3 = _paraxial_third_order_jet(
            local_coefficients,
            magnetic_rigidity,
            curvature_sign,
            gradient_sign,
        )
        A_values[segment] = A
        F2_values[segment] = F2
        F3_values[segment] = F3
        step_count = int(np.ceil(length / step_limit))
        total_steps += step_count
        if total_steps > step_cap:
            raise ValueError("Taylor integration exceeds maximum_steps")
        step = _third_order_rk4_step(A, F2, F3, dA, dF2, dF3, length / step_count)
        local = _identity_third_order(7)
        for _ in range(step_count):
            local = _compose_third_order(step, local)
        embedded = _identity_third_order(parameter_count)
        embedded = (
            local[0],
            local[1],
            local[2],
            np.zeros_like(embedded[3]),
            np.zeros_like(embedded[4]),
            np.zeros_like(embedded[5]),
        )
        indexes = np.asarray(
            [block * lengths.size + segment for block in range(7)], dtype=np.int64
        )
        embedded[3][indexes] = local[3]
        embedded[4][indexes] = local[4]
        embedded[5][indexes] = local[5]
        accumulated = _compose_third_order(embedded, accumulated)

    from .beam import propagate_variational_map

    native = propagate_variational_map(
        lengths,
        A_values,
        F2_values,
        F3_values,
        maximum_order=3,
        maximum_step_m=step_limit,
        maximum_steps=step_cap,
    )
    native_R = np.asarray(native["R"], dtype=float)
    native_T = np.asarray(native["T"], dtype=float)
    native_U = np.asarray(native["U"], dtype=float)
    if (
        not np.allclose(accumulated[0], native_R, rtol=2.0e-12, atol=2.0e-13)
        or not np.allclose(accumulated[1], native_T, rtol=3.0e-11, atol=3.0e-12)
        or not np.allclose(accumulated[2], native_U, rtol=8.0e-10, atol=8.0e-11)
    ):
        raise RuntimeError(
            "forward-AD cubic Taylor propagation disagrees with the native value"
        )
    return ThirdOrderTaylorMap(
        R=native_R,
        T=native_T,
        U=native_U,
        R_jacobian=np.asarray(accumulated[3], dtype=float),
        T_jacobian=np.asarray(accumulated[4], dtype=float),
        U_jacobian=np.asarray(accumulated[5], dtype=float),
        multipole_components=THIRD_ORDER_MULTIPOLE_COMPONENTS,
    )


@dataclass(frozen=True)
class PlanarSecondOrderTaylorMapObjective:
    """Fixed-orbit objective for selected canonical ``R`` and ``T`` entries."""

    orbit: PlanarDesignOrbit
    target_R: np.ndarray
    target_T: np.ndarray
    R_band: np.ndarray | float
    T_band: np.ndarray | float
    normal_dipole_band: np.ndarray | float
    R_entries: tuple[tuple[int, int], ...] = _ALL_R_ENTRIES
    T_entries: tuple[tuple[int, int, int], ...] = _ALL_SYMMETRIC_T_ENTRIES
    curvature_sign: float = 1.0
    gradient_sign: float = 1.0
    maximum_step_m: float = 1.0e-3

    def __post_init__(self):
        if not isinstance(self.orbit, PlanarDesignOrbit):
            raise TypeError("orbit must be a PlanarDesignOrbit")
        target_R = _finite_array(self.target_R, shape=(6, 6), name="target_R")
        target_T = _finite_array(self.target_T, shape=(6, 6, 6), name="target_T")
        symmetry = float(np.max(np.abs(target_T - np.swapaxes(target_T, 1, 2))))
        symmetry_scale = max(1.0, float(np.max(np.abs(target_T))))
        if symmetry > 1.0e-12 * symmetry_scale:
            raise ValueError("target_T input indices must be symmetric")
        R_band = np.broadcast_to(np.asarray(self.R_band, dtype=float), (6, 6)).copy()
        T_band = np.broadcast_to(np.asarray(self.T_band, dtype=float), (6, 6, 6)).copy()
        count = len(self.orbit.segment_lengths)
        dipole_band = np.broadcast_to(
            np.asarray(self.normal_dipole_band, dtype=float), (count,)
        ).copy()
        R_entries = _checked_entries(self.R_entries, rank=2, name="R_entries")
        T_entries = _checked_entries(self.T_entries, rank=3, name="T_entries")
        if (
            not np.all(np.isfinite(R_band))
            or np.any(R_band <= 0.0)
            or not np.all(np.isfinite(T_band))
            or np.any(T_band <= 0.0)
            or not np.all(np.isfinite(dipole_band))
            or np.any(dipole_band <= 0.0)
            or not np.isfinite(float(self.maximum_step_m))
            or float(self.maximum_step_m) <= 0.0
        ):
            raise ValueError("Taylor-map bands and integration step are invalid")
        object.__setattr__(self, "target_R", target_R.copy())
        object.__setattr__(self, "target_T", target_T.copy())
        object.__setattr__(self, "R_band", R_band)
        object.__setattr__(self, "T_band", T_band)
        object.__setattr__(self, "normal_dipole_band", dipole_band)
        object.__setattr__(self, "R_entries", R_entries)
        object.__setattr__(self, "T_entries", T_entries)

    @property
    def raw_field_response_size(self) -> int:
        return 5 * len(self.orbit.segment_lengths)

    @property
    def derivative_backend(self) -> str:
        return "forward-mode-rk4-taylor-ad"

    @property
    def required_normal_dipole(self) -> np.ndarray:
        return (
            self.orbit.magnetic_rigidity
            * self.orbit.signed_curvature
            / self.curvature_sign
        )

    @property
    def response_target(self) -> np.ndarray:
        return np.r_[
            self.required_normal_dipole,
            [self.target_R[index] for index in self.R_entries],
            [self.target_T[index] for index in self.T_entries],
        ]

    @property
    def response_band(self) -> np.ndarray:
        return np.r_[
            self.normal_dipole_band,
            [self.R_band[index] for index in self.R_entries],
            [self.T_band[index] for index in self.T_entries],
        ]

    @property
    def response_slices(self) -> tuple[tuple[str, slice], ...]:
        count = len(self.orbit.segment_lengths)
        R_stop = count + len(self.R_entries)
        return (
            ("normal_dipole", slice(0, count)),
            ("R", slice(count, R_stop)),
            ("T", slice(R_stop, R_stop + len(self.T_entries))),
        )

    def evaluate_taylor_map(self, field_response) -> SecondOrderTaylorMap:
        return second_order_taylor_map_from_multipoles(
            field_response,
            self.orbit.segment_lengths,
            self.orbit.magnetic_rigidity,
            curvature_sign=self.curvature_sign,
            gradient_sign=self.gradient_sign,
            maximum_step_m=self.maximum_step_m,
        )

    def transform(self, field_response) -> np.ndarray:
        values = _unpack_multipoles(field_response, len(self.orbit.segment_lengths))
        transfer = self.evaluate_taylor_map(values.reshape(-1))
        return np.r_[
            values[0],
            [transfer.R[index] for index in self.R_entries],
            [transfer.T[index] for index in self.T_entries],
        ]

    def transform_jacobian(self, field_response) -> np.ndarray:
        transfer = self.evaluate_taylor_map(field_response)
        count = len(self.orbit.segment_lengths)
        dipole = np.zeros((count, self.raw_field_response_size))
        dipole[np.arange(count), np.arange(count)] = 1.0
        R_rows = np.asarray(
            [transfer.R_jacobian[(slice(None),) + index] for index in self.R_entries]
        )
        T_rows = np.asarray(
            [transfer.T_jacobian[(slice(None),) + index] for index in self.T_entries]
        )
        return np.vstack((dipole, R_rows, T_rows))


@dataclass(frozen=True)
class PlanarThirdOrderTaylorMapObjective:
    """Fixed-orbit objective for selected canonical ``R``, ``T``, and ``U``."""

    orbit: PlanarDesignOrbit
    target_R: np.ndarray
    target_T: np.ndarray
    target_U: np.ndarray
    R_band: np.ndarray | float
    T_band: np.ndarray | float
    U_band: np.ndarray | float
    normal_dipole_band: np.ndarray | float
    R_entries: tuple[tuple[int, int], ...] = _ALL_R_ENTRIES
    T_entries: tuple[tuple[int, int, int], ...] = _ALL_SYMMETRIC_T_ENTRIES
    U_entries: tuple[tuple[int, int, int, int], ...] = _ALL_SYMMETRIC_U_ENTRIES
    curvature_sign: float = 1.0
    gradient_sign: float = 1.0
    maximum_step_m: float = 1.0e-3

    def __post_init__(self):
        if not isinstance(self.orbit, PlanarDesignOrbit):
            raise TypeError("orbit must be a PlanarDesignOrbit")
        target_R = _finite_array(self.target_R, shape=(6, 6), name="target_R")
        target_T = _finite_array(self.target_T, shape=(6, 6, 6), name="target_T")
        target_U = _finite_array(self.target_U, shape=(6, 6, 6, 6), name="target_U")
        T_defect = float(np.max(np.abs(target_T - np.swapaxes(target_T, 1, 2))))
        T_scale = max(1.0, float(np.max(np.abs(target_T))))
        U_defect = max(
            float(np.max(np.abs(target_U - np.transpose(target_U, axes))))
            for axes in (
                (0, 1, 3, 2),
                (0, 2, 1, 3),
                (0, 3, 2, 1),
            )
        )
        U_scale = max(1.0, float(np.max(np.abs(target_U))))
        if T_defect > 1.0e-12 * T_scale:
            raise ValueError("target_T input indices must be symmetric")
        if U_defect > 1.0e-12 * U_scale:
            raise ValueError("target_U input indices must be symmetric")
        R_band = np.broadcast_to(np.asarray(self.R_band, dtype=float), (6, 6)).copy()
        T_band = np.broadcast_to(np.asarray(self.T_band, dtype=float), (6, 6, 6)).copy()
        U_band = np.broadcast_to(
            np.asarray(self.U_band, dtype=float), (6, 6, 6, 6)
        ).copy()
        count = len(self.orbit.segment_lengths)
        dipole_band = np.broadcast_to(
            np.asarray(self.normal_dipole_band, dtype=float), (count,)
        ).copy()
        R_entries = _checked_entries(self.R_entries, rank=2, name="R_entries")
        T_entries = _checked_entries(self.T_entries, rank=3, name="T_entries")
        U_entries = _checked_entries(self.U_entries, rank=4, name="U_entries")
        if (
            not np.all(np.isfinite(R_band))
            or np.any(R_band <= 0.0)
            or not np.all(np.isfinite(T_band))
            or np.any(T_band <= 0.0)
            or not np.all(np.isfinite(U_band))
            or np.any(U_band <= 0.0)
            or not np.all(np.isfinite(dipole_band))
            or np.any(dipole_band <= 0.0)
            or not np.isfinite(float(self.maximum_step_m))
            or float(self.maximum_step_m) <= 0.0
        ):
            raise ValueError("Taylor-map bands and integration step are invalid")
        object.__setattr__(self, "target_R", target_R.copy())
        object.__setattr__(self, "target_T", target_T.copy())
        object.__setattr__(self, "target_U", target_U.copy())
        object.__setattr__(self, "R_band", R_band)
        object.__setattr__(self, "T_band", T_band)
        object.__setattr__(self, "U_band", U_band)
        object.__setattr__(self, "normal_dipole_band", dipole_band)
        object.__setattr__(self, "R_entries", R_entries)
        object.__setattr__(self, "T_entries", T_entries)
        object.__setattr__(self, "U_entries", U_entries)

    @property
    def raw_field_response_size(self) -> int:
        return 7 * len(self.orbit.segment_lengths)

    @property
    def derivative_backend(self) -> str:
        return "forward-mode-rk4-taylor-ad"

    @property
    def required_normal_dipole(self) -> np.ndarray:
        return (
            self.orbit.magnetic_rigidity
            * self.orbit.signed_curvature
            / self.curvature_sign
        )

    @property
    def response_target(self) -> np.ndarray:
        return np.r_[
            self.required_normal_dipole,
            [self.target_R[index] for index in self.R_entries],
            [self.target_T[index] for index in self.T_entries],
            [self.target_U[index] for index in self.U_entries],
        ]

    @property
    def response_band(self) -> np.ndarray:
        return np.r_[
            self.normal_dipole_band,
            [self.R_band[index] for index in self.R_entries],
            [self.T_band[index] for index in self.T_entries],
            [self.U_band[index] for index in self.U_entries],
        ]

    @property
    def response_slices(self) -> tuple[tuple[str, slice], ...]:
        count = len(self.orbit.segment_lengths)
        R_stop = count + len(self.R_entries)
        T_stop = R_stop + len(self.T_entries)
        return (
            ("normal_dipole", slice(0, count)),
            ("R", slice(count, R_stop)),
            ("T", slice(R_stop, T_stop)),
            ("U", slice(T_stop, T_stop + len(self.U_entries))),
        )

    def evaluate_taylor_map(self, field_response) -> ThirdOrderTaylorMap:
        return third_order_taylor_map_from_multipoles(
            field_response,
            self.orbit.segment_lengths,
            self.orbit.magnetic_rigidity,
            curvature_sign=self.curvature_sign,
            gradient_sign=self.gradient_sign,
            maximum_step_m=self.maximum_step_m,
        )

    def transform(self, field_response) -> np.ndarray:
        values = _unpack_third_order_multipoles(
            field_response, len(self.orbit.segment_lengths)
        )
        transfer = self.evaluate_taylor_map(values.reshape(-1))
        return np.r_[
            values[0],
            [transfer.R[index] for index in self.R_entries],
            [transfer.T[index] for index in self.T_entries],
            [transfer.U[index] for index in self.U_entries],
        ]

    def transform_jacobian(self, field_response) -> np.ndarray:
        transfer = self.evaluate_taylor_map(field_response)
        count = len(self.orbit.segment_lengths)
        dipole = np.zeros((count, self.raw_field_response_size))
        dipole[np.arange(count), np.arange(count)] = 1.0
        R_rows = np.asarray(
            [transfer.R_jacobian[(slice(None),) + index] for index in self.R_entries]
        )
        T_rows = np.asarray(
            [transfer.T_jacobian[(slice(None),) + index] for index in self.T_entries]
        )
        U_rows = np.asarray(
            [transfer.U_jacobian[(slice(None),) + index] for index in self.U_entries]
        )
        return np.vstack((dipole, R_rows, T_rows, U_rows))


@dataclass(frozen=True)
class TaylorMapReachabilityCertificate:
    """Local TSVD certificate for the reachable part of a map target."""

    numerical_rank: int
    singular_values: np.ndarray
    parameter_step: np.ndarray
    predicted_response: np.ndarray
    normalized_target_residual: np.ndarray
    normalized_unreachable_residual: np.ndarray
    max_unreachable_band_ratio: float
    relative_unreachable_residual: float
    component_max_unreachable_band_ratios: tuple[tuple[str, float], ...]
    linearized_reachable: bool
    relative_tolerance: float


def certify_taylor_map_reachability(
    objective,
    current_multipole_response,
    *,
    field_basis=None,
    relative_tolerance=1.0e-10,
    acceptance_band_ratio=1.0,
) -> TaylorMapReachabilityCertificate:
    """Project a target error onto the local AD-reachable response subspace.

    The certificate is deliberately local: it proves what the current
    linearization and declared field basis can or cannot correct.  It does not
    claim global reachability of a strongly nonlinear target.
    """
    required = ("response_target", "response_band", "response_slices")
    if any(not hasattr(objective, name) for name in required):
        raise TypeError("objective must expose Taylor-map response metadata")
    tolerance = float(relative_tolerance)
    acceptance = float(acceptance_band_ratio)
    if (
        not np.isfinite(tolerance)
        or tolerance < 0.0
        or not np.isfinite(acceptance)
        or acceptance < 0.0
    ):
        raise ValueError("reachability tolerances must be finite and nonnegative")
    automatic = differentiate_transfer_matrix_field_response(
        objective,
        current_multipole_response,
        field_basis=field_basis,
    )
    target = _finite_array(objective.response_target, name="target response").reshape(
        -1
    )
    band = _finite_array(objective.response_band, name="response band").reshape(-1)
    if target.shape != automatic.design_response.shape or band.shape != target.shape:
        raise RuntimeError("objective target metadata is inconsistent")
    normalized = (target - automatic.design_response) / band
    normalized_jacobian = automatic.directional_jacobian / band[:, None]
    left, singular_values, right = np.linalg.svd(
        normalized_jacobian, full_matrices=False
    )
    threshold = tolerance * singular_values[0] if singular_values.size else 0.0
    rank = int(np.count_nonzero(singular_values > threshold))
    if rank:
        reachable_part = left[:, :rank] @ (left[:, :rank].T @ normalized)
        parameter_step = right[:rank].T @ (
            (left[:, :rank].T @ normalized) / singular_values[:rank]
        )
    else:
        reachable_part = np.zeros_like(normalized)
        parameter_step = np.zeros(automatic.field_basis.shape[1])
    unreachable = normalized - reachable_part
    predicted = automatic.design_response + band * reachable_part
    max_ratio = float(np.max(np.abs(unreachable), initial=0.0))
    relative = float(
        np.linalg.norm(unreachable) / max(np.linalg.norm(normalized), 1.0e-300)
    )
    component_ratios = tuple(
        (
            str(name),
            float(np.max(np.abs(unreachable[part]), initial=0.0)),
        )
        for name, part in objective.response_slices
    )
    return TaylorMapReachabilityCertificate(
        numerical_rank=rank,
        singular_values=np.asarray(singular_values, dtype=float),
        parameter_step=np.asarray(parameter_step, dtype=float),
        predicted_response=np.asarray(predicted, dtype=float),
        normalized_target_residual=np.asarray(normalized, dtype=float),
        normalized_unreachable_residual=np.asarray(unreachable, dtype=float),
        max_unreachable_band_ratio=max_ratio,
        relative_unreachable_residual=relative,
        component_max_unreachable_band_ratios=component_ratios,
        linearized_reachable=bool(max_ratio <= acceptance),
        relative_tolerance=tolerance,
    )


def _kan_parameters_and_entry_jacobian(block):
    """Return KAN coordinates and their analytic derivative by block entry."""
    value = _finite_array(block, shape=(2, 2), name="transverse map block")
    first, second = value[0]
    third, fourth = value[1]
    column_norm_squared = float(first * first + third * third)
    if column_norm_squared <= np.finfo(float).tiny:
        raise ValueError("transverse map block has a degenerate first column")
    shear = float(
        (first * second + third * fourth) / column_norm_squared
    )
    parameters = np.asarray(
        [
            np.arctan2(-third, first),
            0.5 * np.log(column_norm_squared),
            shear,
        ],
        dtype=float,
    )
    entry_jacobian = np.asarray(
        [
            [third / column_norm_squared, 0.0,
             -first / column_norm_squared, 0.0],
            [first / column_norm_squared, 0.0,
             third / column_norm_squared, 0.0],
            [
                (second - 2.0 * shear * first) / column_norm_squared,
                first / column_norm_squared,
                (fourth - 2.0 * shear * third) / column_norm_squared,
                third / column_norm_squared,
            ],
        ],
        dtype=float,
    )
    return parameters, entry_jacobian


def _canonical_linear_symplectic_residual(matrix):
    value = _finite_array(matrix, shape=(6, 6), name="first-order map")
    poisson = np.zeros((6, 6), dtype=float)
    for coordinate in (0, 2, 4):
        poisson[coordinate, coordinate + 1] = 1.0
        poisson[coordinate + 1, coordinate] = -1.0
    return float(np.max(np.abs(value.T @ poisson @ value - poisson), initial=0.0))


_DECOUPLED_OFF_BLOCK_ENTRIES = tuple(
    (row, column)
    for row in range(6)
    for column in range(6)
    if row // 2 != column // 2
)


def _default_normal_quadrupole_basis(segment_count):
    count = int(segment_count)
    basis = np.zeros((len(SECOND_ORDER_MULTIPOLE_COMPONENTS) * count, count))
    basis[count : 2 * count] = np.eye(count)
    return basis


def _decoupled_first_order_differential(transfer, field_basis):
    matrix = _finite_array(transfer.R, shape=(6, 6), name="first-order map")
    raw_jacobian = _finite_array(
        transfer.R_jacobian, name="first-order map Jacobian"
    )
    basis = _finite_array(field_basis, name="field_basis")
    if (
        raw_jacobian.ndim != 3
        or raw_jacobian.shape[1:] != (6, 6)
        or basis.ndim != 2
        or basis.shape[0] != raw_jacobian.shape[0]
    ):
        raise ValueError("first-order Jacobian and field_basis shapes do not match")

    parameters = []
    directional_rows = []
    for rows, entries in (
        ((0, 1), ((0, 0), (0, 1), (1, 0), (1, 1))),
        ((2, 3), ((2, 2), (2, 3), (3, 2), (3, 3))),
    ):
        block_parameters, entry_jacobian = _kan_parameters_and_entry_jacobian(
            matrix[np.ix_(rows, rows)]
        )
        entry_field_jacobian = np.asarray(
            [raw_jacobian[:, row, column] for row, column in entries]
        )
        parameters.extend(block_parameters)
        directional_rows.append(entry_jacobian @ entry_field_jacobian @ basis)

    parameters.append(float(matrix[4, 5]))
    directional_rows.append((raw_jacobian[:, 4, 5] @ basis)[None, :])
    off_block_values = np.asarray(
        [matrix[index] for index in _DECOUPLED_OFF_BLOCK_ENTRIES], dtype=float
    )
    off_block_jacobian = np.asarray(
        [raw_jacobian[:, index[0], index[1]] @ basis
         for index in _DECOUPLED_OFF_BLOCK_ENTRIES],
        dtype=float,
    )
    return (
        np.asarray(parameters, dtype=float),
        np.vstack(directional_rows),
        off_block_values,
        off_block_jacobian,
        _canonical_linear_symplectic_residual(matrix),
    )


@dataclass(frozen=True)
class DecoupledFirstOrderReachabilityCertificate:
    """Local AD rank certificate in ``KAN_x + KAN_y + R56`` coordinates."""

    multipole_response: np.ndarray
    transfer_matrix: np.ndarray
    parameters: np.ndarray
    field_basis: np.ndarray
    directional_parameter_jacobian: np.ndarray
    transverse_singular_values: np.ndarray
    response_singular_values: np.ndarray
    transverse_numerical_rank: int
    response_numerical_rank: int
    full_transverse_rank: bool
    decoupling_residual: float
    symplectic_residual: float
    relative_tolerance: float


def certify_decoupled_first_order_reachability(
    multipole_response,
    segment_lengths,
    magnetic_rigidity,
    *,
    field_basis=None,
    relative_tolerance=1.0e-9,
    decoupling_tolerance=1.0e-10,
    maximum_step_m=1.0e-3,
    maximum_steps=1_000_000,
) -> DecoupledFirstOrderReachabilityCertificate:
    """Certify the six transverse ``Sp(2,R) x Sp(2,R)`` directions locally.

    The default field basis varies one normal-quadrupole coefficient per
    segment and leaves dipole, skew, and nonlinear multipoles unchanged.
    ``full_transverse_rank`` is true only when all six independent transverse
    KAN directions survive TSVD at the requested relative tolerance.
    """
    lengths = _finite_array(segment_lengths, name="segment_lengths").reshape(-1)
    if lengths.size == 0 or np.any(lengths <= 0.0):
        raise ValueError("segment_lengths must be non-empty and positive")
    response = _finite_array(
        multipole_response, name="multipole_response"
    ).reshape(-1)
    expected = len(SECOND_ORDER_MULTIPOLE_COMPONENTS) * lengths.size
    if response.shape != (expected,):
        raise ValueError("multipole_response must contain five blocks per segment")
    basis = (
        _default_normal_quadrupole_basis(lengths.size)
        if field_basis is None
        else _finite_array(field_basis, name="field_basis")
    )
    if basis.ndim != 2 or basis.shape[0] != response.size:
        raise ValueError("field_basis must have one row per multipole response")
    tolerance = float(relative_tolerance)
    decoupling_tolerance = float(decoupling_tolerance)
    if (
        not np.isfinite(tolerance)
        or tolerance < 0.0
        or not np.isfinite(decoupling_tolerance)
        or decoupling_tolerance < 0.0
    ):
        raise ValueError("rank and decoupling tolerances must be nonnegative")

    transfer = second_order_taylor_map_from_multipoles(
        response,
        lengths,
        magnetic_rigidity,
        maximum_step_m=maximum_step_m,
        maximum_steps=maximum_steps,
    )
    parameters, jacobian, off_block, _, symplectic = (
        _decoupled_first_order_differential(transfer, basis)
    )
    decoupling = float(np.max(np.abs(off_block), initial=0.0))
    scale = max(1.0, float(np.max(np.abs(transfer.R))))
    if decoupling > decoupling_tolerance * scale:
        raise ValueError(
            "current first-order map is not x/y/longitudinal decoupled; "
            f"residual={decoupling:.3e}"
        )
    transverse_singular = np.linalg.svd(jacobian[:6], compute_uv=False)
    response_singular = np.linalg.svd(jacobian, compute_uv=False)
    transverse_threshold = (
        tolerance * transverse_singular[0] if transverse_singular.size else 0.0
    )
    response_threshold = (
        tolerance * response_singular[0] if response_singular.size else 0.0
    )
    transverse_rank = int(np.count_nonzero(
        transverse_singular > transverse_threshold
    ))
    response_rank = int(np.count_nonzero(response_singular > response_threshold))
    return DecoupledFirstOrderReachabilityCertificate(
        multipole_response=response.copy(),
        transfer_matrix=np.asarray(transfer.R, dtype=float).copy(),
        parameters=parameters.copy(),
        field_basis=np.asarray(basis, dtype=float).copy(),
        directional_parameter_jacobian=jacobian.copy(),
        transverse_singular_values=np.asarray(
            transverse_singular, dtype=float
        ),
        response_singular_values=np.asarray(response_singular, dtype=float),
        transverse_numerical_rank=transverse_rank,
        response_numerical_rank=response_rank,
        full_transverse_rank=bool(transverse_rank == 6),
        decoupling_residual=decoupling,
        symplectic_residual=symplectic,
        relative_tolerance=tolerance,
    )


@dataclass(frozen=True)
class DecoupledFirstOrderContinuationStage:
    """One accepted or rejected rank-monitored homotopy attempt."""

    attempt: int
    start_progress: float
    requested_progress: float
    accepted: bool
    maximum_normalized_residual: float
    decoupling_residual: float
    transverse_numerical_rank: int
    response_numerical_rank: int
    transverse_singular_values: np.ndarray
    achieved_parameters: np.ndarray
    function_evaluations: int
    optimizer_status: str


@dataclass(frozen=True)
class DecoupledFirstOrderContinuationResult:
    """Auditable rank-6 continuation result for a decoupled first-order map."""

    target: DecoupledFirstOrderTarget
    initial_multipole_response: np.ndarray
    achieved_multipole_response: np.ndarray
    field_basis: np.ndarray
    basis_coefficients: np.ndarray
    initial_transfer_matrix: np.ndarray
    achieved_transfer_matrix: np.ndarray
    initial_parameters: np.ndarray
    target_parameters: np.ndarray
    achieved_parameters: np.ndarray
    stages: tuple[DecoupledFirstOrderContinuationStage, ...]
    converged: bool
    status: str
    final_maximum_normalized_parameter_error: float
    final_matrix_error: float
    final_decoupling_residual: float
    final_symplectic_residual: float
    final_transverse_rank: int


def solve_decoupled_first_order_continuation(
    initial_multipole_response,
    segment_lengths,
    magnetic_rigidity,
    target,
    *,
    field_basis=None,
    parameter_band=1.0,
    relative_rank_tolerance=1.0e-9,
    decoupling_tolerance=1.0e-9,
    maximum_target_step=5.0e-4,
    minimum_target_step=1.0e-6,
    stage_tolerance=5.0e-8,
    final_tolerance=1.0e-7,
    maximum_stage_evaluations=80,
    maximum_stages=512,
    maximum_step_m=1.0e-3,
    maximum_steps=1_000_000,
) -> DecoupledFirstOrderContinuationResult:
    """Follow a symplectic KAN homotopy while retaining transverse rank six.

    This is deliberately a *local-to-global attempt*, not an unconditional
    controllability claim.  The starting lattice must already have rank six;
    a drift-only rank-three seed is returned as a loud non-converged result.
    Every accepted stage satisfies the nonlinear KAN/off-block residual and
    retains all six transverse TSVD modes.  A target outside the declared
    field basis, including an unsupported ``R56`` change, terminates with an
    auditable non-converged result.
    """
    from scipy.optimize import least_squares

    lengths = _finite_array(segment_lengths, name="segment_lengths").reshape(-1)
    if lengths.size == 0 or np.any(lengths <= 0.0):
        raise ValueError("segment_lengths must be non-empty and positive")
    initial = _finite_array(
        initial_multipole_response, name="initial_multipole_response"
    ).reshape(-1)
    expected = len(SECOND_ORDER_MULTIPOLE_COMPONENTS) * lengths.size
    if initial.shape != (expected,):
        raise ValueError("initial multipoles must contain five blocks per segment")
    target_value = (
        target
        if isinstance(target, DecoupledFirstOrderTarget)
        else DecoupledFirstOrderTarget.from_matrix(target)
    )
    basis = (
        _default_normal_quadrupole_basis(lengths.size)
        if field_basis is None
        else _finite_array(field_basis, name="field_basis")
    )
    if basis.ndim != 2 or basis.shape[0] != initial.size or basis.shape[1] == 0:
        raise ValueError("field_basis must have shape (5*n_segment,n_control)")
    bands = np.broadcast_to(
        np.asarray(parameter_band, dtype=float), (7,)
    ).copy()
    rank_tolerance = float(relative_rank_tolerance)
    decoupling_tolerance = float(decoupling_tolerance)
    maximum_target_step = float(maximum_target_step)
    minimum_target_step = float(minimum_target_step)
    stage_tolerance = float(stage_tolerance)
    final_tolerance = float(final_tolerance)
    maximum_stage_evaluations = int(maximum_stage_evaluations)
    maximum_stages = int(maximum_stages)
    if (
        not np.all(np.isfinite(bands))
        or np.any(bands <= 0.0)
        or not np.isfinite(rank_tolerance)
        or rank_tolerance < 0.0
        or not np.isfinite(decoupling_tolerance)
        or decoupling_tolerance <= 0.0
        or not 0.0 < minimum_target_step <= maximum_target_step
        or not np.isfinite(maximum_target_step)
        or not np.isfinite(stage_tolerance)
        or stage_tolerance <= 0.0
        or not np.isfinite(final_tolerance)
        or final_tolerance <= 0.0
        or maximum_stage_evaluations < 1
        or maximum_stages < 1
    ):
        raise ValueError("continuation tolerances, steps, and bands are invalid")

    cache = {}

    def evaluate(coefficients):
        coefficients = np.asarray(coefficients, dtype=float)
        key = coefficients.tobytes()
        if cache.get("key") == key:
            return cache["value"]
        response = initial + basis @ coefficients
        transfer = second_order_taylor_map_from_multipoles(
            response,
            lengths,
            magnetic_rigidity,
            maximum_step_m=maximum_step_m,
            maximum_steps=maximum_steps,
        )
        differential = _decoupled_first_order_differential(transfer, basis)
        value = (response, transfer) + differential
        cache["key"] = key
        cache["value"] = value
        return value

    coefficients = np.zeros(basis.shape[1], dtype=float)
    initial_data = evaluate(coefficients)
    initial_parameters = initial_data[2].copy()
    initial_off_block = initial_data[4]
    initial_scale = max(1.0, float(np.max(np.abs(initial_data[1].R))))
    if float(np.max(np.abs(initial_off_block), initial=0.0)) > (
        decoupling_tolerance * initial_scale
    ):
        raise ValueError("continuation requires an initially decoupled map")
    # This also validates the magnetostatic longitudinal block and both
    # transverse determinants at the starting point.
    DecoupledFirstOrderTarget.from_matrix(
        initial_data[1].R, tolerance=max(decoupling_tolerance, 1.0e-10)
    )
    target_parameters = target_value.parameters.copy()
    for rotation_index in (0, 3):
        difference = target_parameters[rotation_index] - initial_parameters[
            rotation_index
        ]
        target_parameters[rotation_index] = initial_parameters[rotation_index] + (
            np.arctan2(np.sin(difference), np.cos(difference))
        )
    stages = []

    def ranks(parameter_jacobian):
        transverse = np.linalg.svd(
            parameter_jacobian[:6] / bands[:6, None], compute_uv=False
        )
        response = np.linalg.svd(
            parameter_jacobian / bands[:, None], compute_uv=False
        )
        transverse_threshold = (
            rank_tolerance * transverse[0] if transverse.size else 0.0
        )
        response_threshold = rank_tolerance * response[0] if response.size else 0.0
        return (
            int(np.count_nonzero(transverse > transverse_threshold)),
            int(np.count_nonzero(response > response_threshold)),
            np.asarray(transverse, dtype=float),
        )

    def parameter_difference(parameters, wanted):
        difference = np.asarray(parameters - wanted, dtype=float)
        for rotation_index in (0, 3):
            difference[rotation_index] = np.arctan2(
                np.sin(difference[rotation_index]),
                np.cos(difference[rotation_index]),
            )
        return difference

    def finish(status, requested_converged):
        data = evaluate(coefficients)
        response, transfer, parameters, jacobian, off_block, _, symplectic = data
        normalized_error = parameter_difference(
            parameters, target_parameters
        ) / bands
        parameter_error = float(np.max(np.abs(normalized_error), initial=0.0))
        matrix_error = float(
            np.max(np.abs(transfer.R - target_value.matrix), initial=0.0)
        )
        decoupling = float(np.max(np.abs(off_block), initial=0.0))
        transverse_rank, _, _ = ranks(jacobian)
        converged = bool(
            requested_converged
            and parameter_error <= final_tolerance
            and decoupling <= decoupling_tolerance
            and transverse_rank == 6
        )
        if requested_converged and not converged:
            status = "continuation reached the endpoint but failed the final gate"
        return DecoupledFirstOrderContinuationResult(
            target=target_value,
            initial_multipole_response=initial.copy(),
            achieved_multipole_response=np.asarray(response, dtype=float).copy(),
            field_basis=np.asarray(basis, dtype=float).copy(),
            basis_coefficients=np.asarray(coefficients, dtype=float).copy(),
            initial_transfer_matrix=np.asarray(initial_data[1].R, dtype=float).copy(),
            achieved_transfer_matrix=np.asarray(transfer.R, dtype=float).copy(),
            initial_parameters=initial_parameters.copy(),
            target_parameters=target_parameters.copy(),
            achieved_parameters=np.asarray(parameters, dtype=float).copy(),
            stages=tuple(stages),
            converged=converged,
            status=status,
            final_maximum_normalized_parameter_error=parameter_error,
            final_matrix_error=matrix_error,
            final_decoupling_residual=decoupling,
            final_symplectic_residual=float(symplectic),
            final_transverse_rank=transverse_rank,
        )

    initial_transverse_rank, _, _ = ranks(initial_data[3])
    if initial_transverse_rank < 6:
        return finish(
            "initial lattice has transverse rank "
            f"{initial_transverse_rank}; provide a non-singular rank-6 seed",
            False,
        )

    total_change = float(np.max(np.abs(
        (target_parameters - initial_parameters) / bands
    ), initial=0.0))
    if total_change <= final_tolerance:
        return finish("target already reached by the rank-6 seed", True)
    progress = 0.0
    progress_step = min(1.0, maximum_target_step / total_change)
    minimum_progress_step = min(1.0, minimum_target_step / total_change)

    for attempt in range(1, maximum_stages + 1):
        if progress >= 1.0 - 8.0 * np.finfo(float).eps:
            return finish("rank-6 KAN continuation converged", True)
        requested_progress = min(1.0, progress + progress_step)
        desired = initial_parameters + requested_progress * (
            target_parameters - initial_parameters
        )

        def residual(trial_coefficients, wanted=desired):
            data = evaluate(trial_coefficients)
            parameter_residual = parameter_difference(data[2], wanted) / bands
            off_block_residual = data[4] / decoupling_tolerance
            return np.r_[parameter_residual, off_block_residual]

        def residual_jacobian(trial_coefficients):
            data = evaluate(trial_coefficients)
            return np.vstack(
                (data[3] / bands[:, None], data[5] / decoupling_tolerance)
            )

        trial = least_squares(
            residual,
            coefficients,
            jac=residual_jacobian,
            method="trf",
            x_scale="jac",
            ftol=1.0e-12,
            xtol=1.0e-12,
            gtol=1.0e-12,
            max_nfev=maximum_stage_evaluations,
        )
        trial_data = evaluate(trial.x)
        maximum_residual = float(
            np.max(np.abs(parameter_difference(
                trial_data[2], desired
            ) / bands), initial=0.0)
        )
        trial_decoupling = float(
            np.max(np.abs(trial_data[4]), initial=0.0)
        )
        transverse_rank, response_rank, transverse_singular = ranks(trial_data[3])
        accepted = bool(
            maximum_residual <= stage_tolerance
            and trial_decoupling <= decoupling_tolerance
            and transverse_rank == 6
        )
        stages.append(
            DecoupledFirstOrderContinuationStage(
                attempt=attempt,
                start_progress=progress,
                requested_progress=requested_progress,
                accepted=accepted,
                maximum_normalized_residual=maximum_residual,
                decoupling_residual=trial_decoupling,
                transverse_numerical_rank=transverse_rank,
                response_numerical_rank=response_rank,
                transverse_singular_values=transverse_singular.copy(),
                achieved_parameters=np.asarray(trial_data[2], dtype=float).copy(),
                function_evaluations=int(trial.nfev),
                optimizer_status=str(trial.message),
            )
        )
        if accepted:
            coefficients = np.asarray(trial.x, dtype=float).copy()
            progress = requested_progress
            progress_step = min(1.0 - progress, 1.5 * progress_step)
        else:
            progress_step *= 0.5
            if progress_step < minimum_progress_step:
                return finish(
                    "rank-6 continuation could not accept the minimum target step",
                    False,
                )

    return finish("rank-6 continuation exceeded maximum_stages", False)


@dataclass(frozen=True)
class SecondOrderTaylorMaterialInversePipelineResult:
    """Auditable multipole -> ``R/T`` -> material screening chain."""

    objective: PlanarSecondOrderTaylorMapObjective
    multipole_distribution: np.ndarray
    realized_R: np.ndarray
    realized_T: np.ndarray
    R_difference: np.ndarray
    T_difference: np.ndarray
    normalized_R_difference: np.ndarray
    normalized_T_difference: np.ndarray
    automatic_differentiation: TransferMatrixAutomaticDifferentiation
    multipole_correction: TransferMatrixFieldCorrection
    material_selection: TSVDElementCandidateSelection
    proposed_R: np.ndarray
    proposed_T: np.ndarray
    proposed_exact_max_band_ratio: float
    stage_order: tuple[str, ...] = (
        "normal-skew-multipole-distribution",
        "forward-ad-second-order-taylor-map",
        "target-R-T-difference",
        "tsvd-minimax-multipole-correction",
        "aca-thin-qr-tsvd-material-inverse",
        "native-exact-R-T-gate",
    )


def differentiate_second_order_taylor_field_response(
    objective: PlanarSecondOrderTaylorMapObjective,
    current_multipole_response,
    *,
    field_basis=None,
) -> TransferMatrixAutomaticDifferentiation:
    """Return the full/directional AD Jacobian of selected ``R/T`` rows."""
    if not isinstance(objective, PlanarSecondOrderTaylorMapObjective):
        raise TypeError("objective must be PlanarSecondOrderTaylorMapObjective")
    current = _finite_array(
        current_multipole_response, name="current_multipole_response"
    ).reshape(-1)
    if current.shape != (objective.raw_field_response_size,):
        raise ValueError("current multipoles do not match the objective")
    basis = (
        np.eye(current.size)
        if field_basis is None
        else _finite_array(field_basis, name="field_basis")
    )
    if basis.ndim != 2 or basis.shape[0] != current.size:
        raise ValueError("field_basis must have one row per raw multipole response")
    design = objective.transform(current)
    jacobian = objective.transform_jacobian(current)
    return TransferMatrixAutomaticDifferentiation(
        backend="forward-mode-rk4-taylor-ad",
        current_field_response=current.copy(),
        field_basis=np.asarray(basis, dtype=float).copy(),
        design_response=np.asarray(design, dtype=float),
        full_jacobian=np.asarray(jacobian, dtype=float),
        directional_jacobian=np.asarray(jacobian @ basis, dtype=float),
    )


def differentiate_third_order_taylor_field_response(
    objective: PlanarThirdOrderTaylorMapObjective,
    current_multipole_response,
    *,
    field_basis=None,
) -> TransferMatrixAutomaticDifferentiation:
    """Return the full/directional AD Jacobian of selected ``R/T/U`` rows."""
    if not isinstance(objective, PlanarThirdOrderTaylorMapObjective):
        raise TypeError("objective must be PlanarThirdOrderTaylorMapObjective")
    return differentiate_transfer_matrix_field_response(
        objective,
        current_multipole_response,
        field_basis=field_basis,
    )


def run_second_order_taylor_material_inverse_pipeline(
    objective: PlanarSecondOrderTaylorMapObjective,
    current_multipole_response,
    *,
    candidate_elements,
    candidate_multipole_response_delta,
    candidate_volumes,
    volume_budget,
    field_basis=None,
    field_inverse_relative_tolerance=1.0e-3,
    field_inverse_maximum_step_scale=1.0,
    field_inverse_line_search_steps=8,
    material_relative_tolerance=1.0e-3,
    material_improvement_capture=0.9,
    ratio_tolerance=1.0e-12,
    active_elements=None,
    predecessor_elements=None,
    candidate_volume_changes=None,
    candidate_material_active=None,
    candidate_exclusion_groups=None,
    maximum_changed_volume=None,
    maximum_changed_elements=None,
    candidate_secondary_cost=None,
) -> SecondOrderTaylorMaterialInversePipelineResult:
    """Run the explicit second-order map-to-binary-material inverse chain."""
    current = _finite_array(
        current_multipole_response, name="current_multipole_response"
    ).reshape(-1)
    automatic = differentiate_second_order_taylor_field_response(
        objective, current, field_basis=field_basis
    )
    correction = solve_transfer_matrix_field_correction(
        objective,
        current,
        field_basis=automatic.field_basis,
        relative_tolerance=field_inverse_relative_tolerance,
        maximum_step_scale=field_inverse_maximum_step_scale,
        line_search_steps=field_inverse_line_search_steps,
    )
    realized = objective.evaluate_taylor_map(current)
    selection = select_tsvd_element_candidates(
        current_response=current,
        response_target=correction.target_field_response,
        response_band=correction.field_response_band,
        candidate_elements=candidate_elements,
        candidate_response_delta=candidate_multipole_response_delta,
        candidate_volumes=candidate_volumes,
        volume_budget=volume_budget,
        active_elements=active_elements,
        predecessor_elements=predecessor_elements,
        relative_tolerance=material_relative_tolerance,
        improvement_capture=material_improvement_capture,
        ratio_tolerance=ratio_tolerance,
        candidate_volume_changes=candidate_volume_changes,
        candidate_material_active=candidate_material_active,
        candidate_exclusion_groups=candidate_exclusion_groups,
        maximum_changed_volume=maximum_changed_volume,
        maximum_changed_elements=maximum_changed_elements,
        candidate_secondary_cost=candidate_secondary_cost,
    )
    proposed = objective.evaluate_taylor_map(selection.predicted_response)
    proposed_design = objective.transform(selection.predicted_response)
    proposed_ratio = float(
        np.max(
            np.abs(
                (proposed_design - objective.response_target) / objective.response_band
            )
        )
    )
    return SecondOrderTaylorMaterialInversePipelineResult(
        objective=objective,
        multipole_distribution=current.copy(),
        realized_R=realized.R.copy(),
        realized_T=realized.T.copy(),
        R_difference=objective.target_R - realized.R,
        T_difference=objective.target_T - realized.T,
        normalized_R_difference=(objective.target_R - realized.R) / objective.R_band,
        normalized_T_difference=(objective.target_T - realized.T) / objective.T_band,
        automatic_differentiation=automatic,
        multipole_correction=correction,
        material_selection=selection,
        proposed_R=proposed.R.copy(),
        proposed_T=proposed.T.copy(),
        proposed_exact_max_band_ratio=proposed_ratio,
    )


@dataclass(frozen=True)
class ThirdOrderTaylorMaterialInversePipelineResult:
    """Auditable multipole -> ``R/T/U`` -> material screening chain."""

    objective: PlanarThirdOrderTaylorMapObjective
    multipole_distribution: np.ndarray
    realized_R: np.ndarray
    realized_T: np.ndarray
    realized_U: np.ndarray
    R_difference: np.ndarray
    T_difference: np.ndarray
    U_difference: np.ndarray
    normalized_R_difference: np.ndarray
    normalized_T_difference: np.ndarray
    normalized_U_difference: np.ndarray
    reachability: TaylorMapReachabilityCertificate
    automatic_differentiation: TransferMatrixAutomaticDifferentiation
    multipole_correction: TransferMatrixFieldCorrection
    material_selection: TSVDElementCandidateSelection
    proposed_R: np.ndarray
    proposed_T: np.ndarray
    proposed_U: np.ndarray
    proposed_exact_max_band_ratio: float
    stage_order: tuple[str, ...] = (
        "normal-skew-multipole-distribution-through-octupole",
        "forward-ad-third-order-taylor-map",
        "target-R-T-U-difference",
        "tsvd-reachability-certificate",
        "tsvd-minimax-multipole-correction",
        "aca-thin-qr-tsvd-material-inverse",
        "native-exact-R-T-U-gate",
    )


def run_third_order_taylor_material_inverse_pipeline(
    objective: PlanarThirdOrderTaylorMapObjective,
    current_multipole_response,
    *,
    candidate_elements,
    candidate_multipole_response_delta,
    candidate_volumes,
    volume_budget,
    field_basis=None,
    reachability_relative_tolerance=1.0e-10,
    reachability_acceptance_band_ratio=1.0,
    field_inverse_relative_tolerance=1.0e-3,
    field_inverse_maximum_step_scale=1.0,
    field_inverse_line_search_steps=8,
    material_relative_tolerance=1.0e-3,
    material_improvement_capture=0.9,
    ratio_tolerance=1.0e-12,
    active_elements=None,
    predecessor_elements=None,
    candidate_volume_changes=None,
    candidate_material_active=None,
    candidate_exclusion_groups=None,
    maximum_changed_volume=None,
    maximum_changed_elements=None,
    candidate_secondary_cost=None,
) -> ThirdOrderTaylorMaterialInversePipelineResult:
    """Run the explicit cubic map-to-binary-material inverse chain."""
    if not isinstance(objective, PlanarThirdOrderTaylorMapObjective):
        raise TypeError("objective must be PlanarThirdOrderTaylorMapObjective")
    current = _finite_array(
        current_multipole_response, name="current_multipole_response"
    ).reshape(-1)
    automatic = differentiate_third_order_taylor_field_response(
        objective, current, field_basis=field_basis
    )
    reachability = certify_taylor_map_reachability(
        objective,
        current,
        field_basis=automatic.field_basis,
        relative_tolerance=reachability_relative_tolerance,
        acceptance_band_ratio=reachability_acceptance_band_ratio,
    )
    correction = solve_transfer_matrix_field_correction(
        objective,
        current,
        field_basis=automatic.field_basis,
        relative_tolerance=field_inverse_relative_tolerance,
        maximum_step_scale=field_inverse_maximum_step_scale,
        line_search_steps=field_inverse_line_search_steps,
    )
    realized = objective.evaluate_taylor_map(current)
    selection = select_tsvd_element_candidates(
        current_response=current,
        response_target=correction.target_field_response,
        response_band=correction.field_response_band,
        candidate_elements=candidate_elements,
        candidate_response_delta=candidate_multipole_response_delta,
        candidate_volumes=candidate_volumes,
        volume_budget=volume_budget,
        active_elements=active_elements,
        predecessor_elements=predecessor_elements,
        relative_tolerance=material_relative_tolerance,
        improvement_capture=material_improvement_capture,
        ratio_tolerance=ratio_tolerance,
        candidate_volume_changes=candidate_volume_changes,
        candidate_material_active=candidate_material_active,
        candidate_exclusion_groups=candidate_exclusion_groups,
        maximum_changed_volume=maximum_changed_volume,
        maximum_changed_elements=maximum_changed_elements,
        candidate_secondary_cost=candidate_secondary_cost,
    )
    proposed = objective.evaluate_taylor_map(selection.predicted_response)
    proposed_design = objective.transform(selection.predicted_response)
    proposed_ratio = float(
        np.max(
            np.abs(
                (proposed_design - objective.response_target) / objective.response_band
            )
        )
    )
    return ThirdOrderTaylorMaterialInversePipelineResult(
        objective=objective,
        multipole_distribution=current.copy(),
        realized_R=realized.R.copy(),
        realized_T=realized.T.copy(),
        realized_U=realized.U.copy(),
        R_difference=objective.target_R - realized.R,
        T_difference=objective.target_T - realized.T,
        U_difference=objective.target_U - realized.U,
        normalized_R_difference=(objective.target_R - realized.R) / objective.R_band,
        normalized_T_difference=(objective.target_T - realized.T) / objective.T_band,
        normalized_U_difference=(objective.target_U - realized.U) / objective.U_band,
        reachability=reachability,
        automatic_differentiation=automatic,
        multipole_correction=correction,
        material_selection=selection,
        proposed_R=proposed.R.copy(),
        proposed_T=proposed.T.copy(),
        proposed_U=proposed.U.copy(),
        proposed_exact_max_band_ratio=proposed_ratio,
    )


@dataclass(frozen=True)
class SecondOrderTaylorTopologyResult:
    """Whole-element HDiv-MMM result scored by an exact native ``R/T`` map."""

    objective: PlanarSecondOrderTaylorMapObjective
    generation: HDivMMMGenerationResult
    realized_multipole_response: np.ndarray
    realized_R: np.ndarray
    realized_T: np.ndarray
    normal_dipole_max_band_ratio: float
    taylor_map_max_band_ratio: float
    topology: GrowthTopologyReport

    @property
    def active_elements(self) -> np.ndarray:
        return self.generation.active_elements

    @property
    def converged(self) -> bool:
        return bool(
            self.normal_dipole_max_band_ratio <= 1.0
            and self.taylor_map_max_band_ratio <= 1.0
        )


def optimize_hdiv_mmm_magnet_from_second_order_taylor_map(
    objective: PlanarSecondOrderTaylorMapObjective,
    *,
    charge_gram,
    fes,
    inv_chi,
    rhs,
    multipole_response_matrix,
    active_elements,
    element_volumes,
    volume_max,
    incident_multipole_response=None,
    **generation_options,
) -> SecondOrderTaylorTopologyResult:
    """Optimize binary material directly against selected native ``R/T`` rows.

    Candidate screening contracts the forward-AD ``R/T`` Jacobian with the
    raw HDiv multipole rows before ACA--thin-QR--TSVD.  Every committed binary
    state is fully re-solved, transformed by the native C++ variational map,
    and scored in the original engineering bands.
    """
    if not isinstance(objective, PlanarSecondOrderTaylorMapObjective):
        raise TypeError("objective must be PlanarSecondOrderTaylorMapObjective")
    response_matrix = _finite_array(
        multipole_response_matrix, name="multipole_response_matrix"
    )
    expected = (objective.raw_field_response_size, int(fes.ndof))
    if response_matrix.ndim != 2 or response_matrix.shape != expected:
        raise ValueError(
            "multipole_response_matrix must have shape (5*n_orbit_segment,fes.ndof)"
        )
    incident = (
        np.zeros(objective.raw_field_response_size)
        if incident_multipole_response is None
        else _finite_array(
            incident_multipole_response, name="incident_multipole_response"
        ).reshape(-1)
    )
    if incident.shape != (objective.raw_field_response_size,):
        raise ValueError("incident multipoles must match the raw response size")
    reserved = {
        "response_matrix",
        "response_target",
        "response_band",
        "response_transform",
        "response_transform_jacobian",
        "incident_response",
    }
    overlap = reserved.intersection(generation_options)
    if overlap:
        raise TypeError(
            "generation_options cannot override the Taylor-map contract: "
            + ", ".join(sorted(overlap))
        )
    generation = grow_hdiv_mmm_by_superposition(
        charge_gram=charge_gram,
        fes=fes,
        inv_chi=inv_chi,
        rhs=rhs,
        response_matrix=response_matrix,
        active_elements=active_elements,
        element_volumes=element_volumes,
        response_target=objective.response_target,
        response_band=objective.response_band,
        volume_max=volume_max,
        incident_response=incident,
        response_transform=objective.transform,
        response_transform_jacobian=objective.transform_jacobian,
        **generation_options,
    )
    raw = np.asarray(generation.response, dtype=float)
    transfer = objective.evaluate_taylor_map(raw)
    count = len(objective.orbit.segment_lengths)
    dipole_ratio = float(
        np.max(
            np.abs(
                (raw[:count] - objective.required_normal_dipole)
                / objective.normal_dipole_band
            )
        )
    )
    R_ratio = max(
        abs((transfer.R[index] - objective.target_R[index]) / objective.R_band[index])
        for index in objective.R_entries
    )
    T_ratio = max(
        abs((transfer.T[index] - objective.target_T[index]) / objective.T_band[index])
        for index in objective.T_entries
    )
    topology = ngsolve_growth_topology(fes.mesh, generation.active_elements)
    return SecondOrderTaylorTopologyResult(
        objective=objective,
        generation=generation,
        realized_multipole_response=raw.copy(),
        realized_R=transfer.R.copy(),
        realized_T=transfer.T.copy(),
        normal_dipole_max_band_ratio=dipole_ratio,
        taylor_map_max_band_ratio=float(max(R_ratio, T_ratio)),
        topology=topology,
    )


@dataclass(frozen=True)
class ThirdOrderTaylorTopologyResult:
    """Whole-element HDiv-MMM result scored by an exact native ``R/T/U`` map."""

    objective: PlanarThirdOrderTaylorMapObjective
    generation: HDivMMMGenerationResult
    realized_multipole_response: np.ndarray
    realized_R: np.ndarray
    realized_T: np.ndarray
    realized_U: np.ndarray
    normal_dipole_max_band_ratio: float
    taylor_map_max_band_ratio: float
    topology: GrowthTopologyReport

    @property
    def active_elements(self) -> np.ndarray:
        return self.generation.active_elements

    @property
    def converged(self) -> bool:
        return bool(
            self.normal_dipole_max_band_ratio <= 1.0
            and self.taylor_map_max_band_ratio <= 1.0
        )


def optimize_hdiv_mmm_magnet_from_third_order_taylor_map(
    objective: PlanarThirdOrderTaylorMapObjective,
    *,
    charge_gram,
    fes,
    inv_chi,
    rhs,
    multipole_response_matrix,
    active_elements,
    element_volumes,
    volume_max,
    incident_multipole_response=None,
    **generation_options,
) -> ThirdOrderTaylorTopologyResult:
    """Optimize binary material directly against selected native ``R/T/U`` rows."""
    if not isinstance(objective, PlanarThirdOrderTaylorMapObjective):
        raise TypeError("objective must be PlanarThirdOrderTaylorMapObjective")
    response_matrix = _finite_array(
        multipole_response_matrix, name="multipole_response_matrix"
    )
    expected = (objective.raw_field_response_size, int(fes.ndof))
    if response_matrix.ndim != 2 or response_matrix.shape != expected:
        raise ValueError(
            "multipole_response_matrix must have shape (7*n_orbit_segment,fes.ndof)"
        )
    incident = (
        np.zeros(objective.raw_field_response_size)
        if incident_multipole_response is None
        else _finite_array(
            incident_multipole_response, name="incident_multipole_response"
        ).reshape(-1)
    )
    if incident.shape != (objective.raw_field_response_size,):
        raise ValueError("incident multipoles must match the raw response size")
    reserved = {
        "response_matrix",
        "response_target",
        "response_band",
        "response_transform",
        "response_transform_jacobian",
        "incident_response",
    }
    overlap = reserved.intersection(generation_options)
    if overlap:
        raise TypeError(
            "generation_options cannot override the Taylor-map contract: "
            + ", ".join(sorted(overlap))
        )
    generation = grow_hdiv_mmm_by_superposition(
        charge_gram=charge_gram,
        fes=fes,
        inv_chi=inv_chi,
        rhs=rhs,
        response_matrix=response_matrix,
        active_elements=active_elements,
        element_volumes=element_volumes,
        response_target=objective.response_target,
        response_band=objective.response_band,
        volume_max=volume_max,
        incident_response=incident,
        response_transform=objective.transform,
        response_transform_jacobian=objective.transform_jacobian,
        **generation_options,
    )
    raw = np.asarray(generation.response, dtype=float)
    transfer = objective.evaluate_taylor_map(raw)
    count = len(objective.orbit.segment_lengths)
    dipole_ratio = float(
        np.max(
            np.abs(
                (raw[:count] - objective.required_normal_dipole)
                / objective.normal_dipole_band
            )
        )
    )
    R_ratio = max(
        abs((transfer.R[index] - objective.target_R[index]) / objective.R_band[index])
        for index in objective.R_entries
    )
    T_ratio = max(
        abs((transfer.T[index] - objective.target_T[index]) / objective.T_band[index])
        for index in objective.T_entries
    )
    U_ratio = max(
        abs((transfer.U[index] - objective.target_U[index]) / objective.U_band[index])
        for index in objective.U_entries
    )
    topology = ngsolve_growth_topology(fes.mesh, generation.active_elements)
    return ThirdOrderTaylorTopologyResult(
        objective=objective,
        generation=generation,
        realized_multipole_response=raw.copy(),
        realized_R=transfer.R.copy(),
        realized_T=transfer.T.copy(),
        realized_U=transfer.U.copy(),
        normal_dipole_max_band_ratio=dipole_ratio,
        taylor_map_max_band_ratio=float(max(R_ratio, T_ratio, U_ratio)),
        topology=topology,
    )


__all__ = [
    "SECOND_ORDER_MULTIPOLE_COMPONENTS",
    "THIRD_ORDER_MULTIPOLE_COMPONENTS",
    "DecoupledFirstOrderContinuationResult",
    "DecoupledFirstOrderContinuationStage",
    "DecoupledFirstOrderReachabilityCertificate",
    "DecoupledFirstOrderTarget",
    "PlanarSecondOrderTaylorMapObjective",
    "PlanarThirdOrderTaylorMapObjective",
    "SecondOrderTaylorMap",
    "SecondOrderTaylorMaterialInversePipelineResult",
    "SecondOrderTaylorTopologyResult",
    "Symplectic2x2KAN",
    "TaylorMapReachabilityCertificate",
    "ThirdOrderTaylorMap",
    "ThirdOrderTaylorMaterialInversePipelineResult",
    "ThirdOrderTaylorTopologyResult",
    "build_planar_orbit_cubic_multipole_response_matrix",
    "build_planar_orbit_multipole_response_matrix",
    "certify_decoupled_first_order_reachability",
    "certify_taylor_map_reachability",
    "differentiate_second_order_taylor_field_response",
    "differentiate_third_order_taylor_field_response",
    "optimize_hdiv_mmm_magnet_from_second_order_taylor_map",
    "optimize_hdiv_mmm_magnet_from_third_order_taylor_map",
    "planar_orbit_cubic_multipole_observations",
    "planar_orbit_multipole_observations",
    "run_second_order_taylor_material_inverse_pipeline",
    "run_third_order_taylor_material_inverse_pipeline",
    "second_order_taylor_map_from_multipoles",
    "solve_decoupled_first_order_continuation",
    "third_order_taylor_map_from_multipoles",
]
