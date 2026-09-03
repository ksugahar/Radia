"""Common observation tube for independent accelerator-magnet field engines.

The module compares physical, gauge-invariant magnetic flux density sampled by
HDiv-MMM, reduced-A, Omega-reduced-Omega, or another independent solver.  It
also samples vector potential when supplied, but deliberately does not compare
raw A values because different gauges are not an error.  A comparisons belong
after projection to the agreed design-orbit gauge or through curl(A)/tracking.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np

ArrayFieldEvaluator = Callable[[np.ndarray], np.ndarray]


def _finite_array(name: str, value, *, ndim: int, last_dim: int | None = None):
    result = np.asarray(value, dtype=float)
    if result.ndim != ndim:
        raise ValueError(f"{name} must have {ndim} dimensions")
    if last_dim is not None and result.shape[-1] != last_dim:
        raise ValueError(f"{name} must end in dimension {last_dim}")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be finite")
    result = np.array(result, copy=True)
    result.setflags(write=False)
    return result


@dataclass(frozen=True)
class CurvilinearObservationTube:
    """One shared design-orbit frame and set of transverse observation points.

    Local component order is ``(x, y, s)``: normal, binormal, tangent.  The
    frame must be right-handed, ``tangent x normal = binormal``.
    """

    station_s: np.ndarray
    center: np.ndarray
    tangent: np.ndarray
    normal: np.ndarray
    binormal: np.ndarray
    transverse_offsets: np.ndarray

    def __post_init__(self):
        station_s = _finite_array("station_s", self.station_s, ndim=1)
        if station_s.size == 0:
            raise ValueError("station_s must contain at least one station")
        if station_s.size > 1 and np.any(np.diff(station_s) <= 0.0):
            raise ValueError("station_s must be strictly increasing")
        count = station_s.size
        arrays = {}
        for name in ("center", "tangent", "normal", "binormal"):
            value = _finite_array(name, getattr(self, name), ndim=2, last_dim=3)
            if value.shape[0] != count:
                raise ValueError(f"{name} station count must match station_s")
            arrays[name] = value
        offsets = _finite_array(
            "transverse_offsets", self.transverse_offsets, ndim=2, last_dim=2)
        if offsets.shape[0] == 0:
            raise ValueError("transverse_offsets must contain at least one point")

        tangent = arrays["tangent"]
        normal = arrays["normal"]
        binormal = arrays["binormal"]
        tolerance = 2.0e-10
        for name, value in (("tangent", tangent), ("normal", normal),
                            ("binormal", binormal)):
            error = np.max(np.abs(np.linalg.norm(value, axis=1) - 1.0))
            if error > tolerance:
                raise ValueError(f"{name} must be unit length; max error={error:.3e}")
        orthogonality = max(
            np.max(np.abs(np.einsum("ij,ij->i", tangent, normal))),
            np.max(np.abs(np.einsum("ij,ij->i", tangent, binormal))),
            np.max(np.abs(np.einsum("ij,ij->i", normal, binormal))),
        )
        if orthogonality > tolerance:
            raise ValueError(
                f"observation frame must be orthogonal; max dot={orthogonality:.3e}")
        handedness = np.max(
            np.linalg.norm(np.cross(tangent, normal) - binormal, axis=1))
        if handedness > tolerance:
            raise ValueError(
                "observation frame must satisfy tangent x normal = binormal; "
                f"max error={handedness:.3e}")

        object.__setattr__(self, "station_s", station_s)
        object.__setattr__(self, "transverse_offsets", offsets)
        for name, value in arrays.items():
            object.__setattr__(self, name, value)

    @property
    def station_count(self) -> int:
        return int(self.station_s.size)

    @property
    def transverse_point_count(self) -> int:
        return int(self.transverse_offsets.shape[0])

    def global_points(self) -> np.ndarray:
        """Return points with shape ``(station, transverse point, xyz)``."""
        x = self.transverse_offsets[:, 0]
        y = self.transverse_offsets[:, 1]
        return (
            self.center[:, None, :]
            + x[None, :, None] * self.normal[:, None, :]
            + y[None, :, None] * self.binormal[:, None, :]
        )

    def vectors_to_local(self, global_vectors) -> np.ndarray:
        """Convert a sampled vector array to local ``(x,y,s)`` components."""
        values = np.asarray(global_vectors, dtype=float)
        expected = (self.station_count, self.transverse_point_count, 3)
        if values.shape != expected:
            raise ValueError(f"global_vectors must have shape {expected}")
        return np.stack((
            np.einsum("spj,sj->sp", values, self.normal),
            np.einsum("spj,sj->sp", values, self.binormal),
            np.einsum("spj,sj->sp", values, self.tangent),
        ), axis=-1)


def circular_transverse_offsets(radius_m: float, sample_count: int = 48) -> np.ndarray:
    """Return one centered transverse circle for station multipole analysis."""
    radius = float(radius_m)
    sample_count = int(sample_count)
    if not np.isfinite(radius) or radius <= 0.0:
        raise ValueError("radius_m must be positive and finite")
    if sample_count < 4:
        raise ValueError("sample_count must be at least four")
    angles = np.linspace(0.0, 2.0 * np.pi, sample_count, endpoint=False)
    return np.column_stack((radius * np.cos(angles), radius * np.sin(angles)))


def rectangular_transverse_offsets(
    x_offsets_m: Sequence[float], y_offsets_m: Sequence[float]
) -> np.ndarray:
    """Return the Cartesian product of physical x and y aperture offsets."""
    x = np.asarray(x_offsets_m, dtype=float).reshape(-1)
    y = np.asarray(y_offsets_m, dtype=float).reshape(-1)
    if x.size == 0 or y.size == 0 or not np.all(np.isfinite(x)) or not np.all(
        np.isfinite(y)
    ):
        raise ValueError("x_offsets_m and y_offsets_m must be finite and non-empty")
    xx, yy = np.meshgrid(x, y, indexing="ij")
    return np.column_stack((xx.ravel(), yy.ravel()))


@dataclass(frozen=True)
class MagnetFieldEngine:
    """Named batch evaluators supplied by one independent field formulation."""

    name: str
    magnetic_flux_density: ArrayFieldEvaluator
    vector_potential: ArrayFieldEvaluator | None = None

    def __post_init__(self):
        if not str(self.name).strip():
            raise ValueError("field-engine name must be non-empty")
        if not callable(self.magnetic_flux_density):
            raise TypeError("magnetic_flux_density must be callable")
        if self.vector_potential is not None and not callable(self.vector_potential):
            raise TypeError("vector_potential must be callable when supplied")


@dataclass(frozen=True)
class MagnetFieldSample:
    """One engine sampled on one immutable observation tube."""

    engine_name: str
    tube: CurvilinearObservationTube
    b_global: np.ndarray
    b_local: np.ndarray
    a_global: np.ndarray | None = None
    a_local: np.ndarray | None = None

    def __post_init__(self):
        if not str(self.engine_name).strip():
            raise ValueError("engine_name must be non-empty")
        if not isinstance(self.tube, CurvilinearObservationTube):
            raise TypeError("tube must be a CurvilinearObservationTube")
        expected = (
            self.tube.station_count,
            self.tube.transverse_point_count,
            3,
        )

        def checked(name, value):
            result = _finite_array(name, value, ndim=3, last_dim=3)
            if result.shape != expected:
                raise ValueError(f"{name} must have shape {expected}")
            return result

        b_global = checked("b_global", self.b_global)
        b_local = checked("b_local", self.b_local)
        if (self.a_global is None) != (self.a_local is None):
            raise ValueError("a_global and a_local must be supplied together")
        a_global = None
        a_local = None
        if self.a_global is not None:
            a_global = checked("a_global", self.a_global)
            a_local = checked("a_local", self.a_local)
        object.__setattr__(self, "b_global", b_global)
        object.__setattr__(self, "b_local", b_local)
        object.__setattr__(self, "a_global", a_global)
        object.__setattr__(self, "a_local", a_local)


def _evaluate_vector_field(
    name: str, evaluator: ArrayFieldEvaluator, points: np.ndarray
) -> np.ndarray:
    flat = np.ascontiguousarray(points.reshape(-1, 3))
    values = np.asarray(evaluator(flat), dtype=float)
    if values.shape != flat.shape:
        raise ValueError(
            f"{name} evaluator must return shape {flat.shape}, got {values.shape}")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} evaluator returned non-finite values")
    return values.reshape(points.shape)


def sample_field_engine(
    engine: MagnetFieldEngine, tube: CurvilinearObservationTube
) -> MagnetFieldSample:
    """Evaluate one engine once on the complete common observation tube."""
    points = tube.global_points()
    b_global = _evaluate_vector_field(
        f"{engine.name}.B", engine.magnetic_flux_density, points)
    b_local = tube.vectors_to_local(b_global)
    a_global = None
    a_local = None
    if engine.vector_potential is not None:
        a_global = _evaluate_vector_field(
            f"{engine.name}.A", engine.vector_potential, points)
        a_local = tube.vectors_to_local(a_global)
    return MagnetFieldSample(
        engine_name=engine.name,
        tube=tube,
        b_global=b_global,
        b_local=b_local,
        a_global=a_global,
        a_local=a_local,
    )


def _same_tube(left: CurvilinearObservationTube, right: CurvilinearObservationTube):
    for name in (
        "station_s", "center", "tangent", "normal", "binormal",
        "transverse_offsets",
    ):
        if not np.array_equal(getattr(left, name), getattr(right, name)):
            return False
    return True


def compare_magnetic_flux_density(
    samples: Sequence[MagnetFieldSample], *, scale_floor_t: float = 1.0e-15
) -> dict:
    """Return gauge-invariant pairwise B residuals on one common tube."""
    values = tuple(samples)
    if len(values) < 2:
        raise ValueError("at least two field samples are required")
    names = [sample.engine_name for sample in values]
    if len(set(names)) != len(names):
        raise ValueError("field-engine names must be unique")
    reference_tube = values[0].tube
    if any(not _same_tube(reference_tube, sample.tube) for sample in values[1:]):
        raise ValueError("all field samples must use the same observation tube")
    scale_floor = float(scale_floor_t)
    if not np.isfinite(scale_floor) or scale_floor <= 0.0:
        raise ValueError("scale_floor_t must be positive and finite")

    pairs = []
    for left_index, left in enumerate(values[:-1]):
        for right in values[left_index + 1:]:
            difference = left.b_local - right.b_local
            point_norm = np.linalg.norm(difference, axis=-1)
            rms_difference = float(np.sqrt(np.mean(point_norm ** 2)))
            symmetric_scale = float(np.sqrt(0.5 * np.mean(
                np.sum(left.b_local ** 2, axis=-1)
                + np.sum(right.b_local ** 2, axis=-1))))
            maximum_scale = float(max(
                np.max(np.linalg.norm(left.b_local, axis=-1)),
                np.max(np.linalg.norm(right.b_local, axis=-1)),
                scale_floor,
            ))
            if reference_tube.station_count > 1:
                left_integral = np.trapezoid(
                    left.b_local, reference_tube.station_s, axis=0)
                right_integral = np.trapezoid(
                    right.b_local, reference_tube.station_s, axis=0)
            else:
                shape = (reference_tube.transverse_point_count, 3)
                left_integral = np.zeros(shape)
                right_integral = np.zeros(shape)
            integral_difference = np.linalg.norm(
                left_integral - right_integral, axis=-1)
            integral_scale = max(
                float(np.max(np.linalg.norm(left_integral, axis=-1))),
                float(np.max(np.linalg.norm(right_integral, axis=-1))),
                scale_floor,
            )
            pairs.append({
                "left": left.engine_name,
                "right": right.engine_name,
                "maximum_vector_error_t": float(np.max(point_norm)),
                "rms_vector_error_t": rms_difference,
                "relative_rms_error": rms_difference / max(
                    symmetric_scale, scale_floor),
                "relative_maximum_error": float(np.max(point_norm)) / maximum_scale,
                "component_maximum_error_t": np.max(
                    np.abs(difference), axis=(0, 1)).tolist(),
                "maximum_integrated_error_t_m": float(
                    np.max(integral_difference)),
                "relative_integrated_error": float(
                    np.max(integral_difference)) / integral_scale,
            })
    return {
        "engine_names": names,
        "station_count": reference_tube.station_count,
        "transverse_point_count": reference_tube.transverse_point_count,
        "raw_vector_potential_compared": False,
        "pairs": pairs,
    }


def longitudinal_reversal_symmetry(
    station_s: np.ndarray,
    profile: np.ndarray,
    *,
    scale_floor: float = 1.0e-30,
) -> dict:
    """Measure the even-in-longitudinal-coordinate defect of a field profile.

    Accelerator magnets with identical entrance and exit geometry must have an
    even main-multipole profile about their longitudinal centre.  This check is
    intentionally performed on a derived physical observable rather than on
    mesh coordinates, so it catches a field evaluator that lost an exact image
    operation or a full-domain discretization that developed a large spurious
    odd component.
    """

    stations = np.asarray(station_s, dtype=float)
    values = np.asarray(profile)
    if stations.ndim != 1 or values.ndim != 1 or stations.size != values.size:
        raise ValueError("station_s and profile must be equal-length 1-D arrays")
    if stations.size < 2:
        raise ValueError("at least two longitudinal stations are required")
    if not np.isfinite(stations).all() or not np.isfinite(values).all():
        raise ValueError("station_s and profile must be finite")
    if not np.allclose(stations, -stations[::-1], rtol=1.0e-12, atol=1.0e-14):
        raise ValueError("station_s must be symmetric about zero")
    floor = float(scale_floor)
    if not np.isfinite(floor) or floor <= 0.0:
        raise ValueError("scale_floor must be positive and finite")

    reversed_values = values[::-1]
    difference = values - reversed_values
    odd = 0.5 * difference
    even = 0.5 * (values + reversed_values)
    rms_scale = max(float(np.sqrt(np.mean(np.abs(values) ** 2))), floor)
    maximum_scale = max(float(np.max(np.abs(values))), floor)
    return {
        "relative_rms_defect": float(np.sqrt(np.mean(np.abs(difference) ** 2)))
        / rms_scale,
        "relative_maximum_defect": float(np.max(np.abs(difference)))
        / maximum_scale,
        "odd_to_even_l2": float(np.linalg.norm(odd))
        / max(float(np.linalg.norm(even)), floor),
        "maximum_absolute_defect": float(np.max(np.abs(difference))),
    }


def project_straight_quadrupole_symmetry(
    points, evaluator: ArrayFieldEvaluator, *, longitudinal_axis: int = 0
) -> np.ndarray:
    """Project a straight normal-quadrupole field onto its exact symmetry.

    The projection averages the four rotations about the longitudinal axis
    with the quadrupole character ``(-1)^k`` and the two longitudinal reversal
    samples.  Transverse B is even through the magnet centre while the
    longitudinal fringe component is odd.  This removes mesh-induced dipole,
    skew, and longitudinal-asymmetry contamination without changing the exact
    quadrupole component.
    """

    query = np.asarray(points, dtype=float)
    if query.ndim != 2 or query.shape[1] != 3 or not np.all(np.isfinite(query)):
        raise ValueError("points must be one finite (n, 3) array")
    axis = int(longitudinal_axis)
    if axis not in (0, 1, 2):
        raise ValueError("longitudinal_axis must be 0, 1, or 2")
    if not callable(evaluator):
        raise TypeError("evaluator must be callable")
    transverse = [index for index in range(3) if index != axis]
    reversal = np.eye(3)
    reversal[axis, axis] = -1.0
    projected = np.zeros_like(query)
    for reverse in (False, True):
        parity = reversal if reverse else np.eye(3)
        for quarter_turn in range(4):
            angle = 0.5 * np.pi * quarter_turn
            cosine, sine = np.cos(angle), np.sin(angle)
            rotation = np.eye(3)
            i, j = transverse
            rotation[i, i] = cosine
            rotation[i, j] = -sine
            rotation[j, i] = sine
            rotation[j, j] = cosine
            transformed_points = query @ rotation.T
            if reverse:
                transformed_points = transformed_points @ reversal
            sampled = np.asarray(evaluator(transformed_points), dtype=float)
            if sampled.shape != query.shape or not np.all(np.isfinite(sampled)):
                raise ValueError("evaluator must return one finite (n, 3) array")
            # Column-vector form: P * chi_k * R^T * B(R P x).
            mapped = ((-1.0) ** quarter_turn) * sampled @ rotation
            if reverse:
                mapped = mapped @ parity
            projected += mapped
    return projected / 8.0


def transverse_multipole_spectrum(
    sample: MagnetFieldSample, maximum_order: int = 6
) -> np.ndarray:
    """Fit ``By+i*Bx = sum C_n (x+i*y)^(n-1)`` at every station.

    The observation offsets must be one centered constant-radius circle.  The
    returned array has shape ``(station_count, maximum_order)``.
    """
    maximum_order = int(maximum_order)
    if maximum_order < 1:
        raise ValueError("maximum_order must be at least one")
    offsets = sample.tube.transverse_offsets
    if offsets.shape[0] < 2 * maximum_order:
        raise ValueError("the transverse circle needs at least twice maximum_order points")
    radii = np.linalg.norm(offsets, axis=1)
    if radii[0] <= 0.0 or not np.allclose(
        radii, radii[0], rtol=1.0e-10, atol=1.0e-14
    ):
        raise ValueError("transverse offsets must form one centered circle")
    coordinate = offsets[:, 0] + 1j * offsets[:, 1]
    transverse = sample.b_local[:, :, 1] + 1j * sample.b_local[:, :, 0]
    basis = np.column_stack([
        coordinate ** power for power in range(maximum_order)
    ])
    coefficients, _, rank, _ = np.linalg.lstsq(
        basis, transverse.T, rcond=None)
    if rank != maximum_order:
        raise ValueError("transverse offsets do not define a full-rank fit")
    return coefficients.T


def compare_integrated_multipole_rows(
    reference_rows: Sequence[dict],
    candidate_rows: Sequence[dict],
    main_order: int,
) -> dict:
    """Compare two integrated-multipole tables on the same reference circle.

    The main harmonic is compared by relative complex-coefficient error.  The
    remaining harmonics are reported as normal/skew accelerator-unit
    differences; no universal pass band is assigned because acceptable higher
    harmonics are magnet- and lattice-specific.
    """

    main_order = int(main_order)
    if main_order < 1:
        raise ValueError("main_order must be at least one")

    def indexed(rows: Sequence[dict], name: str) -> dict[int, dict]:
        result: dict[int, dict] = {}
        for row in rows:
            order = int(row["order"])
            if order in result:
                raise ValueError(f"{name} contains duplicate order {order}")
            result[order] = row
        if not result:
            raise ValueError(f"{name} must not be empty")
        return result

    reference = indexed(reference_rows, "reference_rows")
    candidate = indexed(candidate_rows, "candidate_rows")
    common = sorted(set(reference) & set(candidate))
    if main_order not in common:
        raise ValueError(f"main_order {main_order} is absent from the common rows")

    def coefficient(row: dict) -> complex:
        return complex(
            float(row["integrated_real_t_m_per_m_power"]),
            float(row["integrated_imag_t_m_per_m_power"]),
        )

    reference_main = coefficient(reference[main_order])
    candidate_main = coefficient(candidate[main_order])
    main_relative_error = (
        abs(candidate_main - reference_main) / abs(reference_main)
        if abs(reference_main) > 1.0e-30 else None
    )
    harmonics = []
    for order in common:
        ref = reference[order]
        cand = candidate[order]
        harmonics.append({
            "order": order,
            "reference_normal_units": float(ref["normal_units_at_reference_radius"]),
            "candidate_normal_units": float(cand["normal_units_at_reference_radius"]),
            "normal_units_difference": float(
                cand["normal_units_at_reference_radius"]
                - ref["normal_units_at_reference_radius"]),
            "reference_skew_units": float(ref["skew_units_at_reference_radius"]),
            "candidate_skew_units": float(cand["skew_units_at_reference_radius"]),
            "skew_units_difference": float(
                cand["skew_units_at_reference_radius"]
                - ref["skew_units_at_reference_radius"]),
        })
    return {
        "main_order": main_order,
        "main_reference_coefficient": [reference_main.real, reference_main.imag],
        "main_candidate_coefficient": [candidate_main.real, candidate_main.imag],
        "main_relative_error": main_relative_error,
        "harmonics_at_reference_radius": harmonics,
    }


def radial_field_index(radius_m, main_field_t) -> tuple[np.ndarray, np.ndarray]:
    """Return centered ``d log|B| / d log r`` on a positive radial line."""
    radius = np.asarray(radius_m, dtype=float).reshape(-1)
    field = np.asarray(main_field_t, dtype=float).reshape(-1)
    if radius.size != field.size or radius.size < 3:
        raise ValueError("radius_m and main_field_t need the same length >= 3")
    if (not np.all(np.isfinite(radius)) or not np.all(np.isfinite(field))
            or np.any(radius <= 0.0) or np.any(np.diff(radius) <= 0.0)
            or np.any(field == 0.0)):
        raise ValueError(
            "radius must be finite, positive, increasing and field must be finite/nonzero")
    log_radius = np.log(radius)
    log_field = np.log(np.abs(field))
    index = ((log_field[2:] - log_field[:-2])
             / (log_radius[2:] - log_radius[:-2]))
    return radius[1:-1], index
