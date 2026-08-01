"""Time-domain charged-particle tracking over solved field callables.

Radia owns the electromagnetic field.  SciPy owns adaptive time integration.
This module supplies the small SI-unit adapter needed for closed trajectories,
electrostatic acceleration, and stop-surface accounting that cannot be
represented by the longitudinal-coordinate Xsuite bridge.

Use :mod:`radia.xsuite_bridge` for accelerator beamline tracking.  Use this
module when physical time, rather than a monotone beamline coordinate, is the
independent variable.
"""
from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np


FieldMap = Callable[[object, object, object], tuple[object, object, object]]
SPEED_OF_LIGHT_M_S = 299_792_458.0


@dataclass(frozen=True)
class ParticleSpecies:
    """Particle charge and rest mass in SI units."""

    charge_c: float
    mass_kg: float

    def __post_init__(self) -> None:
        if not np.isfinite((self.charge_c, self.mass_kg)).all():
            raise ValueError("charge_c and mass_kg must be finite")
        if self.charge_c == 0.0:
            raise ValueError("charge_c must be nonzero")
        if self.mass_kg <= 0.0:
            raise ValueError("mass_kg must be positive")


@dataclass(frozen=True)
class TrackingBox:
    """Axis-aligned stop box in metres."""

    minimum_m: tuple[float, float, float]
    maximum_m: tuple[float, float, float]

    def __post_init__(self) -> None:
        if len(self.minimum_m) != 3 or len(self.maximum_m) != 3:
            raise ValueError("tracking box bounds must have three components")
        if not np.isfinite((*self.minimum_m, *self.maximum_m)).all():
            raise ValueError("tracking box bounds must be finite")
        if any(lo >= hi for lo, hi in zip(self.minimum_m, self.maximum_m)):
            raise ValueError("each tracking box minimum must be below its maximum")


@dataclass(frozen=True)
class TrackingPlane:
    """Terminal plane used to record an exact trajectory crossing.

    ``normal`` points toward the accepted crossing direction when ``direction``
    is ``+1`` and away from it when ``direction`` is ``-1``.  Use ``0`` to
    accept a crossing in either direction.  Plane vectors need not be
    normalized; normalization is performed internally.
    """

    point_m: tuple[float, float, float]
    normal: tuple[float, float, float]
    direction: int = 0

    def __post_init__(self) -> None:
        if len(self.point_m) != 3 or len(self.normal) != 3:
            raise ValueError("tracking plane point and normal must have three components")
        if not np.isfinite((*self.point_m, *self.normal)).all():
            raise ValueError("tracking plane point and normal must be finite")
        if np.linalg.norm(self.normal) == 0.0:
            raise ValueError("tracking plane normal must be nonzero")
        if self.direction not in (-1, 0, 1):
            raise ValueError("tracking plane direction must be -1, 0, or +1")

    @property
    def unit_normal(self) -> np.ndarray:
        normal = np.asarray(self.normal, dtype=float)
        return normal / np.linalg.norm(normal)


def speed_from_kinetic_voltage(
    species: ParticleSpecies,
    kinetic_voltage_v: float,
    *,
    relativistic: bool,
) -> float:
    """Convert an accelerating-voltage magnitude to speed.

    The kinetic energy is ``abs(q) * voltage``.  A signed voltage is rejected
    because the launch direction is represented separately.
    """
    voltage = float(kinetic_voltage_v)
    if not np.isfinite(voltage) or voltage < 0.0:
        raise ValueError("kinetic_voltage_v must be finite and nonnegative")
    energy_j = abs(species.charge_c) * voltage
    if not relativistic:
        return float(np.sqrt(2.0 * energy_j / species.mass_kg))
    rest_energy_j = species.mass_kg * SPEED_OF_LIGHT_M_S**2
    gamma = 1.0 + energy_j / rest_energy_j
    return float(SPEED_OF_LIGHT_M_S * np.sqrt(1.0 - 1.0 / gamma**2))


def velocity_from_kinetic_voltage(
    species: ParticleSpecies,
    kinetic_voltage_v: float,
    direction: Sequence[float],
    *,
    relativistic: bool,
) -> np.ndarray:
    """Return a launch velocity from voltage magnitude and direction."""
    unit = np.asarray(direction, dtype=float)
    if unit.shape != (3,) or not np.all(np.isfinite(unit)):
        raise ValueError("direction must contain three finite components")
    norm = float(np.linalg.norm(unit))
    if norm == 0.0:
        raise ValueError("direction must be nonzero")
    return unit / norm * speed_from_kinetic_voltage(
        species, kinetic_voltage_v, relativistic=relativistic
    )


def uniform_magnetic_trajectory(
    species: ParticleSpecies,
    kinetic_voltage_v: float,
    times_s: object,
    *,
    magnetic_flux_density_t: float,
) -> dict[str, np.ndarray]:
    """Closed-form nonrelativistic orbit for ``v(0)=+x`` and ``B=+z``."""
    times = np.asarray(times_s, dtype=float)
    if times.ndim != 1 or not np.all(np.isfinite(times)):
        raise ValueError("times_s must be a finite one-dimensional array")
    bz = float(magnetic_flux_density_t)
    if not np.isfinite(bz) or bz == 0.0:
        raise ValueError("magnetic_flux_density_t must be finite and nonzero")
    speed = speed_from_kinetic_voltage(
        species, kinetic_voltage_v, relativistic=False
    )
    omega = -species.charge_c * bz / species.mass_kg
    phase = omega * times
    zeros = np.zeros_like(times)
    position = np.column_stack(
        (
            speed * np.sin(phase) / omega,
            speed * (1.0 - np.cos(phase)) / omega,
            zeros,
        )
    )
    velocity = np.column_stack(
        (speed * np.cos(phase), speed * np.sin(phase), zeros)
    )
    return {"position_m": position, "velocity_m_s": velocity}


def _field_vector(field: FieldMap | None, position: np.ndarray) -> np.ndarray:
    if field is None:
        return np.zeros(3, dtype=float)
    values = np.asarray(field(*position), dtype=float)
    if values.shape != (3,) or not np.all(np.isfinite(values)):
        raise ValueError("field callable must return three finite components")
    return values


def _kinetic_energy_j(
    normalized_momentum: np.ndarray,
    species: ParticleSpecies,
    relativistic: bool,
) -> float:
    if relativistic:
        gamma = float(
            np.sqrt(1.0 + np.dot(normalized_momentum, normalized_momentum))
        )
        return (gamma - 1.0) * species.mass_kg * SPEED_OF_LIGHT_M_S**2
    return float(
        0.5
        * species.mass_kg
        * SPEED_OF_LIGHT_M_S**2
        * np.dot(normalized_momentum, normalized_momentum)
    )


def track_lorentz_ivp(
    species: ParticleSpecies,
    initial_position_m: Sequence[float],
    initial_velocity_m_s: Sequence[float],
    times_s: Sequence[float],
    *,
    electric_field_v_m: FieldMap | None = None,
    magnetic_flux_density_t: FieldMap | None = None,
    relativistic: bool = False,
    stop_box: TrackingBox | None = None,
    stop_plane: TrackingPlane | None = None,
    rtol: float = 1.0e-11,
    atol: float = 1.0e-14,
) -> dict[str, Any]:
    """Integrate the Lorentz equation on an explicit physical-time grid.

    Momentum is the integrated state, so the relativistic mode uses the same
    equation ``dp/dt = q(E + v x B)`` without a fragile acceleration formula.
    A stop box or plane terminates at the first crossed surface and reports the
    exact event state.  In the nonrelativistic case this is exactly
    ``dr/dt = v, dv/dt = (Q/m) v cross B + (Q/m) E``; therefore a legacy
    velocity-state implementation's scalar ``q`` is ``Q/m``, not charge alone.
    """
    try:
        from scipy.integrate import solve_ivp
    except ImportError as exc:
        raise ImportError(
            "time-domain particle tracking requires the Radia 'beam' extra"
        ) from exc

    position0 = np.asarray(initial_position_m, dtype=float)
    velocity0 = np.asarray(initial_velocity_m_s, dtype=float)
    times = np.asarray(times_s, dtype=float)
    if position0.shape != (3,) or velocity0.shape != (3,):
        raise ValueError("initial position and velocity must have three components")
    if not np.all(np.isfinite(position0)) or not np.all(np.isfinite(velocity0)):
        raise ValueError("initial position and velocity must be finite")
    if times.ndim != 1 or times.size < 2 or not np.all(np.isfinite(times)):
        raise ValueError("times_s must contain at least two finite values")
    if np.any(np.diff(times) <= 0.0):
        raise ValueError("times_s must be strictly increasing")
    speed0 = float(np.linalg.norm(velocity0))
    if relativistic and speed0 >= SPEED_OF_LIGHT_M_S:
        raise ValueError("relativistic initial speed must be below c")
    gamma0 = (
        1.0 / np.sqrt(1.0 - (speed0 / SPEED_OF_LIGHT_M_S) ** 2)
        if relativistic
        else 1.0
    )
    normalized_momentum0 = gamma0 * velocity0 / SPEED_OF_LIGHT_M_S
    y0 = np.concatenate((position0, normalized_momentum0))

    def rhs(_time: float, state: np.ndarray) -> np.ndarray:
        position = state[:3]
        normalized_momentum = state[3:]
        if relativistic:
            gamma = np.sqrt(
                1.0 + np.dot(normalized_momentum, normalized_momentum)
            )
        else:
            gamma = 1.0
        velocity = SPEED_OF_LIGHT_M_S * normalized_momentum / gamma
        electric = _field_vector(electric_field_v_m, position)
        magnetic = _field_vector(magnetic_flux_density_t, position)
        force = species.charge_c * (electric + np.cross(velocity, magnetic))
        normalized_momentum_rate = force / (
            species.mass_kg * SPEED_OF_LIGHT_M_S
        )
        return np.concatenate((velocity, normalized_momentum_rate))

    events = None
    event_names: list[str] = []
    if stop_box is not None:
        events = []
        for axis, axis_name in enumerate(("x", "y", "z")):
            for side, bound in (
                ("minimum", stop_box.minimum_m[axis]),
                ("maximum", stop_box.maximum_m[axis]),
            ):
                sign = 1.0 if side == "minimum" else -1.0

                def face_event(
                    _time: float,
                    state: np.ndarray,
                    *,
                    axis: int = axis,
                    bound: float = bound,
                    sign: float = sign,
                ) -> float:
                    return sign * (state[axis] - bound)

                face_event.terminal = True
                face_event.direction = -1.0
                events.append(face_event)
                event_names.append(f"{axis_name}_{side}")
    if stop_plane is not None:
        if events is None:
            events = []
        plane_point = np.asarray(stop_plane.point_m, dtype=float)
        plane_normal = stop_plane.unit_normal

        def plane_event(_time: float, state: np.ndarray) -> float:
            return float(np.dot(state[:3] - plane_point, plane_normal))

        plane_event.terminal = True
        plane_event.direction = float(stop_plane.direction)
        events.append(plane_event)
        event_names.append("plane")

    solution = solve_ivp(
        rhs,
        (float(times[0]), float(times[-1])),
        y0,
        method="DOP853",
        t_eval=times,
        events=events,
        rtol=float(rtol),
        atol=float(atol),
    )
    positions = solution.y[:3].T
    normalized_momenta = solution.y[3:].T
    if relativistic:
        gamma = np.sqrt(1.0 + np.sum(normalized_momenta**2, axis=1))
    else:
        gamma = np.ones(normalized_momenta.shape[0])
    velocities = SPEED_OF_LIGHT_M_S * normalized_momenta / gamma[:, None]
    energies = np.array(
        [
            _kinetic_energy_j(row, species, relativistic)
            for row in normalized_momenta
        ]
    )
    stop_event = None
    if events is not None:
        for name, event_times, event_states in zip(
            event_names, solution.t_events, solution.y_events
        ):
            if event_times.size:
                event_normalized_momentum = event_states[0, 3:]
                if relativistic:
                    event_gamma = np.sqrt(
                        1.0 + np.dot(
                            event_normalized_momentum, event_normalized_momentum
                        )
                    )
                else:
                    event_gamma = 1.0
                event_velocity = (
                    SPEED_OF_LIGHT_M_S
                    * event_normalized_momentum
                    / event_gamma
                )
                stop_event = {
                    "face": name,
                    "time_s": float(event_times[0]),
                    "position_m": event_states[0, :3].tolist(),
                    "velocity_m_s": event_velocity.tolist(),
                }
                break
    return {
        "schema": "radia-time-domain-lorentz-track/v1",
        "backend": "scipy.integrate.solve_ivp:DOP853",
        "units": {
            "position": "m",
            "velocity": "m/s",
            "time": "s",
            "electric_field": "V/m",
            "magnetic_flux_density": "T",
            "kinetic_energy": "J",
        },
        "relativistic": bool(relativistic),
        "success": bool(solution.success),
        "message": solution.message,
        "time_s": solution.t.tolist(),
        "position_m": positions.tolist(),
        "velocity_m_s": velocities.tolist(),
        "kinetic_energy_j": energies.tolist(),
        "maximum_relative_kinetic_energy_drift": float(
            np.max(np.abs(energies - energies[0])) / max(abs(energies[0]), 1.0e-300)
        ),
        "stop_event": stop_event,
        "function_evaluation_count": int(solution.nfev),
    }


def _momentum_scaled_velocity(
    velocity_m_s: Sequence[float],
    momentum_scale: float,
    *,
    relativistic: bool,
) -> np.ndarray:
    """Scale momentum while retaining the launch direction."""
    velocity = np.asarray(velocity_m_s, dtype=float)
    scale = float(momentum_scale)
    if velocity.shape != (3,) or not np.all(np.isfinite(velocity)):
        raise ValueError("initial_velocity_m_s must contain three finite components")
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError("momentum scale must be finite and positive")
    if not relativistic:
        return scale * velocity
    speed = float(np.linalg.norm(velocity))
    if speed >= SPEED_OF_LIGHT_M_S:
        raise ValueError("relativistic initial speed must be below c")
    beta_gamma = (
        velocity / SPEED_OF_LIGHT_M_S
        / np.sqrt(1.0 - (speed / SPEED_OF_LIGHT_M_S) ** 2)
    )
    scaled = scale * beta_gamma
    gamma = np.sqrt(1.0 + np.dot(scaled, scaled))
    return SPEED_OF_LIGHT_M_S * scaled / gamma


def track_two_momentum_exit_dispersion(
    species: ParticleSpecies,
    initial_position_m: Sequence[float],
    nominal_velocity_m_s: Sequence[float],
    times_s: Sequence[float],
    exit_plane: TrackingPlane,
    *,
    relative_momentum_offset: float,
    transverse_direction: Sequence[float],
    electric_field_v_m: FieldMap | None = None,
    magnetic_flux_density_t: FieldMap | None = None,
    relativistic: bool = False,
    stop_box: TrackingBox | None = None,
    rtol: float = 1.0e-11,
    atol: float = 1.0e-14,
) -> dict[str, Any]:
    """Track ``p0*(1 +/- delta)`` and compare their exit-plane positions.

    The returned scalar ``eta_m`` is the symmetric two-particle dispersion
    criterion ``dot(r_plus-r_minus, e_transverse)/(2*delta)``.  It is a direct
    physical acceptance metric: both particles arrive at the same transverse
    exit position when ``eta_m == 0``.  It is not used as a finite-difference
    substitute for an analytic topology or shape sensitivity.
    """
    delta = float(relative_momentum_offset)
    if not np.isfinite(delta) or not 0.0 < delta < 1.0:
        raise ValueError("relative_momentum_offset must lie strictly between 0 and 1")
    transverse = np.asarray(transverse_direction, dtype=float)
    if transverse.shape != (3,) or not np.all(np.isfinite(transverse)):
        raise ValueError("transverse_direction must contain three finite components")
    # Remove any component normal to the exit plane so eta is intrinsic to it.
    normal = exit_plane.unit_normal
    transverse = transverse - np.dot(transverse, normal) * normal
    transverse_norm = float(np.linalg.norm(transverse))
    if transverse_norm == 0.0:
        raise ValueError("transverse_direction must not be parallel to the exit plane normal")
    transverse /= transverse_norm

    common = dict(
        electric_field_v_m=electric_field_v_m,
        magnetic_flux_density_t=magnetic_flux_density_t,
        relativistic=relativistic,
        stop_box=stop_box,
        stop_plane=exit_plane,
        rtol=rtol,
        atol=atol,
    )
    tracks = {}
    for name, scale in (("minus", 1.0 - delta), ("plus", 1.0 + delta)):
        velocity = _momentum_scaled_velocity(
            nominal_velocity_m_s, scale, relativistic=relativistic
        )
        track = track_lorentz_ivp(
            species, initial_position_m, velocity, times_s, **common
        )
        event = track["stop_event"]
        if not track["success"] or event is None or event["face"] != "plane":
            raise RuntimeError(
                f"{name} momentum trajectory did not reach the exit plane: "
                f"{track['message']}"
            )
        tracks[name] = track

    minus_position = np.asarray(
        tracks["minus"]["stop_event"]["position_m"], dtype=float
    )
    plus_position = np.asarray(
        tracks["plus"]["stop_event"]["position_m"], dtype=float
    )
    difference = plus_position - minus_position
    dispersion_vector = difference / (2.0 * delta)
    return {
        "schema": "radia-two-momentum-exit-dispersion/v1",
        "relative_momentum_offset": delta,
        "minus_position_m": minus_position.tolist(),
        "plus_position_m": plus_position.tolist(),
        "position_difference_m": difference.tolist(),
        "dispersion_vector_m": dispersion_vector.tolist(),
        "transverse_direction": transverse.tolist(),
        "eta_m": float(np.dot(dispersion_vector, transverse)),
        "coincident_exit_error_m": float(np.linalg.norm(difference)),
        "minus_track": tracks["minus"],
        "plus_track": tracks["plus"],
    }
