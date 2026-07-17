"""Energy-consistent angle-periodic motor reduced-order models.

The generalized currents are ordered as physical phase currents followed by
internal eddy-current coordinates.  The latter may be HCurl Eddy Bubble,
conductor-cycle bridge, SIBC/CLN, or other passive reduced coordinates.  A
single periodic flux law

    lambda(theta, q) = L(theta) q + psi_pm(theta) + psi_hys

drives both the voltage equation and the mechanical torque.  This makes the
motion voltage, reluctance torque, PM torque, and eddy-current reaction obey
the same power balance instead of being coupled by unrelated lookup tables.

The module is NumPy-only.  NGSolve/HDiv hysteresis is connected lazily through
``HDivHysteresisPort`` so that Simulink/FMI wrappers can use the core motor ROM
without importing a finite-element runtime.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Protocol

import numpy as np


_TWO_PI = 2.0 * np.pi


def _finite_real(value, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} contains non-finite values")
    return result


def _vector(value, size: int, name: str) -> np.ndarray:
    result = _finite_real(value, name).reshape(-1)
    if result.size != size:
        raise ValueError(f"{name} must contain {size} values")
    return result


def _square(value, size: int, name: str) -> np.ndarray:
    result = _finite_real(value, name)
    if result.shape != (size, size):
        raise ValueError(f"{name} must have shape ({size}, {size})")
    return result


@dataclass(frozen=True, slots=True)
class MotorPortContract:
    """Stable electromechanical port ordering for C/Simulink/FMI adapters."""

    phase_names: tuple[str, ...]
    eddy_names: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        phase = tuple(str(name) for name in self.phase_names)
        eddy = tuple(str(name) for name in self.eddy_names)
        if not phase:
            raise ValueError("phase_names must not be empty")
        names = phase + eddy
        if any(not name for name in names):
            raise ValueError("port names must not be empty")
        if len(set(names)) != len(names):
            raise ValueError("phase and eddy port names must be unique")
        object.__setattr__(self, "phase_names", phase)
        object.__setattr__(self, "eddy_names", eddy)

    @property
    def n_phase(self) -> int:
        return len(self.phase_names)

    @property
    def n_eddy(self) -> int:
        return len(self.eddy_names)

    @property
    def n_generalized(self) -> int:
        return self.n_phase + self.n_eddy

    @property
    def generalized_names(self) -> tuple[str, ...]:
        return self.phase_names + self.eddy_names

    def generalized_voltage(self, phase_voltages) -> np.ndarray:
        result = np.zeros(self.n_generalized)
        result[: self.n_phase] = _vector(
            phase_voltages, self.n_phase, "phase_voltages_V"
        )
        return result

    def diagnostics(self) -> dict[str, object]:
        return {
            "phase_names": self.phase_names,
            "eddy_names": self.eddy_names,
            "n_phase": self.n_phase,
            "n_eddy": self.n_eddy,
            "inputs": (
                "phase_voltages_V",
                "load_torque_Nm",
                "ambient_temperature_K",
            ),
            "states": (
                "rotor_angle_rad",
                "rotor_speed_rad_s",
                "generalized_currents_A",
                "temperature_K",
                "hysteresis_state",
            ),
            "outputs": (
                "phase_currents_A",
                "eddy_currents_A",
                "phase_flux_linkage_Wb",
                "electromagnetic_torque_Nm",
                "resistive_loss_W",
                "hysteresis_loss_W",
            ),
        }


@dataclass(frozen=True, slots=True)
class PeriodicAngleTable:
    """Odd-point Fourier interpolation of one real periodic tensor.

    An odd number of samples removes the Nyquist derivative ambiguity.  Skew is
    applied analytically to every Fourier mode by its sinc factor, so a skewed
    machine does not require a second angle sweep.
    """

    angles_rad: np.ndarray
    values: np.ndarray
    period_rad: float = _TWO_PI
    _origin: float = field(init=False, repr=False)
    _modes: np.ndarray = field(init=False, repr=False)
    _coefficients: np.ndarray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        angles = _finite_real(self.angles_rad, "angles_rad").reshape(-1)
        values = _finite_real(self.values, "values")
        period = float(self.period_rad)
        if not np.isfinite(period) or period <= 0.0:
            raise ValueError("period_rad must be positive")
        if angles.size < 3 or angles.size % 2 == 0:
            raise ValueError("PeriodicAngleTable requires an odd number of at least 3 samples")
        if values.shape[0] != angles.size:
            raise ValueError("values first dimension must match angles_rad")
        unwrapped = np.unwrap(angles * (_TWO_PI / period)) * (period / _TWO_PI)
        if np.any(np.diff(unwrapped) <= 0.0):
            raise ValueError("angles_rad must be strictly increasing")
        spacings = np.diff(np.r_[unwrapped, unwrapped[0] + period])
        expected = period / angles.size
        if not np.allclose(spacings, expected, rtol=1.0e-10, atol=1.0e-12 * period):
            raise ValueError("angles_rad must be a uniform periodic grid without a duplicate endpoint")
        coefficients = np.fft.fft(values, axis=0) / angles.size
        modes = np.fft.fftfreq(angles.size, d=1.0 / angles.size)
        object.__setattr__(self, "angles_rad", angles)
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "period_rad", period)
        object.__setattr__(self, "_origin", float(unwrapped[0]))
        object.__setattr__(self, "_modes", modes)
        object.__setattr__(self, "_coefficients", coefficients)

    @property
    def value_shape(self) -> tuple[int, ...]:
        return self.values.shape[1:]

    def evaluate(
        self,
        angle_rad: float,
        *,
        derivative: int = 0,
        skew_span_rad: float = 0.0,
    ) -> np.ndarray:
        derivative = int(derivative)
        if derivative < 0:
            raise ValueError("derivative must be non-negative")
        angle = float(angle_rad)
        skew = float(skew_span_rad)
        if not np.isfinite(angle) or not np.isfinite(skew) or skew < 0.0:
            raise ValueError("angle_rad must be finite and skew_span_rad non-negative")
        phase = _TWO_PI * (angle - self._origin) / self.period_rad
        angular_modes = _TWO_PI * self._modes / self.period_rad
        factor = (1j * angular_modes) ** derivative
        if skew > 0.0:
            factor = factor * np.sinc(angular_modes * skew / _TWO_PI)
        wave = np.exp(1j * phase * self._modes) * factor
        result = np.tensordot(wave, self._coefficients, axes=(0, 0))
        imag = float(np.max(np.abs(np.imag(result)))) if np.size(result) else 0.0
        scale = max(float(np.max(np.abs(result))) if np.size(result) else 0.0, 1.0)
        if imag > 5.0e-11 * scale:
            raise RuntimeError("periodic interpolation produced a non-negligible imaginary part")
        return np.asarray(np.real(result), dtype=float)

    def diagnostics(self) -> dict[str, object]:
        return {
            "samples": int(self.angles_rad.size),
            "period_rad": self.period_rad,
            "value_shape": self.value_shape,
            "highest_harmonic": int((self.angles_rad.size - 1) // 2),
            "interpolation": "periodic-fourier",
        }


@dataclass(frozen=True, slots=True)
class HysteresisEvaluation:
    """One pure trial evaluation from a committed hysteresis state."""

    flux_linkage_Wb: np.ndarray
    torque_Nm: float
    stored_energy_J: float
    dissipated_energy_increment_J: float
    state: object
    diagnostics: dict[str, object] = field(default_factory=dict)


class HysteresisPort(Protocol):
    """Functional hysteresis adapter; trial calls must not mutate state."""

    def evaluate(
        self,
        rotor_angle_rad: float,
        generalized_currents_A: np.ndarray,
        committed_state: object | None,
    ) -> HysteresisEvaluation: ...


@dataclass(frozen=True, slots=True)
class HDivHysteresisState:
    restart_state: object
    B: np.ndarray
    H: np.ndarray
    stored_energy_J: float
    cumulative_dissipation_J: float


class HDivHysteresisPort:
    """Adapter from persistent ``vim.HDivSolver`` history to motor ports.

    ``field_from_currents`` maps physical phase currents and rotor angle to the
    applied 3D H field accepted by ``HDivSolver.SolveHysteresis``.
    ``flux_from_result`` maps the converged HDiv result back to generalized
    flux linkages.  Repeated evaluations from the same committed state are
    safe because the HDiv material protocol advances history only in the
    returned restart state.
    """

    def __init__(
        self,
        solver,
        material,
        ports: MotorPortContract,
        field_from_currents: Callable[[np.ndarray, float], object],
        flux_from_result: Callable[[object, float], np.ndarray],
        *,
        torque_from_result: Callable[[object, float], float] | None = None,
        solve_options: dict[str, object] | None = None,
    ) -> None:
        self.solver = solver
        self.material = material
        self.ports = ports
        self.field_from_currents = field_from_currents
        self.flux_from_result = flux_from_result
        self.torque_from_result = torque_from_result
        self.solve_options = dict(solve_options or {})
        try:
            import ngsolve as ng

            self._element_volumes = np.asarray(
                ng.Integrate(ng.CoefficientFunction(1.0), solver.mesh, element_wise=True),
                dtype=float,
            )
        except Exception:
            self._element_volumes = None

    def _stored_energy(self, result) -> float:
        state = result["state"]
        if result.get("state_layout") != "element" or self._element_volumes is None:
            return float("nan")
        if not hasattr(self.material, "stored_energy"):
            return float("nan")
        density = np.asarray(
            self.material.stored_energy(state["B"], state["material_states"]),
            dtype=float,
        )
        if density.shape != self._element_volumes.shape:
            return float("nan")
        return float(self._element_volumes @ density)

    def evaluate(
        self,
        rotor_angle_rad: float,
        generalized_currents_A: np.ndarray,
        committed_state: HDivHysteresisState | None,
    ) -> HysteresisEvaluation:
        currents = _vector(
            generalized_currents_A,
            self.ports.n_generalized,
            "generalized_currents_A",
        )
        applied = self.field_from_currents(
            currents[: self.ports.n_phase], float(rotor_angle_rad)
        )
        initial = None if committed_state is None else committed_state.restart_state
        result = self.solver.SolveHysteresis(
            [applied],
            material=self.material,
            initial_state=initial,
            **self.solve_options,
        )
        flux = _vector(
            self.flux_from_result(result, float(rotor_angle_rad)),
            self.ports.n_generalized,
            "hysteresis flux_linkage_Wb",
        )
        torque = (
            0.0
            if self.torque_from_result is None
            else float(self.torque_from_result(result, float(rotor_angle_rad)))
        )
        step = result["steps"][-1]
        B = np.asarray(step["B"], dtype=float)
        H = np.asarray(step["H"], dtype=float)
        stored = self._stored_energy(result)
        dissipated = float("nan")
        cumulative = 0.0 if committed_state is None else committed_state.cumulative_dissipation_J
        if (
            committed_state is not None
            and self._element_volumes is not None
            and B.shape == committed_state.B.shape == H.shape == committed_state.H.shape
        ):
            magnetic_work = float(
                np.sum(
                    self._element_volumes[:, None]
                    * 0.5
                    * (committed_state.H + H)
                    * (B - committed_state.B)
                )
            )
            if np.isfinite(stored) and np.isfinite(committed_state.stored_energy_J):
                dissipated = magnetic_work - (stored - committed_state.stored_energy_J)
                cumulative += dissipated
        next_state = HDivHysteresisState(
            restart_state=result["state"],
            B=B.copy(),
            H=H.copy(),
            stored_energy_J=stored,
            cumulative_dissipation_J=cumulative,
        )
        return HysteresisEvaluation(
            flux_linkage_Wb=flux,
            torque_Nm=torque,
            stored_energy_J=stored,
            dissipated_energy_increment_J=dissipated,
            state=next_state,
            diagnostics={
                "operator_build_count": int(self.solver.operator_build_count),
                "state_layout": result.get("state_layout"),
                "permanent_magnet_model": result.get("permanent_magnet_model"),
            },
        )


@dataclass(frozen=True, slots=True)
class MotorROMState:
    time_s: float
    rotor_angle_rad: float
    rotor_speed_rad_s: float
    generalized_currents_A: np.ndarray
    temperature_K: float
    hysteresis_state: object | None = None
    hysteresis_flux_linkage_Wb: np.ndarray | None = None
    hysteresis_stored_energy_J: float = 0.0


@dataclass(frozen=True, slots=True)
class MotorROMInput:
    phase_voltages_V: np.ndarray
    load_torque_Nm: float = 0.0
    ambient_temperature_K: float | None = None


@dataclass(frozen=True, slots=True)
class MotorROMStepOutput:
    phase_currents_A: np.ndarray
    eddy_currents_A: np.ndarray
    phase_flux_linkage_Wb: np.ndarray
    electromagnetic_torque_Nm: float
    torque_components_Nm: dict[str, float]
    speed_voltage_V: np.ndarray
    resistive_loss_W: float
    hysteresis_loss_W: float
    electrical_input_power_W: float
    mechanical_load_power_W: float
    stored_energy_J: float
    energy_balance_residual_W: float
    nonlinear_iterations: int


class AnglePeriodicMotorROM:
    """Passive angle-periodic electromechanical ROM with internal eddy currents."""

    def __init__(
        self,
        ports: MotorPortContract,
        inductance_H: PeriodicAngleTable,
        resistance_ohm: PeriodicAngleTable,
        pm_flux_linkage_Wb: PeriodicAngleTable,
        *,
        inertia_kg_m2: float,
        viscous_friction_Nm_s: float = 0.0,
        motion_flux_gradient_Wb_per_rad: PeriodicAngleTable | None = None,
        cogging_coenergy_J: PeriodicAngleTable | None = None,
        skew_span_rad: float = 0.0,
        end_winding_inductance_H=None,
        end_winding_resistance_ohm=None,
        reference_temperature_K: float = 293.15,
        resistance_temperature_coefficient_per_K=0.0,
        thermal_capacity_J_per_K: float | None = None,
        thermal_conductance_W_per_K: float = 0.0,
        hysteresis_port: HysteresisPort | None = None,
    ) -> None:
        self.ports = ports
        self.inductance_H = inductance_H
        self.resistance_ohm = resistance_ohm
        self.pm_flux_linkage_Wb = pm_flux_linkage_Wb
        self.motion_flux_gradient_Wb_per_rad = motion_flux_gradient_Wb_per_rad
        self.cogging_coenergy_J = cogging_coenergy_J
        self.inertia_kg_m2 = float(inertia_kg_m2)
        self.viscous_friction_Nm_s = float(viscous_friction_Nm_s)
        self.skew_span_rad = float(skew_span_rad)
        self.reference_temperature_K = float(reference_temperature_K)
        self.thermal_capacity_J_per_K = (
            None if thermal_capacity_J_per_K is None else float(thermal_capacity_J_per_K)
        )
        self.thermal_conductance_W_per_K = float(thermal_conductance_W_per_K)
        self.hysteresis_port = hysteresis_port
        n = ports.n_generalized
        self.end_winding_inductance_H = (
            np.zeros((n, n))
            if end_winding_inductance_H is None
            else _square(end_winding_inductance_H, n, "end_winding_inductance_H")
        )
        self.end_winding_resistance_ohm = (
            np.zeros((n, n))
            if end_winding_resistance_ohm is None
            else _square(end_winding_resistance_ohm, n, "end_winding_resistance_ohm")
        )
        alpha = np.asarray(resistance_temperature_coefficient_per_K, dtype=float)
        if alpha.ndim == 0:
            alpha = np.full(n, float(alpha))
        self.resistance_temperature_coefficient_per_K = _vector(
            alpha, n, "resistance_temperature_coefficient_per_K"
        )
        if self.inertia_kg_m2 <= 0.0:
            raise ValueError("inertia_kg_m2 must be positive")
        if self.viscous_friction_Nm_s < 0.0:
            raise ValueError("viscous_friction_Nm_s must be non-negative")
        if self.skew_span_rad < 0.0:
            raise ValueError("skew_span_rad must be non-negative")
        if self.thermal_capacity_J_per_K is not None and self.thermal_capacity_J_per_K <= 0.0:
            raise ValueError("thermal_capacity_J_per_K must be positive")
        if self.thermal_conductance_W_per_K < 0.0:
            raise ValueError("thermal_conductance_W_per_K must be non-negative")
        if inductance_H.value_shape != (n, n):
            raise ValueError(f"inductance_H values must have shape ({n}, {n})")
        if resistance_ohm.value_shape != (n, n):
            raise ValueError(f"resistance_ohm values must have shape ({n}, {n})")
        if pm_flux_linkage_Wb.value_shape != (n,):
            raise ValueError(f"pm_flux_linkage_Wb values must have shape ({n},)")
        if (
            motion_flux_gradient_Wb_per_rad is not None
            and motion_flux_gradient_Wb_per_rad.value_shape != (n,)
        ):
            raise ValueError(
                f"motion_flux_gradient_Wb_per_rad values must have shape ({n},)"
            )
        if cogging_coenergy_J is not None and cogging_coenergy_J.value_shape != ():
            raise ValueError("cogging_coenergy_J values must be scalar at every angle")
        periods = {
            inductance_H.period_rad,
            resistance_ohm.period_rad,
            pm_flux_linkage_Wb.period_rad,
        }
        if motion_flux_gradient_Wb_per_rad is not None:
            periods.add(motion_flux_gradient_Wb_per_rad.period_rad)
        if cogging_coenergy_J is not None:
            periods.add(cogging_coenergy_J.period_rad)
        if len(periods) != 1:
            raise ValueError("all angle tables must use the same period_rad")
        self.period_rad = periods.pop()
        self._validate_passivity()

    def _validate_passivity(self) -> None:
        count = max(65, 4 * self.inductance_H.angles_rad.size + 1)
        for angle in np.linspace(0.0, self.period_rad, count, endpoint=False):
            L = self.inductance(angle)
            R = self.resistance(angle, self.reference_temperature_K)
            if not np.allclose(L, L.T, rtol=1.0e-10, atol=1.0e-12):
                raise ValueError("interpolated inductance must be symmetric")
            if not np.allclose(R, R.T, rtol=1.0e-10, atol=1.0e-12):
                raise ValueError("interpolated resistance must be symmetric")
            if float(np.min(np.linalg.eigvalsh(L))) <= 0.0:
                raise ValueError("interpolated inductance must be positive definite")
            if float(np.min(np.linalg.eigvalsh(R))) < -1.0e-11 * max(
                float(np.linalg.norm(R)), 1.0
            ):
                raise ValueError("interpolated resistance must be positive semidefinite")

    def _eval(self, table: PeriodicAngleTable, angle: float, derivative: int = 0):
        return table.evaluate(
            angle,
            derivative=derivative,
            skew_span_rad=self.skew_span_rad,
        )

    def inductance(self, angle_rad: float, *, derivative: int = 0) -> np.ndarray:
        value = self._eval(self.inductance_H, angle_rad, derivative)
        if derivative == 0:
            value = value + self.end_winding_inductance_H
        return 0.5 * (value + value.T)

    def resistance(self, angle_rad: float, temperature_K: float) -> np.ndarray:
        base = self._eval(self.resistance_ohm, angle_rad) + self.end_winding_resistance_ohm
        factor = 1.0 + self.resistance_temperature_coefficient_per_K * (
            float(temperature_K) - self.reference_temperature_K
        )
        if np.any(factor <= 0.0):
            raise ValueError("temperature law produced non-positive resistance scale")
        scale = np.sqrt(factor)
        return 0.5 * ((scale[:, None] * base * scale[None, :]) +
                      (scale[:, None] * base * scale[None, :]).T)

    def pm_flux(self, angle_rad: float, *, derivative: int = 0) -> np.ndarray:
        return self._eval(self.pm_flux_linkage_Wb, angle_rad, derivative)

    def motion_flux_gradient(self, angle_rad: float) -> np.ndarray:
        if self.motion_flux_gradient_Wb_per_rad is None:
            return np.zeros(self.ports.n_generalized)
        return self._eval(self.motion_flux_gradient_Wb_per_rad, angle_rad)

    def cogging_coenergy(self, angle_rad: float, *, derivative: int = 0) -> float:
        if self.cogging_coenergy_J is None:
            return 0.0
        return float(self._eval(self.cogging_coenergy_J, angle_rad, derivative))

    def magnetic_energy(self, angle_rad: float, currents) -> float:
        q = _vector(currents, self.ports.n_generalized, "currents")
        return 0.5 * float(q @ self.inductance(angle_rad) @ q)

    def coenergy(self, angle_rad: float, currents) -> float:
        q = _vector(currents, self.ports.n_generalized, "currents")
        return (
            self.magnetic_energy(angle_rad, q)
            + float(q @ self.pm_flux(angle_rad))
            + self.cogging_coenergy(angle_rad)
        )

    def torque_components(self, angle_rad: float, currents) -> dict[str, float]:
        q = _vector(currents, self.ports.n_generalized, "currents")
        reluctance = 0.5 * float(q @ self.inductance(angle_rad, derivative=1) @ q)
        permanent_magnet = float(q @ self.pm_flux(angle_rad, derivative=1))
        motional = float(q @ self.motion_flux_gradient(angle_rad))
        cogging = self.cogging_coenergy(angle_rad, derivative=1)
        return {
            "reluctance": reluctance,
            "permanent_magnet": permanent_magnet,
            "motional_lorentz": motional,
            "cogging": cogging,
            "hysteresis": 0.0,
            "total": reluctance + permanent_magnet + motional + cogging,
        }

    def virtual_work_torque(
        self,
        angle_rad: float,
        currents,
        *,
        delta_angle_rad: float = 1.0e-5,
    ) -> float:
        delta = float(delta_angle_rad)
        if delta <= 0.0:
            raise ValueError("delta_angle_rad must be positive")
        q = _vector(currents, self.ports.n_generalized, "currents")
        return (
            self.coenergy(float(angle_rad) + delta, q)
            - self.coenergy(float(angle_rad) - delta, q)
        ) / (2.0 * delta)

    def torque_audit(
        self,
        angle_rad: float,
        currents,
        *,
        maxwell_torque_Nm: float,
        delta_angle_rad: float = 1.0e-5,
    ) -> dict[str, float]:
        components = self.torque_components(angle_rad, currents)
        coenergy = (
            components["reluctance"]
            + components["permanent_magnet"]
            + components["cogging"]
        )
        virtual = self.virtual_work_torque(
            angle_rad, currents, delta_angle_rad=delta_angle_rad
        )
        maxwell = float(maxwell_torque_Nm)
        values = np.array((maxwell, coenergy, virtual))
        scale = max(float(np.max(np.abs(values))), np.finfo(float).tiny)
        return {
            "maxwell_stress_Nm": maxwell,
            "coenergy_derivative_Nm": coenergy,
            "virtual_work_Nm": virtual,
            "relative_spread": float((np.max(values) - np.min(values)) / scale),
        }

    def initial_state(
        self,
        *,
        rotor_angle_rad: float = 0.0,
        rotor_speed_rad_s: float = 0.0,
        phase_currents_A=None,
        eddy_currents_A=None,
        temperature_K: float | None = None,
    ) -> MotorROMState:
        q = np.zeros(self.ports.n_generalized)
        if phase_currents_A is not None:
            q[: self.ports.n_phase] = _vector(
                phase_currents_A, self.ports.n_phase, "phase_currents_A"
            )
        if eddy_currents_A is not None:
            q[self.ports.n_phase :] = _vector(
                eddy_currents_A, self.ports.n_eddy, "eddy_currents_A"
            )
        temperature = (
            self.reference_temperature_K if temperature_K is None else float(temperature_K)
        )
        hysteresis_state = None
        hysteresis_flux = np.zeros(self.ports.n_generalized)
        hysteresis_energy = 0.0
        if self.hysteresis_port is not None:
            evaluated = self.hysteresis_port.evaluate(
                float(rotor_angle_rad), q, None
            )
            hysteresis_state = evaluated.state
            hysteresis_flux = _vector(
                evaluated.flux_linkage_Wb,
                self.ports.n_generalized,
                "hysteresis flux_linkage_Wb",
            )
            hysteresis_energy = float(evaluated.stored_energy_J)
        return MotorROMState(
            time_s=0.0,
            rotor_angle_rad=float(rotor_angle_rad),
            rotor_speed_rad_s=float(rotor_speed_rad_s),
            generalized_currents_A=q,
            temperature_K=temperature,
            hysteresis_state=hysteresis_state,
            hysteresis_flux_linkage_Wb=hysteresis_flux,
            hysteresis_stored_energy_J=hysteresis_energy,
        )

    def step(
        self,
        state: MotorROMState,
        command: MotorROMInput,
        dt_s: float,
        *,
        max_iterations: int = 30,
        tolerance: float = 1.0e-11,
    ) -> tuple[MotorROMState, MotorROMStepOutput]:
        dt = float(dt_s)
        if not np.isfinite(dt) or dt <= 0.0:
            raise ValueError("dt_s must be positive")
        q0 = _vector(
            state.generalized_currents_A,
            self.ports.n_generalized,
            "state.generalized_currents_A",
        )
        u = self.ports.generalized_voltage(command.phase_voltages_V)
        load = float(command.load_torque_Nm)
        ambient = (
            state.temperature_K
            if command.ambient_temperature_K is None
            else float(command.ambient_temperature_K)
        )
        theta0 = float(state.rotor_angle_rad)
        omega0 = float(state.rotor_speed_rad_s)
        hflux0 = (
            np.zeros(self.ports.n_generalized)
            if state.hysteresis_flux_linkage_Wb is None
            else _vector(
                state.hysteresis_flux_linkage_Wb,
                self.ports.n_generalized,
                "state.hysteresis_flux_linkage_Wb",
            )
        )
        lambda0 = self.inductance(theta0) @ q0 + self.pm_flux(theta0) + hflux0
        theta1 = theta0 + dt * omega0
        omega1 = omega0
        q1 = q0.copy()
        temperature1 = float(state.temperature_K)
        h_eval = None
        iterations = 0
        for iterations in range(1, int(max_iterations) + 1):
            theta_mid = 0.5 * (theta0 + theta1)
            omega_mid = 0.5 * (omega0 + omega1)
            temperature_mid = 0.5 * (state.temperature_K + temperature1)
            R_mid = self.resistance(theta_mid, temperature_mid)
            motion = self.motion_flux_gradient(theta_mid)
            hflux1 = hflux0
            htorque = 0.0
            if self.hysteresis_port is not None:
                h_eval = self.hysteresis_port.evaluate(
                    theta1, q1, state.hysteresis_state
                )
                hflux1 = _vector(
                    h_eval.flux_linkage_Wb,
                    self.ports.n_generalized,
                    "hysteresis flux_linkage_Wb",
                )
                htorque = float(h_eval.torque_Nm)
            rhs = (
                lambda0
                + dt * u
                - 0.5 * dt * (R_mid @ q0)
                - dt * omega_mid * motion
                - self.pm_flux(theta1)
                - hflux1
            )
            q_next = np.linalg.solve(
                self.inductance(theta1) + 0.5 * dt * R_mid,
                rhs,
            )
            q_mid = 0.5 * (q0 + q_next)
            components = self.torque_components(theta_mid, q_mid)
            components["hysteresis"] = htorque
            components["total"] += htorque
            torque = components["total"]
            omega_next = omega0 + dt * (
                torque - load - self.viscous_friction_Nm_s * omega_mid
            ) / self.inertia_kg_m2
            theta_next = theta0 + 0.5 * dt * (omega0 + omega_next)
            loss = float(q_mid @ R_mid @ q_mid)
            temperature_next = float(state.temperature_K)
            if self.thermal_capacity_J_per_K is not None:
                cth = self.thermal_capacity_J_per_K
                hth = self.thermal_conductance_W_per_K
                temperature_next = (
                    state.temperature_K + dt * (loss + hth * ambient) / cth
                ) / (1.0 + dt * hth / cth)
            error = max(
                float(np.linalg.norm(q_next - q1)),
                abs(theta_next - theta1),
                abs(omega_next - omega1),
                abs(temperature_next - temperature1),
            )
            q1 = q_next
            theta1 = theta_next
            omega1 = omega_next
            temperature1 = temperature_next
            if error <= tolerance * max(1.0, float(np.linalg.norm(q1)), abs(omega1)):
                break
        else:
            raise RuntimeError(
                f"motor ROM step did not converge in {max_iterations} iterations"
            )

        theta_mid = 0.5 * (theta0 + theta1)
        omega_mid = 0.5 * (omega0 + omega1)
        q_mid = 0.5 * (q0 + q1)
        R_mid = self.resistance(
            theta_mid, 0.5 * (state.temperature_K + temperature1)
        )
        components = self.torque_components(theta_mid, q_mid)
        hysteresis_loss = 0.0
        next_hstate = state.hysteresis_state
        next_hflux = hflux0
        next_henergy = state.hysteresis_stored_energy_J
        if h_eval is not None:
            components["hysteresis"] = float(h_eval.torque_Nm)
            components["total"] += components["hysteresis"]
            next_hstate = h_eval.state
            next_hflux = _vector(
                h_eval.flux_linkage_Wb,
                self.ports.n_generalized,
                "hysteresis flux_linkage_Wb",
            )
            next_henergy = float(h_eval.stored_energy_J)
            if np.isfinite(h_eval.dissipated_energy_increment_J):
                hysteresis_loss = float(h_eval.dissipated_energy_increment_J / dt)
        next_state = MotorROMState(
            time_s=float(state.time_s) + dt,
            rotor_angle_rad=theta1 % self.period_rad,
            rotor_speed_rad_s=omega1,
            generalized_currents_A=q1,
            temperature_K=temperature1,
            hysteresis_state=next_hstate,
            hysteresis_flux_linkage_Wb=next_hflux,
            hysteresis_stored_energy_J=next_henergy,
        )
        resistive_loss = float(q_mid @ R_mid @ q_mid)
        electrical_power = float(q_mid @ u)
        mechanical_load_power = load * omega_mid
        stored0 = (
            self.magnetic_energy(theta0, q0)
            + 0.5 * self.inertia_kg_m2 * omega0**2
            - self.cogging_coenergy(theta0)
            + float(state.hysteresis_stored_energy_J)
        )
        stored1 = (
            self.magnetic_energy(theta1, q1)
            + 0.5 * self.inertia_kg_m2 * omega1**2
            - self.cogging_coenergy(theta1)
            + float(next_henergy)
        )
        energy_rate = (stored1 - stored0) / dt
        balance = (
            electrical_power
            - resistive_loss
            - hysteresis_loss
            - mechanical_load_power
            - self.viscous_friction_Nm_s * omega_mid**2
            - energy_rate
        )
        speed_voltage = omega_mid * (
            self.inductance(theta_mid, derivative=1) @ q_mid
            + self.pm_flux(theta_mid, derivative=1)
            + self.motion_flux_gradient(theta_mid)
        )
        flux = (
            self.inductance(theta1) @ q1
            + self.pm_flux(theta1)
            + next_hflux
        )
        output = MotorROMStepOutput(
            phase_currents_A=q1[: self.ports.n_phase].copy(),
            eddy_currents_A=q1[self.ports.n_phase :].copy(),
            phase_flux_linkage_Wb=flux[: self.ports.n_phase].copy(),
            electromagnetic_torque_Nm=components["total"],
            torque_components_Nm=dict(components),
            speed_voltage_V=speed_voltage.copy(),
            resistive_loss_W=resistive_loss,
            hysteresis_loss_W=hysteresis_loss,
            electrical_input_power_W=electrical_power,
            mechanical_load_power_W=mechanical_load_power,
            stored_energy_J=stored1,
            energy_balance_residual_W=float(balance),
            nonlinear_iterations=iterations,
        )
        return next_state, output

    def diagnostics(self) -> dict[str, object]:
        angles = np.linspace(0.0, self.period_rad, 129, endpoint=False)
        lmin = min(float(np.min(np.linalg.eigvalsh(self.inductance(a)))) for a in angles)
        rmin = min(
            float(np.min(np.linalg.eigvalsh(self.resistance(a, self.reference_temperature_K))))
            for a in angles
        )
        return {
            "schema": "radia.motor.angle_periodic_rom.v1",
            "ports": self.ports.diagnostics(),
            "period_rad": self.period_rad,
            "skew_span_rad": self.skew_span_rad,
            "minimum_inductance_eigenvalue_H": lmin,
            "minimum_resistance_eigenvalue_ohm": rmin,
            "passive": bool(lmin > 0.0 and rmin >= -1.0e-11),
            "has_motional_v_cross_b": self.motion_flux_gradient_Wb_per_rad is not None,
            "has_cogging_coenergy": self.cogging_coenergy_J is not None,
            "has_hysteresis": self.hysteresis_port is not None,
            "has_thermal_state": self.thermal_capacity_J_per_K is not None,
            "interpolation": self.inductance_H.diagnostics(),
        }

    def save_npz(self, path) -> str:
        """Write the platform-neutral arrays consumed by C/Simulink/FMI wrappers."""

        path = Path(path)
        np.savez_compressed(
            path,
            schema=np.asarray("radia.motor.angle_periodic_rom.v1"),
            angles_rad=self.inductance_H.angles_rad,
            angle_origin_rad=self.inductance_H.angles_rad[0],
            inductance_H=self.inductance_H.values,
            resistance_ohm=self.resistance_ohm.values,
            pm_flux_linkage_Wb=self.pm_flux_linkage_Wb.values,
            motion_flux_gradient_Wb_per_rad=(
                np.zeros_like(self.pm_flux_linkage_Wb.values)
                if self.motion_flux_gradient_Wb_per_rad is None
                else self.motion_flux_gradient_Wb_per_rad.values
            ),
            cogging_coenergy_J=(
                np.zeros(self.inductance_H.angles_rad.size)
                if self.cogging_coenergy_J is None
                else self.cogging_coenergy_J.values
            ),
            period_rad=self.period_rad,
            skew_span_rad=self.skew_span_rad,
            inertia_kg_m2=self.inertia_kg_m2,
            viscous_friction_Nm_s=self.viscous_friction_Nm_s,
            reference_temperature_K=self.reference_temperature_K,
            thermal_capacity_J_per_K=(
                np.nan
                if self.thermal_capacity_J_per_K is None
                else self.thermal_capacity_J_per_K
            ),
            thermal_conductance_W_per_K=self.thermal_conductance_W_per_K,
            resistance_temperature_coefficient_per_K=(
                self.resistance_temperature_coefficient_per_K
            ),
            end_winding_inductance_H=self.end_winding_inductance_H,
            end_winding_resistance_ohm=self.end_winding_resistance_ohm,
            phase_names=np.asarray(self.ports.phase_names),
            eddy_names=np.asarray(self.ports.eddy_names),
            has_motional_v_cross_b=(
                self.motion_flux_gradient_Wb_per_rad is not None
            ),
            has_cogging_coenergy=(self.cogging_coenergy_J is not None),
            external_hysteresis_required=(self.hysteresis_port is not None),
        )
        return str(path)


def LoadAnglePeriodicMotorROM(path, *, hysteresis_port=None) -> AnglePeriodicMotorROM:
    """Reconstruct a motor ROM from :meth:`AnglePeriodicMotorROM.save_npz`.

    Hysteresis history is intentionally not serialized.  A hysteretic bundle
    requires an explicit functional hysteresis port at load time; omitting it
    is an error rather than a silent linear fallback.
    """

    with np.load(Path(path), allow_pickle=False) as payload:
        schema = str(np.asarray(payload["schema"]).item())
        if schema != "radia.motor.angle_periodic_rom.v1":
            raise ValueError(f"unsupported motor ROM schema: {schema}")
        external_hysteresis_required = bool(
            np.asarray(payload["external_hysteresis_required"]).item()
        )
        if external_hysteresis_required and hysteresis_port is None:
            raise ValueError(
                "motor ROM bundle requires an external functional hysteresis port"
            )
        angles = np.asarray(payload["angles_rad"], dtype=float)
        period = float(np.asarray(payload["period_rad"]).item())
        motion = (
            PeriodicAngleTable(
                angles,
                np.asarray(payload["motion_flux_gradient_Wb_per_rad"], dtype=float),
                period,
            )
            if bool(np.asarray(payload["has_motional_v_cross_b"]).item())
            else None
        )
        cogging = (
            PeriodicAngleTable(
                angles,
                np.asarray(payload["cogging_coenergy_J"], dtype=float),
                period,
            )
            if bool(np.asarray(payload["has_cogging_coenergy"]).item())
            else None
        )
        thermal_capacity = float(
            np.asarray(payload["thermal_capacity_J_per_K"]).item()
        )
        if np.isnan(thermal_capacity):
            thermal_capacity = None
        ports = MotorPortContract(
            tuple(str(value) for value in payload["phase_names"].tolist()),
            tuple(str(value) for value in payload["eddy_names"].tolist()),
        )
        return AnglePeriodicMotorROM(
            ports,
            PeriodicAngleTable(
                angles, np.asarray(payload["inductance_H"], dtype=float), period
            ),
            PeriodicAngleTable(
                angles, np.asarray(payload["resistance_ohm"], dtype=float), period
            ),
            PeriodicAngleTable(
                angles,
                np.asarray(payload["pm_flux_linkage_Wb"], dtype=float),
                period,
            ),
            inertia_kg_m2=float(np.asarray(payload["inertia_kg_m2"]).item()),
            viscous_friction_Nm_s=float(
                np.asarray(payload["viscous_friction_Nm_s"]).item()
            ),
            motion_flux_gradient_Wb_per_rad=motion,
            cogging_coenergy_J=cogging,
            skew_span_rad=float(np.asarray(payload["skew_span_rad"]).item()),
            end_winding_inductance_H=np.asarray(
                payload["end_winding_inductance_H"], dtype=float
            ),
            end_winding_resistance_ohm=np.asarray(
                payload["end_winding_resistance_ohm"], dtype=float
            ),
            reference_temperature_K=float(
                np.asarray(payload["reference_temperature_K"]).item()
            ),
            resistance_temperature_coefficient_per_K=np.asarray(
                payload["resistance_temperature_coefficient_per_K"], dtype=float
            ),
            thermal_capacity_J_per_K=thermal_capacity,
            thermal_conductance_W_per_K=float(
                np.asarray(payload["thermal_conductance_W_per_K"]).item()
            ),
            hysteresis_port=hysteresis_port,
        )


def ProjectRigidRotationMotionalFluxGradient(
    current_basis,
    magnetic_flux_density_T,
    *,
    axis=(0.0, 0.0, 1.0),
    origin=(0.0, 0.0, 0.0),
) -> np.ndarray:
    """Project ``v_per_omega x B`` onto sampled divergence-free current modes.

    The returned generalized voltage coefficient has units Wb/rad.  At angular
    speed ``omega`` the speed-voltage term is ``omega * result`` and its paired
    Lorentz torque is ``q @ result``.  Using the same vector in both equations
    is the discrete electromechanical power identity.
    """

    points = _finite_real(current_basis.points, "current_basis.points")
    weights = _finite_real(current_basis.weights, "current_basis.weights").reshape(-1)
    modes = _finite_real(current_basis.modes, "current_basis.modes")
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("current basis points must have shape (n, 3)")
    if modes.ndim != 3 or modes.shape[1:] != points.shape:
        raise ValueError("current basis modes must have shape (n_modes, n, 3)")
    if weights.shape[0] != points.shape[0]:
        raise ValueError("current basis weights must match points")
    axis_v = _vector(axis, 3, "axis")
    norm = float(np.linalg.norm(axis_v))
    if norm <= 0.0:
        raise ValueError("axis must be non-zero")
    axis_v /= norm
    origin_v = _vector(origin, 3, "origin")
    velocity_per_omega = np.cross(axis_v[None, :], points - origin_v[None, :])
    B = magnetic_flux_density_T(points) if callable(magnetic_flux_density_T) else magnetic_flux_density_T
    B = _finite_real(B, "magnetic_flux_density_T")
    if B.ndim == 1:
        B = np.broadcast_to(_vector(B, 3, "magnetic_flux_density_T"), points.shape)
    if B.shape != points.shape:
        raise ValueError("magnetic_flux_density_T must have shape (n, 3) or (3,)")
    vxB = np.cross(velocity_per_omega, B)
    return np.einsum("n,mnk,nk->m", weights, modes, vxB)


def MotorROMFromHybridVIMSweep(
    angles_rad,
    hybrid_systems,
    ports: MotorPortContract,
    *,
    phase_inductance_H,
    phase_resistance_ohm,
    phase_eddy_mutual_H,
    pm_flux_linkage_Wb,
    inertia_kg_m2: float,
    time_domain_eddy_resistance_ohm=None,
    time_domain_eddy_inductance_H=None,
    motion_flux_gradient_Wb_per_rad=None,
    period_rad: float = _TWO_PI,
    **options,
) -> AnglePeriodicMotorROM:
    """Build a motor ROM from an angle sweep of ``HybridVIMSystem`` objects.

    A nonzero ``surface_mass`` means the frequency-domain system still contains
    a DtN/SIBC ``sqrt(s)`` term.  Such a term must first be realized as a
    positive-real CLN state system and supplied through the explicit
    ``time_domain_eddy_*`` arrays.  It is never silently replaced by a constant
    resistance.
    """

    systems = tuple(hybrid_systems)
    angles = _finite_real(angles_rad, "angles_rad").reshape(-1)
    if len(systems) != angles.size:
        raise ValueError("hybrid_systems must match angles_rad")
    if not systems:
        raise ValueError("hybrid_systems must not be empty")
    ne = systems[0].n_modes
    if ports.n_eddy != ne:
        raise ValueError("ports.eddy_names must match HybridVIMSystem mode count")
    if any(system.n_modes != ne for system in systems):
        raise ValueError("all HybridVIMSystem objects must have the same mode count")
    if any(system.blocks != systems[0].blocks for system in systems):
        raise ValueError("hybrid block ordering must be angle invariant")
    has_sibc = any(np.linalg.norm(system.surface_mass) > 1.0e-14 for system in systems)
    if has_sibc and (
        time_domain_eddy_resistance_ohm is None
        or time_domain_eddy_inductance_H is None
    ):
        raise ValueError(
            "frequency-domain DtN/SIBC surface_mass is active; provide its positive-real "
            "time-domain CLN realization through time_domain_eddy_resistance_ohm and "
            "time_domain_eddy_inductance_H"
        )
    np_ = ports.n_phase
    n = ports.n_generalized

    def angle_array(value, shape, name):
        arr = _finite_real(value, name)
        if arr.shape == shape:
            arr = np.broadcast_to(arr, (angles.size,) + shape).copy()
        if arr.shape != (angles.size,) + shape:
            raise ValueError(f"{name} must have shape {shape} or {(angles.size,) + shape}")
        return arr

    Lpp = angle_array(phase_inductance_H, (np_, np_), "phase_inductance_H")
    Rpp = angle_array(phase_resistance_ohm, (np_, np_), "phase_resistance_ohm")
    Mpe = angle_array(phase_eddy_mutual_H, (np_, ne), "phase_eddy_mutual_H")
    if time_domain_eddy_inductance_H is None:
        Lee = np.stack([system.inductance for system in systems])
    else:
        Lee = angle_array(
            time_domain_eddy_inductance_H,
            (ne, ne),
            "time_domain_eddy_inductance_H",
        )
    if time_domain_eddy_resistance_ohm is None:
        Ree = np.stack([system.resistance for system in systems])
    else:
        Ree = angle_array(
            time_domain_eddy_resistance_ohm,
            (ne, ne),
            "time_domain_eddy_resistance_ohm",
        )
    L = np.zeros((angles.size, n, n))
    R = np.zeros_like(L)
    L[:, :np_, :np_] = Lpp
    L[:, :np_, np_:] = Mpe
    L[:, np_:, :np_] = np.swapaxes(Mpe, 1, 2)
    L[:, np_:, np_:] = Lee
    R[:, :np_, :np_] = Rpp
    R[:, np_:, np_:] = Ree
    psi = angle_array(pm_flux_linkage_Wb, (n,), "pm_flux_linkage_Wb")
    motion_table = None
    if motion_flux_gradient_Wb_per_rad is not None:
        motion = angle_array(
            motion_flux_gradient_Wb_per_rad,
            (n,),
            "motion_flux_gradient_Wb_per_rad",
        )
        motion_table = PeriodicAngleTable(angles, motion, period_rad)
    return AnglePeriodicMotorROM(
        ports,
        PeriodicAngleTable(angles, L, period_rad),
        PeriodicAngleTable(angles, R, period_rad),
        PeriodicAngleTable(angles, psi, period_rad),
        inertia_kg_m2=inertia_kg_m2,
        motion_flux_gradient_Wb_per_rad=motion_table,
        **options,
    )


__all__ = [
    "AnglePeriodicMotorROM",
    "HDivHysteresisPort",
    "HDivHysteresisState",
    "HysteresisEvaluation",
    "HysteresisPort",
    "LoadAnglePeriodicMotorROM",
    "MotorPortContract",
    "MotorROMFromHybridVIMSweep",
    "MotorROMInput",
    "MotorROMState",
    "MotorROMStepOutput",
    "PeriodicAngleTable",
    "ProjectRigidRotationMotionalFluxGradient",
]
