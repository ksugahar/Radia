"""Backward-Euler transient driver for the reduced HDiv/HCurl system.

The harmonic system stores the reciprocal blocks as

    A_m m - K j / mu = b_m
    R j + L dj/dt + K^H dm/dt = b_e.

This module supplies the time-discrete driver around those existing blocks.
It deliberately accepts only real instantaneous surface resistance in the
time domain.  A complex ESIM/SIBC impedance belongs to a convolution
quadrature driver and is rejected here instead of being silently misused.
"""

from __future__ import annotations

import json
from pathlib import Path
import numpy as np

from ._eddy_hybrid import MU0, SurfaceImpedanceGram
from ._nonlinear import _bh_inverse_funcs


def _dense(value, name: str) -> np.ndarray:
    if hasattr(value, "to_dense"):
        value = value.to_dense()
    out = np.asarray(value)
    if out.ndim != 2 or out.shape[0] != out.shape[1]:
        raise ValueError(f"{name} must be a square matrix")
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} contains non-finite values")
    return np.asarray(out, dtype=np.result_type(out, np.complex128))


def _vector(value, size: int, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=np.result_type(value, np.complex128)).reshape(-1)
    if out.size != size:
        raise ValueError(f"{name} must have {size} entries, got {out.size}")
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} contains non-finite values")
    return out


def _resolve_source(source, step: int, time_s: float, previous: dict[str, np.ndarray], size: int, name: str):
    if source is None:
        return np.zeros(size, dtype=complex)
    if callable(source):
        value = source(step, time_s, previous)
    else:
        value = source
    return _vector(value, size, name)


def _surface_value(surface_impedance, step: int, time_s: float, dt: float):
    if surface_impedance is None:
        return 0.0
    if callable(surface_impedance):
        value = surface_impedance(step, time_s, dt)
    elif isinstance(surface_impedance, (list, tuple)):
        if step >= len(surface_impedance):
            raise ValueError("surface_impedance sequence is shorter than the time grid")
        value = surface_impedance[step]
    else:
        value = surface_impedance
    if isinstance(value, SurfaceImpedanceGram):
        matrix = np.asarray(value.matrix, dtype=complex)
        scale = max(float(np.linalg.norm(matrix)), 1.0)
        if np.linalg.norm(matrix.imag) > 1.0e-12 * scale:
            raise ValueError(
                "solve_hdiv_hcurl_transient requires real instantaneous surface "
                "resistance; use a convolution-quadrature driver for complex ESIM/SIBC"
            )
        return value
    array = np.asarray(value)
    scale = max(float(np.linalg.norm(array)), 1.0)
    if np.iscomplexobj(array) and np.linalg.norm(array.imag) > 1.0e-12 * scale:
        raise ValueError(
            "solve_hdiv_hcurl_transient rejects complex surface impedance; "
            "use a convolution-quadrature driver for ESIM/SIBC"
        )
    real = np.asarray(array.real if np.iscomplexobj(array) else array, dtype=float)
    if real.ndim != 0:
        raise ValueError(
            "surface_impedance arrays must be wrapped in SurfaceImpedanceGram"
        )
    if np.any(real < 0.0):
        raise ValueError("surface resistance must be non-negative")
    return real.item() if real.ndim == 0 else real


def _surface_term(eddy_system, surface_value) -> np.ndarray:
    if surface_value is None or np.all(np.asarray(surface_value) == 0.0):
        return np.zeros_like(_dense(eddy_system.resistance, "resistance"))
    zero = _dense(eddy_system.resistance, "resistance")
    return np.asarray(
        eddy_system.impedance(0.0, surface_impedance=surface_value) - zero,
        dtype=complex,
    )


def _real_matrix(value, name: str) -> np.ndarray:
    out = _dense(value, name)
    scale = max(float(np.linalg.norm(out)), 1.0)
    if np.linalg.norm(out.imag) > 1.0e-12 * scale:
        raise ValueError(f"{name} must be real in the local time-domain driver")
    return np.asarray(out.real, dtype=float)


def _real_array2d(value, name: str) -> np.ndarray:
    if hasattr(value, "to_dense"):
        value = value.to_dense()
    out = np.asarray(value)
    if out.ndim != 2:
        raise ValueError(f"{name} must be a two-dimensional matrix")
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} contains non-finite values")
    scale = max(float(np.linalg.norm(out)), 1.0)
    if np.iscomplexobj(out) and np.linalg.norm(out.imag) > 1.0e-12 * scale:
        raise ValueError(f"{name} must be real in the local time-domain driver")
    return np.asarray(out.real if np.iscomplexobj(out) else out, dtype=float)


def _real_vector(value, size: int, name: str) -> np.ndarray:
    out = _vector(value, size, name)
    scale = max(float(np.linalg.norm(out)), 1.0)
    if np.linalg.norm(out.imag) > 1.0e-12 * scale:
        raise ValueError(f"{name} must be real in the local time-domain driver")
    return np.asarray(out.real, dtype=float)


class _ReducedBHConstitutiveLaw:
    """Project an isotropic single-valued B-H law onto sampled HDiv modes."""

    def __init__(self, basis, bh_curve):
        modes = np.asarray(basis.modes)
        weights = np.asarray(basis.weights, dtype=float).reshape(-1)
        scale = max(float(np.linalg.norm(modes)), 1.0)
        if np.iscomplexobj(modes) and np.linalg.norm(modes.imag) > 1.0e-12 * scale:
            raise ValueError("magnetization basis must be real for nonlinear time stepping")
        self.modes = np.asarray(modes.real, dtype=float)
        self.weights = weights
        curve = np.asarray(bh_curve, dtype=float)
        if curve.ndim != 2 or curve.shape[1] != 2 or curve.shape[0] < 2:
            raise ValueError("bh_curve must have shape (n, 2) with n >= 2")
        if not np.all(np.isfinite(curve)):
            raise ValueError("bh_curve contains non-finite values")
        if curve[0, 0] != 0.0 or curve[0, 1] != 0.0:
            raise ValueError("bh_curve must begin at (0, 0)")
        if np.any(np.diff(curve[:, 0]) <= 0.0) or np.any(np.diff(curve[:, 1]) <= 0.0):
            raise ValueError("bh_curve H and B values must be strictly increasing")
        self.curve = np.array(curve, copy=True)
        self._fields, self._coenergy, self.maximum_magnetization = _bh_inverse_funcs(
            curve[:, 0], curve[:, 1]
        )

    @property
    def n_modes(self) -> int:
        return int(self.modes.shape[0])

    def evaluate(self, coefficients):
        coefficients = np.asarray(coefficients, dtype=float).reshape(-1)
        if coefficients.size != self.n_modes:
            raise ValueError(
                f"magnetization coefficients must have {self.n_modes} entries"
            )
        samples = np.einsum("a,aik->ik", coefficients, self.modes)
        magnitude = np.linalg.norm(samples, axis=1)
        nu_sec, nu_diff = self._fields(np.maximum(magnitude, 1.0e-30))
        field = nu_sec[:, np.newaxis] * samples
        residual = np.einsum(
            "aik,ik,i->a", self.modes, field, self.weights
        )

        direction = np.zeros_like(samples)
        nonzero = magnitude > 1.0e-30
        direction[nonzero] = samples[nonzero] / magnitude[nonzero, np.newaxis]
        tangent = np.einsum(
            "aik,bik,i->ab",
            self.modes,
            self.modes,
            self.weights * nu_sec,
        )
        parallel_modes = np.einsum("aik,ik->ai", self.modes, direction)
        tangent += np.einsum(
            "ai,bi,i->ab",
            parallel_modes,
            parallel_modes,
            self.weights * (nu_diff - nu_sec),
        )
        tangent = 0.5 * (tangent + tangent.T)
        coenergy = float(
            np.dot(self.weights, self._coenergy(np.maximum(magnitude, 0.0)))
        )
        return residual, tangent, coenergy, samples


def _resolve_nonlinear_demag(system, demag_operator, constitutive):
    if demag_operator is None:
        reduction = getattr(system, "hdiv_reduction", None)
        demag_operator = None if reduction is None else reduction.demag
    if demag_operator is None:
        raise ValueError(
            "demag_operator is required unless system.hdiv_reduction stores one"
        )
    demag = _real_matrix(demag_operator, "demag_operator")
    expected = (constitutive.n_modes, constitutive.n_modes)
    if demag.shape != expected:
        raise ValueError(f"demag_operator must have shape {expected}")
    hermitian_error = float(np.linalg.norm(demag - demag.T)) / max(
        float(np.linalg.norm(demag)), 1.0e-300
    )
    if hermitian_error > 1.0e-10:
        raise ValueError("demag_operator must be symmetric")
    return 0.5 * (demag + demag.T)


def solve_hdiv_hcurl_nonlinear_transient(
    system,
    times,
    *,
    bh_curve,
    demag_operator=None,
    magnetic_rhs=None,
    eddy_rhs=None,
    initial_magnetization=None,
    initial_eddy=None,
    mu: float = MU0,
    surface_impedance=None,
    residual_tolerance: float = 1.0e-8,
    nonlinear_max_iterations: int = 40,
    line_search_minimum: float = 2.0**-20,
    energy_balance_absolute_tolerance: float = 1.0e-9,
    energy_balance_relative_tolerance: float = 1.0e-8,
    enforce_energy_balance: bool = True,
):
    """Advance a nonlinear B-H HDiv-MMM/HCurl-VIM system.

    Each backward-Euler step solves the magnetic constitutive residual and the
    eddy-current equation in one Newton system.  The magnetic residual is

    ``int H(M).q dV + N m - K j / mu - b_m``

    with the consistent differential-reluctivity tangent projected on the
    sampled HDiv basis.  This is a bulk nonlinear B-H solve; it is distinct
    from the frequency-domain local-ESIM/SIBC outer iteration.
    """

    if callable(system) or not (
        hasattr(system, "coupling") and hasattr(system, "eddy_system")
    ):
        raise TypeError(
            "system must be one fixed CoupledHDivHybridVIMSystem"
        )
    if not np.isfinite(mu) or mu <= 0.0:
        raise ValueError("mu must be positive")
    tolerance = float(residual_tolerance)
    if not np.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("residual_tolerance must be positive")
    max_iterations = int(nonlinear_max_iterations)
    if max_iterations < 1:
        raise ValueError("nonlinear_max_iterations must be >= 1")
    line_search_minimum = float(line_search_minimum)
    if not 0.0 < line_search_minimum <= 1.0:
        raise ValueError("line_search_minimum must be in (0, 1]")
    energy_abs_tolerance = float(energy_balance_absolute_tolerance)
    energy_rel_tolerance = float(energy_balance_relative_tolerance)
    if not np.isfinite(energy_abs_tolerance) or energy_abs_tolerance <= 0.0:
        raise ValueError("energy_balance_absolute_tolerance must be positive")
    if not np.isfinite(energy_rel_tolerance) or energy_rel_tolerance <= 0.0:
        raise ValueError("energy_balance_relative_tolerance must be positive")

    grid = np.asarray(times, dtype=float).reshape(-1)
    if grid.size < 2 or not np.all(np.isfinite(grid)) or np.any(np.diff(grid) <= 0.0):
        raise ValueError("times must contain at least two finite, strictly increasing values")

    constitutive = _ReducedBHConstitutiveLaw(system.magnetization_basis, bh_curve)
    n_m = constitutive.n_modes
    n_e = int(system.n_hcurl_vim_modes)
    if int(system.n_hdiv_modes) != n_m:
        raise ValueError("magnetization basis and system mode counts do not match")
    demag = _resolve_nonlinear_demag(system, demag_operator, constitutive)
    K = _real_array2d(system.coupling, "coupling")
    if K.shape != (n_m, n_e):
        raise ValueError("system coupling shape does not match its mode counts")
    R = _real_matrix(system.eddy_system.resistance, "eddy resistance")
    L = _real_matrix(system.eddy_system.inductance, "eddy inductance")
    if R.shape != (n_e, n_e) or L.shape != (n_e, n_e):
        raise ValueError("eddy system matrices do not match the HCurl mode count")

    m_prev = (
        np.zeros(n_m)
        if initial_magnetization is None
        else _real_vector(initial_magnetization, n_m, "initial_magnetization")
    )
    j_prev = (
        np.zeros(n_e)
        if initial_eddy is None
        else _real_vector(initial_eddy, n_e, "initial_eddy")
    )
    previous = {"magnetization": m_prev.copy(), "eddy": j_prev.copy()}
    _, _, coenergy_prev, _ = constitutive.evaluate(m_prev)
    magnetic_energy_prev = float(
        mu * (coenergy_prev + 0.5 * np.dot(m_prev, demag @ m_prev))
    )
    eddy_energy_prev = float(0.5 * np.dot(j_prev, L @ j_prev))
    stored_energy_prev = magnetic_energy_prev + eddy_energy_prev
    snapshots = [
        {
            "step": 0,
            "time_s": float(grid[0]),
            "dt_s": 0.0,
            "magnetization_coefficients": m_prev.copy(),
            "eddy_coefficients": j_prev.copy(),
            "nonlinear_iterations": 0,
            "residual_relative_norm": 0.0,
            "joule_loss_w": 0.0,
            "magnetic_energy_j": magnetic_energy_prev,
            "eddy_energy_j": eddy_energy_prev,
            "stored_energy_j": stored_energy_prev,
            "energy_balance_residual_w": 0.0,
            "energy_balance_relative_norm": 0.0,
            "energy_balance_mixed_norm": 0.0,
        }
    ]
    states = []

    for step in range(1, grid.size):
        time_s = float(grid[step])
        dt = float(grid[step] - grid[step - 1])
        surface_value = _surface_value(surface_impedance, step - 1, time_s, dt)
        surface_term = _real_matrix(
            _surface_term(system.eddy_system, surface_value),
            "surface resistance",
        )
        R_eff = R + surface_term
        magnetic_source = (
            magnetic_rhs
            if magnetic_rhs is not None
            else getattr(system, "magnetic_rhs", None)
        )
        eddy_source = (
            eddy_rhs if eddy_rhs is not None else getattr(system, "eddy_rhs", None)
        )
        bm = _real_vector(
            _resolve_source(magnetic_source, step, time_s, previous, n_m, "magnetic_rhs"),
            n_m,
            "magnetic_rhs",
        )
        be = _real_vector(
            _resolve_source(eddy_source, step, time_s, previous, n_e, "eddy_rhs"),
            n_e,
            "eddy_rhs",
        )

        m_now = m_prev.copy()
        j_now = j_prev.copy()
        history = []

        def residual_and_tangent(magnetic, eddy):
            c_residual, c_tangent, coenergy, samples = constitutive.evaluate(magnetic)
            magnetic_residual = c_residual + demag @ magnetic - K @ eddy / mu - bm
            eddy_residual = (
                R_eff @ eddy
                + L @ (eddy - j_prev) / dt
                + K.T @ (magnetic - m_prev) / dt
                - be
            )
            magnetic_scale = max(
                float(np.linalg.norm(c_residual)),
                float(np.linalg.norm(demag @ magnetic)),
                float(np.linalg.norm(K @ eddy / mu)),
                float(np.linalg.norm(bm)),
                1.0e-300,
            )
            eddy_scale = max(
                float(np.linalg.norm(R_eff @ eddy)),
                float(np.linalg.norm(L @ (eddy - j_prev) / dt)),
                float(np.linalg.norm(K.T @ (magnetic - m_prev) / dt)),
                float(np.linalg.norm(be)),
                1.0e-300,
            )
            relative = max(
                float(np.linalg.norm(magnetic_residual)) / magnetic_scale,
                float(np.linalg.norm(eddy_residual)) / eddy_scale,
            )
            return (
                np.concatenate((magnetic_residual, eddy_residual)),
                c_tangent,
                coenergy,
                samples,
                relative,
            )

        converged = False
        for iteration in range(1, max_iterations + 1):
            residual, c_tangent, coenergy_now, samples_now, relative = residual_and_tangent(
                m_now, j_now
            )
            history.append(
                {
                    "iteration": iteration,
                    "residual_relative_norm": relative,
                    "line_search_factor": 0.0,
                }
            )
            if relative <= tolerance:
                converged = True
                break
            jacobian = np.block(
                [
                    [c_tangent + demag, -K / mu],
                    [K.T / dt, R_eff + L / dt],
                ]
            )
            delta = np.linalg.solve(jacobian, -residual)
            factor = 1.0
            accepted = False
            while factor >= line_search_minimum:
                candidate_m = m_now + factor * delta[:n_m]
                candidate_j = j_now + factor * delta[n_m:]
                candidate = residual_and_tangent(candidate_m, candidate_j)
                if candidate[4] < relative:
                    m_now = candidate_m
                    j_now = candidate_j
                    accepted = True
                    break
                factor *= 0.5
            history[-1]["line_search_factor"] = factor if accepted else 0.0
            if not accepted:
                raise RuntimeError(
                    f"nonlinear transient step {step} line search failed at "
                    f"relative residual {relative:.3e}"
                )

        if not converged:
            raise RuntimeError(
                f"nonlinear transient step {step} did not converge after "
                f"{max_iterations} iterations"
            )

        residual, _, coenergy_now, samples_now, relative = residual_and_tangent(
            m_now, j_now
        )
        magnetic_energy = float(
            mu * (coenergy_now + 0.5 * np.dot(m_now, demag @ m_now))
        )
        eddy_energy = float(0.5 * np.dot(j_now, L @ j_now))
        stored_energy = magnetic_energy + eddy_energy
        delta_m = m_now - m_prev
        delta_j = j_now - j_prev
        constitutive_now, _, _, _ = constitutive.evaluate(m_now)
        magnetic_be_dissipation = float(
            mu
            * (
                np.dot(delta_m, constitutive_now)
                - (coenergy_now - coenergy_prev)
                + 0.5 * np.dot(delta_m, demag @ delta_m)
            )
            / dt
        )
        eddy_be_dissipation = float(0.5 * np.dot(delta_j, L @ delta_j) / dt)
        backward_euler_dissipation = magnetic_be_dissipation + eddy_be_dissipation
        joule_loss = float(np.dot(j_now, R_eff @ j_now))
        source_power = float(
            mu * np.dot(delta_m / dt, bm) + np.dot(j_now, be)
        )
        balance = float(
            source_power
            - joule_loss
            - (stored_energy - stored_energy_prev) / dt
            - backward_euler_dissipation
        )
        balance_scale = max(
            abs(source_power),
            abs(joule_loss),
            abs((stored_energy - stored_energy_prev) / dt),
            abs(backward_euler_dissipation),
            1.0e-300,
        )
        balance_relative = abs(balance) / balance_scale
        balance_mixed = abs(balance) / (
            energy_abs_tolerance + energy_rel_tolerance * balance_scale
        )
        if enforce_energy_balance and balance_mixed > 1.0:
            raise RuntimeError(
                f"nonlinear transient step {step} energy balance mixed norm "
                f"{balance_mixed:.3e} exceeds one"
            )
        if magnetic_be_dissipation < -1.0e-10 * max(balance_scale, 1.0):
            raise RuntimeError(
                f"nonlinear transient step {step} has negative magnetic "
                "backward-Euler dissipation"
            )

        state = {
            "step": step,
            "time_s": time_s,
            "dt_s": dt,
            "magnetization_coefficients": m_now.copy(),
            "eddy_coefficients": j_now.copy(),
            "magnetization_samples": samples_now.copy(),
            "nonlinear_iterations": len(history),
            "nonlinear_history": history,
            "residual_relative_norm": relative,
            "joule_loss_w": joule_loss,
            "magnetic_energy_j": magnetic_energy,
            "eddy_energy_j": eddy_energy,
            "stored_energy_j": stored_energy,
            "source_power_w": source_power,
            "magnetic_backward_euler_dissipation_w": magnetic_be_dissipation,
            "eddy_backward_euler_dissipation_w": eddy_be_dissipation,
            "backward_euler_dissipation_w": backward_euler_dissipation,
            "energy_balance_residual_w": balance,
            "energy_balance_relative_norm": balance_relative,
            "energy_balance_scale_w": balance_scale,
            "energy_balance_mixed_norm": balance_mixed,
        }
        states.append(state)
        snapshots.append(state.copy())
        m_prev = m_now
        j_prev = j_now
        coenergy_prev = coenergy_now
        stored_energy_prev = stored_energy
        previous = {"magnetization": m_prev.copy(), "eddy": j_prev.copy()}

    return {
        "schema": "cae-ai-lab.radia-vim.hdiv-hcurl-nonlinear-transient.v1",
        "times_s": grid.copy(),
        "bh_curve": np.array(constitutive.curve, copy=True),
        "states": states,
        "snapshots": snapshots,
        "final_magnetization": m_prev.copy(),
        "final_eddy": j_prev.copy(),
        "n_steps": len(states),
        "n_snapshots": len(snapshots),
        "all_steps_converged": True,
        "max_nonlinear_iterations": max(state["nonlinear_iterations"] for state in states),
        "max_residual_relative_norm": max(state["residual_relative_norm"] for state in states),
        "max_abs_energy_balance_residual_w": max(
            abs(state["energy_balance_residual_w"]) for state in states
        ),
        "max_energy_balance_relative_norm": max(
            state["energy_balance_relative_norm"] for state in states
        ),
        "max_energy_balance_mixed_norm": max(
            state["energy_balance_mixed_norm"] for state in states
        ),
        "all_energy_steps_balanced": all(
            state["energy_balance_mixed_norm"] <= 1.0 for state in states
        ),
        "contract": {
            "time_integrator": "backward_euler",
            "bulk_material": "single-valued isotropic B-H in HDiv magnetization form",
            "nonlinear_solver": "fully coupled Newton with residual line search",
            "magnetic_residual": "int H(M).q dV + N m - K j / mu - b_m",
            "eddy_residual": "R j + L dj/dt + K^T dm/dt - b_e",
            "joule_loss": "j^T R_eff j",
            "stored_energy": "mu*(int Wco(M) dV + m^T N m/2) + j^T L j/2",
            "complex_surface_impedance": "rejected; use convolution quadrature",
            "energy_balance_absolute_tolerance_w": energy_abs_tolerance,
            "energy_balance_relative_tolerance": energy_rel_tolerance,
            "enforce_energy_balance": bool(enforce_energy_balance),
        },
    }


def solve_hdiv_hcurl_transient(
    system,
    times,
    *,
    magnetic_operator=None,
    magnetic_rhs=None,
    eddy_rhs=None,
    initial_magnetization=None,
    initial_eddy=None,
    mu: float = MU0,
    surface_impedance=None,
    solver: str = "dense",
    residual_tolerance: float = 1.0e-10,
    energy_balance_absolute_tolerance: float = 1.0e-10,
    energy_balance_relative_tolerance: float = 1.0e-10,
    enforce_energy_balance: bool = True,
):
    """Advance a reduced HDiv-MMM/HCurl-VIM system with backward Euler.

    ``system`` may be a fixed coupled system or a ``(step, time, previous)``
    provider returning the moved/reassembled coupled system.  This is the
    motion hook used by moving-source validation cases.

    ``magnetic_operator`` may be a matrix or ``(step, time, previous)``
    callback returning the tangent/operator for that step.  If omitted, the
    operator stored on the current system is used, which is the production
    path for a moving/reassembled system.  The callback is a stepping hook,
    not an implicit nonlinear Newton solver: a nonlinear material must provide
    a converged tangent and rhs before each call.
    ``magnetic_rhs`` and ``eddy_rhs`` accept the same callback form.

    ``surface_impedance`` is an instantaneous real resistance (or a sequence/
    callback of such values).  Complex ESIM/SIBC values are rejected because
    they require convolution quadrature and cannot be represented by a local
    backward-Euler matrix.

    By default, a step is rejected when its mixed energy-balance norm exceeds
    one.  Set ``enforce_energy_balance=False`` only for diagnostic probes that
    intentionally use a non-Hermitian or otherwise non-energy-consistent
    operator.
    """

    if not callable(system) and (
        not hasattr(system, "coupling") or not hasattr(system, "eddy_system")
    ):
        raise TypeError("system must be a CoupledHDivHybridVIMSystem")
    if not np.isfinite(mu) or mu <= 0.0:
        raise ValueError("mu must be positive")
    if solver != "dense":
        raise ValueError("the transient reduced driver currently supports solver='dense' only")
    tolerance = float(residual_tolerance)
    if not np.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("residual_tolerance must be positive")
    energy_abs_tolerance = float(energy_balance_absolute_tolerance)
    energy_rel_tolerance = float(energy_balance_relative_tolerance)
    if not np.isfinite(energy_abs_tolerance) or energy_abs_tolerance <= 0.0:
        raise ValueError("energy_balance_absolute_tolerance must be positive")
    if not np.isfinite(energy_rel_tolerance) or energy_rel_tolerance <= 0.0:
        raise ValueError("energy_balance_relative_tolerance must be positive")

    grid = np.asarray(times, dtype=float).reshape(-1)
    if grid.size < 2 or not np.all(np.isfinite(grid)) or np.any(np.diff(grid) <= 0.0):
        raise ValueError("times must contain at least two finite, strictly increasing values")

    system_provider = system if callable(system) else None
    initial_system = (
        system_provider(0, float(grid[0]), {"magnetization": np.zeros(0), "eddy": np.zeros(0)})
        if system_provider is not None
        else system
    )
    if not hasattr(initial_system, "n_hdiv_modes") or not hasattr(initial_system, "n_hcurl_vim_modes"):
        raise TypeError("system must be a CoupledHDivHybridVIMSystem or a system provider")
    n_m = int(initial_system.n_hdiv_modes)
    n_e = int(initial_system.n_hcurl_vim_modes)

    m_prev = np.zeros(n_m, dtype=complex) if initial_magnetization is None else _vector(
        initial_magnetization, n_m, "initial_magnetization"
    )
    j_prev = np.zeros(n_e, dtype=complex) if initial_eddy is None else _vector(
        initial_eddy, n_e, "initial_eddy"
    )
    states = []
    previous = {"magnetization": m_prev.copy(), "eddy": j_prev.copy()}
    initial_inductance = _dense(initial_system.eddy_system.inductance, "eddy inductance")
    initial_operator_source = (
        magnetic_operator
        if magnetic_operator is not None
        else getattr(initial_system, "magnetic_operator", None)
    )
    initial_operator = (
        initial_operator_source(0, float(grid[0]), previous)
        if callable(initial_operator_source)
        else initial_operator_source
    )
    if initial_operator is None:
        raise ValueError(
            "magnetic_operator is required unless the system stores one"
        )
    initial_operator = _dense(initial_operator, "magnetic_operator")
    if initial_operator.shape != (n_m, n_m):
        raise ValueError(f"magnetic_operator must have shape {(n_m, n_m)}")
    previous_energy = float(
        0.5 * mu * np.real(np.vdot(m_prev, initial_operator @ m_prev))
        + 0.5 * np.real(np.vdot(j_prev, initial_inductance @ j_prev))
    )
    initial_snapshot = {
        "step": 0,
        "time_s": float(grid[0]),
        "dt_s": 0.0,
        "magnetization_coefficients": m_prev.copy(),
        "eddy_coefficients": j_prev.copy(),
        "residual_relative_norm": 0.0,
        "joule_loss_w": 0.0,
        "magnetic_energy_j": float(
            0.5 * mu * np.real(np.vdot(m_prev, initial_operator @ m_prev))
        ),
        "eddy_energy_j": float(
            0.5 * np.real(np.vdot(j_prev, initial_inductance @ j_prev))
        ),
        "stored_energy_j": previous_energy,
        "source_power_w": 0.0,
        "backward_euler_dissipation_w": 0.0,
        "operator_motion_work_w": 0.0,
        "magnetic_operator_motion_work_w": 0.0,
        "eddy_operator_motion_work_w": 0.0,
        "energy_balance_residual_w": 0.0,
        "energy_balance_relative_norm": 0.0,
        "energy_balance_scale_w": 0.0,
        "energy_balance_mixed_norm": 0.0,
    }
    snapshots = [initial_snapshot]
    # A nonzero initial state represents an already equilibrated source
    # problem (for example a solver's TIME STEP = 0 output).  Use its
    # operator as the first motion-work reference so the first moved interval
    # is included in the discrete energy identity.
    previous_magnetic_operator = initial_operator.copy()
    previous_eddy_inductance = initial_inductance.copy()

    for step in range(1, grid.size):
        time_s = float(grid[step])
        dt = float(grid[step] - grid[step - 1])
        current_surface = _surface_value(surface_impedance, step - 1, time_s, dt)
        current_system = (
            system_provider(step, time_s, previous)
            if system_provider is not None
            else system
        )
        operator_source = (
            magnetic_operator
            if magnetic_operator is not None
            else getattr(current_system, "magnetic_operator", None)
        )
        A = operator_source(step, time_s, previous) if callable(operator_source) else operator_source
        if A is None:
            raise ValueError(
                "magnetic_operator is required unless the current system stores one"
            )
        A = _dense(A, "magnetic_operator")
        if A.shape != (n_m, n_m):
            raise ValueError(f"magnetic_operator must have shape {(n_m, n_m)}")
        magnetic_source = (
            magnetic_rhs
            if magnetic_rhs is not None
            else getattr(current_system, "magnetic_rhs", None)
        )
        eddy_source = (
            eddy_rhs
            if eddy_rhs is not None
            else getattr(current_system, "eddy_rhs", None)
        )
        bm = _resolve_source(magnetic_source, step, time_s, previous, n_m, "magnetic_rhs")
        be = _resolve_source(eddy_source, step, time_s, previous, n_e, "eddy_rhs")
        if (
            int(current_system.n_hdiv_modes) != n_m
            or int(current_system.n_hcurl_vim_modes) != n_e
        ):
            raise ValueError("moving system provider changed the reduced mode count")
        K = np.asarray(current_system.coupling, dtype=complex)
        if K.shape != (n_m, n_e):
            raise ValueError("system coupling shape does not match its mode counts")
        R = _dense(current_system.eddy_system.resistance, "eddy resistance")
        L = _dense(current_system.eddy_system.inductance, "eddy inductance")
        if R.shape != (n_e, n_e) or L.shape != (n_e, n_e):
            raise ValueError("eddy system matrices do not match the HCurl mode count")
        surface_term = _surface_term(current_system.eddy_system, current_surface)
        Z = R + L / dt + surface_term
        upper = -K / mu
        lower = K.conj().T / dt
        operator = np.block([[A, upper], [lower, Z]])
        rhs_e = be + (L @ j_prev + K.conj().T @ m_prev) / dt
        rhs = np.concatenate([bm, rhs_e])
        solution = np.linalg.solve(operator, rhs)
        m_now = solution[:n_m]
        j_now = solution[n_m:]
        residual = operator @ solution - rhs
        residual_relative = float(np.linalg.norm(residual) / max(np.linalg.norm(rhs), 1.0e-300))
        if not np.isfinite(residual_relative) or residual_relative > tolerance:
            raise RuntimeError(
                f"transient step {step} residual {residual_relative:.3e} exceeds "
                f"tolerance {tolerance:.3e}"
            )

        r_eff = R + surface_term
        joule_loss = float(np.real(np.vdot(j_now, r_eff @ j_now)))
        magnetic_energy = float(0.5 * mu * np.real(np.vdot(m_now, A @ m_now)))
        eddy_energy = float(0.5 * np.real(np.vdot(j_now, L @ j_now)))
        stored_energy = magnetic_energy + eddy_energy
        delta_m = (m_now - m_prev) / dt
        backward_euler_dissipation = float(
            0.5 * mu * np.real(np.vdot(m_now - m_prev, A @ (m_now - m_prev))) / dt
            + 0.5 * np.real(np.vdot(j_now - j_prev, L @ (j_now - j_prev))) / dt
        )
        magnetic_motion_work = 0.0
        eddy_motion_work = 0.0
        if previous_magnetic_operator is not None:
            magnetic_motion_work = float(
                0.5 * mu * np.real(
                    np.vdot(m_prev, (A - previous_magnetic_operator) @ m_prev)
                ) / dt
            )
        if previous_eddy_inductance is not None:
            eddy_motion_work = float(
                0.5 * np.real(
                    np.vdot(j_prev, (L - previous_eddy_inductance) @ j_prev)
                ) / dt
            )
        operator_motion_work = magnetic_motion_work + eddy_motion_work
        source_power = float(
            np.real(mu * np.vdot(delta_m, bm) + np.vdot(j_now, be))
        )
        balance = (
            source_power
            - joule_loss
            - (stored_energy - previous_energy) / dt
            - backward_euler_dissipation
            + operator_motion_work
        )
        balance_scale = max(
            abs(source_power),
            abs(joule_loss),
            abs((stored_energy - previous_energy) / dt),
            abs(backward_euler_dissipation),
            abs(operator_motion_work),
            1.0e-300,
        )
        balance_relative = float(abs(balance) / balance_scale)
        balance_mixed = float(
            abs(balance)
            / (energy_abs_tolerance + energy_rel_tolerance * balance_scale)
        )
        if enforce_energy_balance and balance_mixed > 1.0:
            raise RuntimeError(
                f"transient step {step} energy balance mixed norm "
                f"{balance_mixed:.3e} exceeds one"
            )
        states.append(
            {
                "step": step,
                "time_s": time_s,
                "dt_s": dt,
                "magnetization_coefficients": m_now.copy(),
                "eddy_coefficients": j_now.copy(),
                "residual_relative_norm": residual_relative,
                "joule_loss_w": joule_loss,
                "magnetic_energy_j": magnetic_energy,
                "eddy_energy_j": eddy_energy,
                "stored_energy_j": stored_energy,
                "source_power_w": source_power,
                "backward_euler_dissipation_w": backward_euler_dissipation,
                "operator_motion_work_w": operator_motion_work,
                "magnetic_operator_motion_work_w": magnetic_motion_work,
                "eddy_operator_motion_work_w": eddy_motion_work,
                "energy_balance_residual_w": float(balance),
                "energy_balance_relative_norm": balance_relative,
                "energy_balance_scale_w": balance_scale,
                "energy_balance_mixed_norm": balance_mixed,
            }
        )
        snapshots.append(states[-1].copy())
        previous_energy = stored_energy
        previous_magnetic_operator = A.copy()
        previous_eddy_inductance = L.copy()
        m_prev = m_now
        j_prev = j_now
        previous = {"magnetization": m_prev.copy(), "eddy": j_prev.copy()}

    return {
        "schema": "cae-ai-lab.radia-vim.hdiv-hcurl-transient.v1",
        "times_s": grid.copy(),
        "states": states,
        "snapshots": snapshots,
        "final_magnetization": m_prev.copy(),
        "final_eddy": j_prev.copy(),
        "n_steps": len(states),
        "n_snapshots": len(snapshots),
        "all_steps_converged": True,
        "max_residual_relative_norm": max(state["residual_relative_norm"] for state in states),
        "max_abs_energy_balance_residual_w": max(
            abs(state["energy_balance_residual_w"]) for state in states
        ),
        "max_energy_balance_relative_norm": max(
            state["energy_balance_relative_norm"] for state in states
        ),
        "max_energy_balance_mixed_norm": max(
            state["energy_balance_mixed_norm"] for state in states
        ),
        "all_energy_steps_balanced": all(
            state["energy_balance_mixed_norm"] <= 1.0 for state in states
        ),
        "contract": {
            "time_integrator": "backward_euler",
            "magnetic_operator": "per-step linear/tangent operator",
            "eddy_operator": "R + L/dt + instantaneous surface resistance",
            "complex_surface_impedance": "rejected; use convolution quadrature",
            "joule_loss": "j^H R_eff j",
            "energy_identity": (
                "source_power + operator_motion_work = Joule_loss + "
                "d(stored_energy)/dt + backward_Euler_dissipation; "
                "operator_motion_work = magnetic_operator_motion_work + "
                "eddy_operator_motion_work"
            ),
            "energy_balance_gate": (
                "abs(residual) <= energy_balance_absolute_tolerance + "
                "energy_balance_relative_tolerance * balance_scale"
            ),
            "energy_balance_absolute_tolerance_w": energy_abs_tolerance,
            "energy_balance_relative_tolerance": energy_rel_tolerance,
            "enforce_energy_balance": bool(enforce_energy_balance),
        },
    }


def save_transient_artifact(result: dict, path, *, metadata=None) -> None:
    """Write a JSON-safe transient result with optional case metadata."""

    def encode(value):
        if isinstance(value, np.ndarray):
            return encode(value.tolist())
        if isinstance(value, (np.floating, np.integer, np.bool_)):
            return value.item()
        if isinstance(value, complex):
            return {"real": value.real, "imag": value.imag}
        if isinstance(value, dict):
            return {str(key): encode(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [encode(item) for item in value]
        return value

    record = encode(result)
    if metadata:
        record["metadata"] = encode(metadata)
    Path(path).write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")


__all__ = [
    "solve_hdiv_hcurl_transient",
    "solve_hdiv_hcurl_nonlinear_transient",
    "save_transient_artifact",
]
