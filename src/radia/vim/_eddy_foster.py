"""Foster modal form of a reduced HCurl eddy-current VIM model.

The reduced harmonic contract is

    (R + s L + Zs(s) M_surface) c = P u,

with ``R`` positive semidefinite and ``L`` positive definite.  The generalised
eigenproblem ``R v = lambda L v`` with ``V^T L V = I`` diagonalises the
SIBC-free operator: ``V^T (R + s L) V = diag(lambda) + s I``.  In modal
coordinates ``c = V z`` the system is a Foster sum of first-order sections,

    z_k = (V^T P u)_k / (s + lambda_k),

exactly equal to the direct solve, with real poles ``-lambda_k <= 0``.  The
continuous state-space form for the input ``u = -di/dt`` is therefore
diagonal: ``z' = -diag(lambda) z + (V^T P) u``, ``y = (V^T P)^T z``.

A nonzero SIBC block ``Zs(s) M_surface`` is not rational in ``s``; the model
keeps the direct frequency-domain solve for it and refuses state-space export,
as the reduced contract always has.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._eddy_hybrid import (
    HybridVIMSystem,
    _is_matrix_free_operator,
    _min_hermitian_eigenvalue,
    _mode_names,
    _port_rhs_matrix,
    _solve_reduced_linear,
    _surface_impedance_term,
)


def _dense(matrix) -> np.ndarray:
    if _is_matrix_free_operator(matrix):
        return np.asarray(matrix.to_dense())
    return np.asarray(matrix)


def _hermitian_part(matrix: np.ndarray) -> np.ndarray:
    return 0.5 * (matrix + matrix.conj().T)


@dataclass(frozen=True)
class HCurlEddyFosterModel:
    """Reduced HCurl eddy-current model in Foster (modal) form.

    Construct from reduced matrices; the modal decomposition is computed once.
    ``solve`` returns reduced coordinates ``c`` so existing current and force
    operators apply unchanged.
    """

    resistance: np.ndarray
    inductance: np.ndarray
    surface_mass: np.ndarray
    port_rhs: np.ndarray
    basis_names: tuple[str, ...] = ()
    blocks: dict[str, tuple[int, int]] | None = None
    passive_tol: float = 1.0e-10

    def __post_init__(self) -> None:
        from scipy.linalg import eigh

        resistance = _dense(self.resistance)
        inductance = _dense(self.inductance)
        surface_mass = np.asarray(self.surface_mass)
        n = resistance.shape[0] if resistance.ndim == 2 else -1
        if resistance.ndim != 2 or resistance.shape[0] != resistance.shape[1]:
            raise ValueError("resistance must be square")
        if inductance.shape != resistance.shape:
            raise ValueError("inductance must match resistance")
        if surface_mass.shape != resistance.shape:
            raise ValueError("surface_mass must match resistance")
        for name, value in (("resistance", resistance), ("inductance", inductance),
                            ("surface_mass", surface_mass)):
            if not np.all(np.isfinite(value)):
                raise ValueError(f"{name} must contain only finite values")
        for name, value in (("resistance", resistance), ("inductance", inductance)):
            scale = max(float(np.linalg.norm(value)), np.finfo(float).tiny)
            if np.linalg.norm(value - value.conj().T) > 1.0e-10 * scale:
                raise ValueError(f"{name} must be Hermitian for the Foster form")
        resistance = _hermitian_part(resistance)
        inductance = _hermitian_part(inductance)
        ports = _port_rhs_matrix(self.port_rhs, n)
        names = _mode_names(self.basis_names or None, n)
        blocks = {} if self.blocks is None else dict(self.blocks)
        for name, bounds in blocks.items():
            if len(bounds) != 2:
                raise ValueError(f"block {name!r} must contain (start, stop)")
            start, stop = (int(bounds[0]), int(bounds[1]))
            if start < 0 or stop < start or stop > n:
                raise ValueError(f"block {name!r} is out of range")
            blocks[name] = (start, stop)
        l_min = _min_hermitian_eigenvalue(inductance)
        l_max = float(np.max(np.linalg.eigvalsh(inductance).real))
        if not l_min > self.passive_tol * max(l_max, np.finfo(float).tiny):
            raise ValueError("inductance must be positive definite for the Foster form")
        r_min = _min_hermitian_eigenvalue(resistance)
        r_max = float(np.max(np.linalg.eigvalsh(resistance).real, initial=0.0))
        if r_min < -self.passive_tol * max(r_max, np.finfo(float).tiny):
            raise ValueError("resistance must be positive semidefinite for the Foster form")
        rates, vectors = eigh(resistance, inductance)
        rates = np.clip(np.real(rates), 0.0, None)
        object.__setattr__(self, "resistance", np.array(resistance, copy=True))
        object.__setattr__(self, "inductance", np.array(inductance, copy=True))
        object.__setattr__(self, "surface_mass", np.array(surface_mass, copy=True))
        object.__setattr__(self, "port_rhs", np.array(ports, copy=True))
        object.__setattr__(self, "basis_names", names)
        object.__setattr__(self, "blocks", blocks)
        object.__setattr__(self, "_rates", rates)
        object.__setattr__(self, "_vectors", vectors)
        object.__setattr__(self, "_modal_port", vectors.conj().T @ ports)

    @property
    def state_order(self) -> int:
        return int(self.resistance.shape[0])

    @property
    def port_count(self) -> int:
        return int(self.port_rhs.shape[1])

    @property
    def has_sibc_termination(self) -> bool:
        return bool(np.linalg.norm(self.surface_mass) > 1.0e-14)

    @property
    def decay_rates(self) -> np.ndarray:
        """Foster decay rates ``lambda_k >= 0`` (poles at ``-lambda_k``), ascending."""
        return np.array(self._rates, copy=True)

    @property
    def modes(self) -> np.ndarray:
        """``L``-orthonormal modal vectors ``V`` (columns), ``c = V z``."""
        return np.array(self._vectors, copy=True)

    @property
    def modal_port_rhs(self) -> np.ndarray:
        """Modal input matrix ``V^* P``."""
        return np.array(self._modal_port, copy=True)

    def modal_force_operator(self, force_operator) -> np.ndarray:
        """Contract a reduced-coordinate force operator ``K(k, a, b)`` onto the modes."""
        operator = np.asarray(force_operator)
        if operator.shape != (3, self.state_order, self.port_count):
            raise ValueError(
                f"force_operator must have shape (3, {self.state_order}, {self.port_count})")
        return np.einsum("kab,aj->kjb", operator, self._vectors)

    def impedance(self, s, *, surface_impedance=0.0) -> np.ndarray:
        s = complex(s)
        if not (np.isfinite(s.real) and np.isfinite(s.imag)):
            raise ValueError("s must be finite")
        return (
            self.resistance
            + s * self.inductance
            + _surface_impedance_term(self.surface_mass, surface_impedance)
        )

    def solve(self, s, drive, *, surface_impedance=0.0) -> np.ndarray:
        """Solve for reduced current coefficients ``c`` under generalized drive."""

        values = np.asarray(drive)
        vector_drive = values.ndim <= 1
        if values.ndim == 0:
            if self.port_count != 1:
                raise ValueError(
                    f"scalar drive is valid only for one port, not {self.port_count}")
            values = values.reshape(1, 1)
        elif values.ndim == 1:
            values = values[:, np.newaxis]
        if values.ndim != 2 or values.shape[0] != self.port_count:
            raise ValueError(f"drive must have shape ({self.port_count}, n_cases)")
        s = complex(s)
        if not (np.isfinite(s.real) and np.isfinite(s.imag)):
            raise ValueError("s must be finite")
        sibc = _surface_impedance_term(self.surface_mass, surface_impedance)
        if np.linalg.norm(sibc) > 0.0:
            solved = _solve_reduced_linear(
                self.resistance + s * self.inductance + sibc, self.port_rhs @ values)
        else:
            denominator = s + self._rates
            if np.any(np.abs(denominator) == 0.0):
                raise ValueError("s coincides with a Foster pole")
            modal = (self._modal_port @ values) / denominator[:, np.newaxis]
            solved = self._vectors @ modal
        return solved[:, 0] if vector_drive else solved

    def solve_vector_potential_drive(self, s, coil_current, *, surface_impedance=0.0) -> np.ndarray:
        """Solve ``(R+sL+ZsMs)c = -s P i`` for coil current ``i``."""

        return self.solve(s, -complex(s) * np.asarray(coil_current),
                          surface_impedance=surface_impedance)

    def port_admittance(self, s, *, surface_impedance=0.0) -> np.ndarray:
        response = self.solve(s, np.eye(self.port_count), surface_impedance=surface_impedance)
        return self.port_rhs.conj().T @ response

    def derivative_input_state_space(self) -> dict[str, np.ndarray]:
        """Return the diagonal modal ``z_dot=A z+B u``, ``y=C z`` for ``u=-i_dot``."""

        if self.has_sibc_termination:
            raise ValueError("SIBC termination must be rationalized before state-space export")
        a = -np.diag(self._rates)
        b = self._modal_port
        c = self._modal_port.conj().T
        d = np.zeros((self.port_count, self.port_count), dtype=np.result_type(a, b, c))
        return {"A": a, "B": b, "C": c, "D": d}

    def diagnostics(self) -> dict[str, object]:
        rmin = _min_hermitian_eigenvalue(self.resistance)
        lmin = _min_hermitian_eigenvalue(self.inductance)
        smin = _min_hermitian_eigenvalue(self.surface_mass)
        tol = self.passive_tol
        return {
            "state_order": self.state_order,
            "port_count": self.port_count,
            "form": "foster_modal",
            "blocks": {name: [start, stop] for name, (start, stop) in (self.blocks or {}).items()},
            "input_convention": "u=-d(coil_current)/dt for vector-potential ports",
            "has_sibc_termination": self.has_sibc_termination,
            "finite_rl_state_space": not self.has_sibc_termination,
            "min_decay_rate": float(np.min(self._rates)) if self._rates.size else 0.0,
            "max_decay_rate": float(np.max(self._rates)) if self._rates.size else 0.0,
            "min_resistance_eigenvalue": rmin,
            "min_inductance_eigenvalue": lmin,
            "min_surface_mass_eigenvalue": smin,
            "passive": rmin >= -tol and lmin >= -tol and smin >= -tol,
        }


def HCurlEddyFosterModelFromVIM(system: HybridVIMSystem, rhs) -> HCurlEddyFosterModel:
    """Create the Foster modal model from an HCurl Eddy Bubble VIM system."""

    if not isinstance(system, HybridVIMSystem):
        raise TypeError("system must be a HybridVIMSystem")
    return HCurlEddyFosterModel(
        resistance=system.resistance,
        inductance=system.inductance,
        surface_mass=system.surface_mass,
        port_rhs=rhs,
        basis_names=system.basis_names,
        blocks=system.blocks,
    )


__all__ = ["HCurlEddyFosterModel", "HCurlEddyFosterModelFromVIM"]
