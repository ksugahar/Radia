"""Hand PRIMA-reduced port models to LTspice and Simulink.

A PRIMA reduction is a real descriptor port model

    E x' = A x + B u,    y = C x + D u,

whose input ``u`` is either the port currents (impedance orientation, output
port voltages) or the port voltages (admittance orientation, output port
currents).  This module

* builds that model from :class:`radia.prima_hacapk.PRIMAHACApKModel`
  (:func:`from_prima_hacapk`) or from explicit matrices
  (:class:`DescriptorPortModel`),
* switches orientation exactly through the descriptor inverse,
* writes an LTspice ``.subckt`` whose dynamics are ``Laplace=`` controlled
  sources, one real first- or second-order pole section each
  (:func:`write_ltspice_subckt`),
* writes the JSON exchange that ``radia.simulink.loadPrimaPortModel`` turns
  into a ``dss``/``ss`` object for the Simulink LTI System block
  (:func:`export_prima_lti_json`).

The LTspice sections come from a pole-residue expansion of the descriptor
model.  It is reconstructed against the descriptor response and rejected,
never approximated, when the model has unstable or defective poles or a
polynomial part of degree above one.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import scipy.linalg

SCHEMA = "radia.prima.port_model.v1"
RECONSTRUCTION_LIMIT = 1e-6
_QUANTITIES = {"current": "A", "voltage": "V"}
_ORIENTATION_INPUT = {"impedance": "current", "admittance": "voltage"}


def _real_matrix(value, label: str, shape: tuple[int, int]) -> np.ndarray:
    array = np.asarray(value)
    if np.iscomplexobj(array):
        raise ValueError(f"{label} must be real")
    array = np.array(array, dtype=np.float64)
    if array.shape != shape:
        raise ValueError(f"{label} must have shape {shape}, got {array.shape}")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{label} must be finite")
    return array


@dataclass(frozen=True)
class DescriptorPortModel:
    """Real descriptor port model ``E x' = A x + B u``, ``y = C x + D u``.

    ``input_quantity`` is ``"current"`` (impedance orientation) or
    ``"voltage"`` (admittance orientation); the output is the other
    quantity at the same ports, in SI units.
    """

    E: np.ndarray
    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    D: np.ndarray
    input_quantity: str
    port_names: tuple[str, ...]
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.input_quantity not in _QUANTITIES:
            raise ValueError("input_quantity must be 'current' or 'voltage'")
        names = tuple(str(name) for name in self.port_names)
        if not names or len(set(names)) != len(names) or not all(names):
            raise ValueError("port_names must be unique non-empty names")
        n = np.asarray(self.A).shape[0] if np.ndim(self.A) == 2 else -1
        p = len(names)
        object.__setattr__(self, "port_names", names)
        object.__setattr__(self, "E", _real_matrix(self.E, "E", (n, n)))
        object.__setattr__(self, "A", _real_matrix(self.A, "A", (n, n)))
        object.__setattr__(self, "B", _real_matrix(self.B, "B", (n, p)))
        object.__setattr__(self, "C", _real_matrix(self.C, "C", (p, n)))
        object.__setattr__(self, "D", _real_matrix(self.D, "D", (p, p)))
        object.__setattr__(self, "metadata", dict(self.metadata))

    @property
    def orientation(self) -> str:
        return "impedance" if self.input_quantity == "current" else "admittance"

    @property
    def output_quantity(self) -> str:
        return "voltage" if self.input_quantity == "current" else "current"

    @property
    def state_count(self) -> int:
        return self.A.shape[0]

    @property
    def port_count(self) -> int:
        return len(self.port_names)

    def response(self, s: complex) -> np.ndarray:
        """Port transfer matrix ``C (sE - A)^-1 B + D`` at complex frequency ``s``."""
        pencil = complex(s) * self.E - self.A
        try:
            factor = scipy.linalg.lu_factor(pencil, check_finite=False)
        except scipy.linalg.LinAlgError as exc:
            raise ValueError(f"descriptor pencil is singular at s={s}") from exc
        if np.min(np.abs(np.diag(factor[0]))) == 0.0:
            raise ValueError(f"descriptor pencil is singular at s={s}")
        return self.C @ scipy.linalg.lu_solve(factor, self.B.astype(np.complex128)) + self.D

    def inverse(self) -> "DescriptorPortModel":
        """The same ports in the other orientation (exact descriptor inverse).

        States ``[x; u]``: ``E x' = A x + B u`` and ``0 = C x + D u - y`` with
        the old output ``y`` as input and the old input ``u`` as output.
        """
        n, p = self.state_count, self.port_count
        E = np.zeros((n + p, n + p))
        E[:n, :n] = self.E
        A = np.block([[self.A, self.B], [self.C, self.D]])
        B = np.vstack([np.zeros((n, p)), -np.eye(p)])
        C = np.hstack([np.zeros((p, n)), np.eye(p)])
        return DescriptorPortModel(
            E, A, B, C, np.zeros((p, p)), self.output_quantity, self.port_names,
            self.metadata)

    def oriented(self, orientation: str) -> "DescriptorPortModel":
        if orientation not in _ORIENTATION_INPUT:
            raise ValueError("orientation must be 'impedance' or 'admittance'")
        return self if orientation == self.orientation else self.inverse()


def from_prima_hacapk(model, *, port_name: str = "P1") -> DescriptorPortModel:
    """Descriptor port model of a :class:`radia.prima_hacapk.PRIMAHACApKModel`.

    The descriptor branch stores the MNA projection ``(G_q + s E_q) x = b``
    (in ``R_q``/``L_q``) with port voltage ``b^T x`` per unit port current
    (impedance orientation).  The exact
    series-chain branch ``Z = R + sL`` is returned as its admittance
    ``L i' = -R i + v``.
    """
    metadata = {"source": "radia.prima_hacapk.PRIMAHACApKModel",
                "q": int(model.q), "n_full": int(model.n_full)}
    if model.build_stats.get("is_series"):
        if model.q != 1:
            raise ValueError("series-chain PRIMA model must have q == 1")
        return DescriptorPortModel(
            np.asarray(model.L_q, dtype=np.float64), -np.asarray(model.R_q, dtype=np.float64),
            np.ones((1, 1)), np.ones((1, 1)), np.zeros((1, 1)), "voltage",
            (port_name,), metadata)
    b = np.asarray(model.port_vec_q, dtype=np.float64).reshape(-1)
    return DescriptorPortModel(
        model.L_q, -np.asarray(model.R_q, dtype=np.float64), b[:, None], b[None, :],
        np.zeros((1, 1)), "current", (port_name,), metadata)


@dataclass(frozen=True)
class PoleResidueForm:
    """``H(s) = D0 + s D1 + sum_k R_k / (s - p_k)`` with stable simple poles."""

    poles: np.ndarray       # (m,) complex
    residues: np.ndarray    # (m, p, p) complex
    D0: np.ndarray          # (p, p) real
    D1: np.ndarray          # (p, p) real

    def response(self, s: complex) -> np.ndarray:
        s = complex(s)
        return (self.D0 + s * self.D1
                + np.tensordot(1.0 / (s - self.poles), self.residues, axes=1))

    def real_sections(self) -> list[dict]:
        """Real sections: ``{"order": 1, "num": [r], "den": [a]}`` is
        ``r / (s + a)``; ``{"order": 2, "num": [b1, b0], "den": [a1, a0]}`` is
        ``(b1 s + b0) / (s^2 + a1 s + a0)``.  ``num`` entries are (p, p)."""
        sections = []
        for pole, residue in zip(self.poles, self.residues):
            if pole.imag == 0.0:
                sections.append({"order": 1, "num": [residue.real], "den": [-pole.real]})
            elif pole.imag > 0.0:
                sections.append({
                    "order": 2,
                    "num": [2.0 * residue.real, -2.0 * (residue * np.conj(pole)).real],
                    "den": [-2.0 * pole.real, abs(pole) ** 2],
                })
        return sections


def pole_residue(model: DescriptorPortModel) -> PoleResidueForm:
    """Pole-residue expansion of ``model``, checked against its response."""
    alpha_beta, left, right = scipy.linalg.eig(
        model.A, model.E, left=True, right=True, homogeneous_eigvals=True)
    alpha, beta = alpha_beta
    scale = np.abs(alpha) + np.abs(beta)
    finite = np.abs(beta) > 1e-12 * scale
    poles = alpha[finite] / beta[finite]
    vr, vl = right[:, finite], left[:, finite]
    norm = np.einsum("ik,ij,jk->k", vl.conj(), model.E, vr)
    reference = np.linalg.norm(vl, axis=0) * np.linalg.norm(model.E @ vr, axis=0)
    if np.any(np.abs(norm) <= 1e-10 * np.maximum(reference, np.finfo(float).tiny)):
        raise ValueError("descriptor model has a defective or ill-conditioned pole")
    residues = np.einsum("pk,kq->kpq", model.C @ vr, (vl.conj().T @ model.B)) / norm[:, None, None]

    radius = np.abs(poles)
    unstable = poles.real > 1e-9 * np.maximum(radius, 1.0)
    if np.any(unstable):
        raise ValueError(f"descriptor model has unstable poles: {poles[unstable]}")
    poles = _pair_conjugates(poles)

    proper = PoleResidueForm(poles, residues, np.zeros_like(model.D), np.zeros_like(model.D))
    if np.all(finite):
        D0, D1 = model.D.copy(), np.zeros_like(model.D)
    else:
        D0, D1 = _polynomial_part(model, proper, radius)
    form = PoleResidueForm(poles, residues, D0, D1)
    _check_reconstruction(model, form, radius)
    return form


def _polynomial_part(model, proper, radius):
    """``D0 + s D1`` of ``H - proper part``, sampled where it dominates.

    The samples sit three decades above the geometric-mean pole radius, so
    the proper part has decayed and the round-off of ``H`` (relative to
    ``|H|``) cannot swamp ``D0``.  A third sample rejects higher degrees.
    """
    nonzero = radius[radius > 0.0]
    s_ref = 1e3 * (float(np.exp(np.mean(np.log(nonzero)))) if nonzero.size else 1.0)
    samples = [1j * s_ref * factor for factor in (1.0, 2.0, 3.5)]
    exact = [model.response(s) for s in samples]
    P1, P2, P3 = (h - proper.response(s) for h, s in zip(exact, samples))
    D1 = (P2 - P1) / (samples[1] - samples[0])
    D0 = P1 - samples[0] * D1
    size = max(np.max(np.abs(h)) for h in exact)
    if np.max(np.abs(D0 + samples[2] * D1 - P3)) > 1e-9 * size:
        raise ValueError("descriptor model has a polynomial part of degree > 1")
    if np.max(np.abs(D1)) * abs(samples[2]) <= 1e-9 * size:
        D1 = np.zeros_like(D1)  # proper: the slope is round-off
    return D0.real, D1.real


def _pair_conjugates(poles: np.ndarray) -> np.ndarray:
    """Snap numerically real poles to the real axis; require conjugate pairs."""
    poles = poles.copy()
    tiny = 1e-12 * np.maximum(np.abs(poles), 1.0)
    poles[np.abs(poles.imag) <= tiny] = poles[np.abs(poles.imag) <= tiny].real
    upper = np.sort_complex(poles[poles.imag > 0.0])
    lower = np.sort_complex(np.conj(poles[poles.imag < 0.0]))
    if upper.shape != lower.shape or not np.allclose(upper, lower, rtol=1e-9, atol=0.0):
        raise ValueError("descriptor poles do not form conjugate pairs")
    return poles


def _check_reconstruction(model, form, radius):
    if radius.size:
        lo, hi = max(radius.min(), 1e-12) / 30.0, max(radius.max(), 1e-12) * 30.0
    else:
        lo, hi = 1e-3, 1e3
    worst = 0.0
    for omega in np.geomspace(lo, hi, 61):
        exact = model.response(1j * omega)
        error = np.linalg.norm(form.response(1j * omega) - exact)
        worst = max(worst, error / max(np.linalg.norm(exact), np.finfo(float).tiny))
    if worst > RECONSTRUCTION_LIMIT:
        raise ValueError(
            f"pole-residue reconstruction error {worst:.3e} exceeds {RECONSTRUCTION_LIMIT:g}")


def _num(value: float) -> str:
    return repr(float(value))


def _signed(value: float) -> str:
    """``+x`` or ``-x`` so that expressions never contain ``+-``."""
    value = float(value)
    return ("-" if value < 0.0 else "+") + repr(abs(value))


def ltspice_subckt(model: DescriptorPortModel, *, name: str = "PRIMA") -> str:
    """LTspice ``.subckt`` realizing ``model``'s port relation.

    Pins are ``P<k>_p P<k>_n`` per port, in ``port_names`` order.  A circuit
    element fixes the port relation, not the drive, so the subcircuit works
    under current or voltage excitation either way.  It is realized in the
    orientation whose transfer function is proper, preferring the model's
    own: impedance as series E sources controlled by the sensed port
    currents, admittance as parallel G sources controlled by the port
    voltages.  Each real pole section is one ``Laplace=`` source.  A model
    improper in both orientations is refused.
    """
    if not name.isidentifier():
        raise ValueError("subckt name must be an identifier")
    oriented, form = model, pole_residue(model)
    if np.any(form.D1):
        oriented = model.inverse()
        form = pole_residue(oriented)
        if np.any(form.D1):
            raise ValueError("model is improper in both orientations")
    sections = form.real_sections()
    impedance = oriented.orientation == "impedance"
    ports = [f"P{k + 1}" for k in range(oriented.port_count)]
    lines = [
        f"* Radia PRIMA port model, {oriented.orientation} realization: "
        f"{', '.join(f'{p}={n}' for p, n in zip(ports, oriented.port_names))}",
        f"* states={oriented.state_count} poles={form.poles.size} schema={SCHEMA}",
        f".subckt {name} {' '.join(f'{port}_p {port}_n' for port in ports)}",
    ]
    control = []
    for port in ports:
        if impedance:
            # The sensed port current becomes the voltage of node u_<port>.
            lines.append(f"Vsense_{port} {port}_p {port}_s0 0")
            lines.append(f"Hu_{port} u_{port} 0 Vsense_{port} 1")
            control.append(f"u_{port} 0")
        else:
            control.append(f"{port}_p {port}_n")

    for k, port in enumerate(ports):
        terms = []
        for j in range(len(ports)):
            if form.D0[k, j] != 0.0:
                terms.append((f"D{j + 1}", control[j], _num(form.D0[k, j])))
            for m, section in enumerate(sections):
                num = [coefficient[k, j] for coefficient in section["num"]]
                if not any(num):
                    continue
                if section["order"] == 1:
                    expr = f"{_num(num[0])}/(s{_signed(section['den'][0])})"
                else:
                    a1, a0 = section["den"]
                    expr = (f"({_num(num[0])}*s{_signed(num[1])})"
                            f"/(s*s{_signed(a1)}*s{_signed(a0)})")
                terms.append((f"S{m + 1}_{j + 1}", control[j], f"Laplace={expr}"))
        if not terms:
            raise ValueError(f"port {oriented.port_names[k]} has an identically zero response")
        if impedance:
            node = f"{port}_s0"
            for index, (label, ctrl, value) in enumerate(terms):
                nxt = f"{port}_n" if index == len(terms) - 1 else f"{port}_s{index + 1}"
                lines.append(f"E{label}_{port} {node} {nxt} {ctrl} {value}")
                node = nxt
        else:
            for label, ctrl, value in terms:
                lines.append(f"G{label}_{port} {port}_p {port}_n {ctrl} {value}")
    lines.append(f".ends {name}")
    return "\n".join(lines) + "\n"


def write_ltspice_subckt(model: DescriptorPortModel, path, *, name: str = "PRIMA") -> Path:
    """Write :func:`ltspice_subckt` to ``path`` (ASCII; include with ``.include``)."""
    path = Path(path)
    path.write_text(ltspice_subckt(model, name=name), encoding="ascii", newline="\n")
    return path


def _encode(array: np.ndarray) -> dict:
    array = np.asarray(array, dtype=np.float64)
    return {"shape": list(array.shape), "values": array.reshape(-1).tolist()}


def _check_frequencies(model: DescriptorPortModel) -> np.ndarray:
    alpha, beta = scipy.linalg.eigvals(model.A, model.E, homogeneous_eigvals=True)
    finite = np.abs(beta) > 1e-12 * (np.abs(alpha) + np.abs(beta))
    radius = np.abs(alpha[finite] / beta[finite])
    if radius.size:
        lo, hi = max(radius.min(), 1e-12) / 10.0, max(radius.max(), 1e-12) * 10.0
    else:
        lo, hi = 2.0 * math.pi * 1.0, 2.0 * math.pi * 1e6
    return np.geomspace(lo, hi, 9) / (2.0 * math.pi)


def export_prima_lti_json(model: DescriptorPortModel, path, *,
                          orientation: str | None = None,
                          metadata: dict | None = None) -> Path:
    """Write the exchange read by ``radia.simulink.loadPrimaPortModel``.

    Arrays are row-major ``{shape, values}``.  ``check`` holds the port
    response at a few frequencies; the MATLAB loader refuses the file when
    its ``dss`` object does not reproduce them.
    """
    oriented = model.oriented(orientation or model.orientation)
    frequencies = _check_frequencies(oriented)
    response = np.stack([oriented.response(2j * math.pi * f) for f in frequencies])
    payload = {
        "schema": SCHEMA,
        "orientation": oriented.orientation,
        "input_quantity": oriented.input_quantity,
        "output_quantity": oriented.output_quantity,
        "input_unit": _QUANTITIES[oriented.input_quantity],
        "output_unit": _QUANTITIES[oriented.output_quantity],
        "port_names": list(oriented.port_names),
        "state_count": oriented.state_count,
        "port_count": oriented.port_count,
        "arrays": {key: _encode(getattr(oriented, key)) for key in ("E", "A", "B", "C", "D")},
        "check": {
            "frequency_hz": frequencies.tolist(),
            "response_real": _encode(response.real),
            "response_imag": _encode(response.imag),
            "relative_limit": RECONSTRUCTION_LIMIT,
        },
        "metadata": {**oriented.metadata, **(metadata or {})},
    }
    path = Path(path)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, allow_nan=False)
        handle.write("\n")
    return path


__all__ = [
    "SCHEMA",
    "DescriptorPortModel",
    "PoleResidueForm",
    "export_prima_lti_json",
    "from_prima_hacapk",
    "ltspice_subckt",
    "pole_residue",
    "write_ltspice_subckt",
]
