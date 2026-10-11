"""LTspice and MATLAB exchange adapters for CPU-owned descriptor port models."""
from __future__ import annotations
from pathlib import Path
import json
import math
import numpy as np
import scipy.linalg
from radia.prima_export import DescriptorPortModel, pole_residue, SCHEMA, RECONSTRUCTION_LIMIT, _QUANTITIES

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

