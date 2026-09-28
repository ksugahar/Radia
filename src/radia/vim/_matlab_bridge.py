"""MATLAB exchange contract for reduced HCurl Eddy Bubble/CLN models.

NGSolve owns the mesh and finite-element assembly.  This module exports only
the reduced numeric contract consumed by ``radia.simulink`` in MATLAB.  JSON
is used deliberately: it is available in base MATLAB and keeps the exchange
independent of SciPy and MATLAB's version-specific MAT-file writers.
"""

from __future__ import annotations

from pathlib import Path
import json
from typing import Mapping

import numpy as np

from ._eddy_hybrid import HCurlEddyCLNModel


def _real_finite_array(value, name: str) -> np.ndarray:
    array = np.asarray(value)
    if np.iscomplexobj(array) and np.max(np.abs(np.imag(array)), initial=0.0) > 1.0e-13:
        raise ValueError(f"{name} must be real for the MATLAB state-space contract")
    array = np.asarray(np.real(array), dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains non-finite values")
    return array


def _flat_row_major(value, name: str, shape: tuple[int, ...] | None = None) -> dict[str, object]:
    array = _real_finite_array(value, name)
    if shape is not None and tuple(array.shape) != tuple(shape):
        raise ValueError(f"{name} must have shape {shape}, got {array.shape}")
    return {
        "shape": [int(size) for size in array.shape],
        "values": array.ravel(order="C").tolist(),
    }


def _model_payload(
    model: HCurlEddyCLNModel,
    *,
    force_operator=None,
    sample_time_s: float,
    metadata: Mapping[str, object] | None = None,
) -> dict[str, object]:
    if not isinstance(model, HCurlEddyCLNModel):
        raise TypeError("model must be an HCurlEddyCLNModel")
    if model.has_sibc_termination:
        raise ValueError(
            "SIBC termination must be rationalized before MATLAB state-space export"
        )
    n_state = model.state_order
    n_port = model.port_count
    arrays: dict[str, object] = {
        "resistance": _flat_row_major(model.resistance, "resistance", (n_state, n_state)),
        "inductance": _flat_row_major(model.inductance, "inductance", (n_state, n_state)),
        "surface_mass": _flat_row_major(model.surface_mass, "surface_mass", (n_state, n_state)),
        "port_rhs": _flat_row_major(model.port_rhs, "port_rhs", (n_state, n_port)),
    }
    if force_operator is not None:
        arrays["force_operator"] = _flat_row_major(
            force_operator,
            "force_operator",
            (3, n_state, n_port),
        )
    return {
        "state_order": n_state,
        "port_count": n_port,
        "input_convention": "u=-d(coil_current)/dt",
        "force_convention": "F=0.5*real(sum(K(k,a,b)*c(a)*conj(i(b))))",
        "has_sibc_termination": False,
        "basis_names": list(model.basis_names),
        "blocks": {
            name: [int(start), int(stop)]
            for name, (start, stop) in (model.blocks or {}).items()
        },
        "arrays": arrays,
        "metadata": dict(metadata or {}),
        "sample_time_s": sample_time_s,
    }


def ExportHCurlEddyCLNJSON(
    model: HCurlEddyCLNModel,
    filename: str | Path,
    *,
    force_operator=None,
    metadata: Mapping[str, object] | None = None,
    sample_time_s: float = 1.0e-5,
) -> Path:
    """Export a reduced HCurl-VIM/CLN model for MATLAB and Simulink.

    Parameters
    ----------
    model:
        ``HCurlEddyCLNModel`` produced by the NGSolve/Radia path.
    filename:
        JSON destination.  Matrix values are stored in row-major flattened
        arrays with explicit shapes, so MATLAB can reconstruct them without
        depending on JSON nested-array orientation rules.
    force_operator:
        Optional real array with shape ``(3, n_state, n_port)``.  Its
        ``(k,a,b)`` entry is the volume integral of the cross product of
        current mode ``a`` with the conjugate external magnetic field for
        port ``b``.  MATLAB evaluates the physical time-average force as
        ``0.5*real(sum(K .* c .* conj(i)))``.
    metadata:
        JSON-compatible provenance such as frequency and material values.
    sample_time_s:
        Suggested Simulink sample time for the exported state-space model.

    SIBC surface terms are intentionally rejected.  A frequency-dependent
    SIBC DtN term must first be rationalized into additional states before it
    can be represented by this finite-dimensional state-space exchange.
    """

    sample_time_s = float(sample_time_s)
    if not np.isfinite(sample_time_s) or sample_time_s <= 0.0:
        raise ValueError("sample_time_s must be positive and finite")
    model_payload = _model_payload(
        model,
        force_operator=force_operator,
        sample_time_s=sample_time_s,
        metadata=metadata,
    )
    payload = {
        "schema": "radia.hcurl.eddy_cln.exchange.v1",
        **model_payload,
    }
    destination = Path(filename)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(payload, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return destination


def ExportHCurlEddyCLNFamilyJSON(
    snapshots,
    filename: str | Path,
    *,
    sample_time_s: float = 1.0e-5,
    metadata: Mapping[str, object] | None = None,
    shared_state_basis: bool = True,
) -> Path:
    """Export a height-indexed family of common-basis HCurl CLN snapshots.

    ``snapshots`` is an iterable of mappings with ``height_m``, ``model``,
    and optional ``force_operator`` keys.  All snapshots must have identical
    state and port dimensions because the family is intended for moving
    mechanical coupling without a hidden state-coordinate transfer.
    """

    sample_time_s = float(sample_time_s)
    if not np.isfinite(sample_time_s) or sample_time_s <= 0.0:
        raise ValueError("sample_time_s must be positive and finite")
    if not shared_state_basis:
        raise ValueError(
            "a moving family requires a common reduced state basis or explicit transfer maps"
        )

    normalized = []
    for snapshot in snapshots:
        if not isinstance(snapshot, Mapping):
            raise TypeError("each family snapshot must be a mapping")
        if "height_m" not in snapshot or "model" not in snapshot:
            raise ValueError("each family snapshot needs height_m and model")
        height_m = float(snapshot["height_m"])
        if not np.isfinite(height_m):
            raise ValueError("snapshot height_m must be finite")
        model = snapshot["model"]
        model_payload = _model_payload(
            model,
            force_operator=snapshot.get("force_operator"),
            sample_time_s=sample_time_s,
            metadata=snapshot.get("metadata"),
        )
        normalized.append((height_m, model_payload))
    if not normalized:
        raise ValueError("at least one family snapshot is required")
    normalized.sort(key=lambda item: item[0])
    heights = np.asarray([item[0] for item in normalized], dtype=float)
    if np.any(np.diff(heights) <= 0.0):
        raise ValueError("snapshot heights must be strictly increasing")
    reference = normalized[0][1]
    for _, model_payload in normalized[1:]:
        if (
            model_payload["state_order"] != reference["state_order"]
            or model_payload["port_count"] != reference["port_count"]
        ):
            raise ValueError("all family snapshots must share state and port dimensions")

    payload = {
        "schema": "radia.hcurl.eddy_cln.family.v1",
        "shared_state_basis": True,
        "height_unit": "m",
        "sample_time_s": sample_time_s,
        "state_order": reference["state_order"],
        "port_count": reference["port_count"],
        "interpolation_default": "linear",
        "extrapolation_default": "error",
        "metadata": dict(metadata or {}),
        "snapshots": [
            {"height_m": height_m, **model_payload}
            for height_m, model_payload in normalized
        ],
    }
    destination = Path(filename)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(payload, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return destination


def _foster_payload(model, *, force_operator, sample_time_s, metadata, modes=None):
    from ._eddy_foster import HCurlEddyFosterModel

    if not isinstance(model, HCurlEddyFosterModel):
        raise TypeError("model must be an HCurlEddyFosterModel")
    if model.has_sibc_termination:
        raise ValueError(
            "SIBC termination must be rationalized before MATLAB state-space export"
        )
    n_state = model.state_order
    n_port = model.port_count
    vectors = model.modes if modes is None else np.asarray(modes)
    modal_port = vectors.conj().T @ model.port_rhs
    arrays: dict[str, object] = {
        "modal_port_rhs": _flat_row_major(modal_port, "modal_port_rhs", (n_state, n_port)),
    }
    if force_operator is not None:
        operator = np.asarray(force_operator)
        if operator.shape != (3, n_state, n_port):
            raise ValueError(f"force_operator must have shape (3, {n_state}, {n_port})")
        arrays["modal_force_operator"] = _flat_row_major(
            np.einsum("kab,aj->kjb", operator, vectors),
            "modal_force_operator",
            (3, n_state, n_port),
        )
    return {
        "state_order": n_state,
        "port_count": n_port,
        "input_convention": "u=-d(coil_current)/dt",
        "state_convention": "z_dot=-diag(decay_rates)*z+modal_port_rhs*u, y=modal_port_rhs'*z",
        "force_convention": "F=0.5*real(sum(K(k,j,b)*z(j)*conj(i(b)))) in modal coordinates",
        "has_sibc_termination": False,
        "basis_names": list(model.basis_names),
        "blocks": {
            name: [int(start), int(stop)]
            for name, (start, stop) in (model.blocks or {}).items()
        },
        "arrays": arrays,
        "metadata": dict(metadata or {}),
        "sample_time_s": sample_time_s,
    }


def _write_json(payload, filename) -> Path:
    destination = Path(filename)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(payload, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return destination


def ExportHCurlEddyFosterJSON(
    model,
    filename: str | Path,
    *,
    force_operator=None,
    metadata: Mapping[str, object] | None = None,
    sample_time_s: float = 1.0e-5,
) -> Path:
    """Export a reduced HCurl eddy-current model in Foster modal form for MATLAB.

    The payload holds the decay rates ``lambda_k`` (poles ``-lambda_k``), the
    modal input matrix ``V^T P`` and, optionally, the force operator contracted
    onto the modes.  MATLAB integrates the diagonal system
    ``z' = -diag(lambda) z + (V^T P) u`` for ``u = -di/dt`` and evaluates the
    physical time-average force ``0.5*real(sum(K .* z .* conj(i)))``.  SIBC
    surface terms are rejected, as for the reduced exchange.
    """

    sample_time_s = float(sample_time_s)
    if not np.isfinite(sample_time_s) or sample_time_s <= 0.0:
        raise ValueError("sample_time_s must be positive and finite")
    body = _foster_payload(model, force_operator=force_operator,
                           sample_time_s=sample_time_s, metadata=metadata)
    payload = {
        "schema": "radia.hcurl.eddy_foster.exchange.v1",
        "decay_rates": _flat_row_major(model.decay_rates, "decay_rates", (model.state_order,)),
        **body,
    }
    return _write_json(payload, filename)


def ExportHCurlEddyFosterFamilyJSON(
    snapshots,
    filename: str | Path,
    *,
    sample_time_s: float = 1.0e-5,
    metadata: Mapping[str, object] | None = None,
    shared_operator_rtol: float = 1.0e-12,
) -> Path:
    """Export a height-indexed family of Foster models with one shared mode set.

    ``snapshots`` is an iterable of mappings with ``height_m``, ``model`` and
    optional ``force_operator``.  Every snapshot must share the reduced
    ``R``/``L`` pair (a fixed conductor mesh; only the source ports and force
    operators move), so one set of decay rates and modes serves all heights and
    the modal input and force operators interpolate linearly in height.  A
    family whose ``R`` or ``L`` changes with height has no common Foster modes
    and is rejected.
    """

    sample_time_s = float(sample_time_s)
    if not np.isfinite(sample_time_s) or sample_time_s <= 0.0:
        raise ValueError("sample_time_s must be positive and finite")
    items = []
    for snapshot in snapshots:
        if not isinstance(snapshot, Mapping):
            raise TypeError("each family snapshot must be a mapping")
        if "height_m" not in snapshot or "model" not in snapshot:
            raise ValueError("each family snapshot needs height_m and model")
        height_m = float(snapshot["height_m"])
        if not np.isfinite(height_m):
            raise ValueError("snapshot height_m must be finite")
        items.append((height_m, snapshot))
    if not items:
        raise ValueError("at least one family snapshot is required")
    items.sort(key=lambda item: item[0])
    heights = np.asarray([item[0] for item in items], dtype=float)
    if np.any(np.diff(heights) <= 0.0):
        raise ValueError("snapshot heights must be strictly increasing")
    reference = items[0][1]["model"]
    for _, snapshot in items[1:]:
        model = snapshot["model"]
        if model.state_order != reference.state_order or model.port_count != reference.port_count:
            raise ValueError("all family snapshots must share state and port dimensions")
        for name in ("resistance", "inductance", "surface_mass"):
            ref = getattr(reference, name)
            cur = getattr(model, name)
            scale = max(float(np.linalg.norm(ref)), np.finfo(float).tiny)
            if float(np.linalg.norm(cur - ref)) > shared_operator_rtol * scale:
                raise ValueError(
                    f"family snapshots differ in {name}; a moving Foster family needs "
                    "one shared R/L pair (fixed conductor mesh)")
    modes = reference.modes
    payload = {
        "schema": "radia.hcurl.eddy_foster.family.v1",
        "shared_modes": True,
        "height_unit": "m",
        "sample_time_s": sample_time_s,
        "state_order": reference.state_order,
        "port_count": reference.port_count,
        "decay_rates": _flat_row_major(reference.decay_rates, "decay_rates",
                                       (reference.state_order,)),
        "interpolation_default": "linear",
        "extrapolation_default": "error",
        "metadata": dict(metadata or {}),
        "snapshots": [
            {"height_m": height_m,
             **_foster_payload(snapshot["model"], force_operator=snapshot.get("force_operator"),
                               sample_time_s=sample_time_s, metadata=snapshot.get("metadata"),
                               modes=modes)}
            for height_m, snapshot in items
        ],
    }
    return _write_json(payload, filename)


__all__ = [
    "ExportHCurlEddyCLNJSON",
    "ExportHCurlEddyCLNFamilyJSON",
    "ExportHCurlEddyFosterJSON",
    "ExportHCurlEddyFosterFamilyJSON",
]
