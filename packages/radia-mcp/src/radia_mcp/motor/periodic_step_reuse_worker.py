"""Owned worker: run Radia's PeriodicStepReuseSolver on the qualified synthetic ring.

Reads one JSON request on stdin and prints one JSON result on stdout.  The
numerical method is Radia's ``radia.periodic_step_reuse.PeriodicStepReuseSolver``;
this worker only builds the qualified synthetic case, observes setup reuse,
recomputes every step's true residual against its own effective right-hand side,
and writes a replayable artifact.

Scope: the periodic scalar ring with a pitch rotation (optionally signed) that
Radia's surrogate tests qualify for the ``amg`` and ``sparsecholesky`` routes.
It is a linear-solve check of the reuse route.  It is not a motor validation and
does not exercise AMS.
"""

from __future__ import annotations

import hashlib
import json
import math
import numbers
import os
import platform
import sys
import uuid
from contextlib import contextmanager
from pathlib import Path

SCHEMA = "radia-mcp.motor.periodic-step-reuse.v1"
CASE = "synthetic_periodic_ring"
METHODS = ("amg", "sparsecholesky")
# n <= 3000 keeps Radia's dense Cholesky definiteness certificate.
MAX_DOFS = 3000
MAX_STEPS = 200
RHS_KINDS = ("random", "zero")
_KEYS = {"case", "n", "steps", "dt", "pitch", "signed", "method", "rhs", "seed", "output_dir"}

SCOPE = {
    "case": CASE,
    "claim": "linear-solve check of setup reuse and pitch-mapped state on a synthetic ring",
    "motor_validation": False,
    "ams": False,
}


class UnsupportedRuntime(RuntimeError):
    """The installed Radia runtime does not provide the required entry."""


def _integer(request: dict, key: str, low: int, high: int) -> int:
    value = request.get(key)
    if isinstance(value, bool) or not isinstance(value, numbers.Integral) or not low <= value <= high:
        raise ValueError(f"{key} must be an integer in [{low}, {high}], got {value!r}")
    return int(value)


def parse_request(request) -> dict:
    """Validate the request; reject unknown keys and anything outside the qualified case."""
    if not isinstance(request, dict):
        raise ValueError("request must be a JSON object")
    unknown = sorted(set(request) - _KEYS)
    if unknown:
        raise ValueError(f"unknown request keys: {unknown}")
    if request.get("case") != CASE:
        raise ValueError(f"case must be {CASE!r}, got {request.get('case')!r}")
    n = _integer(request, "n", 4, MAX_DOFS)
    steps = _integer(request, "steps", 1, MAX_STEPS)
    pitch = _integer(request, "pitch", 1, n - 1)
    seed = _integer(request, "seed", 0, 2**32 - 1)
    dt = request.get("dt")
    if isinstance(dt, bool) or not isinstance(dt, numbers.Real) or not math.isfinite(dt) or dt <= 0.0:
        raise ValueError(f"dt must be a finite positive number, got {dt!r}")
    signed = request.get("signed")
    if not isinstance(signed, bool):
        raise ValueError(f"signed must be true or false, got {signed!r}")
    method = request.get("method")
    if method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}, got {method!r}")
    rhs = request.get("rhs")
    if rhs not in RHS_KINDS:
        raise ValueError(f"rhs must be one of {RHS_KINDS}, got {rhs!r}")
    output_dir = request.get("output_dir")
    if not isinstance(output_dir, str) or not os.path.isabs(output_dir):
        raise ValueError("output_dir must be an absolute path")
    if not Path(output_dir).is_dir():
        raise ValueError(f"output_dir does not exist: {output_dir}")
    return {"case": CASE, "n": n, "steps": steps, "dt": float(dt), "pitch": pitch,
            "signed": signed, "method": method, "rhs": rhs, "seed": seed,
            "output_dir": output_dir}


def ring_system(n: int):
    """The qualified SPD periodic ring K and lumped conductivity mass M."""
    import numpy as np
    from scipy.sparse import csc_matrix, diags

    off = np.full(n - 1, -1.0)
    K = diags([off, np.full(n, 4.0), off], [-1, 0, 1], shape=(n, n), format="lil")
    K[0, n - 1] = -1.0
    K[n - 1, 0] = -1.0
    M = diags([np.full(n, 0.5)], [0], shape=(n, n), format="csc")
    return csc_matrix(K), M


def pitch_mapping(n: int, pitch: int, signed: bool):
    import numpy as np

    perm = np.roll(np.arange(n), -pitch)
    sign = np.where(np.arange(n) % 2 == 0, 1.0, -1.0) if signed else None
    return perm, sign


@contextmanager
def _observe_setups(module):
    """Count calls into Radia's setup builders while the solver is built and stepped."""
    calls = []
    names = ("_sparsecholesky", "_build_preconditioner")
    originals = {}
    for name in names:
        original = getattr(module, name, None)
        if not callable(original):
            raise UnsupportedRuntime(f"radia.periodic_step_reuse has no {name}")
        originals[name] = original

        def counting(*args, _name=name, _original=original, **kwargs):
            calls.append(_name)
            return _original(*args, **kwargs)

        setattr(module, name, counting)
    try:
        yield calls
    finally:
        for name, original in originals.items():
            setattr(module, name, original)


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_artifact(directory: Path, params: dict, payload: bytes) -> Path:
    """Write one artifact owned by this invocation and expose it only when complete.

    The name carries the request digest plus a random invocation id, so identical
    concurrent requests never share a file.  Bytes go to an exclusively created
    partial file, are flushed to disk, then renamed to the final name; the final
    file is reread and must match the payload before it is reported.
    """
    tag = hashlib.sha256(json.dumps(params, sort_keys=True).encode("utf-8")).hexdigest()[:12]
    stem = f"periodic_step_reuse_{params['method']}_{tag}_{uuid.uuid4().hex}"
    partial = directory / f"{stem}.json.partial"
    final = directory / f"{stem}.json"
    with open(partial, "xb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    try:
        os.rename(partial, final)
    except OSError:
        partial.unlink(missing_ok=True)
        raise
    if final.read_bytes() != payload:
        raise RuntimeError(f"artifact {final} does not hold the bytes that were written")
    return final


def _provenance(module) -> dict:
    import numpy
    import scipy

    import radia

    try:
        import ngsolve
        ngsolve_version = getattr(ngsolve, "__version__", "unknown")
    except ImportError:
        ngsolve_version = None
    source = Path(module.__file__)
    return {
        "radia_version": getattr(radia, "__version__", "unknown"),
        "entry": "radia.periodic_step_reuse.PeriodicStepReuseSolver",
        "entry_source_sha256": _sha256_file(source),
        "ngsolve_version": ngsolve_version,
        "numpy_version": numpy.__version__,
        "scipy_version": scipy.__version__,
        "python": platform.python_version(),
        "platform": platform.platform(),
    }


def run(request) -> dict:
    """Run the qualified case.  Raises ValueError, UnsupportedRuntime or RuntimeError."""
    params = parse_request(request)
    try:
        import numpy as np

        import radia.periodic_step_reuse as module
        from radia._residual_gate import RELATIVE_LIMIT, residual_scale
    except ImportError as exc:
        raise UnsupportedRuntime(f"Radia periodic step reuse is unavailable: {exc}") from None
    if not hasattr(module, "PeriodicStepReuseSolver"):
        raise UnsupportedRuntime("radia.periodic_step_reuse has no PeriodicStepReuseSolver")

    n, dt = params["n"], params["dt"]
    K, M = ring_system(n)
    perm, sign = pitch_mapping(n, params["pitch"], params["signed"])
    # Independent operator and mapping used only for the residual recheck.
    A = (K + M / dt).tocsr()
    mapping_sign = np.ones(n) if sign is None else sign
    rng = np.random.default_rng(params["seed"])

    records = []
    x_prev = np.zeros(n)
    with _observe_setups(module) as setup_calls:
        try:
            solver = module.PeriodicStepReuseSolver(
                K, M, dt=dt, perm=perm, sign=sign, method=params["method"])
        except ImportError as exc:
            raise UnsupportedRuntime(f"solver backend unavailable: {exc}") from None
        builds_after_setup = len(setup_calls)
        for index in range(params["steps"]):
            rhs = rng.standard_normal(n) if params["rhs"] == "random" else np.zeros(n)
            b_eff = rhs + (M @ (mapping_sign * x_prev[perm])) / dt
            x_next = solver.step(rhs, x_prev=x_prev)
            b_norm = float(np.linalg.norm(b_eff))
            relative = float(np.linalg.norm(A @ x_next - b_eff)) / residual_scale(b_norm)
            if not math.isfinite(relative) or relative > RELATIVE_LIMIT:
                raise RuntimeError(
                    f"step {index}: independent true relative residual {relative:.3e} "
                    f"exceeds {RELATIVE_LIMIT:g}")
            records.append({
                "step": index,
                "relative_residual": relative,
                "solver_relative_residual": float(solver.last_relative_residual),
                "b_eff_norm": b_norm,
                "zero_rhs": b_norm == 0.0,
            })
            x_prev = x_next
        builds_total = len(setup_calls)

    observed = {
        "steps": len(records),
        "solver_steps_solved": int(solver.steps_solved),
        "setup_builds_observed": builds_total,
        "setup_builds_during_steps": builds_total - builds_after_setup,
        "setup_count_reported": int(solver.setup_count),
        "reuse_count": len(records) - (builds_total - builds_after_setup),
    }
    if builds_total != 1 or observed["setup_builds_during_steps"] != 0:
        raise RuntimeError(f"setup was not reused: {observed}")
    if observed["solver_steps_solved"] != observed["steps"]:
        raise RuntimeError(f"solver step count disagrees: {observed}")

    residual_max = max(record["relative_residual"] for record in records)
    provenance = _provenance(module)
    artifact_body = {
        "schema": SCHEMA,
        "scope": SCOPE,
        "request": params,
        "residual_limit": RELATIVE_LIMIT,
        "definiteness": solver.definiteness,
        "observed": observed,
        "steps": records,
        "final_state": [float(value) for value in x_prev],
        "provenance": provenance,
    }
    payload = (json.dumps(artifact_body, indent=2, sort_keys=True) + "\n").encode("utf-8")
    path = write_artifact(Path(params["output_dir"]), params, payload)
    return {
        "schema": SCHEMA,
        "status": "ok",
        "scope": SCOPE,
        "request": params,
        "backend": {"method": params["method"], "definiteness_route": solver.definiteness["route"]},
        "residual_limit": RELATIVE_LIMIT,
        "residual_max": residual_max,
        "observed": observed,
        "steps": records,
        "artifact": {"path": str(path), "sha256": _sha256_file(path), "bytes": len(payload)},
        "provenance": provenance,
    }


def main() -> int:
    try:
        request = json.loads(sys.stdin.read())
        result = run(request)
    except UnsupportedRuntime as exc:
        result = {"schema": SCHEMA, "status": "unsupported_runtime", "error": str(exc)}
    except (json.JSONDecodeError, ValueError, TypeError) as exc:
        result = {"schema": SCHEMA, "status": "invalid_input", "error": str(exc)}
    except RuntimeError as exc:
        result = {"schema": SCHEMA, "status": "gate_failed", "error": str(exc)}
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
