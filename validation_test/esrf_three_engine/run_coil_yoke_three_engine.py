"""Three-engine nonlinear validation for coil-driven ESRF Examples 6 and 7.

HDiv-MMM receives its checked iron-only Cubit Q2 mesh.  HCurl reduced-A and
the TOSCA-style mixed total/reduced Omega route receive an independently
checked, conforming physical-air plus Kelvin mesh built from the same iron
STEP.  Every formulation is driven by one mesh-free CoilBuilder source tree.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE_ROOT = HERE.parents[1] / "src"
SOURCE_PACKAGE = SOURCE_ROOT / "radia"
if (SOURCE_PACKAGE / "__init__.py").is_file():
    # A checkout run must exercise the checkout's native extension and Python
    # modules as one unit.  A released compute-node run has no source tree and
    # therefore resolves entirely from the installed wheel; copying a .pyd is
    # never a valid substitute for either route.
    sys.path.insert(0, str(SOURCE_ROOT))

import ngsolve as ng
import numpy as np
import radia as rad

from radia.electromagnet_validation import (
    require_static_electromagnet_three_engine_contract,
)
from radia.esrf_examples import get_esrf_bh_table
from radia.kelvin_identify_ngsolve import detect_kelvin_offset, has_kelvin_identification
from radia.kelvin_solver import MixedOmegaPicardNotConverged

from esrf_coil_yoke import (
    average_observation_field,
    build_radia_coil_source,
    core_selector,
    get_case,
    observation_volume_quadrature,
)


CTYPE_RUNNER_PATH = HERE.parent / "c_type_three_engine" / "run_three_engine.py"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_shared_engines():
    """Reuse the C-type's tested three formulation adapters.

    Geometry/source ownership stays in this runner; the formulation adapters
    remain a single implementation until they are promoted into a public
    validation support module.  Loading is explicit so invoking this script
    never depends on the current working directory.
    """
    if not CTYPE_RUNNER_PATH.is_file():
        raise FileNotFoundError(CTYPE_RUNNER_PATH)
    spec = importlib.util.spec_from_file_location("_radia_c_type_engines", CTYPE_RUNNER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load shared three-engine adapters from {CTYPE_RUNNER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _relative_rms(reference: np.ndarray, candidate: np.ndarray) -> float:
    denominator = float(np.sqrt(np.mean(np.sum(reference * reference, axis=1))))
    if denominator <= 0.0:
        raise RuntimeError("three-engine comparison has a zero reference field")
    return float(
        np.sqrt(np.mean(np.sum((candidate - reference) ** 2, axis=1))) / denominator
    )


CHECKPOINT_SCHEMA = "radia.validation.esrf-coil-yoke-checkpoint.v3"
LEGACY_CHECKPOINT_SCHEMAS = ("radia.validation.esrf-coil-yoke-checkpoint.v2",)
STATE_SCHEMA = "radia.validation.esrf-coil-yoke-picard-state.v1"
# A v2 checkpoint carried the iteration cap inside its identity.  The cap does not
# change a converged solution, so it is provenance, not identity, from v3 on.
LEGACY_CAP_KEY = "nonlinear_maximum_iterations"


def _checkpoint_contract(**values: object) -> dict[str, object]:
    return dict(values)


def _legacy_contract(payload: dict[str, object]) -> dict[str, object]:
    """Strip the iteration cap a v2 checkpoint stored inside its contract."""
    contract = dict(payload.get("contract", {}))
    contract.pop(LEGACY_CAP_KEY, None)
    settings = contract.get("engine_settings")
    if isinstance(settings, dict):
        settings = dict(settings)
        settings.pop(LEGACY_CAP_KEY, None)
        contract["engine_settings"] = settings
    return contract


def _is_converged_result(diagnostics: dict[str, object]) -> bool:
    if not diagnostics.get("nonlinear", True):
        return True
    return bool(dict(diagnostics.get("nonlinear_stats") or {}).get("converged", False))


def _read_checkpoint(path: Path, contract: dict[str, object]):
    """Return a converged result whose identity matches ``contract``, or None.

    A stored contract that differs, an unknown schema, or a non-converged solve
    all raise: a resumed run must never silently continue from a different
    problem or from a partial iterate presented as a result.
    """
    if not path.is_file():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    schema = payload.get("schema")
    if schema == CHECKPOINT_SCHEMA:
        stored = payload.get("contract")
    elif schema in LEGACY_CHECKPOINT_SCHEMAS:
        stored = _legacy_contract(payload)
    else:
        raise RuntimeError(f"unsupported checkpoint schema {schema!r}: remove {path}")
    if stored != contract:
        raise RuntimeError(f"checkpoint contract changed: remove {path}")
    diagnostics = dict(payload["diagnostics"])
    if not _is_converged_result(diagnostics):
        raise RuntimeError(
            f"checkpoint holds a non-converged solve and is not a result: remove {path}")
    return np.asarray(payload["field_T"], dtype=float), diagnostics


def _write_checkpoint(
    path: Path, contract: dict[str, object], field: np.ndarray,
    diagnostics: dict[str, object], provenance: dict[str, object],
) -> None:
    if not _is_converged_result(diagnostics):
        raise RuntimeError(
            "refusing to write a non-converged solve as a result checkpoint")
    path.write_text(
        json.dumps(
            {
                "schema": CHECKPOINT_SCHEMA,
                "contract": contract,
                "provenance": provenance,
                "field_T": np.asarray(field, dtype=float).tolist(),
                "diagnostics": diagnostics,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def _read_state(path: Path, contract: dict[str, object]):
    """Return the persisted partial Picard state for ``contract``, or None."""
    if not path.is_file():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != STATE_SCHEMA or payload.get("contract") != contract:
        raise RuntimeError(f"Picard state contract changed: remove {path}")
    return dict(payload["state"])


def _write_state(
    path: Path, contract: dict[str, object], state: dict[str, object],
    provenance: dict[str, object],
) -> None:
    """Persist a non-converged engine's per-element state, explicitly partial."""
    path.write_text(
        json.dumps(
            {
                "schema": STATE_SCHEMA,
                "contract": contract,
                "provenance": provenance,
                "converged": False,
                "state": state,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def _picard_state(name: str, stats: dict[str, object], material_key: str) -> dict[str, object]:
    return {
        "engine": name,
        "element_numbers": list(stats["element_numbers"]),
        material_key: list(stats[material_key]),
        "iterations": int(stats.get("iterations", 0)),
        "resumed_iterations": int(stats.get("resumed_iterations", 0)),
        "relative_B_change": stats.get("relative_B_change", stats.get("final_relative_change")),
        "contraction_rate_estimate": stats.get("contraction_rate_estimate"),
        "history": list(stats.get("history", ())),
    }


def _pairwise(fields: dict[str, np.ndarray], selector: np.ndarray) -> dict[str, object]:
    rows: dict[str, object] = {}
    names = tuple(fields)
    for index, left in enumerate(names):
        for right in names[index + 1 :]:
            left_values = fields[left][selector]
            right_values = fields[right][selector]
            rows[f"{left}__vs__{right}"] = {
                "relative_rms": _relative_rms(left_values, right_values),
                "maximum_absolute_difference_T": float(
                    np.max(np.linalg.norm(left_values - right_values, axis=1))
                ),
            }
    return rows


def _require_mesh_contract(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    contract = payload.get("contract", {})
    if (
        payload.get("cubit_returncode") != 0
        or not contract.get("finite_outer_air_box_forbidden", False)
        or set(payload.get("materials", ())) != {"iron", "air", "kelvin"}
        or not {"iron_air_interface", "kelvin_int", "kelvin_ext"}
        <= set(payload.get("boundaries", ()))
        or "GND" not in set(payload.get("bbboundaries", ()))
        or int(payload.get("identification_count", 0)) <= 0
    ):
        raise RuntimeError(f"invalid coil-yoke FEM mesh contract: {path}")
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=(6, 7), type=int, required=True)
    parser.add_argument("--assets-dir", type=Path, required=True)
    parser.add_argument("--fem-mesh", type=Path, required=True)
    parser.add_argument("--fem-mesh-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--hdiv-order", choices=(1, 2), type=int, default=2)
    parser.add_argument("--fem-order", type=int, default=2)
    parser.add_argument("--hdiv-gram-eps", type=float, default=1.0e-12)
    parser.add_argument(
        "--hdiv-image",
        default=None,
        help="Explicit HDiv IMA mirror contract, for example --hdiv-image=-x-y for one ESRF quadrupole pole.",
    )
    parser.add_argument("--reduced-a-solver", choices=("direct", "bddc", "ams", "auto"), default="direct")
    parser.add_argument("--reduced-a-relaxation", type=float, default=0.1)
    parser.add_argument("--nonlinear-tolerance", type=float, default=2.0e-5)
    parser.add_argument("--nonlinear-maximum-iterations", type=int, default=80)
    parser.add_argument(
        "--mixed-nonlinear-maximum-iterations",
        type=int,
        default=None,
        help=(
            "Optional Picard cap for mixed total/reduced Omega only. "
            "HDiv and reduced-A retain --nonlinear-maximum-iterations."
        ),
    )
    parser.add_argument(
        "--mixed-relaxation",
        type=float,
        default=0.3,
        help="Damped-Picard relaxation of the mixed total/reduced Omega material update.",
    )
    parser.add_argument(
        "--mixed-anderson-depth",
        type=int,
        default=0,
        help=(
            "Constrained Anderson mixing depth for the mixed total/reduced Omega Picard "
            "loop (0 = the plain damped Picard the earlier runs used).  Opt-in: on a "
            "shielded-knee probe the safeguarded mixing rejected a third of its steps "
            "without beating plain Picard, so it is a per-case choice, and it is part of "
            "the mixed checkpoint identity."
        ),
    )
    parser.add_argument(
        "--reduced-a-anderson-depth",
        type=int,
        default=0,
        help=(
            "Constrained Anderson mixing depth for the reduced-A Picard loop (a smooth "
            "saturating cube converged in 14 instead of 24 iterations at depth 2).  The "
            "default 0 keeps the reduced-A checkpoint identity of earlier runs; a positive "
            "depth enters the identity."
        ),
    )
    parser.add_argument("--source-trace-tolerance", type=float, default=0.05)
    parser.add_argument("--relative-rms-tolerance", type=float, default=0.03)
    parser.add_argument("--observation-half-width", type=float, default=2.0e-5)
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--preflight", action="store_true")
    options = parser.parse_args(argv)
    if options.fem_order < 1:
        raise ValueError("--fem-order must be positive")
    if options.nonlinear_maximum_iterations < 1:
        raise ValueError("--nonlinear-maximum-iterations must be positive")
    if (
        options.mixed_nonlinear_maximum_iterations is not None
        and options.mixed_nonlinear_maximum_iterations < 1
    ):
        raise ValueError("--mixed-nonlinear-maximum-iterations must be positive")
    if not 0.0 < options.hdiv_gram_eps < 1.0:
        raise ValueError("--hdiv-gram-eps must lie in (0, 1)")
    if not 0.0 < options.reduced_a_relaxation <= 1.0:
        raise ValueError("--reduced-a-relaxation must lie in (0, 1]")
    if not 0.0 < options.mixed_relaxation <= 1.0:
        raise ValueError("--mixed-relaxation must lie in (0, 1]")
    if options.mixed_anderson_depth < 0 or options.reduced_a_anderson_depth < 0:
        raise ValueError("Anderson depths must be non-negative")
    if not 0.0 < options.source_trace_tolerance < 1.0:
        raise ValueError("--source-trace-tolerance must lie in (0, 1)")
    if (
        not np.isfinite(options.observation_half_width)
        or options.observation_half_width <= 0.0
    ):
        raise ValueError("--observation-half-width must be finite and positive")
    if options.threads > 0:
        ng.SetNumThreads(options.threads)

    case = get_case(options.case)
    assets_dir = options.assets_dir.resolve()
    iron_mesh_path = assets_dir / "model.vol"
    fem_mesh_path = options.fem_mesh.resolve()
    for path in (iron_mesh_path, fem_mesh_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    fem_mesh_report = _require_mesh_contract(options.fem_mesh_report.resolve())
    iron_mesh = ng.Mesh(str(iron_mesh_path))
    fem_mesh = ng.Mesh(str(fem_mesh_path))
    if not has_kelvin_identification(fem_mesh):
        raise RuntimeError("FEM mesh has no Kelvin point identification")
    points, field_points = observation_volume_quadrature(
        case.number, half_width_m=options.observation_half_width
    )
    for point in field_points:
        if not fem_mesh(*map(float, point)):
            raise RuntimeError(
                f"observation quadrature point lies outside FEM mesh: {point.tolist()}"
            )
    rad.UtiDelAll()
    coil, coil_manifest = build_radia_coil_source(case.number)
    source_field = average_observation_field(
        np.asarray(rad.Fld(coil, "b", field_points), dtype=float), len(points)
    )
    if not np.isfinite(source_field).all():
        raise RuntimeError("CoilBuilder source field is non-finite at an observation point")
    kelvin_center = tuple(float(value) for value in detect_kelvin_offset(fem_mesh))
    engines = _load_shared_engines()
    bh_table = get_esrf_bh_table(case.number)
    mixed_nonlinear_maximum_iterations = (
        int(options.nonlinear_maximum_iterations)
        if options.mixed_nonlinear_maximum_iterations is None
        else int(options.mixed_nonlinear_maximum_iterations)
    )
    output = options.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    common = _checkpoint_contract(
        case=int(case.number),
        iron_mesh_sha256=_sha256(iron_mesh_path),
        fem_mesh_sha256=_sha256(fem_mesh_path),
        fem_mesh_report_sha256=_sha256(options.fem_mesh_report.resolve()),
        hdiv_order=int(options.hdiv_order),
        fem_order=int(options.fem_order),
        hdiv_gram_eps=float(options.hdiv_gram_eps),
        hdiv_image=options.hdiv_image,
        nonlinear_tolerance=float(options.nonlinear_tolerance),
        observation_points_m=points.tolist(),
        observation_half_width_m=float(options.observation_half_width),
        observation_quadrature="tensor_gauss_2x2x2",
    )
    fields: dict[str, np.ndarray] = {}
    diagnostics: dict[str, dict[str, object]] = {}

    def run_hdiv(_state):
        return engines.solve_hdiv(
            iron_mesh,
            coil,
            bh_table,
            nonlinear=True,
            order=options.hdiv_order,
            gram_eps=options.hdiv_gram_eps,
            nonlinear_tolerance=options.nonlinear_tolerance,
            nonlinear_maximum_iterations=options.nonlinear_maximum_iterations,
            points=field_points,
            image=options.hdiv_image,
        )

    def run_reduced_a(state):
        return engines.solve_reduced_a(
            fem_mesh,
            coil,
            bh_table,
            nonlinear=True,
            order=options.fem_order,
            linear_solver=options.reduced_a_solver,
            relax=options.reduced_a_relaxation,
            nonlinear_tolerance=options.nonlinear_tolerance,
            nonlinear_maximum_iterations=options.nonlinear_maximum_iterations,
            nonlinear_verbose=False,
            kelvin_center=kelvin_center,
            kelvin_radius=case.kelvin_radius_m,
            points=field_points,
            anderson_depth=options.reduced_a_anderson_depth,
            nu_initial=None if state is None else np.asarray(state["nu_elements"], dtype=float),
            observation_points=field_points,
        )

    def run_mixed(state):
        return engines.solve_omega(
            fem_mesh,
            coil,
            bh_table,
            nonlinear=True,
            order=options.fem_order,
            nonlinear_tolerance=options.nonlinear_tolerance,
            nonlinear_maximum_iterations=mixed_nonlinear_maximum_iterations,
            nonlinear_verbose=False,
            kelvin_center=kelvin_center,
            kelvin_radius=case.kelvin_radius_m,
            points=field_points,
            source_trace_tolerance=options.source_trace_tolerance,
            relaxation=options.mixed_relaxation,
            anderson_depth=options.mixed_anderson_depth,
            mu_r_initial=(1000.0 if state is None
                          else np.asarray(state["mu_r_elements"], dtype=float)),
            observation_points=field_points,
        )

    # (engine, solve(state), per-element state key or None)
    solver_specs = (
        ("hdiv_mmm", run_hdiv, None),
        ("reduced_a", run_reduced_a, "nu_elements"),
        ("mixed_total_reduced_omega", run_mixed, "mu_r_elements"),
    )
    reduced_a_settings = {
        "linear_solver": options.reduced_a_solver,
        "relaxation": float(options.reduced_a_relaxation),
    }
    if options.reduced_a_anderson_depth:
        reduced_a_settings["anderson_depth"] = int(options.reduced_a_anderson_depth)
    engine_settings = {
        "hdiv_mmm": {
            "linear_tolerance": 1.0e-8,
            "linear_maximum_iterations": 12000,
            "preconditioner": "auto",
        },
        "reduced_a": reduced_a_settings,
        "mixed_total_reduced_omega": {
            "source_potential_contract": "total_hodge",
            "source_trace_tolerance": float(options.source_trace_tolerance),
            "relaxation": float(options.mixed_relaxation),
            "anderson_depth": int(options.mixed_anderson_depth),
        },
    }
    iteration_caps = {
        "hdiv_mmm": int(options.nonlinear_maximum_iterations),
        "reduced_a": int(options.nonlinear_maximum_iterations),
        "mixed_total_reduced_omega": int(mixed_nonlinear_maximum_iterations),
    }

    def provenance(name: str) -> dict[str, object]:
        return {
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "machine": platform.node(),
            "nonlinear_maximum_iterations": iteration_caps[name],
            "engine_settings": engine_settings[name],
        }
    if options.preflight:
        payload = {
            "schema": "radia.validation.esrf-coil-yoke-preflight.v2",
            "passed": True,
            "case": int(case.number),
            "source": coil_manifest,
            "source_field_rms_T": float(np.sqrt(np.mean(np.sum(source_field * source_field, axis=1)))),
            "hdiv_image": options.hdiv_image,
            "iron_mesh_sha256": _sha256(iron_mesh_path),
            "fem_mesh_sha256": _sha256(fem_mesh_path),
            "fem_mesh_contract": fem_mesh_report,
            "kelvin_center_m": list(kelvin_center),
            "observation_points_m": points.tolist(),
        }
        output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        return 0
    for name, solve, state_key in solver_specs:
        checkpoint = output.with_suffix(f".{name}.checkpoint.json")
        state_path = output.with_suffix(f".{name}.state.json")
        contract = _checkpoint_contract(
            **common, engine=name, engine_settings=engine_settings[name]
        )
        resumed = _read_checkpoint(checkpoint, contract) if options.resume else None
        if resumed is not None:
            fields[name], diagnostics[name] = resumed
            diagnostics[name]["resumed_from_checkpoint"] = True
            continue
        # A partial state left by an earlier non-converged run of the SAME
        # problem is a warm start, never a result.
        state = (_read_state(state_path, contract)
                 if options.resume and state_key is not None else None)
        try:
            field_samples, diagnostics[name] = solve(state)
        except MixedOmegaPicardNotConverged as exc:
            if state_key is None:
                raise
            stats = dict(exc.state["nonlinear_stats"])
            if state is not None:
                stats["resumed_iterations"] = (
                    int(state.get("resumed_iterations", 0)) + int(state["iterations"]))
            _write_state(state_path, contract, _picard_state(name, stats, state_key),
                         provenance(name))
            raise RuntimeError(
                f"{name} did not converge; its per-element state was saved to "
                f"{state_path} for a warm restart with --resume"
            ) from exc
        stats = dict(diagnostics[name].get("nonlinear_stats") or {})
        if state is not None:
            stats["resumed_iterations"] = (
                int(state.get("resumed_iterations", 0)) + int(state["iterations"]))
            stats["warm_start_state"] = str(state_path)
            diagnostics[name]["nonlinear_stats"] = stats
        if not _is_converged_result(diagnostics[name]):
            if state_key is not None and state_key in stats:
                _write_state(state_path, contract, _picard_state(name, stats, state_key),
                             provenance(name))
                raise RuntimeError(
                    f"{name} did not converge; its per-element state was saved to "
                    f"{state_path} for a warm restart with --resume"
                )
            raise RuntimeError(f"{name} did not converge")
        fields[name] = average_observation_field(field_samples, len(points))
        diagnostics[name]["observation_contract"] = {
            "kind": "volume_average",
            "quadrature": "tensor_gauss_2x2x2",
            "half_width_m": float(options.observation_half_width),
            "sample_count_per_centre": 8,
        }
        _write_checkpoint(checkpoint, contract, fields[name], diagnostics[name],
                          provenance(name))
        if state_path.is_file():
            state_path.unlink()
    formulation_contract = require_static_electromagnet_three_engine_contract(diagnostics)
    all_points = np.ones(len(points), dtype=bool)
    core = core_selector(case.number, points)
    raw_pairs = _pairwise(fields, all_points)
    core_pairs = _pairwise(fields, core)
    maximum_core_relative_rms = max(
        float(row["relative_rms"]) for row in core_pairs.values()
    )
    nonlinear_converged = all(
        bool(row.get("nonlinear_stats", {}).get("converged", False))
        for row in diagnostics.values()
    )
    passed = nonlinear_converged and maximum_core_relative_rms <= options.relative_rms_tolerance
    result = {
        "schema": "radia.validation.esrf-coil-yoke-three-engine.v2",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "passed": bool(passed),
        "machine": platform.node(),
        "python": sys.version,
        "case": int(case.number),
        "formulation_contract": formulation_contract,
        "shared_input_contract": {
            "iron_cad": "same Cubit STEP authority for HDiv and FEM meshes",
            "coil_source": coil_manifest,
            "bh_table": "shared monotone PCHIP B(H)",
            "finite_outer_air_box_forbidden": True,
            "hdiv_air_mesh_forbidden": True,
            "hdiv_image": options.hdiv_image,
            "fem_kelvin_mesh_shared": True,
            "mixed_h1_source_potential": "total_volume_hodge_with_harmonic_cut",
        },
        "iron_mesh": str(iron_mesh_path),
        "iron_mesh_sha256": _sha256(iron_mesh_path),
        "fem_mesh": str(fem_mesh_path),
        "fem_mesh_sha256": _sha256(fem_mesh_path),
        "fem_mesh_contract": fem_mesh_report,
        "kelvin_center_m": list(kelvin_center),
        "observation_points_m": points.tolist(),
        "observation_contract": {
            "kind": "volume_average",
            "quadrature": "tensor_gauss_2x2x2",
            "half_width_m": float(options.observation_half_width),
            "sample_count_per_centre": 8,
        },
        "core_point_count": int(np.count_nonzero(core)),
        "engines": diagnostics,
        "fields_T": {name: value.tolist() for name, value in fields.items()},
        "pairwise_raw": raw_pairs,
        "pairwise_core": core_pairs,
        "maximum_core_pairwise_relative_rms": maximum_core_relative_rms,
        "relative_rms_tolerance": float(options.relative_rms_tolerance),
        "nonlinear_converged": nonlinear_converged,
        "provenance": {
            "nonlinear_maximum_iterations": iteration_caps,
            "engine_settings": engine_settings,
            "checkpoint_schema": CHECKPOINT_SCHEMA,
        },
    }
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if not passed:
        raise RuntimeError(f"three-engine comparison did not pass; see {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
