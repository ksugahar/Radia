"""Executable quality contracts for the standalone MATLAB Optuna package.

The official MathWorks MATLAB MCP Server remains the execution owner.  These
helpers inspect Radia's distribution evidence, construct commands for that
official server, and validate release evidence without proxying upstream
Optuna MCP tools.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
from collections import Counter
from collections.abc import Mapping
from pathlib import Path
from typing import Any

_NOTICE_TOKENS = (
    "Copyright (c) 2018 Preferred Networks, Inc.",
    "Copyright (c) 2025 Preferred Networks, Inc.",
    (
        "Optuna, the Optuna logo and any related marks are trademarks of "
        "Preferred Networks, Inc."
    ),
    "independent, unofficial project",
)
_TABLE_NAMES = {
    "TrialTable",
    "ObjectiveTable",
    "ParamTable",
    "IntermediateTable",
    "UserAttrTable",
    "ConstraintTable",
    "SamplerStateTable",
}
_ORACLE_SCOPE_FILES = {
    "all": ("test_optuna",),
    "shared": (
        "test_optuna_core_parity",
        "test_optuna_pruner_parity",
        "test_optuna_upstream_oracle",
    ),
    "upstream_mcp": ("test_optuna_mcp_oracle",),
    "matlab_integration": (
        "test_optuna_nsgaii_joint",
        "test_optuna_reliability",
        "test_optuna_sampler_wrappers",
        "test_optuna_simulink_block",
        "test_optuna_storage_bridge",
        "test_optuna_table",
    ),
    "samplers": (
        "test_optuna_nsgaii_joint",
        "test_optuna_sampler_wrappers",
        "test_optuna_upstream_oracle",
    ),
    "storage": (
        "test_optuna_reliability",
        "test_optuna_storage_bridge",
        "test_optuna_table",
        "test_optuna_upstream_oracle",
    ),
    "simulink": ("test_optuna_simulink_block", "test_optuna_table"),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def _sha256_lf(path: Path) -> str:
    """Hash generated text with the repository's platform-neutral LF form."""
    content = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(content).hexdigest().upper()


def _read_json(path: Path, errors: list[str], label: str) -> dict[str, Any]:
    if not path.is_file():
        errors.append(f"missing {label}: {path}")
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        errors.append(f"invalid {label}: {error}")
        return {}
    if not isinstance(value, dict):
        errors.append(f"invalid {label}: expected a JSON object")
        return {}
    return value


def _repository_candidates() -> list[Path]:
    seeds = [Path.cwd(), Path(__file__).resolve()]
    result: list[Path] = []
    for seed in seeds:
        for candidate in (seed, *seed.parents):
            if candidate not in result:
                result.append(candidate)
    return result


def _repository_layout(root: Path) -> dict[str, Path | str] | None:
    matlab_root = root / "matlab"
    manifest = root / "packages" / "radia-optuna" / "src" / "radia_optuna" / "manifest.json"
    if not (matlab_root / "optuna49_api_coverage.json").is_file():
        return None
    if not manifest.is_file():
        return None
    return {
        "source_kind": "repository",
        "root": root,
        "matlab_root": matlab_root,
        "manifest": manifest,
        "notice": root / "packages" / "radia-optuna" / "THIRD_PARTY_NOTICES.md",
        "oracle": root / "tests" / "matlab" / "fixtures" / "optuna49_oracle.json",
        "mcp_oracle": root
        / "tests"
        / "matlab"
        / "fixtures"
        / "optuna49_mcp_oracle.json",
        "inventory": root
        / "tests"
        / "matlab"
        / "fixtures"
        / "optuna49_public_api.json",
        "test_manifest": root
        / "tests"
        / "matlab"
        / "fixtures"
        / "optuna_test_manifest.json",
    }


def _installed_layout(package_root: Path) -> dict[str, Path | str] | None:
    matlab_root = package_root / "matlab"
    manifest = package_root / "manifest.json"
    if not manifest.is_file() or not (matlab_root / "optuna49_api_coverage.json").is_file():
        return None
    return {
        "source_kind": "installed-distribution",
        "root": package_root,
        "matlab_root": matlab_root,
        "manifest": manifest,
        "notice": matlab_root / "THIRD_PARTY_NOTICES.md",
        "oracle": Path(),
        "mcp_oracle": Path(),
        "inventory": Path(),
        "test_manifest": Path(),
    }


def _resolve_layout(location: str = "") -> dict[str, Path | str]:
    if location:
        supplied = Path(location).expanduser().resolve()
        candidates = [supplied]
        if supplied.is_file():
            candidates.insert(0, supplied.parent)
        for candidate in candidates:
            layout = _repository_layout(candidate)
            if layout is not None:
                return layout
            if candidate.name.lower() == "matlab":
                layout = _repository_layout(candidate.parent)
                if layout is not None:
                    return layout
            layout = _installed_layout(candidate)
            if layout is not None:
                return layout
            if candidate.name.lower() == "matlab":
                layout = _installed_layout(candidate.parent)
                if layout is not None:
                    return layout
        raise FileNotFoundError(
            f"{supplied} is neither a Radia repository nor an installed radia-optuna package"
        )

    for candidate in _repository_candidates():
        layout = _repository_layout(candidate)
        if layout is not None:
            return layout

    spec = importlib.util.find_spec("radia_optuna")
    if spec is not None:
        locations = list(spec.submodule_search_locations or [])
        if spec.origin:
            locations.append(str(Path(spec.origin).parent))
        for value in locations:
            layout = _installed_layout(Path(value).resolve())
            if layout is not None:
                return layout
    raise FileNotFoundError(
        "radia-optuna was not found in a repository or installed Python environment"
    )


def _matlab_top_level_inventory(root: Path) -> tuple[list[str], list[str]]:
    files = sorted((root / "+radia" / "+optuna").glob("*.m"))
    classes = sorted(
        path.stem
        for path in files
        if "classdef" in path.read_text(encoding="utf-8", errors="ignore")
    )
    functions = sorted(path.stem for path in files if path.stem not in classes)
    return classes, functions


def _optuna_source_commands(source_path: Path) -> list[str]:
    """Read the standalone Optuna command table from the shared MEX source."""
    if not source_path.is_file():
        return []
    source = source_path.read_text(encoding="utf-8", errors="ignore")
    match = re.search(
        r"mxArray\*\s+Commands\(\)\s*\{\s*"
        r"#ifdef\s+RADIA_OPTUNA_MEX_ONLY\s*"
        r"static\s+const\s+char\*\s+names\[\]\s*=\s*\{(.*?)\};\s*"
        r"#else",
        source,
        flags=re.DOTALL,
    )
    if match is None:
        return []
    return re.findall(r'"([^"]+)"', match.group(1))


def matlab_optuna_health(distribution_path: str = "") -> dict[str, Any]:
    """Inspect manifest, public API oracle coverage, MEX, and notices."""
    errors: list[str] = []
    try:
        layout = _resolve_layout(distribution_path)
    except (OSError, FileNotFoundError) as error:
        return {
            "schema": "radia-mcp.matlab-optuna-health/v1",
            "ok": False,
            "status": "error",
            "errors": [str(error)],
        }

    root = Path(layout["root"])
    matlab_root = Path(layout["matlab_root"])
    manifest_path = Path(layout["manifest"])
    coverage_path = matlab_root / "optuna49_api_coverage.json"
    compatibility_path = matlab_root / "optuna_upstream_compatibility.json"
    notice_path = Path(layout["notice"])
    mex_path = matlab_root / "optuna_mex.mexw64"
    source_path = root / "src" / "matlab" / "radia_mex.cpp"

    manifest = _read_json(manifest_path, errors, "distribution manifest")
    coverage = _read_json(coverage_path, errors, "API coverage")
    compatibility = _read_json(
        compatibility_path, errors, "compatibility contract"
    )

    matlab_files = sorted((matlab_root / "+radia" / "+optuna").rglob("*.m"))
    classes, functions = _matlab_top_level_inventory(matlab_root)
    expected_files = int(manifest.get("matlab_file_count", -1))
    if len(matlab_files) != expected_files:
        errors.append(
            f"MATLAB file inventory mismatch: {len(matlab_files)} != {expected_files}"
        )

    entries = coverage.get("entries", [])
    if not isinstance(entries, list):
        entries = []
        errors.append("API coverage entries must be a list")
    public_api = {
        "entry_count": int(coverage.get("surface_entry_count", -1)),
        "present_count": int(coverage.get("surface_present_count", -1)),
        "missing_count": int(coverage.get("surface_missing_count", -1)),
        "verified_count": int(coverage.get("oracle_verified_count", -1)),
        "asserted_count": int(coverage.get("oracle_asserted_count", -1)),
        "partial_count": int(coverage.get("oracle_partial_count", -1)),
        "unmapped_count": int(coverage.get("oracle_unmapped_count", -1)),
        "required_count": int(coverage.get("required_entry_count", -1)),
        "required_present_count": int(
            coverage.get("required_present_count", -1)
        ),
        "required_verified_count": int(
            coverage.get("required_oracle_mapped_count", -1)
        ),
        "required_asserted_count": int(
            coverage.get("required_oracle_asserted_count", -1)
        ),
        "required_unmapped_count": int(
            coverage.get("required_oracle_unmapped_count", -1)
        ),
        "complete": bool(coverage.get("full_compatibility_complete", False)),
    }
    if public_api["entry_count"] != len(entries):
        errors.append("API coverage entry count differs from the entry list")
    if not (
        public_api["complete"]
        and public_api["present_count"] == public_api["entry_count"]
        and public_api["verified_count"] + public_api["asserted_count"]
        == public_api["entry_count"]
        and public_api["missing_count"] == 0
        and public_api["partial_count"] == 0
        and public_api["unmapped_count"] == 0
        and public_api["required_present_count"] == public_api["required_count"]
        and public_api["required_verified_count"] == public_api["required_count"]
        and public_api["required_asserted_count"] == 0
        and public_api["required_unmapped_count"] == 0
    ):
        errors.append(
            "required public API scope is not completely evidence-mapped"
        )
    closure = compatibility.get("complete_api_closure", {})
    if not isinstance(closure, dict) or not closure.get("complete", False):
        errors.append("compatibility contract does not declare complete API closure")

    notice_text = (
        notice_path.read_text(encoding="utf-8", errors="ignore")
        if notice_path.is_file()
        else ""
    )
    missing_notices = [token for token in _NOTICE_TOKENS if token not in notice_text]
    if missing_notices:
        errors.append("upstream license/trademark notices are incomplete")

    simulink_entries = [
        str(item) for item in manifest.get("simulink_entry_points", [])
    ]
    missing_simulink = [
        relative for relative in simulink_entries if not (matlab_root / relative).is_file()
    ]
    if missing_simulink:
        errors.append("missing standalone Simulink entries: " + ", ".join(missing_simulink))
    if manifest.get("simulink_standalone") is not True:
        errors.append("distribution manifest does not declare standalone Simulink")
    if manifest.get("radia_runtime_required") is not False:
        errors.append("radia-optuna unexpectedly requires the Radia runtime")
    source_commands = (
        _optuna_source_commands(source_path)
        if layout["source_kind"] == "repository"
        else []
    )
    expected_commands = int(manifest.get("native_command_count", -1))
    if layout["source_kind"] == "repository":
        if not source_commands:
            errors.append(f"missing standalone Optuna MEX source contract: {source_path}")
        elif len(source_commands) != expected_commands:
            errors.append(
                "Optuna MEX source command inventory mismatch: "
                f"{len(source_commands)} != {expected_commands}"
            )
    elif not mex_path.is_file():
        errors.append(f"missing native gateway: {mex_path}")

    oracle_path = Path(layout["oracle"])
    recorded_oracle_sha = str(coverage.get("upstream_oracle_sha256", "")).upper()
    actual_oracle_sha = _sha256_lf(oracle_path) if oracle_path.is_file() else None
    if actual_oracle_sha is not None and actual_oracle_sha != recorded_oracle_sha:
        errors.append("upstream oracle SHA256 differs from the coverage record")

    inventory_path = Path(layout["inventory"])
    recorded_inventory_sha = str(
        coverage.get("upstream_inventory_sha256", "")
    ).upper()
    actual_inventory_sha = (
        _sha256_lf(inventory_path) if inventory_path.is_file() else None
    )
    if actual_inventory_sha is not None and actual_inventory_sha != recorded_inventory_sha:
        errors.append("public API inventory SHA256 differs from the coverage record")

    test_manifest_path = Path(layout["test_manifest"])
    test_manifest: dict[str, Any] = {}
    if test_manifest_path.is_file():
        test_manifest = _read_json(test_manifest_path, errors, "test policy manifest")
    test_entries = test_manifest.get("entries", [])
    if not isinstance(test_entries, list):
        test_entries = []
        errors.append("test policy entries must be a list")
    classifications = Counter(
        str(entry.get("classification", ""))
        for entry in test_entries
        if isinstance(entry, dict)
    )

    return {
        "schema": "radia-mcp.matlab-optuna-health/v1",
        "ok": not errors,
        "status": "ready" if not errors else "error",
        "source_kind": layout["source_kind"],
        "root": str(root),
        "matlab_root": str(matlab_root),
        "distribution": {
            "name": "radia-optuna",
            "version": manifest.get("version"),
            "upstream_version": manifest.get("oracle_version"),
            "manifest": str(manifest_path),
            "matlab_file_count": len(matlab_files),
            "expected_matlab_file_count": expected_files,
            "top_level_class_count": len(classes),
            "top_level_function_count": len(functions),
            "top_level_classes": classes,
            "top_level_functions": functions,
            "radia_runtime_required": manifest.get("radia_runtime_required"),
        },
        "public_api": public_api,
        "oracle": {
            "coverage": str(coverage_path),
            "coverage_sha256": _sha256(coverage_path) if coverage_path.is_file() else None,
            "fixture": str(oracle_path) if oracle_path.is_file() else None,
            "recorded_fixture_sha256": recorded_oracle_sha,
            "actual_fixture_sha256": actual_oracle_sha,
            "inventory": str(inventory_path) if inventory_path.is_file() else None,
            "recorded_inventory_sha256": recorded_inventory_sha,
            "actual_inventory_sha256": actual_inventory_sha,
        },
        "native": {
            "gateway": manifest.get("native_gateway"),
            "command_count": expected_commands,
            "path": str(mex_path),
            "present": mex_path.is_file(),
            "runtime_ready": mex_path.is_file(),
            "build_required": (
                layout["source_kind"] == "repository" and not mex_path.is_file()
            ),
            "size_bytes": mex_path.stat().st_size if mex_path.is_file() else None,
            "sha256": _sha256(mex_path) if mex_path.is_file() else None,
            "source_contract": (
                str(source_path) if layout["source_kind"] == "repository" else None
            ),
            "source_contract_present": bool(source_commands),
            "source_command_count": (
                len(source_commands)
                if layout["source_kind"] == "repository"
                else None
            ),
            "source_commands": source_commands,
            "source_sha256": (
                _sha256(source_path) if source_commands else None
            ),
        },
        "simulink": {
            "standalone": manifest.get("simulink_standalone"),
            "entry_count": len(simulink_entries),
            "entries": simulink_entries,
            "missing_entries": missing_simulink,
        },
        "test_policy": {
            "available": test_manifest_path.is_file(),
            "manifest": str(test_manifest_path) if test_manifest_path.is_file() else None,
            "entry_count": len(test_entries),
            "classifications": dict(sorted(classifications.items())),
        },
        "stewardship": {
            "notice": str(notice_path),
            "complete": not missing_notices,
            "missing_tokens": missing_notices,
            "independent_unofficial": True,
        },
        "errors": errors,
    }


def _matlab_literal(value: str) -> str:
    return "'" + value.replace("'", "''").replace("\\", "/") + "'"


def _matlab_cellstr(values: tuple[str, ...]) -> str:
    return "{" + ",".join(_matlab_literal(value) for value in values) + "}"


def matlab_optuna_oracle_plan(
    scope: str = "all", repository_path: str = "", output_path: str = ""
) -> dict[str, Any]:
    """Build an official-MATLAB-MCP-ready differential test plan."""
    key = str(scope or "all").strip().lower().replace("-", "_")
    if key not in _ORACLE_SCOPE_FILES:
        raise ValueError(
            f"unknown scope {scope!r}; expected {sorted(_ORACLE_SCOPE_FILES)}"
        )
    health = matlab_optuna_health(repository_path)
    if not health.get("ok") or health.get("source_kind") != "repository":
        return {
            "schema": "radia-mcp.matlab-optuna-oracle-plan/v1",
            "ok": False,
            "status": "error",
            "scope": key,
            "errors": health.get("errors", [])
            or ["oracle execution requires a Radia repository checkout"],
        }

    root = Path(str(health["root"]))
    prefixes = _ORACLE_SCOPE_FILES[key]
    prefix_code = _matlab_cellstr(prefixes)
    root_code = _matlab_literal(str(root))
    output = Path(output_path or rf"C:\temp\radia-optuna-{key}-summary.json")
    output_code = _matlab_literal(str(output))
    matlab_code = (
        f"repoRoot=string({root_code}); cd(repoRoot); "
        "addpath(fullfile(repoRoot,'matlab')); "
        "suite=testsuite(fullfile(repoRoot,'tests','matlab'),"
        "'IncludeSubfolders',false); names=string({suite.Name}); "
        f"prefixes=string({prefix_code}); selected=false(size(names)); "
        "for prefix=prefixes, selected=selected|startsWith(names,prefix); end; "
        "suite=suite(selected); results=run(suite); "
        "summary=struct('schema','radia.optuna.oracle-test-summary.v1',"
        "'scope'," + _matlab_literal(key) + ","
        "'total',numel(results),'passed',sum([results.Passed]),"
        "'failed',sum([results.Failed]),"
        "'incomplete',sum([results.Incomplete])); "
        "encoded=jsonencode(summary,PrettyPrint=true); disp(encoded); "
        f"fileId=fopen({output_code},'w'); "
        "assert(fileId>=0,'radia:optuna:OracleSummaryOutput',"
        "'Could not open oracle summary output'); "
        "cleanupFile=onCleanup(@()fclose(fileId)); "
        "fprintf(fileId,'%s\\n',encoded); clear cleanupFile; "
        "assert(all([results.Passed]),"
        "'radia:optuna:OraclePlanFailed','Optuna oracle plan failed');"
    )

    fixture_commands: list[str] = []
    if key in {"all", "shared", "samplers", "storage"}:
        fixture_commands.extend(
            [
                "python tests/matlab/fixtures/generate_optuna49_api_inventory.py",
                "python tests/matlab/fixtures/generate_optuna49_oracle.py",
                "python tests/matlab/fixtures/generate_optuna49_api_coverage.py",
                "python tests/matlab/fixtures/generate_optuna_test_manifest.py",
            ]
        )
    if key in {"all", "upstream_mcp"}:
        fixture_commands.append(
            "python tests/matlab/fixtures/generate_optuna49_mcp_oracle.py"
        )

    return {
        "schema": "radia-mcp.matlab-optuna-oracle-plan/v1",
        "ok": True,
        "status": "ready",
        "scope": key,
        "runtime_owner": "MathWorks official MATLAB MCP Server",
        "execute_with": "evaluate_matlab_code",
        "repository_root": str(root),
        "test_prefixes": list(prefixes),
        "expected_test_count": health["test_policy"]["entry_count"]
        if key == "all"
        else None,
        "matlab_code": matlab_code,
        "output": str(output),
        "fixture_regeneration": {
            "runtime": "pinned direct Python and official optuna-mcp stdio",
            "working_directory": str(root),
            "commands": fixture_commands,
            "byte_stability_rule": (
                "Run each selected generator twice and require unchanged SHA256."
            ),
        },
        "oracle": health["oracle"],
        "ownership": {
            "shared_study_trial_tools": "optuna/optuna-mcp",
            "seeded_numeric_oracle": "direct optuna==4.9.0",
            "matlab_execution": "MathWorks official MATLAB MCP Server",
            "matlab_difference_contract": "radia-mcp.matlab",
        },
    }


def _python_constant(source: str, name: str) -> float:
    match = re.search(rf"(?m)^{re.escape(name)}\s*=\s*([0-9.]+)$", source)
    if match is None:
        raise RuntimeError(f"benchmark constant {name} is missing")
    return float(match.group(1))


def matlab_optuna_benchmark_plan(
    repository_path: str = "", output_directory: str = r"C:\temp\radia-optuna-benchmark"
) -> dict[str, Any]:
    """Build a same-host cold/warm Python-versus-MATLAB benchmark plan."""
    health = matlab_optuna_health(repository_path)
    if not health.get("ok") or health.get("source_kind") != "repository":
        return {
            "schema": "radia-mcp.matlab-optuna-benchmark-plan/v1",
            "ok": False,
            "status": "error",
            "errors": health.get("errors", [])
            or ["benchmarking requires a Radia repository checkout"],
        }
    root = Path(str(health["root"]))
    benchmark_root = root / "validation_test" / "optimization"
    python_script = benchmark_root / "benchmark_optuna49_python.py"
    matlab_script = benchmark_root / "benchmark_matlab_optuna49.m"
    cold_script = benchmark_root / "benchmark_optuna_mex_cold_start.ps1"
    missing = [str(path) for path in (python_script, matlab_script, cold_script) if not path.is_file()]
    if missing:
        return {
            "schema": "radia-mcp.matlab-optuna-benchmark-plan/v1",
            "ok": False,
            "status": "error",
            "errors": ["missing benchmark scripts: " + ", ".join(missing)],
        }

    source = python_script.read_text(encoding="utf-8")
    trials = int(_python_constant(source, "TRIALS"))
    repeats = int(_python_constant(source, "REPEATS"))
    warmups = int(_python_constant(source, "WARMUP_REPEATS"))
    sampler_seeds = [
        int(value)
        for value in re.findall(r"TPESampler\(\s*seed=(\d+)", source)
    ]
    if len(sampler_seeds) != 2:
        raise RuntimeError(
            "benchmark must define exactly the scalar and grouped TPE seeds"
        )
    scalar_seed, group_seed = sampler_seeds
    output = Path(output_directory).expanduser()
    python_output = output / "optuna49_python.json"
    matlab_output = output / "optuna49_matlab.json"
    cold_output = output / "optuna_mex_cold.json"
    root_literal = _matlab_literal(str(root))
    matlab_output_literal = _matlab_literal(str(matlab_output))
    expected_native_commands = int(health["native"]["command_count"])
    matlab_code = (
        f"repoRoot=string({root_literal}); cd(repoRoot); "
        "addpath(fullfile(repoRoot,'matlab')); "
        "addpath(fullfile(repoRoot,'validation_test','optimization')); "
        f"result=benchmark_matlab_optuna49(string({matlab_output_literal})); "
        "nativeStatus=radia.optuna.nativeStatus(); "
        "assert(nativeStatus.mex_available); "
        "assert(result.versions.optuna_mex_command_count=="
        f"{expected_native_commands});"
    )
    return {
        "schema": "radia-mcp.matlab-optuna-benchmark-plan/v1",
        "ok": True,
        "status": "ready",
        "repository_root": str(root),
        "output_directory": str(output),
        "same_host_required": True,
        "idle_host_required": True,
        "settings": {
            "trials": trials,
            "total_repeats": repeats,
            "warmup_repeats": warmups,
            "reported_repeats": repeats - warmups,
            "scalar_seed": scalar_seed,
            "grouped_conditional_seed": group_seed,
        },
        "cold_start": {
            "scope": "first optuna_mex api.info call; MATLAB startup excluded",
            "powershell_command": (
                f"pwsh -File \"{cold_script}\" -Repeats 7 -Output \"{cold_output}\""
            ),
            "output": str(cold_output),
        },
        "python_warmed": {
            "runtime": "direct optuna==4.9.0",
            "command": f'python "{python_script}" --output "{python_output}"',
            "output": str(python_output),
        },
        "matlab_warmed": {
            "runtime_owner": "MathWorks official MATLAB MCP Server",
            "execute_with": "evaluate_matlab_code",
            "matlab_code": matlab_code,
            "output": str(matlab_output),
        },
        "acceptance": {
            "checksum_absolute_tolerance": 1e-12,
            "same_settings": True,
            "same_host": True,
            "max_matlab_to_python_warmed_time_ratio": 1.0,
            "workloads": ["scalar", "grouped_conditional", "trials_dataframe"],
            "startup_measured_separately": True,
        },
    }


def _mapping(value: Any, label: str, errors: list[str]) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        errors.append(f"{label} must be an object")
        return {}
    return value


def _positive_float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if result > 0 else None


def _integer(value: Any, default: int = -1) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def matlab_optuna_release_gate(
    evidence_json: str | Mapping[str, Any],
    repository_path: str = "",
    max_warmed_time_ratio: float = 1.0,
) -> dict[str, Any]:
    """Validate installed-wheel, MATLAB, resume, and performance evidence."""
    if max_warmed_time_ratio <= 0:
        raise ValueError("max_warmed_time_ratio must be positive")
    if isinstance(evidence_json, str):
        evidence = json.loads(evidence_json)
    else:
        evidence = evidence_json
    if not isinstance(evidence, Mapping):
        raise TypeError("evidence_json must encode a JSON object")

    health = matlab_optuna_health(repository_path)
    errors = list(health.get("errors", []))
    expected_version = health.get("distribution", {}).get("version")
    expected_files = health.get("distribution", {}).get("expected_matlab_file_count")
    expected_commands = health.get("native", {}).get("command_count")
    expected_simulink = health.get("simulink", {}).get("entry_count")

    required = {
        "wheel_verification",
        "doctor",
        "matlab_tests",
        "simulink_e2e",
        "table_resume",
        "performance",
    }
    missing_sections = sorted(required.difference(evidence))
    if missing_sections:
        errors.append("missing evidence sections: " + ", ".join(missing_sections))

    wheel = _mapping(evidence.get("wheel_verification"), "wheel_verification", errors)
    if wheel:
        if wheel.get("schema") != "radia-optuna.wheel-verification.v1" or wheel.get("ok") is not True:
            errors.append("wheel verifier did not pass")
        if wheel.get("version") != expected_version:
            errors.append("wheel version differs from the distribution manifest")
        if wheel.get("matlab_file_count") != expected_files:
            errors.append("wheel MATLAB file count differs from the manifest")
        if wheel.get("native_command_count") != expected_commands:
            errors.append("wheel native command count differs from the manifest")
        if wheel.get("simulink_entry_count") != expected_simulink:
            errors.append("wheel Simulink entry count differs from the manifest")
        if wheel.get("source_fidelity_verified") is not True:
            errors.append("wheel payload was not matched to the checked source tree")
        if _integer(wheel.get("source_fidelity_file_count")) < _integer(
            expected_files
        ):
            errors.append("wheel source-fidelity inventory is incomplete")

    doctor = _mapping(evidence.get("doctor"), "doctor", errors)
    if doctor:
        if doctor.get("schema") != "radia-optuna.doctor.v1" or doctor.get("ok") is not True:
            errors.append("installed-wheel doctor did not pass")
        if doctor.get("version") != expected_version:
            errors.append("installed-wheel version differs from the manifest")
        if doctor.get("matlab_file_count") != expected_files:
            errors.append("installed-wheel MATLAB file count differs from the manifest")
        if doctor.get("native_command_count") != expected_commands:
            errors.append("installed-wheel native command count differs from the manifest")
        if doctor.get("simulink_standalone") is not True:
            errors.append("installed-wheel doctor did not find standalone Simulink")
        if doctor.get("simulink_entry_count") != expected_simulink:
            errors.append("installed-wheel Simulink entry count differs from the manifest")
        if doctor.get("upstream_notices_complete") is not True:
            errors.append("installed-wheel upstream notices are incomplete")
        doctor_sha = str(doctor.get("mex_sha256", "")).casefold()
        if re.fullmatch(r"[0-9a-f]{64}", doctor_sha) is None:
            errors.append("installed-wheel doctor has no valid MEX SHA256")
        source_binary_sha = health.get("native", {}).get("sha256")
        if source_binary_sha and doctor_sha != str(source_binary_sha).casefold():
            errors.append("installed-wheel MEX differs from the source-build binary")

    tests = _mapping(evidence.get("matlab_tests"), "matlab_tests", errors)
    if tests:
        if tests.get("schema") != "radia.optuna.oracle-test-summary.v1":
            errors.append("MATLAB test summary has an invalid schema")
        expected_tests = _integer(
            health.get("test_policy", {}).get("entry_count"), default=0
        )
        total = _integer(tests.get("total"))
        if expected_tests and total != expected_tests:
            errors.append(f"MATLAB test count differs from policy: {total} != {expected_tests}")
        if _integer(tests.get("passed")) != total:
            errors.append("not every MATLAB Optuna test passed")
        if _integer(tests.get("failed")) != 0 or _integer(
            tests.get("incomplete")
        ) != 0:
            errors.append("MATLAB Optuna tests contain failed or incomplete results")

    simulink = _mapping(evidence.get("simulink_e2e"), "simulink_e2e", errors)
    if simulink:
        if simulink.get("schema") != "radia-optuna.standalone-simulink-test.v1" or simulink.get("ok") is not True:
            errors.append("standalone installed-wheel Simulink E2E did not pass")
        if _integer(simulink.get("complete_trials"), default=0) < 1:
            errors.append("Simulink E2E has no completed trial")
        if _integer(simulink.get("failed_trials"), default=0) < 1:
            errors.append("Simulink E2E did not exercise typed failure telemetry")
        if _integer(simulink.get("simulink_block_trials"), default=0) < 1:
            errors.append("Simulink optimization block produced no persisted trials")

    resume = _mapping(evidence.get("table_resume"), "table_resume", errors)
    if resume:
        if resume.get("schema") != "radia-optuna.table-resume-test.v1" or resume.get("ok") is not True:
            errors.append("table/MAT resume evidence did not pass")
        source_trials = _integer(resume.get("source_trial_count"))
        restored_trials = _integer(resume.get("restored_trial_count"), default=-2)
        if source_trials < 1 or restored_trials != source_trials:
            errors.append("table/MAT resume did not restore every source trial")
        tables = {str(name) for name in resume.get("table_names", [])}
        missing_tables = sorted(_TABLE_NAMES.difference(tables))
        if missing_tables:
            errors.append("table/MAT resume is missing tables: " + ", ".join(missing_tables))

    performance = _mapping(evidence.get("performance"), "performance", errors)
    ratios: dict[str, float] = {}
    if performance:
        python_result = _mapping(performance.get("python"), "performance.python", errors)
        matlab_result = _mapping(performance.get("matlab"), "performance.matlab", errors)
        cold = _mapping(performance.get("mex_cold"), "performance.mex_cold", errors)
        for result, runtime in ((python_result, "python-upstream"), (matlab_result, "matlab")):
            if result and (
                result.get("schema") != "radia.validation.optuna49-performance-runtime.v1"
                or result.get("runtime") != runtime
            ):
                errors.append(f"invalid {runtime} performance evidence")
        if python_result and matlab_result:
            python_host = str(python_result.get("host", "")).casefold()
            matlab_host = str(matlab_result.get("host", "")).casefold()
            if not python_host or python_host != matlab_host:
                errors.append("Python and MATLAB performance evidence used different hosts")
            if python_result.get("settings") != matlab_result.get("settings"):
                errors.append("Python and MATLAB performance settings differ")
            python_versions = _mapping(
                python_result.get("versions"), "performance.python.versions", errors
            )
            matlab_versions = _mapping(
                matlab_result.get("versions"), "performance.matlab.versions", errors
            )
            if python_versions.get("optuna") != health.get("distribution", {}).get(
                "upstream_version"
            ):
                errors.append("Python performance evidence is not pinned to the oracle")
            if matlab_versions.get("optuna_mex_command_count") != expected_commands:
                errors.append("MATLAB performance MEX inventory differs from the manifest")
            for workload in ("scalar", "grouped_conditional", "trials_dataframe"):
                py_work = _mapping(python_result.get(workload), f"python.{workload}", errors)
                ml_work = _mapping(matlab_result.get(workload), f"matlab.{workload}", errors)
                py_time = _positive_float(py_work.get("median_warmed_seconds"))
                ml_time = _positive_float(ml_work.get("median_warmed_seconds"))
                if py_time is None or ml_time is None:
                    errors.append(f"{workload} has invalid warmed timing")
                    continue
                ratio = ml_time / py_time
                ratios[workload] = ratio
                if ratio > max_warmed_time_ratio:
                    errors.append(
                        f"{workload} MATLAB/Python warmed time ratio {ratio:.6g} "
                        f"exceeds {max_warmed_time_ratio:.6g}"
                    )
                if workload != "trials_dataframe":
                    try:
                        checksum_error = abs(
                            float(py_work.get("checksum")) - float(ml_work.get("checksum"))
                        )
                    except (TypeError, ValueError):
                        errors.append(f"{workload} checksum is missing")
                    else:
                        if checksum_error > 1e-12:
                            errors.append(f"{workload} seeded checksum differs")
                else:
                    if py_work.get("rows") != ml_work.get("rows"):
                        errors.append("trials_dataframe row count differs")
                    if py_work.get("columns") != ml_work.get("columns"):
                        errors.append("trials_dataframe column count differs")
        if cold and (
            cold.get("schema") != "radia.validation.optuna-mex-first-call.v1"
            or _positive_float(cold.get("median_first_call_seconds")) is None
        ):
            errors.append("invalid optuna_mex cold-start evidence")
        if cold:
            cold_host = str(cold.get("host", "")).casefold()
            if python_result and cold_host != str(
                python_result.get("host", "")
            ).casefold():
                errors.append("cold MEX and warmed performance used different hosts")
            cold_sha = str(cold.get("binary_sha256", "")).casefold()
            expected_sha = str(doctor.get("mex_sha256", "")).casefold()
            if not cold_sha or cold_sha != expected_sha:
                errors.append("cold-start MEX binary differs from the installed wheel")

    return {
        "schema": "radia-mcp.matlab-optuna-release-gate/v1",
        "ok": not errors,
        "status": "ready" if not errors else "error",
        "distribution_version": expected_version,
        "upstream_version": health.get("distribution", {}).get("upstream_version"),
        "public_api": health.get("public_api"),
        "performance": {
            "max_warmed_time_ratio": max_warmed_time_ratio,
            "observed_matlab_to_python_ratios": ratios,
        },
        "ownership": {
            "execution": "MathWorks official MATLAB MCP Server",
            "shared_optuna_tools": "optuna/optuna-mcp",
            "release_gate": "radia-mcp.matlab",
        },
        "errors": errors,
    }


__all__ = [
    "matlab_optuna_benchmark_plan",
    "matlab_optuna_health",
    "matlab_optuna_oracle_plan",
    "matlab_optuna_release_gate",
]
