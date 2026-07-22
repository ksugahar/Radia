"""Headless application runner used by Radia Simulink analysis blocks.

The Simulink block is an operating surface, not a second solver.  This module
loads an application ``DesignSpec``, builds its validated ``calc_*.py`` command,
and records a stable artifact bundle for MATLAB to consume.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import importlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time
from typing import Any, Mapping, Sequence


CONFIG_SCHEMA = "radia.simulink.application_config.v1"
RESULT_SCHEMA = "radia.simulink.application_run.v1"


@dataclass(frozen=True)
class ApplicationDefinition:
    spec_module: str
    spec_class: str
    primary_keys: tuple[str, ...]


APPLICATIONS: dict[str, ApplicationDefinition] = {
    "em": ApplicationDefinition(
        "radia.em_design",
        "EMDesignSpec",
        ("B_origin_mag_T", "Bz_T", "B_peak_T", "M_avg", "rms"),
    ),
    "pcb": ApplicationDefinition(
        "radia.pcb_design",
        "PCBDesignSpec",
        ("inductance_H", "L_H", "resistance_ohm", "impedance_ohm", "n_conductors"),
    ),
    "motor": ApplicationDefinition(
        "radia.motor_design",
        "MotorDesignSpec",
        ("mean_torque_Nm", "torque_Nm", "peak_torque_Nm", "total_loss_W", "loss_W"),
    ),
    "streamfunction": ApplicationDefinition(
        "radia.streamfunction_design",
        "StreamFunctionDesignSpec",
        ("rms", "wire_homogeneity", "peak_J", "inductance_H"),
    ),
    "ih": ApplicationDefinition(
        "radia.ih_design",
        "IHDesignSpec",
        ("P_wp_W", "L_total_nH", "L_coil_nH", "L_nH", "T_mean_C"),
    ),
}


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _command_line(command: Sequence[str]) -> str:
    if sys.platform == "win32":
        return subprocess.list2cmdline(list(command))
    return shlex.join(command)


def _load_config(
    path: Path,
) -> tuple[dict[str, Any], str | None, Path | None, str | None]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise ValueError("The Simulink application config must be a JSON object.")

    schema = payload.get("schema")
    if schema is not None and schema != CONFIG_SCHEMA:
        raise ValueError(f"Unsupported config schema: {schema!r}")

    if "settings" in payload:
        settings = payload["settings"]
        if not isinstance(settings, dict):
            raise ValueError("config.settings must be a JSON object.")
    else:
        settings = {
            key: value
            for key, value in payload.items()
            if key not in {
                "schema",
                "application",
                "primary_key",
                "working_directory",
            }
        }

    primary_key = payload.get("primary_key")
    if primary_key is not None and not isinstance(primary_key, str):
        raise ValueError("config.primary_key must be a string when provided.")
    working_directory = payload.get("working_directory")
    cwd = Path(working_directory) if working_directory else None
    declared_application = payload.get("application")
    if declared_application is not None and not isinstance(declared_application, str):
        raise ValueError("config.application must be a string when provided.")
    return settings, primary_key, cwd, declared_application


def _load_spec(application: str, settings: Mapping[str, Any]) -> tuple[Any, ApplicationDefinition]:
    try:
        definition = APPLICATIONS[application]
    except KeyError as exc:
        supported = ", ".join(sorted(APPLICATIONS))
        raise ValueError(
            f"Unknown application {application!r}; choose one of {supported}."
        ) from exc
    module = importlib.import_module(definition.spec_module)
    spec_type = getattr(module, definition.spec_class)
    return spec_type(**dict(settings)), definition


def _replace_output_path(command: Sequence[str], output_path: Path) -> list[str]:
    result = [str(part) for part in command]
    try:
        index = result.index("--output")
    except ValueError as exc:
        raise ValueError(
            "The application command does not expose the required --output contract."
        ) from exc
    if index + 1 >= len(result):
        raise ValueError("The application command has --output without a path.")
    result[index + 1] = str(output_path)
    return result


def _replace_optional_output_path(
    command: Sequence[str],
    option: str,
    output_path: Path,
) -> tuple[list[str], bool]:
    """Redirect an optional runner-owned artifact path into the run directory."""

    result = [str(part) for part in command]
    try:
        index = result.index(option)
    except ValueError:
        return result, False
    if index + 1 >= len(result):
        raise ValueError(f"The application command has {option} without a path.")
    result[index + 1] = str(output_path)
    return result, True


def _is_gmsh_v41(path: Path) -> bool:
    """Return whether *path* begins with the repository's ASCII MSH v4.1 header."""

    try:
        header = path.read_bytes()[:256].replace(b"\r\n", b"\n")
    except OSError:
        return False
    return header.startswith(b"$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")


def _find_key(value: Any, key: str) -> Any:
    if isinstance(value, dict):
        if key in value:
            return value[key]
        for child in value.values():
            found = _find_key(child, key)
            if found is not None:
                return found
    elif isinstance(value, list):
        for child in value:
            found = _find_key(child, key)
            if found is not None:
                return found
    return None


def _find_path(value: Any, dotted_key: str) -> Any:
    current = value
    for key in dotted_key.split("."):
        if not isinstance(current, dict) or key not in current:
            return None
        current = current[key]
    return current


def _numeric_scalar(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    result = float(value)
    return result if math.isfinite(result) else None


def _primary_value(
    solver_payload: Any,
    requested_key: str | None,
    candidates: Sequence[str],
) -> tuple[str | None, float | None]:
    keys = (requested_key,) if requested_key else tuple(candidates)
    for key in keys:
        if not key:
            continue
        value = _find_path(solver_payload, key) if "." in key else _find_key(solver_payload, key)
        scalar = _numeric_scalar(value)
        if scalar is not None:
            return key, scalar
    return None, None


def _radia_version() -> str:
    try:
        return importlib.metadata.version("radia")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def run_application(
    application: str,
    config_path: str | Path,
    run_dir: str | Path,
    *,
    timeout_s: float = 3600.0,
) -> dict[str, Any]:
    """Run one DesignSpec-backed application and always write ``result.json``."""

    config_path = Path(config_path).resolve()
    run_dir = Path(run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    log_path = run_dir / "run.log"
    command_path = run_dir / "command.txt"
    solver_result_path = run_dir / "solver_result.json"
    result_path = run_dir / "result.json"
    gmsh_output_path = run_dir / f"{application}_fields.msh"

    started_at = _utc_now()
    started = time.monotonic()
    status = "failed"
    returncode: int | None = None
    command: list[str] = []
    error: str | None = None
    primary_key: str | None = None
    primary_value: float | None = None
    gmsh_requested = False

    with log_path.open("w", encoding="utf-8") as log:
        try:
            for stale_path in (command_path, solver_result_path, result_path):
                stale_path.unlink(missing_ok=True)
            for stale_gmsh in run_dir.rglob("*.msh"):
                stale_gmsh.unlink()
            if not config_path.is_file():
                raise FileNotFoundError(f"Configuration JSON not found: {config_path}")
            settings, requested_primary, cwd, declared_application = _load_config(config_path)
            if declared_application is not None and declared_application != application:
                raise ValueError(
                    f"Configuration targets {declared_application!r}, not {application!r}."
                )
            if cwd is not None and not cwd.is_dir():
                raise FileNotFoundError(f"Working directory not found: {cwd}")
            spec, definition = _load_spec(application, settings)
            missing = (
                list(spec.missing_required_inputs())
                if hasattr(spec, "missing_required_inputs")
                else []
            )
            if missing:
                raise ValueError(
                    "Missing required inputs: "
                    + ", ".join(str(item) for item in missing)
                )

            command = _replace_output_path(
                spec.build_command(python=sys.executable),
                solver_result_path,
            )
            command, gmsh_requested = _replace_optional_output_path(
                command,
                "--msh-output",
                gmsh_output_path,
            )
            command_text = _command_line(command)
            command_path.write_text(command_text + "\n", encoding="utf-8")
            log.write(command_text + "\n\n")
            log.flush()
            try:
                completed = subprocess.run(
                    command,
                    cwd=str(cwd) if cwd else None,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    text=True,
                    timeout=None if timeout_s <= 0 else timeout_s,
                    check=False,
                )
                returncode = completed.returncode
                status = "passed" if returncode == 0 else "failed"
            except subprocess.TimeoutExpired:
                status = "timeout"
                error = f"Application exceeded timeout_s={timeout_s:g}."
                log.write("\n" + error + "\n")

            if solver_result_path.is_file():
                try:
                    solver_payload = json.loads(solver_result_path.read_text(encoding="utf-8-sig"))
                    primary_key, primary_value = _primary_value(
                        solver_payload,
                        requested_primary,
                        definition.primary_keys,
                    )
                except Exception as exc:
                    if status == "passed":
                        status = "failed"
                        error = f"Could not read solver_result.json: {exc}"
                    log.write(f"\nCould not read solver_result.json: {exc}\n")
            elif status == "passed":
                status = "failed"
                error = "The application exited successfully without writing solver_result.json."
                log.write("\n" + error + "\n")

            if status == "passed" and gmsh_requested:
                if not gmsh_output_path.is_file():
                    status = "failed"
                    error = (
                        "The spatial application exited successfully without writing "
                        f"the required GMSH artifact: {gmsh_output_path}"
                    )
                    log.write("\n" + error + "\n")
                elif not _is_gmsh_v41(gmsh_output_path):
                    status = "failed"
                    error = (
                        "The spatial application wrote a GMSH artifact that is not "
                        f"ASCII .msh v4.1: {gmsh_output_path}"
                    )
                    log.write("\n" + error + "\n")
        except Exception as exc:
            error = str(exc)
            log.write(error + "\n")

    elapsed_s = time.monotonic() - started
    gmsh_artifacts = sorted(
        str(path.resolve()) for path in run_dir.rglob("*.msh") if path.is_file()
    )
    payload: dict[str, Any] = {
        "radia_result": {
            "schema": RESULT_SCHEMA,
            "application": application,
            "backend": "python-headless-cli",
            "status": status,
            "returncode": returncode,
            "error": error,
            "executed_at_utc": started_at,
            "completed_at_utc": _utc_now(),
            "elapsed_s": round(elapsed_s, 6),
            "timeout_s": timeout_s,
            "runtime_radia_version": _radia_version(),
            "runtime_python": platform.python_version(),
            "runtime_platform": platform.platform(),
            "config": str(config_path),
            "config_sha256": (
                hashlib.sha256(config_path.read_bytes()).hexdigest()
                if config_path.is_file()
                else None
            ),
            "run_dir": str(run_dir),
            "log": str(log_path),
            "command_file": str(command_path) if command_path.is_file() else None,
            "solver_result": str(solver_result_path) if solver_result_path.is_file() else None,
            "command": command,
            "command_line": _command_line(command) if command else "",
            "primary": {"key": primary_key, "value": primary_value},
            "artifacts": {
                "gmsh_policy": "required" if gmsh_requested else "not-applicable",
                "gmsh_format": "msh-v4.1" if gmsh_requested else None,
                "gmsh_primary": (
                    str(gmsh_output_path) if gmsh_output_path.is_file() else None
                ),
                "gmsh": gmsh_artifacts,
            },
        }
    }
    result_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return payload


def build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--application", required=True, choices=sorted(APPLICATIONS))
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--timeout", type=float, default=3600.0)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_argparser().parse_args(argv)
    payload = run_application(
        args.application,
        args.config,
        args.run_dir,
        timeout_s=args.timeout,
    )
    result_path = Path(args.run_dir).resolve() / "result.json"
    print(result_path)
    return 0 if payload["radia_result"]["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
