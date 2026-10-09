"""Generate real-solve IH Simulink operator validation evidence."""

from __future__ import annotations

import json
import hashlib
import importlib.metadata
import os
import platform
import socket
import subprocess
import tempfile
import time
from datetime import UTC, datetime
from pathlib import Path

import ngsolve as ng
import radia
import test_ih_operator_physics_golden as ih_golden

OUTPUT = Path(__file__).with_name("ih_operator_physics_results.json")
ROOT = Path(__file__).resolve().parents[2]


def _git(*args: str) -> str:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, check=True, capture_output=True, text=True,
    )
    return result.stdout.strip()


def _source_identity() -> tuple[str, bool]:
    try:
        return (
            _git("rev-parse", "HEAD"),
            bool(_git("status", "--porcelain", "--untracked-files=no")),
        )
    except (OSError, subprocess.CalledProcessError):
        commit = os.environ.get("RADIA_SOURCE_COMMIT", "").strip()
        dirty = os.environ.get("RADIA_SOURCE_DIRTY", "").strip().lower()
        if not commit or dirty not in {"0", "1", "false", "true"}:
            raise RuntimeError(
                "source checkout has no readable Git identity; set "
                "RADIA_SOURCE_COMMIT and RADIA_SOURCE_DIRTY explicitly"
            )
        return commit, dirty in {"1", "true"}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    started = time.perf_counter()
    source_commit, source_dirty = _source_identity()
    native_module = Path(radia.__file__).resolve().with_name("_radia_pybind.pyd")
    native_build_commit = os.environ.get("RADIA_NATIVE_BUILD_COMMIT", "").strip()
    if not native_build_commit:
        raise RuntimeError("RADIA_NATIVE_BUILD_COMMIT is required for native provenance")
    temp_root = Path(os.environ.get("RADIA_TEMP_DIR", r"C:\temp"))
    temp_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="radia-ih-golden-", dir=temp_root) as directory:
        with ng.TaskManager():
            assembled = ih_golden._assemble(Path(directory))
            direct_power = ih_golden._direct_power(assembled)
            metrics = ih_golden._physics_metrics(assembled["config"])
        config = assembled["config"]

    power = metrics["electromagnetic_power_W"]
    power_mismatch = abs(power - direct_power) / abs(direct_power)
    observed = {
        **metrics,
        "direct_calc_inductance_power_W": direct_power,
        "assembler_to_direct_power_relative_difference": power_mismatch,
        "operator_basis": config["operator_basis"],
        "surrogate": config["surrogate"],
        "n_eddy_unknown": config["n_eddy_unknown"],
        "n_heat": config["n_heat"],
    }
    thresholds = {
        "electromagnetic_power_W_minimum": ih_golden.P_WP_BAND_W[0],
        "electromagnetic_power_W_maximum": ih_golden.P_WP_BAND_W[1],
        "assembler_to_direct_power_relative_difference_maximum": 1e-9,
        "relative_power_error_maximum": 1e-6,
        "thermal_mass_relative_error_maximum": 1e-10,
        "stiffness_nullspace_relative_error_maximum": 1e-10,
        "convection_area_relative_error_maximum": 1e-10,
        "temperature_weight_row_sum_relative_error_maximum": 1e-10,
    }
    mass_reference = ih_golden.RHO_CP * ih_golden.BOX_VOLUME_M3
    stiffness_scale = metrics["stiffness_maximum_absolute_entry"]
    passed = (
        ih_golden.P_WP_BAND_W[0] <= power <= ih_golden.P_WP_BAND_W[1]
        and power_mismatch
        <= thresholds["assembler_to_direct_power_relative_difference_maximum"]
        and metrics["relative_power_error"]
        < thresholds["relative_power_error_maximum"]
        and abs(metrics["mass_sum_J_per_K"] - mass_reference) / mass_reference
        <= thresholds["thermal_mass_relative_error_maximum"]
        and metrics["stiffness_constant_nullspace_maximum"]
        <= thresholds["stiffness_nullspace_relative_error_maximum"] * stiffness_scale
        and abs(metrics["convection_sum_m2"] - ih_golden.BOX_AREA_M2)
        / ih_golden.BOX_AREA_M2
        <= thresholds["convection_area_relative_error_maximum"]
        and metrics["temperature_weight_row_sum_maximum_relative_error"]
        <= thresholds["temperature_weight_row_sum_relative_error_maximum"]
        and metrics["temperature_weight_minimum_J_per_K"] > 0.0
        and config["operator_basis"] == "exact-single-current-linear-response"
        and config["surrogate"] is False
        and config["n_eddy_unknown"] == 1
        and config["eddy_matrix_real"] == [1.0]
        and config["eddy_matrix_imag"] == [0.0]
        and config["eddy_rhs_real"] == [1.0]
        and config["eddy_rhs_imag"] == [0.0]
    )
    result = {
        "schema": "radia.validation.ih-operator-physics.v1",
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "provenance": {
            "command": (
                "python validation_test/induction_heating/"
                "generate_ih_operator_physics_results.py"
            ),
            "python": platform.python_version(),
            "ngsolve": getattr(ng, "__version__", "unknown"),
            "platform_class": platform.system(),
            "commit": source_commit,
            "dirty": source_dirty,
            "source_sha256": {
                path: _sha256(ROOT / path)
                for path in (
                    "src/radia/simulink/ih_operator_assembly.py",
                    "validation_test/induction_heating/test_ih_operator_physics_golden.py",
                    "validation_test/induction_heating/generate_ih_operator_physics_results.py",
                    "packages/cubit-mesh-export/src/cubit_mesh_export/check.py",
                )
            },
            "native_runtime": {
                "distribution_version": importlib.metadata.version("radia"),
                "build_commit": native_build_commit,
                "module_sha256": _sha256(native_module),
                "wheel_sha256": os.environ.get("RADIA_NATIVE_WHEEL_SHA256", "").strip()
                or None,
            },
        },
        "configuration": {
            "frequency_hz": 7000.0,
            "workpiece_conductivity_S_per_m": 5.0e6,
            "workpiece_relative_permeability": 100.0,
            "coil_conductivity_S_per_m": 5.8e7,
            "coupling_mode": "weak",
            "bem_backend": "intree-dense",
            "workpiece_box_volume_m3": ih_golden.BOX_VOLUME_M3,
            "workpiece_box_area_m2": ih_golden.BOX_AREA_M2,
        },
        "runtime_seconds": time.perf_counter() - started,
        "thresholds": thresholds,
        "observed": observed,
        "passed": passed,
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)
    print("PASS" if passed else "FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
