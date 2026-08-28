from __future__ import annotations

import importlib.util
import inspect
import json
import re
import struct
import tomllib
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PACKAGE_ROOT.parents[1]
VERIFIER_PATH = PACKAGE_ROOT / "verify_wheel.py"
VERIFIER_SPEC = importlib.util.spec_from_file_location(
    "radia_optuna_wheel_verifier", VERIFIER_PATH
)
assert VERIFIER_SPEC is not None and VERIFIER_SPEC.loader is not None
wheel_verifier = importlib.util.module_from_spec(VERIFIER_SPEC)
VERIFIER_SPEC.loader.exec_module(wheel_verifier)


def test_independent_version_and_radia_extras_are_synchronized():
    root_project = tomllib.loads(
        (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    )["project"]
    package_project = tomllib.loads(
        (PACKAGE_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    )["project"]
    version = package_project["version"]
    assert version != root_project["version"]
    assert root_project["optional-dependencies"]["optuna"] == [
        f"radia-optuna=={version}"
    ]
    assert root_project["optional-dependencies"]["optuna-upstream"] == [
        f"radia-optuna[upstream]=={version}"
    ]
    init_source = (
        PACKAGE_ROOT / "src" / "radia_optuna" / "__init__.py"
    ).read_text(encoding="utf-8")
    assert f'__version__ = "{version}"' in init_source


def test_only_declared_adapters_call_other_radia_matlab_namespaces():
    manifest = json.loads(
        (PACKAGE_ROOT / "src" / "radia_optuna" / "manifest.json").read_text(
            encoding="utf-8"
        )
    )
    source = REPO_ROOT / "matlab"
    observed: set[str] = set()
    pattern = re.compile(r"\bradia\.(?!optuna\b)")
    for path in (source / "+radia" / "+optuna").rglob("*.m"):
        if pattern.search(path.read_text(encoding="utf-8", errors="ignore")):
            observed.add(path.relative_to(source).as_posix())
    assert observed == set(manifest["radia_integration_adapters"])
    assert manifest["matlab_file_count"] == len(
        list((source / "+radia" / "+optuna").rglob("*.m"))
    )
    assert manifest["simulink_standalone"] is True
    for relative in manifest["simulink_entry_points"]:
        assert (source / relative).is_file()
    assert set(manifest["simulink_entry_points"]) >= {
        "+radia/+simulink/buildOptunaBlock.m",
        "+radia/+simulink/optunaSFunction.m",
        "+radia/+simulink/optunaRuntimeStore.m",
        "+radia/+simulink/addOptunaMonitor.m",
        "radia_optuna_sfun.m",
    }


def test_staging_refuses_a_partial_native_distribution():
    setup_source = (PACKAGE_ROOT / "setup.py").read_text(encoding="utf-8")
    assert "Build.ps1 -OptunaMexOnly" in setup_source
    assert "optuna_mex.mexw64" in setup_source
    assert "optuna_upstream_compatibility.json" in setup_source
    assert "optuna49_api_coverage.json" in setup_source
    assert "THIRD_PARTY_NOTICES.md" in setup_source
    assert 'return "py3", "none", "win_amd64"' in setup_source
    assert '"bdist_wheel": bdist_wheel' in setup_source


def test_wheel_verifier_rejects_solver_boundary_leaks():
    verifier = (PACKAGE_ROOT / "verify_wheel.py").read_text(encoding="utf-8")
    for forbidden in ("radia_mex", "ngsolve", "netgen", "mkl_", "radia_pybind"):
        assert forbidden in verifier
    assert "py3-none-win_amd64" in verifier
    assert "optuna\\s*==\\s*4\\.9\\.0" in verifier
    assert "THIRD_PARTY_NOTICES.md" in verifier
    assert "artifact_integrity_verified" in verifier
    assert "source_fidelity_verified" in verifier
    assert "wheel payload differs from the checked monorepo source" in verifier


def test_wheel_verifier_defaults_to_source_fidelity_and_exposes_artifact_only():
    parameter = inspect.signature(wheel_verifier.verify).parameters["source_fidelity"]
    source = VERIFIER_PATH.read_text(encoding="utf-8")
    assert parameter.default is True
    assert '"--artifact-only"' in source
    assert "source_fidelity=not args.artifact_only" in source


def test_wheel_source_comparison_normalizes_text_line_endings():
    member = "radia_optuna/matlab/optuna_upstream_compatibility.json"
    assert wheel_verifier._normalized_payload(member, b"{\r\n}\r\n") == b"{\n}\n"
    assert wheel_verifier._normalized_payload(member, b"{\n}\n") == b"{\n}\n"


def _synthetic_pe(coff_timestamp: int, debug_timestamp: int) -> bytes:
    payload = bytearray(0x400)
    payload[:2] = b"MZ"
    struct.pack_into("<I", payload, 0x3C, 0x80)
    payload[0x80:0x84] = b"PE\0\0"
    coff = 0x84
    struct.pack_into("<H", payload, coff + 2, 1)
    struct.pack_into("<I", payload, coff + 4, coff_timestamp)
    struct.pack_into("<H", payload, coff + 16, 0xF0)
    optional = coff + 20
    struct.pack_into("<H", payload, optional, 0x20B)
    data_directory = optional + 112
    struct.pack_into("<I", payload, optional + 108, 16)
    struct.pack_into("<II", payload, data_directory + 6 * 8, 0x1000, 28)
    section = optional + 0xF0
    struct.pack_into("<IIII", payload, section + 8, 0x100, 0x1000, 0x100, 0x200)
    struct.pack_into("<I", payload, 0x204, debug_timestamp)
    payload[0x220] = 0x5A
    return bytes(payload)


def test_wheel_source_comparison_ignores_only_pe_build_timestamps():
    first = _synthetic_pe(1, 2)
    second = _synthetic_pe(3, 4)
    assert wheel_verifier._normalize_pe_timestamps(first) == (
        wheel_verifier._normalize_pe_timestamps(second)
    )
    tampered = bytearray(second)
    tampered[0x220] ^= 0x01
    assert wheel_verifier._normalize_pe_timestamps(first) != (
        wheel_verifier._normalize_pe_timestamps(bytes(tampered))
    )


def test_upstream_notices_and_trademark_attribution_are_checked():
    notices = (PACKAGE_ROOT / "THIRD_PARTY_NOTICES.md").read_text(
        encoding="utf-8"
    )
    assert "Copyright (c) 2018 Preferred Networks, Inc." in notices
    assert "Copyright (c) 2025 Preferred Networks, Inc." in notices
    assert (
        "Optuna, the Optuna logo and any related marks are trademarks of "
        "Preferred Networks, Inc."
    ) in notices
    assert "independent, unofficial project" in notices
    assert "not affiliated with, sponsored by, or endorsed by" in notices
    assert "does not use the Optuna logo" in notices
