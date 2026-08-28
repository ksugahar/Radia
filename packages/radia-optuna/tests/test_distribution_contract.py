from __future__ import annotations

import json
import re
import tomllib
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PACKAGE_ROOT.parents[1]


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
    assert "source_fidelity_verified" in verifier
    assert "wheel payload differs from the checked monorepo source" in verifier


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
