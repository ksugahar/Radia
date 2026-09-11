from __future__ import annotations

import importlib.util
import hashlib
import ast
import json
import re
import tomllib
from pathlib import Path
import sys

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PACKAGE_ROOT.parents[1]
sys.path.insert(0, str(PACKAGE_ROOT / "src"))


VERIFY_SPEC = importlib.util.spec_from_file_location(
    "radia_optuna_verify_wheel", PACKAGE_ROOT / "verify_wheel.py"
)
assert VERIFY_SPEC is not None and VERIFY_SPEC.loader is not None
verify_wheel = importlib.util.module_from_spec(VERIFY_SPEC)
VERIFY_SPEC.loader.exec_module(verify_wheel)


def test_api_coverage_is_reproducible_and_names_real_matlab_packages(monkeypatch, tmp_path):
    path = REPO_ROOT / "tests/matlab/fixtures/generate_optuna50_api_coverage.py"
    spec = importlib.util.spec_from_file_location("optuna_coverage_audit", path)
    assert spec is not None and spec.loader is not None
    generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(generator)
    committed = json.loads(generator.DESTINATION.read_text(encoding="utf-8"))
    assert committed["upstream_oracle_sha256"].lower() == hashlib.sha256(
        generator.ORACLE_PATH.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    ).hexdigest()
    regenerated = json.loads(json.dumps(generator.build_coverage()))
    assert committed == regenerated, "Regenerate Optuna API coverage."
    serialized = json.dumps(generator.build_coverage(), indent=2, sort_keys=True) + "\n"
    assert generator.DESTINATION.read_bytes().replace(b"\r\n", b"\n") == serialized.encode()
    unrelated = ast.parse('server.stop(0); other.params; names = "stop"')
    assert not generator._qualified_optuna_references(unrelated)
    actual = ast.parse('optuna.study.Study.stop(study)')
    assert "optuna.study.Study.stop" in generator._qualified_optuna_references(actual)
    assert not generator._oracle_sections_for("optuna.study.Study.stop")
    assert not committed["full_compatibility_complete"]
    entries = [e for e in committed["entries"] if e["kind"] != "module"]
    for entry in entries:
        name = entry["matlab_name"]
        if name is None:
            continue
        parts = name.split(".")
        if entry["kind"].startswith("class-"):
            parts = parts[:-1]
        source = REPO_ROOT / "matlab"
        for package in parts[:-1]:
            source /= "+" + package
        assert (source / (parts[-1] + ".m")).is_file(), name
    # Public package scan ignores private helpers and recognizes @Class layout.
    private = tmp_path / "private"
    private.mkdir()
    (private / "hidden.m").write_text("function hidden()\nend\n", encoding="utf-8")
    legacy = tmp_path / "@Legacy"
    legacy.mkdir()
    (legacy / "Legacy.m").write_text("classdef Legacy\nend\n", encoding="utf-8")
    (legacy / "step.m").write_text("function step(obj)\nend\n", encoding="utf-8")
    monkeypatch.setattr(generator, "MATLAB_DIRECTORY", tmp_path)
    assert generator._matlab_qualified_names() == {"Legacy": "radia.optuna.Legacy"}
    names, members = generator._matlab_surface()
    assert names == {"Legacy"}
    assert "step" in members["Legacy"]


def test_test_manifest_matches_its_generator():
    path = REPO_ROOT / "tests/matlab/fixtures/generate_optuna_test_manifest.py"
    spec = importlib.util.spec_from_file_location("optuna_test_manifest_audit", path)
    assert spec is not None and spec.loader is not None
    generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(generator)
    expected = json.dumps(generator.build_manifest(), indent=2, sort_keys=True) + "\n"
    actual = path.with_name("optuna_test_manifest.json").read_bytes().replace(b"\r\n", b"\n")
    assert actual == expected.encode()


def test_constructor_defaults_are_audited_against_the_pinned_oracle():
    """Every upstream constructor default is compared, or declared, or open.

    The API inventory records classes and their public members but not
    __init__, so before this audit a changed constructor default was invisible
    to the coverage gate while changing every seeded result. Optuna 5.0 is
    that case: TPESampler moved multivariate from False to None and
    constant_liar from False to True, and only the second had a named test.
    """
    oracle = json.loads(
        (REPO_ROOT / "tests/matlab/fixtures/optuna50_oracle.json").read_text(
            encoding="utf-8"
        )
    )
    record = oracle["constructor_defaults"]
    recorded = sum(
        len(parameters)
        for classes in record["modules"].values()
        for parameters in classes.values()
    )
    assert recorded == record["parameter_count"]
    assert record["parameter_count"] > 0

    coverage = json.loads(
        (REPO_ROOT / "matlab/optuna50_api_coverage.json").read_text(encoding="utf-8")
    )
    audit = coverage["constructor_default_audit"]

    # Every audited parameter lands in exactly one bucket.
    accounted = (
        audit["matched_count"]
        + audit["declared_equivalent_count"]
        + audit["declared_not_implemented_count"]
        + audit["declared_divergent_count"]
    )
    assert accounted > 0
    assert len(audit["declared_equivalences"]) == audit["declared_equivalent_count"]
    assert (
        len(audit["declared_not_implemented"])
        == audit["declared_not_implemented_count"]
    )
    assert len(audit["declared_divergences"]) == audit["declared_divergent_count"]

    # Every declaration states a reason; an empty one is not a declaration.
    for table in (
        "declared_equivalences",
        "declared_not_implemented",
        "declared_divergences",
    ):
        for name, reason in audit[table].items():
            assert reason.strip(), f"{table}[{name}] declares no reason"

    # A divergence is a recorded defect, not an equivalence. A new one must
    # fail here instead of being absorbed into a table; fixing one only
    # shrinks the set, and the generator then rejects the stale declaration.
    assert set(audit["declared_divergences"]) <= {
        "FanovaImportanceEvaluator.seed",
        "MeanDecreaseImpurityImportanceEvaluator.seed",
        "Terminator.improvement_evaluator",
    }


def test_matlab_path_names_the_layout_it_resolved():
    import radia_optuna
    from radia_optuna import cli

    resolved = radia_optuna.layout()
    assert resolved in {"wheel", "checkout"}
    root = radia_optuna.matlab_path()
    assert (root / "+radia" / "+optuna").is_dir()
    assert radia_optuna.mex_path().name == "optuna_mex.mexw64"
    payload = cli._doctor_payload()
    assert payload["layout"] == resolved
    assert payload["matlab_file_count"] == payload["expected_matlab_file_count"]


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
    assert "optuna50_api_coverage.json" in setup_source
    assert "THIRD_PARTY_NOTICES.md" in setup_source
    assert 'return "py3", "none", "win_amd64"' in setup_source
    assert '"bdist_wheel": bdist_wheel' in setup_source


def test_focused_optuna_build_keeps_python_and_pybind_out_of_the_build_contract():
    build = (REPO_ROOT / "Build.ps1").read_text(encoding="utf-8")
    optuna_guard = build.index("if (-not $OptunaMexOnly) {")
    pybind_probe = build.index("import pybind11; print(pybind11.get_cmake_dir())")
    assert optuna_guard < pybind_probe
    assert "$PythonCMakeArgs = if ($OptunaMexOnly)" in build
    assert "$NINJA_COMMAND = Get-Command ninja" in build
    assert "-InstallToSitePackages does not apply" in build
    assert '-m pytest "$PROJECT_DIR\\packages\\radia-optuna\\tests" -q' in build


def test_wheel_verifier_rejects_solver_boundary_leaks():
    verifier = (PACKAGE_ROOT / "verify_wheel.py").read_text(encoding="utf-8")
    for forbidden in ("radia_mex", "ngsolve", "netgen", "mkl_", "radia_pybind"):
        assert forbidden in verifier
    assert "py3-none-win_amd64" in verifier
    assert "optuna\\s*==\\s*5\\.0\\.0" in verifier
    assert "THIRD_PARTY_NOTICES.md" in verifier
    assert "source_fidelity_verified" in verifier
    assert "wheel payload differs from the checked monorepo source" in verifier


def test_release_candidate_fidelity_normalizes_text_but_not_semantics():
    member = "radia_optuna/matlab/optuna_upstream_compatibility.json"
    assert verify_wheel._payload_matches(
        member,
        b'{\r\n  "ok": true\r\n}\r\n',
        b'{\n  "ok": true\n}\n',
        release_candidate=True,
    )
    assert not verify_wheel._payload_matches(
        member,
        b'{\r\n  "ok": false\r\n}\r\n',
        b'{\n  "ok": true\n}\n',
        release_candidate=True,
    )
    assert verify_wheel._payload_matches(
        "radia_optuna/matlab/optuna_mex.mexw64",
        b"different-ci-build",
        b"local-build",
        release_candidate=True,
    )
    assert verify_wheel._payload_matches(
        "radia_optuna/matlab/LICENSE",
        b"license text\r\n",
        b"license text\n",
        release_candidate=True,
    )
    mex_member = "radia_optuna/matlab/optuna_mex.mexw64"
    assert not verify_wheel._requires_source_payload(
        mex_member, release_candidate=True
    )
    assert verify_wheel._requires_source_payload(
        mex_member, release_candidate=False
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
