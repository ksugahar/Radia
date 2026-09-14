import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from validation_test.cubit import run_vol_multi_geometry_validation as runner


def _invoke(monkeypatch, tmp_path, child_behavior):
    plugin_dir = tmp_path / "plugin"
    plugin_dir.mkdir()
    ccm_path = plugin_dir / "cubit_mesh_export.ccm"
    ccm_path.write_bytes(b"tested-ccm")
    results_path = tmp_path / "custom" / "evidence.json"
    log_path = tmp_path / "run.log"

    def fake_run(command, *, cwd, env, **kwargs):
        return child_behavior(command, cwd, env, ccm_path, results_path)

    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    monkeypatch.setattr(
        runner.sys,
        "argv",
        [
            "runner",
            "--plugin-dir", str(plugin_dir),
            "--log", str(log_path),
            "--results", str(results_path),
        ],
    )
    return runner.main(), results_path, log_path


def _fresh_evidence(env, ccm_path):
    diagnostic = "Same-material internal surface 6 is explicitly labelled 'intentional_cut'"
    return {
        "run_id": env["CUBIT_VALIDATION_RUN_ID"],
        "all_passed": True,
        "native_provenance": {
            "payloads": {
                "cubit_mesh_export.ccm": {
                    "sha256": runner._sha256(ccm_path),
                    "size": ccm_path.stat().st_size,
                }
            }
        },
        "labelled_same_material_internal_surface": {
            "expected_diagnostic": diagnostic,
        },
    }


def test_custom_result_path_requires_fresh_run_and_attests_loaded_ccm(monkeypatch, tmp_path):
    def child(command, cwd, env, ccm_path, expected_results):
        assert Path(env["CUBIT_VALIDATION_RESULTS"]) == expected_results.resolve()
        evidence = _fresh_evidence(env, ccm_path)
        expected_results.parent.mkdir(parents=True, exist_ok=True)
        expected_results.write_text(json.dumps(evidence), encoding="utf-8")
        stdout = (
            f"Loading Plugin: '{ccm_path.resolve()}'\n"
            f"{evidence['labelled_same_material_internal_surface']['expected_diagnostic']}\n"
        )
        return SimpleNamespace(returncode=0, stdout=stdout)

    code, results_path, log_path = _invoke(monkeypatch, tmp_path, child)
    assert code == 0
    result = json.loads(results_path.read_text(encoding="utf-8"))
    assert result["loaded_ccm"]["observed_in_stdout"] is True
    assert result["loaded_ccm"]["matches_manifest"] is True
    assert result["loaded_ccm"]["stable_during_run"] is True
    assert result["labelled_same_material_internal_surface"]["diagnostic_matched_in_stdout"] is True
    assert log_path.is_file()


@pytest.mark.parametrize("returncode", [1, 7])
def test_child_nonzero_fails_even_if_result_exists(monkeypatch, tmp_path, returncode):
    def child(command, cwd, env, ccm_path, expected_results):
        expected_results.parent.mkdir(parents=True, exist_ok=True)
        expected_results.write_text(json.dumps(_fresh_evidence(env, ccm_path)), encoding="utf-8")
        return SimpleNamespace(returncode=returncode, stdout="child failed\n")

    code, _, _ = _invoke(monkeypatch, tmp_path, child)
    assert code == 1


def test_missing_result_fails(monkeypatch, tmp_path):
    def child(command, cwd, env, ccm_path, expected_results):
        return SimpleNamespace(returncode=0, stdout="no result\n")

    code, results_path, _ = _invoke(monkeypatch, tmp_path, child)
    assert code == 1
    assert not results_path.exists()


def test_stale_result_run_id_fails(monkeypatch, tmp_path):
    def child(command, cwd, env, ccm_path, expected_results):
        evidence = _fresh_evidence(env, ccm_path)
        evidence["run_id"] = "stale-run"
        expected_results.parent.mkdir(parents=True, exist_ok=True)
        expected_results.write_text(json.dumps(evidence), encoding="utf-8")
        return SimpleNamespace(returncode=0, stdout="stale result\n")

    code, _, _ = _invoke(monkeypatch, tmp_path, child)
    assert code == 1


def test_unrelated_export_failure_diagnostic_fails(monkeypatch, tmp_path):
    def child(command, cwd, env, ccm_path, expected_results):
        evidence = _fresh_evidence(env, ccm_path)
        expected_results.write_text(json.dumps(evidence), encoding="utf-8")
        stdout = f"Loading Plugin: '{ccm_path.resolve()}'\nERROR: unrelated failure\n"
        return SimpleNamespace(returncode=0, stdout=stdout)

    code, _, _ = _invoke(monkeypatch, tmp_path, child)
    assert code == 1


def test_ccm_changed_during_child_run_fails(monkeypatch, tmp_path):
    def child(command, cwd, env, ccm_path, expected_results):
        evidence = _fresh_evidence(env, ccm_path)
        expected_results.write_text(json.dumps(evidence), encoding="utf-8")
        ccm_path.write_bytes(b"replaced-during-run")
        stdout = (
            f"Loading Plugin: '{ccm_path.resolve()}'\n"
            f"{evidence['labelled_same_material_internal_surface']['expected_diagnostic']}\n"
        )
        return SimpleNamespace(returncode=0, stdout=stdout)

    code, results_path, _ = _invoke(monkeypatch, tmp_path, child)
    assert code == 1
    result = json.loads(results_path.read_text(encoding="utf-8"))
    assert result["loaded_ccm"]["stable_during_run"] is False
