"""Static native provenance routing contracts run without native builds."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_provenance_helper_and_beam_contracts_are_in_ci_scope():
    import yaml

    helper = "tools/native_build_provenance.ps1"
    for name, expected in (("sparsesolv.yml", True), ("radia-optuna.yml", False)):
        workflow = yaml.safe_load((ROOT / ".github/workflows" / name).read_text())
        triggers = workflow.get("on", workflow.get(True))
        for event in ("push", "pull_request"):
            # -OptunaMexOnly returns before Build.ps1 loads the helper.
            assert (helper in triggers[event]["paths"]) is expected
    manifest = json.loads((ROOT / "tests/test_tier_manifest.json").read_text())
    assert manifest["profiles"]["native-build-contracts"]["paths"] == [
        "tests/test_native_build_provenance.py"
    ]
    workflow = yaml.safe_load((ROOT / ".github/workflows/sparsesolv.yml").read_text())
    triggers = workflow.get("on", workflow.get(True))
    for event in ("push", "pull_request"):
        assert "tests/test_native_build_provenance.py" in triggers[event]["paths"]
    steps = workflow["jobs"]["ams-regression"]["steps"]
    contract = next(step for step in steps if step.get("name") == "Verify native build provenance contracts")
    assert "--profile native-build-contracts" in contract["run"]
    assert "exit $LASTEXITCODE" in contract["run"]
    rules = manifest["impact_rules"]
    for name in (
        "benchmark_beam_transfer_mex_python.py",
        "benchmark_beam_transfer_mex_matlab.m",
        "beam_transfer_benchmark_case.json",
    ):
        assert "tests/test_beam_transfer_mex_pybind_validation.py" in rules[
            f"validation_test/accelerator/{name}"
        ]

