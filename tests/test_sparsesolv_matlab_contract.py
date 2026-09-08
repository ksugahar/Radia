"""Keep the formerly omitted compiled module in the parity inventory."""
import json
import yaml
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_sparsesolv_has_native_and_explicit_fallback_owners():
    manifest = json.loads((ROOT / "matlab/python_api_parity_manifest.json").read_text())
    entry = next(e for e in manifest["binary_extensions"] if e["python"] == "sparsesolv_ngsolve.pyd")
    assert entry["status"] == "focused-native-commands"
    assert (ROOT / entry["fallback"]).is_file()
    source = (ROOT / "src/matlab/radia_mex.cpp").read_text()
    for command in entry["native_commands"]:
        assert f'command == "{command}"' in source
    for owner in entry["matlab"].split("; "):
        assert (ROOT / owner).is_file()


def test_ams_setup_uses_ngsolve_gradient_without_taskmanager():
    source = (ROOT / "src/matlab/radia_mex.cpp").read_text()
    body = source.split("void SparseSolvAMS(")[1].split("void SparseSolvIC(")[0]
    assert "hc->CreateGradient()" in body
    assert "RegionTaskManager" not in body
    assert "a.fespace != space.fespace" in body


def test_matlab_lane_runs_engine_and_retains_json_on_mdx():
    workflow = yaml.safe_load((ROOT/".github/workflows/sparsesolv.yml").read_text())
    job = workflow["jobs"]["ams-regression"]
    assert "mdx" in job["runs-on"]
    step = next(s for s in job["steps"] if s.get("name") == "Build and verify native MATLAB parity")
    assert step["if"] == "steps.matlab-impact.outputs.required == 'true'"
    assert "-MatlabMexOnly" in step["run"]
    assert "run_sparsesolv_parity.py --output" in step["run"]
    assert "sparsesolv-matlab.json" in job["steps"][-1]["with"]["path"]
