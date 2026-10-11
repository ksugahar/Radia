from pathlib import Path
import importlib.util

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "check_matlab_python_parity", ROOT / "tools" / "check_matlab_python_parity.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_every_python_module_has_an_explicit_cpu_ownership_classification():
    result = MODULE.audit(ROOT)
    assert result["ok"], result["errors"]
    assert result["classified_file_count"] == result["python_file_count"] > 0
    assert all(item["classification"] == "cpu-python-only" and not item["matlab"]
               for item in result["assignments"])
    assert not any(item["python"].startswith(("simulink/", "ltspice/"))
                   for item in result["assignments"])
    assert all(item["classification"] == "cpu-python-only" for item in result["binary_extensions"])
    assert [item["id"] for item in result["gpu_features"]] == ["csr-pcg"]


def test_audit_rejects_reintroduced_cpu_mex_and_missing_gpu_source(tmp_path):
    import json
    target = tmp_path / "packages/radia-gpu/matlab/python_api_parity_manifest.json"
    target.parent.mkdir(parents=True)
    manifest = json.loads((ROOT / target.relative_to(tmp_path)).read_text())
    manifest["cpu_mex_supported"] = True
    manifest["gpu_features"][0]["source"] = "../outside.cu"
    target.write_text(json.dumps(manifest))
    (tmp_path / "src/radia").mkdir(parents=True)
    (tmp_path / "src/radia/demo.py").write_text("pass")
    result = MODULE.audit(tmp_path)
    assert not result["ok"]
    assert any("CPU MEX" in error for error in result["errors"])
    assert any("escaping source" in error for error in result["errors"])
