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
    assert [item["id"] for item in result["gpu_features"]] == [
        "csr-pcg", "ams-experimental", "native-hcurl-newton-ams", "reusable-linear-hcurl-ams",
        "reusable-gpu-ic0", "laplace-fmm-operators", "bounded-hdiv-mmm",
        "normalized-leaf-resident-boundary", "explicit-ng-host-hdiv-bh",
    ]


def test_audit_rejects_reintroduced_cpu_mex_and_missing_gpu_source(tmp_path):
    import json
    target = tmp_path / "packages/radia-simulink/matlab/python_api_parity_manifest.json"
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


def test_shared_radia_gpu_source_is_allowed_but_arbitrary_sources_are_rejected(tmp_path):
    """Guards: source relocation preserves GPU ownership and the MATLAB boundary.

    Assumes: path ownership does not establish numerical parity.
    """
    import json
    manifest_path = tmp_path / "packages/radia-simulink/matlab/python_api_parity_manifest.json"
    manifest_path.parent.mkdir(parents=True)
    manifest = json.loads((ROOT / manifest_path.relative_to(tmp_path)).read_text())
    feature = dict(manifest["gpu_features"][-1])
    manifest["gpu_features"] = [feature]
    (tmp_path / "src/radia").mkdir(parents=True)
    (tmp_path / "src/radia/demo.py").write_text("pass")
    for field in ("matlab", "validation"):
        target = tmp_path / feature[field]
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("independent path fixture")
    for source, allowed in (("src/core/gpu/hdiv/shared.cuh", True),
                            ("src/radia/not_gpu.cpp", False),
                            ("packages/radia-simulink/not_native.cpp", False)):
        target = tmp_path / source
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("independent source fixture")
        feature["source"] = source
        manifest_path.write_text(json.dumps(manifest))
        result = MODULE.audit(tmp_path)
        assert result["ok"] is allowed, (source, result["errors"])
    feature["source"] = "src/core/gpu/hdiv/shared.cuh"
    feature["matlab"] = "src/radia/not_gpu.cpp"
    manifest_path.write_text(json.dumps(manifest))
    assert any("matlab is not GPU-owned" in error for error in MODULE.audit(tmp_path)["errors"])


def test_hdiv_frontend_parity_is_pending_with_an_executable_comparator():
    """Guards: ownership moves cannot silently claim same-kernel MEX parity.

    Assumes: source contracts do not execute either GPU frontend.
    """
    import json
    manifest = json.loads((ROOT / "packages/radia-simulink/matlab/python_api_parity_manifest.json").read_text())
    for feature in manifest["gpu_features"]:
        if feature["id"] not in {"normalized-leaf-resident-boundary", "explicit-ng-host-hdiv-bh"}:
            continue
        contract = feature["same_kernel_bitwise_parity"]
        assert contract["status"] == "PENDING", contract
        assert contract["case"] == feature["id"]
        assert (ROOT / contract["protocol"]).is_file()
        spec = importlib.util.spec_from_file_location("bitwise_parity", ROOT / contract["test_reference"])
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        assert callable(module.compare)


def test_bitwise_comparator_rejects_conflicting_compiler_overrides():
    """Guards: matching receipts cannot hide contradictory compiler controls.

    Assumes: this is source-only validation, not actual compiler/GPU qualification.
    """
    import pytest
    spec = importlib.util.spec_from_file_location(
        "bitwise_parity_flags", ROOT / "validation_test/hmatrix_cuda/check_pybind_mex_bitwise.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    flags = ["--fmad=false", "--ftz=false", "--prec-div=true", "--prec-sqrt=true"]
    module.validate_deterministic_flags(flags)
    module.validate_deterministic_flags(list(reversed(flags)))
    for override in ("--fmad=true", "--ftz=true", "--prec-div=false",
                     "--prec-sqrt=false", "--use_fast_math", "-fmad=true", flags[0]):
        with pytest.raises(ValueError, match="overrides forbidden"):
            module.validate_deterministic_flags(flags + [override])
    with pytest.raises(ValueError):
        module.validate_deterministic_flags(flags[:-1])
