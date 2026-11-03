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
        "reusable-gpu-ic0", "laplace-fmm-operators", "bounded-hdiv-mmm", "device-hdiv-bh",
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


def _bitwise_module():
    spec = importlib.util.spec_from_file_location(
        "bitwise_parity_contract", ROOT / "validation_test/hmatrix_cuda/check_pybind_mex_bitwise.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_hdiv_frontend_parity_status_is_bound_to_its_evidence():
    """Guards: same-kernel parity is claimed only for the executed case, with hashes for that case.

    Assumes: source contracts execute neither GPU frontend, so they check the record, not the run.
    """
    import json
    module = _bitwise_module()
    manifest = json.loads((ROOT / "packages/radia-simulink/matlab/python_api_parity_manifest.json").read_text())
    contracts = [feature["same_kernel_bitwise_parity"] for feature in manifest["gpu_features"]
                 if "same_kernel_bitwise_parity" in feature]
    assert contracts, manifest["gpu_features"]
    for contract in contracts:
        assert module.check_contract(contract) in ("PENDING", "PASS"), contract
        assert (ROOT / contract["protocol"]).is_file(), contract


def test_bitwise_contract_refuses_a_pass_beyond_its_evidence():
    """Guards: a PASS cannot name other outputs or another case, or omit a hash.

    Assumes: the evidence field list is the comparator's own.
    """
    import pytest
    module = _bitwise_module()
    case = "device-flat-tet"
    evidence = {key: "a"*(40 if key == "source_sha" else 64) for key in module.EVIDENCE_HASHES}
    evidence.update(case=case, outputs=sorted(module.CASES[case]), cuda_compiler="release 12.8", architecture="sm_75")
    good = dict(status="PASS", case=case, evidence=evidence)
    assert module.check_contract(good) == "PASS"
    assert module.check_contract(dict(status="PENDING", case=case)) == "PENDING"
    bad = [dict(status="PENDING", case=case, evidence=evidence), dict(status="PASS", case=case),
           dict(good, status="UNKNOWN"), dict(good, case="device-hex"),
           dict(good, evidence=dict(evidence, outputs=["apply"])),
           dict(good, evidence=dict(evidence, source_sha="short")),
           dict(good, evidence=dict(evidence, input_sha256="G"*64)),
           dict(good, evidence={key: value for key, value in evidence.items() if key != "cuda_compiler"})]
    for contract in bad:
        with pytest.raises(ValueError):
            module.check_contract(contract)


def test_bitwise_comparator_requires_one_kernel_object_and_the_loaded_binaries():
    """Guards: matching outputs cannot pass across different kernels, inputs or loaded binaries.

    Assumes: this is source-only validation, not actual compiler/GPU qualification.
    """
    import pytest
    spec = importlib.util.spec_from_file_location(
        "bitwise_parity_receipts", ROOT / "validation_test/hmatrix_cuda/check_pybind_mex_bitwise.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    shared = dict(source_sha="a"*40, kernel_object_sha256="b"*64, cuda_compiler="release 12.8, V12.8.61",
                  python_binary_sha256="c"*64, mex_binary_sha256="d"*64, input_sha256="e"*64,
                  device_uuid="GPU-0", architecture="sm_75", case="device-flat-tet",
                  output_names=sorted(module.CASES["device-flat-tet"]))
    python = dict(shared, frontend="pybind11", artifact_sha256="c"*64)
    mex = dict(shared, frontend="mex", artifact_sha256="d"*64)
    assert module.check_receipts(python, mex) == module.CASES["device-flat-tet"]
    for side, key, value, message in (
            ("mex", "kernel_object_sha256", "f"*64, "kernel_object_sha256"),
            ("mex", "input_sha256", "f"*64, "input_sha256"),
            ("mex", "frontend", "mexcuda", "pybind11 and mex"),
            ("python", "artifact_sha256", "f"*64, "loaded frontend binary"),
            ("mex", "artifact_sha256", "c"*64, "loaded frontend binary")):
        a, b = dict(python), dict(mex)
        (a if side == "python" else b)[key] = value
        with pytest.raises(ValueError, match=message):
            module.check_receipts(a, b)
    with pytest.raises(ValueError, match="Expected outputs"):
        module.check_receipts(dict(python, output_names=["apply"]), dict(mex, output_names=["apply"]))
