"""Guards: the device charge-Gram route takes primitives only and reports what it uploads.

Assumes: no GPU, CUDA context or native artifact. Static source scans and
host-side mapping checks; CPU/GPU Gram agreement is a separate validation run.
"""
import importlib.util
from pathlib import Path
import re

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
HDIV = ROOT / "src/core/gpu/hdiv"
_spec = importlib.util.spec_from_file_location("hdiv_device_gram_contract", ROOT / "src/radia/hdiv_gpu.py")
hdiv_gpu = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(hdiv_gpu)


def staged(far=True, n_el=1):
    inputs = dict(cell_verts=np.arange(12.), face_verts=np.arange(9.), n_el=n_el,
                  charge_host=[0, 0], charge_kind=[0, 1], charge_expo=[0, 0, 0, 1, 0, 0],
                  ref_tet_pts=np.full(6, .25), ref_tet_w=[.1, .2],
                  ref_tri_pts=np.full(2, 1/3), ref_tri_w=[.5],
                  eps=1e-10, leaf=64, eta=2.)
    if far:
        inputs.update(ref_tet_pts_lo=np.full(3, .25), ref_tet_w_lo=[1/6],
                      ref_tri_pts_lo=np.full(2, 1/3), ref_tri_w_lo=[.5], ho_far_factor=2.)
    return inputs


def tables(inputs):
    """Staged charge tables consistent with ``inputs`` (one cell, one face monomial)."""
    expo = np.asarray(inputs["charge_expo"], dtype=np.int32)
    return (dict(cell_class=np.zeros(1, np.int32), cell_dofs=np.arange(6, dtype=np.int32)[None],
                 cell_tables=np.ones((1, 1, 6)), cell_exponents=expo[None, :3],
                 face_class=np.zeros(1, np.int32), face_dofs=np.arange(3, dtype=np.int32)[None],
                 face_tables=np.ones((1, 1, 3)), face_exponents=expo[None, 3:5]),
            (np.asarray(inputs["cell_verts"], dtype=np.float64), np.asarray(inputs["face_verts"], dtype=np.float64)))


def gram(inputs):
    charge, geometry = tables(inputs)
    return hdiv_gpu._device_gram_arguments(charge, inputs, geometry)


def test_staged_inputs_map_to_primitive_rows_with_explicit_aca_parameters():
    arguments = gram(staged())
    assert arguments["compression"] == "aca" and arguments["aca"] == hdiv_gpu._GRAM_ACA, arguments["aca"]
    assert arguments["aca"] is not hdiv_gpu._GRAM_ACA, "callers must not share the module constant"
    assert arguments["tet_rule"].shape == (2, 4) and arguments["tri_rule"].shape == (1, 4)
    assert arguments["tri_rule"][0, 2] == 0 and arguments["tet_low_rule"][0, 3] == pytest.approx(1/6)
    assert arguments["far_factor"] == 2.
    assert {k: str(v.dtype) for k, v in arguments.items() if k.startswith("charge_")} == {
        "charge_host": "int32", "charge_kind": "int32", "charge_exponents": "int32"}
    assert not {"eps", "leaf", "eta", "leaves", "coefficients"} & set(arguments), sorted(arguments)


def test_missing_low_rules_disable_the_far_split():
    arguments = gram(staged(far=False))
    assert arguments["far_factor"] == 1e30 and arguments["tet_low_rule"].shape == (0, 4)


@pytest.mark.parametrize("change,error", [
    (dict(n_el=2), ValueError), (dict(image_masks=[1]), NotImplementedError)])
def test_inconsistent_or_image_payloads_fail_before_the_device(change, error):
    inputs = staged()
    inputs.update(change)
    with pytest.raises(error):
        gram(inputs)


def ledger(arguments, *extra, drop=()):
    rows = [dict(payload=p, direction=d, bytes=8 if b is None else b)
            for (p, d), b in hdiv_gpu._expected_device_ledger(arguments).items() if p not in drop]
    return dict(route="device-flat-tet", charges=int(arguments["charge_host"].size),
                compression=arguments["compression"], hmatrix={k: 0 for k in hdiv_gpu._ACA_COUNTERS},
                transfer_ledger=rows + [dict(payload=p, direction=d, bytes=b) for p, d, b in extra])


def test_assembly_ledger_requires_the_exact_primitive_rows():
    arguments = gram(staged())
    assert hdiv_gpu._checked_device_assembly(ledger(arguments), arguments)["charges"] == 2
    flags_only = dict(route="device-flat-tet", charges=2, compression="aca",
                      hmatrix={k: 0 for k in hdiv_gpu._ACA_COUNTERS}, transfer_ledger=[
        dict(payload="prepare_flag", direction="D2H", bytes=4),
        dict(payload="energy_flag", direction="D2H", bytes=4)])
    duplicate = ledger(arguments, ("tet_rule", "H2D", 64))
    resized = ledger(arguments)
    resized["transfer_ledger"][0]["bytes"] += 8
    forged = [flags_only, duplicate, resized, ledger(arguments, drop=("resident_plan",)),
              ledger(arguments, ("host_hmatrix_coefficients", "H2D", 8)),
              ledger(arguments, ("cell_vertices_m", "D2H", 96)), dict(ledger(arguments), route="host-leaves"),
              ledger(arguments, ("charge_map_csr_host_built", "H2D", 64)), ledger(arguments, drop=("charge_cell_tables",)),
              dict(ledger(arguments), compression="dense"), dict(ledger(arguments), hmatrix={}),
              dict(ledger(arguments), hmatrix=dict({k: 0 for k in hdiv_gpu._ACA_COUNTERS}, rank_sum=-1)),
              ledger(arguments, drop=("aca_permutation",)), ledger(arguments, ("aca_lowrank_factors", "D2H", 64)),
              ledger(arguments, ("aca_leaf_status", "D2H", 0)),
              ledger(arguments, ("aca_leaf_status", "D2H", 8), ("aca_leaf_status", "D2H", 8))]
    for bad in ("8", -8, True, 0):
        row = ledger(arguments)
        next(r for r in row["transfer_ledger"] if r["payload"] == "resident_plan")["bytes"] = bad
        forged.append(row)
    for case in forged:
        with pytest.raises(RuntimeError):
            hdiv_gpu._checked_device_assembly(case, arguments)


def test_aca_metadata_rows_are_optional_but_positive():
    arguments = gram(staged())
    assert hdiv_gpu._checked_device_assembly(ledger(arguments, *[(p, d, 8) for p, d in hdiv_gpu._ACA_OPTIONAL_ROWS]),
                                             arguments)["compression"] == "aca"
    assert ("aca_permutation", "H2D") in hdiv_gpu._expected_device_ledger(arguments)
    assert hdiv_gpu._expected_device_ledger(arguments)[("aca_permutation", "H2D")] == 8*2


def test_disabled_far_split_has_no_low_rule_rows():
    arguments = gram(staged(far=False))
    payloads = {p for p, _ in hdiv_gpu._expected_device_ledger(arguments)}
    assert not {"tet_low_rule", "tri_low_rule"} & payloads, sorted(payloads)
    hdiv_gpu._checked_device_assembly(ledger(arguments), arguments)
    with pytest.raises(RuntimeError, match="rejected"):
        hdiv_gpu._checked_device_assembly(ledger(arguments, ("tet_low_rule", "H2D", 0)), arguments)


@pytest.mark.parametrize("expo", [
    [2**31 - 1, 2**31 - 1, 0], [2**31 - 1, 0, 0], [-1, 0, 0], [0, -2**31, 0], [1, 1, 0], [0, 0, 2]])
def test_cell_exponents_are_bounded_before_any_sum(expo):
    inputs = staged()
    inputs["charge_expo"] = expo + [0, 0, 0]
    with pytest.raises(ValueError, match="cell degree <= 1"):
        gram(inputs)


@pytest.mark.parametrize("expo,ok", [
    ([2, 0, 0], True), ([1, 1, 0], True), ([0, 0, 1], False), ([2, 1, 0], False), ([3, 0, 0], False)])
def test_face_exponents_accept_the_degree_two_border_only(expo, ok):
    inputs = staged()
    inputs["charge_expo"] = [1, 0, 0] + expo
    if ok:
        assert gram(inputs)["charge_exponents"][3:].tolist() == expo
    else:
        with pytest.raises(ValueError, match="face degree <= 2"):
            gram(inputs)


def test_native_guards_bound_components_before_the_degree_sum():
    source = (HDIV / "hdiv_flat_gram.cuh").read_text(encoding="utf-8")
    kernel = source[source.index("__global__ void prepareFlatCharges"):source.index("__global__ void assembleFlatGram")]
    assert kernel.index("q.expo[2]>limit") < kernel.index("q.degree=q.expo[0]+q.expo[1]+q.expo[2]"), "device sums before bounding"
    validate = source[source.index("void validate()const"):source.index("__device__ inline double flatPow")]
    assert "expo[3*a+c]>limit" in validate, "host FlatBlock::validate lacks per-component exponent bounds"


def test_device_gram_source_downloads_only_validity_flags():
    source = (HDIV / "hdiv_flat_gram.cuh").read_text(encoding="utf-8")
    downloads = re.findall(r"cudaMemcpy\w*\(([^;]*)DeviceToHost", source)
    assert len(downloads) == 1 and "sizeof(int)" in downloads[0], downloads
    assert "cudaMemcpyToSymbol" not in source and "cudaMemcpyHostToDevice" not in source


def test_device_aca_source_downloads_only_integer_metadata():
    source = (HDIV / "hdiv_flat_aca.cuh").read_text(encoding="utf-8")
    downloads = re.findall(r"cudaMemcpy\w*\(([^;]*)DeviceToHost", source)
    assert len(downloads) == 2 and all("sizeof(int)" in d for d in downloads), downloads
    assert {d.split(",")[1] for d in downloads} == {"active.data", "report.data"}, downloads
    assert re.search(r"Buffer<int> active;", source) and re.search(r"Buffer<int> report;", source)
    report = source[source.index("__global__ void acaReport"):source.index("struct AcaDenseSpan")]
    assert "frobenius2" not in report and "maxAbs" not in report, "ACA report leaks norms"
    assert "cudaMemcpyHostToDevice" not in source and "cudaMemcpyToSymbol" not in source


def test_flat_binding_accepts_no_matrix_payload():
    source = (HDIV / "hdiv_resident_python.cpp").read_text(encoding="utf-8")
    start = source.index('def_static("from_flat_charges"')
    block = source[start:source.index('.def("assembly"', start)]
    names = re.findall(r'py::arg\("([a-z_]+)"\)', block)
    assert names[:10] == ["cell_vertices", "face_vertices", "charge_host", "charge_kind",
                          "charge_exponents", "tet_rule", "tri_rule", "tet_low_rule",
                          "tri_low_rule", "far_factor"], names
    assert not {"leaves", "permutation", "scale", "a1", "a2"} & set(names), names
