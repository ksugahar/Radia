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


class Map:
    shape = (2, 7)


def test_staged_inputs_map_to_primitive_rows_without_hmatrix_parameters():
    arguments = hdiv_gpu._device_gram_arguments(Map, staged())
    assert arguments["tet_rule"].shape == (2, 4) and arguments["tri_rule"].shape == (1, 4)
    assert arguments["tri_rule"][0, 2] == 0 and arguments["tet_low_rule"][0, 3] == pytest.approx(1/6)
    assert arguments["far_factor"] == 2.
    assert {k: str(v.dtype) for k, v in arguments.items() if k.startswith("charge_")} == {
        "charge_host": "int32", "charge_kind": "int32", "charge_exponents": "int32"}
    assert not {"eps", "leaf", "eta", "leaves", "coefficients"} & set(arguments), sorted(arguments)


def test_missing_low_rules_disable_the_far_split():
    arguments = hdiv_gpu._device_gram_arguments(Map, staged(far=False))
    assert arguments["far_factor"] == 1e30 and arguments["tet_low_rule"].shape == (0, 4)


@pytest.mark.parametrize("change,error", [
    (dict(n_el=2), ValueError), (dict(image_masks=[1]), NotImplementedError)])
def test_inconsistent_or_image_payloads_fail_before_the_device(change, error):
    inputs = staged()
    inputs.update(change)
    with pytest.raises(error):
        hdiv_gpu._device_gram_arguments(Map, inputs)


def ledger(*extra):
    rows = [(name, "H2D") for name in hdiv_gpu._DEVICE_GRAM_UPLOADS]
    rows += [("resident_plan", "H2D"), ("charge_map_csr_host_built", "H2D"),
             ("prepare_flag", "D2H"), ("energy_flag", "D2H"), *extra]
    return dict(route="device-flat-tet", charges=4,
                transfer_ledger=[dict(payload=p, direction=d, bytes=8) for p, d in rows])


def test_assembly_ledger_rejects_host_coefficients_and_gram_downloads():
    assert hdiv_gpu._checked_device_assembly(ledger(), 4)["charges"] == 4
    for forged in (ledger(("host_hmatrix_coefficients", "H2D")), ledger(("cell_vertices_m", "D2H")),
                   dict(ledger(), route="host-leaves")):
        with pytest.raises(RuntimeError):
            hdiv_gpu._checked_device_assembly(forged, 4)
    missing = ledger()
    missing["transfer_ledger"] = [r for r in missing["transfer_ledger"] if r["payload"] != "energy_flag"]
    with pytest.raises(RuntimeError, match="ledger"):
        hdiv_gpu._checked_device_assembly(missing, 4)


def test_device_gram_source_downloads_only_validity_flags():
    source = (HDIV / "hdiv_flat_gram.cuh").read_text(encoding="utf-8")
    downloads = re.findall(r"cudaMemcpy\w*\(([^;]*)DeviceToHost", source)
    assert len(downloads) == 1 and "sizeof(int)" in downloads[0], downloads
    assert "cudaMemcpyToSymbol" not in source and "cudaMemcpyHostToDevice" not in source


def test_flat_binding_accepts_no_matrix_payload():
    source = (HDIV / "hdiv_resident_python.cpp").read_text(encoding="utf-8")
    block = source[source.index('def_static("from_flat_charges"'):source.index('.def("assembly"')]
    names = re.findall(r'py::arg\("([a-z_]+)"\)', block)
    assert names[:10] == ["cell_vertices", "face_vertices", "charge_host", "charge_kind",
                          "charge_exponents", "tet_rule", "tri_rule", "tet_low_rule",
                          "tri_low_rule", "far_factor"], names
    assert not {"leaves", "permutation", "scale", "a1", "a2"} & set(names), names
