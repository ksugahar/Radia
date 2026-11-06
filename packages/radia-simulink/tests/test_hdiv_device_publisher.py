"""Guards: the MATLAB device-setup publisher hands the MEX owned flat C-order
float64/int32 columns, keeps charge and evaluation numbering, sends empty ACA
for dense compression and passes the explicit backend through.
Assumes: Radia owns FE setup (``radia.hdiv_gpu.prepare_hdiv_device_setup``).
"""
import importlib.util
from pathlib import Path
import sys
import types

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[3]


def load(relative):
    spec = importlib.util.spec_from_file_location("hdiv_matlab_publisher", ROOT/relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _setup(compression):
    i32 = lambda *v: np.array(v, dtype=np.int32)
    gram = dict(cell_vertices=np.arange(12.).reshape(4, 3), face_vertices=np.zeros((0, 3)),
                charge_host=i32(0, 0), charge_kind=i32(0, 0), charge_exponents=i32(0, 0, 0, 1, 0, 0),
                tet_rule=np.ones((2, 4)), tri_rule=np.ones((1, 4)), tet_low_rule=np.zeros((0, 4)),
                tri_low_rule=np.zeros((0, 4)), cell_class=i32(0), cell_dofs=i32(3, 1, 0, 2),
                cell_tables=np.arange(6.).reshape(2, 3), cell_exponents=i32(0, 0, 0, 1, 0, 0), face_class=i32(),
                face_dofs=i32(), face_tables=np.zeros(0), face_exponents=np.zeros((0, 2), dtype=np.int32),
                far_factor=1e30, compression=compression,
                aca=dict(eps=1e-8, eta=2., leaf_size=32, max_rank=64, rank_cap="fail"))
    tables = dict(rule=np.ones((3, 4)), reference_shapes=np.arange(8.).reshape(2, 4), vertices=np.zeros((4, 3)),
                  element_class=i32(0), element_dofs=i32(3, 1, 0, 2), element_points=i32(0))
    return dict(schema="radia.hdiv-device-setup.v1", backend="gpu", family="affine-tet", order=2, dofs=4,
                mesh_sha256="mesh-hash", label_contract_sha256="contract-hash", evaluation_tables=dict(tables=[tables]),
                charge_gram=gram, material_budget={}, rhs_operator=np.arange(12.).reshape(4, 3),
                mean_operator=np.arange(12.).reshape(3, 4))


@pytest.mark.parametrize("compression", ["aca", "dense"])
def test_device_publisher_owns_flat_columns_and_numbering(tmp_path, monkeypatch, compression):
    module = load("packages/radia-simulink/src/radia_simulink/hdiv.py")
    data = _setup(compression)
    calls = []

    def prepare(path, *, backend, order, **checked):
        if backend != "gpu":
            raise ValueError("explicit backend")
        calls.append(order)
        return data
    monkeypatch.setitem(sys.modules, "radia.hdiv_gpu", types.SimpleNamespace(prepare_hdiv_device_setup=prepare))
    result = module.prepare_hdiv_device_setup(tmp_path/"mesh.vol", backend="gpu", label_contract="labels.json",
                                              mesh_check_receipt="receipt.json", mesh_check_receipt_sha256="a"*64,
                                              mesh_check_provenance={}, order=2)
    assert calls == [2] and result["schema"] == "radia.hdiv-device-setup.v1" and result["dofs"] == 4.
    charges, evaluation = result["charges"], result["evaluation"][0]
    for name, value in list(charges.items()) + list(evaluation.items()):
        if isinstance(value, np.ndarray) and name != "aca":
            assert value.ndim == 2 and value.shape[1] == 1 and value.dtype in (np.float64, np.int32), name
    np.testing.assert_array_equal(charges["cell_dofs"].ravel(), [3, 1, 0, 2])
    np.testing.assert_array_equal(evaluation["element_dofs"].ravel(), [3, 1, 0, 2])
    np.testing.assert_array_equal(charges["cell_tables"].ravel(), np.arange(6.))
    np.testing.assert_array_equal(result["rhs_operator"], data["rhs_operator"])
    assert charges["dofs"] == 4. and charges["compression"] == compression
    if compression == "dense":
        assert charges["aca"].size == 0
    else:
        assert charges["aca"] == dict(eps=1e-8, eta=2., leaf_size=32., max_rank=64., rank_cap="fail")
    charges["cell_tables"][0, 0] = 99
    assert data["charge_gram"]["cell_tables"][0, 0] == 0.  # Publisher owns values.
    data["charge_gram"]["cell_dofs"] = data["charge_gram"]["cell_dofs"].astype(np.int64)
    with pytest.raises(ValueError, match="int32"):
        module.prepare_hdiv_device_setup(tmp_path/"mesh.vol", backend="gpu", label_contract="labels.json",
                                         mesh_check_receipt="receipt.json", mesh_check_receipt_sha256="a"*64,
                                         mesh_check_provenance={})
    with pytest.raises(ValueError, match="explicit backend"):
        module.prepare_hdiv_device_setup(tmp_path/"mesh.vol", backend="cpu", label_contract="labels.json",
                                         mesh_check_receipt="receipt.json", mesh_check_receipt_sha256="a"*64,
                                         mesh_check_provenance={})


@pytest.mark.parametrize("compression", ["aca", "dense"])
def test_setup_file_keeps_tables_and_rejects_unknown_requests(tmp_path, monkeypatch, compression):
    """Guards: the process boundary keeps every table, dtype and scalar and rejects unknown fields."""
    h5py = pytest.importorskip("h5py")
    module = load("packages/radia-simulink/src/radia_simulink/hdiv.py")
    data = _setup(compression)
    aca = dict(eps=1e-8, eta=2, leaf_size=32, max_rank=64, rank_cap="fail") if compression == "aca" else None
    seen = {}

    def prepare(path, **checked):
        seen.update(checked, path=path)
        return data
    monkeypatch.setitem(sys.modules, "radia.hdiv_gpu", types.SimpleNamespace(prepare_hdiv_device_setup=prepare))
    request = dict(mesh_path="mesh.vol", backend="gpu", label_contract="labels.json", mesh_check_receipt="receipt.json",
                   mesh_check_receipt_sha256="a"*64, mesh_check_provenance={}, order=2, curve_order=1,
                   gram_compression=compression)
    if aca:
        request["gram_aca"] = aca
    module.write_hdiv_device_setup(request, tmp_path/"setup.h5")
    assert seen["path"] == "mesh.vol" and seen["order"] == 2
    if aca:
        assert isinstance(seen["gram_aca"]["eta"], float) and isinstance(seen["gram_aca"]["leaf_size"], int)
    expected = module.prepare_hdiv_device_setup("mesh.vol", **{k: v for k, v in request.items() if k != "mesh_path"})
    with h5py.File(tmp_path/"setup.h5", "r") as f:
        assert f.attrs["schema"] == "radia.hdiv-device-setup.v1" and f.attrs["index_base"] == 0 and f.attrs["units"] == "SI"
        np.testing.assert_array_equal(f["rhs_operator"], expected["rhs_operator"])
        for group, tables in ((f["charges"], expected["charges"]), (f["evaluation/group1"], expected["evaluation"][0])):
            for key, value in tables.items():
                if isinstance(value, np.ndarray) and key != "aca":
                    assert group[key].dtype == value.dtype, key
                    np.testing.assert_array_equal(group[key], value.ravel())
        assert f["charges"].attrs["compression"] == compression and ("aca" in f["charges"]) == bool(aca)
    with pytest.raises(ValueError, match="request fields"):
        module.write_hdiv_device_setup(dict(request, stray=1), tmp_path/"other.h5")
