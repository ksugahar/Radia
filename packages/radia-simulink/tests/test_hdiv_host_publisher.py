"""Guards: MATLAB publisher keeps numbering/owned arrays; Assumes: Radia owns FE setup."""
import importlib.util
from pathlib import Path
import sys,types
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[3]
def load(relative):
    spec=importlib.util.spec_from_file_location("hdiv_matlab_publisher",ROOT/relative)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def test_host_publisher_preserves_factors_numbering_and_checked_identity(tmp_path, monkeypatch):
    from scipy.sparse import eye, csr_matrix
    module = load("packages/radia-simulink/src/radia_simulink/hdiv.py")
    mesh, labels = tmp_path / "mesh.vol", tmp_path / "labels.json"
    mesh.write_bytes(b"mesh bytes"); labels.write_text("{}")
    leaf = dict(row=0, column=0, rows=2, columns=2, rank=0,
                a1=np.array([3., 1., 1., 4.]), a2=np.empty(0))
    snapshot = dict(dimension=2, leaves=[leaf], permutation=np.array([1, 0]), scale=np.array([2., 3.]))
    data = dict(dofs=2, snapshot=snapshot, charge_map=eye(2, format="csr"),
        mass=eye(2, format="csr"), rhs_operator=np.ones((2, 3)), mean_operator=np.ones((3, 2)),
        evaluation=csr_matrix(np.ones((3, 2))), weights=np.ones(1), quadrature_order=6)
    def prepare(path, *, backend, label_contract, order, maximum_dofs):
        if backend != "gpu":raise ValueError("explicit backend")
        assert maximum_dofs == 4096
        return dict(data, mesh_sha256="mesh-hash", label_contract_sha256="contract-hash")
    monkeypatch.setitem(sys.modules, "radia.hdiv_gpu", types.SimpleNamespace(prepare_hdiv_host_setup=prepare))
    result = module.prepare_hdiv_host_setup(mesh, backend="gpu", label_contract=labels, order=2)
    assert result["schema"] == "radia.hdiv-host-setup.v1"
    assert result["order"] == 2.
    np.testing.assert_array_equal(result["operator"]["permutation"].ravel(), [2., 1.])
    np.testing.assert_array_equal(result["operator"]["scale"].ravel(), [2., 3.])
    np.testing.assert_array_equal(result["operator"]["leaves"][0]["a1"].ravel(), leaf["a1"])
    assert result["operator"]["leaves"][0]["row"] == 1.
    result["operator"]["leaves"][0]["a1"][0] = 99
    assert leaf["a1"][0] == 3.  # Publisher owns values, never mutates CPU factors.
    with pytest.raises(ValueError, match="explicit backend"):
        module.prepare_hdiv_host_setup(mesh, backend="cpu", label_contract=labels)
