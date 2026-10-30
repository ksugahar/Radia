"""Guards: explicit backend, owned snapshots and workflow lifetime.

Assumes: mocked native owners do not establish numerical GPU validation.
"""
import importlib.util
import json
from pathlib import Path
import sys
import threading
import types

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def load(relative):
    spec = importlib.util.spec_from_file_location("workflow_under_test", ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


API = "src/radia/hdiv_gpu.py"


def test_explicit_backend_before_native_import():
    module = load(API)
    for backend in (None, "cpu", "auto"):
        with pytest.raises(ValueError, match="explicit backend"):
            module.PreparedHDivNonlinear("missing.vol", {}, backend=backend, label_contract="missing.json", mesh_check_receipt="missing.receipt", mesh_check_receipt_sha256="a"*64, mesh_check_provenance={})


def test_solve_reuses_native_owner_and_close_releases_after_error():
    module = load(API)
    owner = module.PreparedHDivNonlinear.__new__(module.PreparedHDivNonlinear)
    owner._lock = threading.RLock()
    owner._closed = False
    owner.dofs = 3
    owner.order = 1
    owner.mesh_sha256 = "fixture"
    owner._rhs = np.eye(3)
    owner._mass = types.SimpleNamespace(indptr=np.arange(4), indices=np.arange(3), data=np.ones(3))
    calls = []

    class Native:
        def solve_nonlinear_material(self, material, *args):
            calls.append((material, args))
            return dict(x=args[3].copy(), relative_residual=0.)

        def close(self):
            calls.append("operator_closed")

    class Material:
        def close(self):
            calls.append("material_closed")
            raise RuntimeError("material close failure")

    owner._operator, owner._material = Native(), Material()
    for field in ([1., 2., 3.], [3., 2., 1.]):
        result = owner.solve(field)
        np.testing.assert_array_equal(result["x"], field)
        assert result["backend"] == "gpu"
    assert calls[0][0] is calls[1][0]
    with pytest.raises(ValueError, match="three finite"):
        owner.solve([float("nan"), 0, 0])
    with pytest.raises(RuntimeError, match="material close failure"):
        owner.close()
    assert calls[-2:] == ["material_closed", "operator_closed"]
    with pytest.raises(RuntimeError, match="closed"):
        owner.solve([0, 0, 0])


@pytest.mark.parametrize("options", [dict(preconditioner="mass-cholesky"), dict(tolerance=True),
    dict(tolerance=1e-5), dict(degree=65), dict(degree=1.5), dict(lower_fraction=1.),
    dict(max_iterations=False), dict(linear_tolerance=float("nan"))])
def test_invalid_solve_options_before_native_work(options):
    module = load(API)
    owner = module.PreparedHDivNonlinear.__new__(module.PreparedHDivNonlinear)
    with pytest.raises(ValueError):
        owner.solve([0., 0., 0.], **options)


def test_snapshot_uses_exact_checked_bytes_after_caller_mutation(tmp_path, monkeypatch):
    import hashlib
    module = load(API)
    mesh = tmp_path / "mesh.vol"; mesh.write_bytes(b"original mesh")
    contract = tmp_path / "labels.json"; contract.write_text("{}")
    receipt = tmp_path / "receipt.json"; receipt.write_text("{}")
    def verify(mesh_bytes, contract_bytes, receipt_bytes, *pins):
        assert mesh_bytes == b"original mesh"
        mesh.write_bytes(b"changed original")
        contract.write_text("changed contract")
    monkeypatch.setattr(module, "_verify_mesh_receipt", verify)
    snapshots = []
    def initialize(self, path, law, order):
        snapshots.append(path)
        assert path.read_bytes() == b"original mesh"
        assert self.mesh_sha256 == hashlib.sha256(b"original mesh").hexdigest()
    monkeypatch.setattr(module.PreparedHDivNonlinear, "_initialize", initialize)
    module.PreparedHDivNonlinear(mesh, {}, backend="gpu", label_contract=contract,
        mesh_check_receipt=receipt, mesh_check_receipt_sha256="a"*64, mesh_check_provenance={})
    assert not snapshots[0].exists()
