"""Workflow dispatch/lifetime contracts; native physics has separate validation."""
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


API = "packages/radia-simulink/src/radia_simulink/hdiv.py"
MCP = "packages/radia-mcp/src/radia_mcp/radia_ngsolve/hdiv_gpu.py"


def test_explicit_backend_before_native_import():
    module = load(API)
    for backend in (None, "cpu", "auto"):
        with pytest.raises(ValueError, match="explicit backend"):
            module.PreparedHDivNonlinear("missing.vol", {}, backend=backend, label_contract="missing.json")


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


@pytest.mark.parametrize("fail", [False, True])
@pytest.mark.parametrize("table_form", [False, True])
@pytest.mark.parametrize("serialization_failure", [False, True])
def test_mcp_reuse_failure_evidence_and_release(tmp_path, monkeypatch, fail, table_form, serialization_failure):
    module = load(MCP)
    events = []

    class Region:
        def __enter__(self): events.append("region_enter")
        def __exit__(self, *_): events.append("region_exit")

    class Prepared:
        def __init__(self, *args, **kwargs): events.append("prepare")
        @classmethod
        def from_bh_table(cls, *args, **kwargs):
            events.append("table")
            return cls(*args, **kwargs)
        def __enter__(self): return self
        def __exit__(self, *_): events.append("close")
        def stats(self): return dict(closed=True)
        def solve(self, field, **kwargs):
            events.append("solve")
            if fail and events.count("solve") == 2:
                raise RuntimeError("native convergence failure")
            if serialization_failure and events.count("solve") == 2:
                return dict(x=np.array(field), unsupported=object())
            return dict(x=np.array(field), relative_residual=np.float64(0.), nested={"iterations": np.int64(3)})

    monkeypatch.setitem(sys.modules, "ngsolve", types.SimpleNamespace(TaskManager=Region))
    monkeypatch.setitem(sys.modules, "radia_simulink.hdiv", types.SimpleNamespace(PreparedHDivNonlinear=Prepared))
    law = tmp_path / "law.json"
    law.write_text(json.dumps(dict(schema="radia.bh-table.v1", table=[[0., 0.], [1., .1]])) if table_form else "{}")
    result = tmp_path / "result.json"
    args = ("mesh.vol", "labels.json", str(law), [[1., 2., 3.], [3., 2., 1.]], str(result), "gpu")
    if fail or serialization_failure:
        with pytest.raises((RuntimeError, AttributeError)):
            module.hdiv_nonlinear_gpu_run(*args)
    else:
        assert module.hdiv_nonlinear_gpu_run(*args)["closed"]
    saved = json.loads(result.read_text())
    assert saved["accepted"] is (not (fail or serialization_failure))
    assert saved["independent_cpu_certified"] is False
    assert len(saved["cases"]) == (1 if fail or serialization_failure else 2)
    assert events == ["region_enter"] + (["table"] if table_form else []) + ["prepare", "solve", "solve", "close", "region_exit"]
    assert saved["preconditioner"] == "reference-diagonal"
    assert len(saved["material_sha256"]) == 64
    with pytest.raises(FileExistsError):
        module.hdiv_nonlinear_gpu_run(*args)


@pytest.mark.parametrize("options", [dict(preconditioner="mass-cholesky"), dict(tolerance=True),
    dict(tolerance=1e-5), dict(degree=65), dict(degree=1.5), dict(lower_fraction=1.),
    dict(max_iterations=False), dict(linear_tolerance=float("nan"))])
def test_invalid_solve_options_before_native_work(options):
    module = load(API)
    owner = module.PreparedHDivNonlinear.__new__(module.PreparedHDivNonlinear)
    with pytest.raises(ValueError):
        owner.solve([0., 0., 0.], **options)


def test_checked_mesh_snapshot_identity_and_stderr(tmp_path, monkeypatch):
    import hashlib
    module = load(API)
    mesh = tmp_path / "mesh.vol"
    mesh.write_bytes(b"original mesh")
    contract = tmp_path / "labels.json"
    contract.write_text("{}")
    snapshots = []
    def check(command, **kwargs):
        snapshot = Path(command[3])
        snapshots.append(snapshot)
        assert snapshot.read_bytes() == b"original mesh"
        mesh.write_bytes(b"changed original")
        return types.SimpleNamespace(returncode=0, stdout='{"passed":true}', stderr="")
    def initialize(self, path, law, order):
        assert path == snapshots[-1]
        assert path.read_bytes() == b"original mesh"
        assert self.mesh_sha256 == hashlib.sha256(b"original mesh").hexdigest()
    monkeypatch.setattr(module.subprocess, "run", check)
    monkeypatch.setattr(module.PreparedHDivNonlinear, "_initialize", initialize)
    module.PreparedHDivNonlinear(mesh, {}, backend="gpu", label_contract=contract)
    assert not snapshots[0].exists()
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k:
        types.SimpleNamespace(returncode=1, stdout="", stderr="No module named cubit_mesh_export"))
    with pytest.raises(RuntimeError, match="cubit-mesh-export required.*No module"):
        module.PreparedHDivNonlinear(mesh, {}, backend="gpu", label_contract=contract)


def test_host_publisher_preserves_factors_numbering_and_checked_identity(tmp_path, monkeypatch):
    from scipy.sparse import eye, csr_matrix
    module = load(API)
    mesh, labels = tmp_path / "mesh.vol", tmp_path / "labels.json"
    mesh.write_bytes(b"mesh bytes"); labels.write_text("{}")
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k:
        types.SimpleNamespace(returncode=0, stdout='{"passed":true}', stderr=""))
    leaf = dict(row=0, column=0, rows=2, columns=2, rank=0,
                a1=np.array([3., 1., 1., 4.]), a2=np.empty(0))
    snapshot = dict(dimension=2, leaves=[leaf], permutation=np.array([1, 0]), scale=np.array([2., 3.]))
    data = dict(dofs=2, snapshot=snapshot, charge_map=eye(2, format="csr"),
        mass=eye(2, format="csr"), rhs_operator=np.ones((2, 3)), mean_operator=np.ones((3, 2)),
        evaluation=csr_matrix(np.ones((3, 2))), weights=np.ones(1), quadrature_order=6)
    monkeypatch.setattr(module, "_prepare_host_components", lambda path, order, **kwargs: data)
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


@pytest.mark.parametrize("compressed", [False, True])
def test_real_label_checker_is_filename_neutral_and_preserves_compression(tmp_path,monkeypatch,compressed):
    import gzip,hashlib,os
    module=load(API)
    contents=(ROOT/"tests/fixtures/hdiv/hex_pair.vol").read_bytes()
    if compressed:contents=gzip.compress(contents,mtime=0)
    path=tmp_path/("caller-name.vol.gz" if compressed else "caller-name.vol");path.write_bytes(contents)
    contract=tmp_path/"caller-labels.json"
    contract.write_text(json.dumps(dict(schema="radia.vol-label-contract.v1",
        required=dict(materials=["body"],boundaries=["exterior"]),
        allowed=dict(materials=["body"],boundaries=["exterior"]))))
    monkeypatch.setenv("PYTHONPATH",str(ROOT/"packages/cubit-mesh-export/src")+os.pathsep+os.environ.get("PYTHONPATH",""))
    def initialize(self,snapshot,law,order):
        assert snapshot.name == ("mesh.vol.gz" if compressed else "mesh.vol")
        assert snapshot.read_bytes() == contents
        assert self.mesh_sha256 == hashlib.sha256(contents).hexdigest()
    monkeypatch.setattr(module.PreparedHDivNonlinear,"_initialize",initialize)
    module.PreparedHDivNonlinear(path,{},backend="gpu",label_contract=contract)
