"""Guards: thin HDiv MCP dispatch/lifetime/evidence; Assumes: Radia owns solver math."""
from pathlib import Path
import importlib.util
import json,sys,types
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]
MCP="src/radia_mcp/radia_ngsolve/hdiv_gpu.py"
def load(relative):
    spec=importlib.util.spec_from_file_location("mcp_hdiv_workflow_test",ROOT/relative)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

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
    monkeypatch.setitem(sys.modules, "radia.hdiv_gpu", types.SimpleNamespace(PreparedHDivNonlinear=Prepared))
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
