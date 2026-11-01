"""Regression lock for the solved prolate magnetization response on the named fixture.

Guards: changes to radia.vim.Solve affecting the computed M_avg[2]; operator,
source-RHS and graded-P1 charge-Gram selector changes (the selector gate-off
mutation was detected at this bound); loss of the magnetic_body / exterior labels.
Assumes: frozen named 1338-element mesh 27c37786 checked against its versioned
label contract, mu_r=1000, H_ext=(0,0,1) A/m, CPU, NGSolve/library threads=1.
Heavy numerical tier: one recorded run took about 760 s; not a fast surrogate.
The stored source/native hashes are provenance, not a current-code identity gate.
This is past-output regression, not an independent analytic validation; the
retired unnamed 76c lock is history only (see baseline provenance).
"""
import hashlib
import importlib.util
import json
from pathlib import Path

import ngsolve as ng
from threadpoolctl import threadpool_limits


def test_prolate_solved_magnetization_lock():
    root = Path(__file__).resolve().parents[2]
    helper_spec = importlib.util.spec_from_file_location("prolate_regression_helper", root / "tools/regression_lock.py")
    helper = importlib.util.module_from_spec(helper_spec)
    helper_spec.loader.exec_module(helper)
    fixtures = Path(__file__).parent / "data" / "prolate_named"
    baseline = helper.load_baseline(fixtures / "baseline.json")
    conditions = baseline["conditions"]
    mesh = fixtures / "prolate_named.vol"
    contract_path = fixtures / "prolate-label-contract.json"
    actual_mesh = hashlib.sha256(mesh.read_bytes()).hexdigest()
    assert actual_mesh == conditions["mesh_sha256"], f"mesh sha256: actual={actual_mesh}; expected={conditions['mesh_sha256']}"
    actual_contract = hashlib.sha256(contract_path.read_bytes()).hexdigest()
    assert actual_contract == conditions["label_contract_sha256"], f"contract sha256: actual={actual_contract}; expected={conditions['label_contract_sha256']}"
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    assert contract["required"] == conditions["labels"], f"contract labels: actual={contract['required']}; expected={conditions['labels']}"
    spec = importlib.util.spec_from_file_location("prolate_validation_case", root / "validation_test/feec/test_hdiv_vim_ellipsoid.py")
    case = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(case)
    expected_labels = (tuple(conditions["labels"]["materials"]), tuple(conditions["labels"]["boundaries"]))
    assert case.PROLATE_LABELS == expected_labels, f"helper labels: actual={case.PROLATE_LABELS}; expected={expected_labels}"
    previous_threads = ng.GetNumThreads()
    ng.SetNumThreads(1)
    try:
        with threadpool_limits(limits=1):
            output = case.run_prolate_case(mesh_path=mesh)
    finally:
        ng.SetNumThreads(previous_threads)
    helper.assert_regression({"M_avg_z": float(output["M_avg"][2])}, baseline, conditions=conditions)
