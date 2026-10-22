"""Proposed surrogate test for the solved prolate magnetization response.

Guards: changes to radia.vim.Solve affecting the computed M_avg[2]; operator
and source-RHS changes; both production mutations were numerically detected.
Assumes: adoption record ID and convergence level 1 correspondence are pending;
no adopted-record claim is made yet. Frozen 1338-element mesh, mu_r=1000, H_ext=(0,0,1) A/m, CPU,
NGSolve/library threads=1. Baseline source/native identities are provenance,
not restrictions on future implementations. No scale reduction; heavy 32.5 s.
The nonzero solved response depends on demag and excitation, not an input copy.
The frozen fixture 76c79f3f differs slightly in coordinates from the heavy
validation mesh eb732dbe; the heavy analytic D check (0.342%) establishes the
same production route, not an analytic oracle on this fixture. The stored
source/native hashes are provenance, not a current-code identity gate.
This is past-output regression, not an independent analytic validation.
"""
import hashlib
import importlib.util
from pathlib import Path

import ngsolve as ng
from threadpoolctl import threadpool_limits
from tools.regression_lock import assert_regression, load_baseline


def test_prolate_solved_magnetization_lock():
    root = Path(__file__).resolve().parents[2]
    fixtures = Path(__file__).parent / "data" / "prolate_mavg_z"
    baseline = load_baseline(fixtures / "baseline.json")
    mesh = fixtures / "prolate.vol"
    assert hashlib.sha256(mesh.read_bytes()).hexdigest() == baseline["conditions"]["mesh_sha256"], "frozen mesh identity differs"
    spec = importlib.util.spec_from_file_location("prolate_validation_case", root / "validation_test/feec/test_hdiv_vim_ellipsoid.py")
    case = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(case)
    previous_threads = ng.GetNumThreads()
    ng.SetNumThreads(1)
    try:
        with threadpool_limits(limits=1):
            output = case.run_prolate_case(maxh=0.6, mesh_path=mesh)
    finally:
        ng.SetNumThreads(previous_threads)
    assert_regression({"M_avg_z": float(output["M_avg"][2])}, baseline, conditions=baseline["conditions"])
