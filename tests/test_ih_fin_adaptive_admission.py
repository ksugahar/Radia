"""Budget and provenance admission of the adaptive fin reference.

Guards: ``adaptive_sequence`` checks the actual deadline before each solve and before accepting
a checkpoint, and the actual DOF count before a solve; ``load_sequences`` admits a saved sequence
only with this source, script and runtime, the configured settings, one row per configured
(order, exterior) and a convergence that re-derives from its levels.
Assumes: a deterministic stand-in solver whose quantities converge as the DOFs double; no FEM runs.
"""

from __future__ import annotations

import copy
import importlib.util
import json
import sys
import types
from argparse import Namespace
from pathlib import Path

import pytest

ng = pytest.importorskip("ngsolve")

SCRIPT = Path(__file__).resolve().parents[1]/"validation_test"/"induction_heating"/"axisym_volume_peec_fin.py"
TOL = {k: 0.005 for k in ("R", "partial_L", "root_heat", "tip_heat")}


@pytest.fixture()
def fin(monkeypatch):
    spec = importlib.util.spec_from_file_location("axisym_volume_peec_fin_under_test", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class FakeMesh:
    def __init__(self, ndof):
        self.ndof = ndof

    @property
    def ne(self):
        return self.ndof

    def Elements(self, vt):
        return [types.SimpleNamespace(nr=i) for i in range(4)]

    def SetRefinementFlag(self, el, flag):
        pass

    def Refine(self):
        self.ndof *= 2


def _value(ndof):
    return 1.0 + 100.0/ndof          # converges as the DOFs double


class Clock:
    """Wall clock that advances by ``step`` per call, plus a jump at the n-th solve."""

    def __init__(self, step=1.0, jump_at_solve=None, jump=0.0):
        self.t, self.step, self.solves, self.jump_at_solve, self.jump = 0.0, step, 0, jump_at_solve, jump

    def __call__(self):
        self.t += self.step
        return self.t

    def solved(self):
        self.solves += 1
        if self.solves == self.jump_at_solve:
            self.t += self.jump


@pytest.fixture()
def stand_in(fin, monkeypatch):
    calls = {"solves": 0}
    clock = Clock()

    def build_section_ring(radius_mm, **kw):
        return FakeMesh(1000), 1.0, 0.0

    monkeypatch.setitem(sys.modules, "aphi_beak_section_sweep",
                        types.SimpleNamespace(build_section_ring=build_section_ring))
    monkeypatch.setattr(fin, "level_ndof", lambda mesh, order: mesh.ndof)

    def solve_level(mesh, *a, **kw):
        calls["solves"] += 1
        clock.solved()
        v = _value(mesh.ndof)
        row = dict(order=1, ndof=mesh.ndof, elements=mesh.ndof, seconds=1.0, R_ohm=v, L_H=v,
                   root_heat_W_per_A2=v, tip_heat_W_per_A2=v, terminal_vs_energy_R=0.0)
        return row, None

    monkeypatch.setattr(fin, "solve_level", solve_level)
    monkeypatch.setattr(fin, "zz_indicator", lambda res, order: __import__("numpy").ones(4))
    monkeypatch.setattr(fin, "mark_whole_mesh",
                        lambda mesh, eta2, theta: (__import__("numpy").ones(4, dtype=bool), 2, 2))
    return calls, clock


def _run(fin, clock, max_ndof=1e9, budget_s=1e6, max_levels=40):
    return fin.adaptive_sequence(200.0, 150e3, 5.8e7, 1, 1.5, 0.05, 20.0, 20.0, 0.5, TOL,
                                 max_ndof, max_levels, budget_s, clock=clock)


def test_valid_sequence_converges_within_budget(fin, stand_in):
    calls, clock = stand_in
    seq = _run(fin, clock)
    assert seq["stop"] == "converged", seq["stop"]
    assert seq["converged"] is True
    assert fin.certify_sequence(seq, _expected(seq)) == []


def test_initially_expired_deadline_solves_nothing(fin, stand_in):
    calls, clock = stand_in
    seq = _run(fin, clock, budget_s=0.5)
    assert (calls["solves"], seq["stop"], seq["converged"]) == (0, "budget: deadline", False)


def test_deadline_passed_while_counting_dofs_solves_nothing(fin, stand_in, monkeypatch):
    calls, clock = stand_in

    def slow_count(mesh, order):
        clock.t += 20.0                              # counting takes longer than the budget
        return mesh.ndof

    monkeypatch.setattr(fin, "level_ndof", slow_count)
    seq = _run(fin, clock, budget_s=10.0)
    assert (calls["solves"], seq["stop"], seq["converged"]) == (0, "budget: deadline", False)


@pytest.mark.parametrize("cap", [0, -1, True, 1.5])
def test_nonpositive_or_non_integer_level_cap_is_refused_before_any_solve(fin, stand_in, cap):
    calls, clock = stand_in
    with pytest.raises(ValueError, match="max_levels"):
        _run(fin, clock, max_levels=cap)
    assert calls["solves"] == 0


def test_dof_overflow_is_refused_before_the_solve(fin, stand_in):
    calls, clock = stand_in
    seq = _run(fin, clock, max_ndof=10)
    assert calls["solves"] == 0, calls
    assert seq["stop"].startswith("budget: next level 1000 DOFs"), seq["stop"]
    assert seq["converged"] is False


def test_dof_overflow_after_refinement_stops_unresolved(fin, stand_in):
    calls, clock = stand_in
    seq = _run(fin, clock, max_ndof=3000)
    assert max(r["ndof"] for r in seq["levels"]) <= 3000
    assert (seq["stop"].startswith("budget: next level 4000 DOFs"), seq["converged"]) == (True, False)


def test_checkpoint_finished_past_the_deadline_cannot_converge(fin, stand_in):
    calls, clock = stand_in
    n = len(_run(fin, Clock())["levels"])              # the converging solve is the n-th
    clock.jump_at_solve, clock.jump = clock.solves + n, 1e7
    seq = _run(fin, clock)
    last = seq["levels"][-1]
    assert (seq["stop"], seq["converged"], last["within_budget"], last["checkpoint"]) \
        == ("budget: deadline", False, False, False), (seq["stop"], last)


# ---------------------------------------------------------------- saved-sequence admission

def _expected(seq):
    return dict(start_maxh_mm=0.05, air_maxh_mm=20.0, air_divisor=20.0, theta=0.5, checkpoint_dof_growth=2.0,
                tolerance=TOL, max_ndof=1e9, max_levels=40, budget_s=1e6)


def _args():
    return Namespace(ref_order=[1, 2], ref_exterior_factor=[1.5, 3.0], ref_theta=0.5, ref_max_ndof=1e9,
                     ref_max_levels=40, ref_hours=1e6/3600, radius_mm=200.0, frequency_hz=150e3, sigma=5.8e7)


PROVENANCE = dict(natives="n", version="v", ring_solver="r", ring_solver_sha256="a", axifem="x",
                  axifem_sha256="b")


@pytest.fixture()
def saved(fin, stand_in, monkeypatch, tmp_path):
    monkeypatch.setattr(fin, "git_state", lambda: ("head", False))
    base = _run(fin, Clock())
    docs = {}
    for o in (1, 2):
        for f in (1.5, 3.0):
            seq = copy.deepcopy(base)
            seq.update(order=o, exterior_factor=f)
            docs[(o, f)] = dict(schema=fin.SEQUENCE_SCHEMA, git_head="head", git_dirty=False,
                                script_sha256=fin.script_sha256(),
                                reference_runtime=dict(PROVENANCE), radius_mm=200.0, frequency_hz=150e3,
                                sigma_S_per_m=5.8e7, sequence=seq)

    def write(mutate=None, extra=None):
        out = []
        for i, (key, doc) in enumerate(sorted(docs.items())):
            d = copy.deepcopy(doc)
            if mutate and i == 0:
                mutate(d)
            p = tmp_path/f"seq{i}.json"
            p.write_text(json.dumps(d), encoding="utf-8")
            out.append(p)
        if extra:
            p = tmp_path/"extra.json"
            d = copy.deepcopy(docs[extra])
            p.write_text(json.dumps(d), encoding="utf-8")
            out.append(p)
        return out

    settings = {1: (0.05, 20.0, 20.0), 2: (0.05, 20.0, 20.0)}
    load = lambda paths: fin.load_sequences(paths, _args(), PROVENANCE, settings, TOL)
    return write, load


def test_valid_saved_sequences_are_admitted(saved):
    write, load = saved
    assert len(load(write())) == 4


@pytest.mark.parametrize("name, mutate, reason", [
    ("script", lambda d: d.update(script_sha256="0"*64), "script_sha256"),
    ("native", lambda d: d["reference_runtime"].update(axifem_sha256="c"), "reference_runtime"),
    ("source", lambda d: d.update(git_head="other"), "git_head"),
    ("dirty", lambda d: d.update(git_dirty=True), "git_dirty"),
    ("schema", lambda d: d.update(schema="v0"), "schema"),
    ("theta", lambda d: d["sequence"].update(theta=0.01), "theta"),
    ("max_ndof", lambda d: d["sequence"].update(max_ndof=1), "max_ndof"),
    ("budget", lambda d: d["sequence"].update(budget_s=1.0), "budget_s"),
    ("level_dofs", lambda d: d["sequence"]["levels"][-1].update(ndof=int(2e9)), "DOFs above"),
    ("late_level", lambda d: d["sequence"]["levels"][-1].update(within_budget=False), "wall-clock"),
    ("not_converged", lambda d: d["sequence"].update(stop="budget: deadline", converged=False), "stop"),
    ("claimed_change", lambda d: d["sequence"]["levels"][-1].update(R_ohm=2.0), "last checkpoints differ"),
    ("sequence_over_budget", lambda d: d["sequence"].update(elapsed_s=1e6 + 1), "above the budget"),
    ("sequence_before_last", lambda d: d["sequence"].update(elapsed_s=0.0), "before its last level"),
    ("sequence_elapsed_nan", lambda d: d["sequence"].update(elapsed_s=float("nan")), "sequence elapsed_s"),
    ("level_elapsed_nan", lambda d: d["sequence"]["levels"][-1].update(elapsed_s=float("nan")), "elapsed_s"),
    ("level_elapsed_negative", lambda d: d["sequence"]["levels"][0].update(elapsed_s=-1.0), "elapsed_s"),
    ("level_elapsed_decreasing", lambda d: d["sequence"]["levels"][1].update(elapsed_s=0.0), "decrease"),
    ("ndof_nan", lambda d: d["sequence"]["levels"][-1].update(ndof=float("nan")), "not a positive integer"),
    ("ndof_float", lambda d: d["sequence"]["levels"][-1].update(ndof=1000.5), "not a positive integer"),
    ("ndof_bool", lambda d: d["sequence"]["levels"][0].update(ndof=True), "not a positive integer"),
    ("ndof_zero", lambda d: d["sequence"]["levels"][0].update(ndof=0), "not a positive integer"),
    ("level_number_bool", lambda d: d["sequence"]["levels"][1].update(level=True), "not an integer"),
    ("quantity_nan", lambda d: d["sequence"]["levels"][-1].update(tip_heat_W_per_A2=float("nan")), "quantities"),
    ("no_levels", lambda d: d["sequence"].update(levels=[]), "no recorded levels"),
])
def test_invalid_saved_sequence_is_refused(saved, name, mutate, reason):
    write, load = saved
    with pytest.raises(ValueError, match="saved sequences refused") as err:
        load(write(mutate))
    assert reason in str(err.value), str(err.value)


def test_duplicate_saved_sequence_is_refused(saved):
    write, load = saved
    with pytest.raises(ValueError, match="given twice"):
        load(write(extra=(1, 1.5)))


def test_missing_saved_sequence_is_refused(saved):
    write, load = saved
    with pytest.raises(ValueError, match="do not cover"):
        load(write()[:3])


def test_checkpoint_on_the_deadline_with_a_late_final_sample_does_not_converge(fin, stand_in):
    calls, clock = stand_in
    probe = _run(fin, Clock())
    last = probe["levels"][-1]["elapsed_s"]          # the converging checkpoint, sampled exactly at the deadline
    seq = _run(fin, clock, budget_s=last)
    assert seq["levels"][-1]["within_budget"] is True and seq["levels"][-1]["checkpoint"] is True
    assert seq["elapsed_s"] > seq["budget_s"], (seq["elapsed_s"], seq["budget_s"])
    assert (seq["converged"], seq["stop"]) == (False, "budget: deadline after the last checkpoint")
    assert fin.certify_sequence(seq, {**_expected(seq), "budget_s": last}) != []


def _main_with_live_sequences(fin, monkeypatch, tmp_path, edit):
    base = _run(fin, Clock())
    base["tolerance"] = {k: 0.1*0.05 for k in TOL}         # main(): one tenth of the 5 % limits
    v = _value(base["levels"][-1]["ndof"])
    edit(base)
    monkeypatch.setattr(fin, "_branch_ring_solver", lambda: (None, dict(PROVENANCE)))
    monkeypatch.setattr(fin, "adaptive_sequence",
                        lambda radius, f, s, order, factor, *a, **kw: {**copy.deepcopy(base), "order": order,
                                                                        "exterior_factor": factor})
    monkeypatch.setattr(fin, "git_state", lambda: ("head", False))
    monkeypatch.setattr(fin, "_peec_runtime", lambda: {})
    monkeypatch.setattr(fin, "section_mesh", lambda *a: (None, None, None))
    monkeypatch.setattr(fin, "solve_peec", lambda *a: dict(R_ohm=v, L_H=v, root_heat_W_per_A2=v,
                                                           tip_heat_W_per_A2=v))
    out = tmp_path/"report.json"
    monkeypatch.setattr(sys, "argv", ["fin", "--ref-start-maxh-mm", "0.05", "0.05", "--ref-air-mm", "20", "20",
                                      "--ref-air-divisor", "20", "20", "--ref-max-ndof", "1e9",
                                      "--ref-max-levels", "40", "--ref-hours", str(1e6/3600),
                                      "--peec-boundary-h-mm", "0.06", "--output", str(out)])
    fin.main()
    return json.loads(out.read_text(encoding="utf-8"))


def test_live_sequences_are_certified_like_saved_ones(fin, stand_in, monkeypatch, tmp_path):
    report = _main_with_live_sequences(fin, monkeypatch, tmp_path, lambda seq: None)
    assert report["reference_converged"] is True, report["reference_certification"]
    assert not any(report["reference_certification"].values())


def test_live_sequence_over_budget_leaves_the_reference_unresolved(fin, stand_in, monkeypatch, tmp_path):
    def over_budget(seq):
        seq["elapsed_s"] = seq["budget_s"] + 1.0      # still claims converged=True
    report = _main_with_live_sequences(fin, monkeypatch, tmp_path, over_budget)
    assert (report["reference_converged"], report["reference_resolved"]) == (False, False)
    assert all("above the budget" in " ".join(r) for r in report["reference_certification"].values()), \
        report["reference_certification"]


def test_live_sequence_without_a_checkpoint_records_an_unresolved_verdict(fin, stand_in, monkeypatch, tmp_path):
    def no_levels(seq):
        seq.update(levels=[], converged=False, stop="budget: next level 1000 DOFs above 1")
    report = _main_with_live_sequences(fin, monkeypatch, tmp_path, no_levels)
    assert (report["reference_converged"], report["reference_resolved"]) == (False, False)
    assert report["verdict"]["status"].startswith("NOT_ADOPTED"), report["verdict"]
    assert report["reference_sequences_without_checkpoint"] == ["o1_f1.5", "o1_f3.0", "o2_f1.5", "o2_f3.0"]
