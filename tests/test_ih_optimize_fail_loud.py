"""Failure boundaries for the sequential IH design optimizer."""

import math

import pytest

from radia.coil_spec import CoilSpec
from radia.ih_optimize import IHOptimizer


SPEC = CoilSpec(r_wire=1.0e-3)


class Context:
    def __init__(self, outcome):
        self.outcome = outcome

    def evaluate(self, _coil, **_kwargs):
        if isinstance(self.outcome, BaseException):
            raise self.outcome
        return self.outcome


def optimizer(monkeypatch, outcome, objective_fn=None):
    monkeypatch.setattr("radia.ih_optimize.build_coil_from_spec", lambda _spec: object())
    return IHOptimizer(
        Context(outcome), lambda _history: SPEC, objective_fn=objective_fn
    )


@pytest.mark.parametrize(
    "error",
    [ValueError("invalid candidate"), RuntimeError("linear solve failed")],
)
def test_recoverable_candidate_failure_is_recorded(monkeypatch, error):
    opt = optimizer(monkeypatch, error)
    assert opt.run(1, verbose=False) is None
    assert type(error).__name__ in opt.history[0].failure


def test_evaluator_programming_error_propagates(monkeypatch):
    opt = optimizer(monkeypatch, KeyError("missing implementation field"))
    with pytest.raises(KeyError, match="missing implementation field"):
        opt.run(1, verbose=False)
    assert opt.history == []


def test_objective_callback_error_propagates(monkeypatch):
    def broken_objective(_metrics):
        raise RuntimeError("objective bug")

    opt = optimizer(monkeypatch, {"P_total": 1.0}, broken_objective)
    with pytest.raises(RuntimeError, match="objective bug"):
        opt.run(1, verbose=False)
    assert opt.history == []


@pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf])
def test_nonfinite_objective_is_a_failed_trial(monkeypatch, value):
    metrics = {"P_total": 2.0}
    opt = optimizer(monkeypatch, metrics, lambda _metrics: value)
    assert opt.run(1, verbose=False) is None
    trial = opt.history[0]
    assert trial.metrics is metrics
    assert "non-finite" in trial.failure
    assert trial.objective == -math.inf


def test_summary_handles_all_trials_failed(monkeypatch):
    opt = optimizer(monkeypatch, ValueError("bad geometry"))
    opt.run(2, verbose=False)
    assert opt.summary() == (
        "trials: 2 total, 0 ok, 2 failed\n"
        "best  : none (all trials failed)"
    )
