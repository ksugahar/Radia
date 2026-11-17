"""Guards: the dual_fem_reference_estimate MCP tool is a thin adapter that
passes its inputs to radia.dual_fem_delta.dual_bracket_reference unchanged and
propagates its ValueError instead of returning a verdict.
Assumes: the radia numerical entry owns all input validation."""
import sys
import types

import pytest

from radia_mcp.radia_ngsolve import server


def _fake_radia(monkeypatch, entry):
    package = types.ModuleType("radia")
    package.__path__ = []
    module = types.ModuleType("radia.dual_fem_delta")
    module.dual_bracket_reference = entry
    monkeypatch.setitem(sys.modules, "radia", package)
    monkeypatch.setitem(sys.modules, "radia.dual_fem_delta", module)


def test_adapter_forwards_inputs_unchanged(monkeypatch):
    seen = {}

    def entry(levels, coarse=None, convergence_fraction=0.1):
        seen.update(levels=levels, coarse=coarse, fraction=convergence_fraction)
        return {"ok": True}

    _fake_radia(monkeypatch, entry)
    levels = [{"q_lower": 1.0, "q_upper": 2.0}, {"q_lower": 1.4, "q_upper": 1.6}]
    coarse = {"q_lower": 12.0, "q_upper": 8.0}
    assert server.dual_fem_reference_estimate(levels, coarse=coarse, convergence_fraction=0.2) == {"ok": True}
    assert seen == {"levels": levels, "coarse": coarse, "fraction": 0.2}


def test_adapter_propagates_entry_value_error(monkeypatch):
    def entry(levels, coarse=None, convergence_fraction=0.1):
        raise ValueError("coarse needs q_lower <= q_upper")

    _fake_radia(monkeypatch, entry)
    with pytest.raises(ValueError, match="coarse needs q_lower <= q_upper"):
        server.dual_fem_reference_estimate(
            [{"q_lower": 9.9, "q_upper": 10.1}] * 2, coarse={"q_lower": 12.0, "q_upper": 8.0})
