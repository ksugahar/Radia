"""Fail-loud propagation of unconverged ESIM cell solves."""

from types import SimpleNamespace

import numpy as np
import pytest

from radia.esim_cell_problem import ESITable, require_esim_converged
from radia.esim_multiport import ESIMMultiportSolver


def test_require_converged_preserves_cell_diagnostics():
    result = {"converged": False, "iterations": 7,
              "relative_change": 2.5e-3, "Z": 1.0 + 2.0j}
    with pytest.raises(RuntimeError, match=r"iterations=7.*relative_change=0.0025"):
        require_esim_converged(result, "panel 4")
    accepted = {"converged": True, "iterations": 2, "Z": 3.0 + 4.0j}
    assert require_esim_converged(accepted) is accepted


def test_esi_table_rejects_unconverged_generated_rows():
    table = np.array([
        [1.0, 1.0, 2.0, 3.0, 4.0, 100.0, 1.0],
        [2.0, 1.1, 2.1, 3.1, 4.1, 110.0, 0.0],
    ])
    with pytest.raises(RuntimeError, match=r"rows \[1\].*H0=\[2.0\]"):
        ESITable(table_data=table)


def test_legacy_five_column_table_remains_loadable():
    table = np.array([
        [1.0, 1.0, 2.0, 3.0, 4.0],
        [2.0, 1.1, 2.1, 3.1, 4.1],
        [3.0, 1.2, 2.2, 3.2, 4.2],
        [4.0, 1.3, 2.3, 3.3, 4.3],
    ])
    loaded = ESITable(table_data=table)
    assert loaded.table.shape == (4, 5)


class _FakeCell:
    def __init__(self, converged):
        self.converged = converged

    def solve(self, _field, **_kwargs):
        return {"converged": self.converged, "iterations": 3,
                "relative_change": 0.1, "Z": 1.0 + 1.0j}


def test_multiport_rejects_the_specific_unconverged_port():
    solver = ESIMMultiportSolver.__new__(ESIMMultiportSolver)
    solver.ports = [SimpleNamespace(name="rolling"),
                    SimpleNamespace(name="transverse")]
    solver._solvers = [_FakeCell(True), _FakeCell(False)]
    with pytest.raises(RuntimeError, match=r"port 1 \(transverse\).*"):
        solver.solve([10.0, 20.0])
