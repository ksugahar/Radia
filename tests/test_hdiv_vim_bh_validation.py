"""Public semantic validation for the origin-anchored nonlinear BH route."""

import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve.meshes import MakeStructured3DMesh

from radia import vim


def _mesh():
    return MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1)


@pytest.mark.parametrize(
    "table, message",
    [
        ([[-1.0, 0.0], [0.0, 1.0]], "H=0"),
        ([[0.0, 0.1], [1.0, 1.0]], "B=0"),
        ([[0.0, 0.0], [1.0, 1.0], [1.0, 1.1]], "strictly increasing"),
        ([[0.0, 0.0], [1.0, 1.0], [2.0, 0.9]], "non-decreasing"),
    ],
)
def test_soft_iron_bh_table_rejects_noncanonical_curves(table, message):
    with pytest.raises(ValueError, match=message):
        vim.Solve(
            _mesh(),
            bh_table=table,
            H_ext=ng.CoefficientFunction((0.0, 0.0, 0.0)),
        )
