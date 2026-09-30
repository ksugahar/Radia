"""The MCP direct gate must not accept a large residual through a huge solution."""
import numpy as np
import pytest
from radia_mcp.radia_ngsolve._direct import check_residual

class Identity:
    def CSR(self):
        return np.ones(1), np.zeros(1, dtype=int), np.array([0, 1])


@pytest.mark.parametrize('solution,residual,accepted', [
    (1., 1e-10, True), (1e12, 1e-7, True), (1e12, 2e-6, False), (float('inf'), 0., False),
    (1., float('nan'), False),
])
def test_strict_true_residual_gate(solution, residual, accepted):
    x, r, b = (np.array([v]) for v in (solution, residual, 1.))
    if accepted:
        assert check_residual(Identity(), r, x, b, np.array([True]), 'test') == residual
    else:
        with pytest.raises(RuntimeError, match='residual'):
            check_residual(Identity(), r, x, b, np.array([True]), 'test')
