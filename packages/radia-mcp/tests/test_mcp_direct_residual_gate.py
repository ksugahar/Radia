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

def test_gate_operates_without_radia_in_fresh_process():
    import os
    from pathlib import Path
    import subprocess
    import sys
    source = Path(__file__).resolve().parents[1] / 'src'
    code = r'''
import importlib.abc, sys
class NoRadia(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "radia" or fullname.startswith("radia."):
            raise AssertionError("Standalone MCP gate imported Radia: " + fullname)
sys.meta_path.insert(0, NoRadia())
import numpy as np
from radia_mcp.radia_ngsolve._direct import check_residual
from radia_mcp.radia_ngsolve._vendor.residual_gate import RELATIVE_LIMIT
assert RELATIVE_LIMIT == 1e-6
class Matrix:
    def CSR(self):
        return np.ones(1), np.zeros(1, dtype=int), np.array([0, 1])
free = np.array([True])
b = np.array([1e-30]); x = np.ones(1)
assert check_residual(Matrix(), np.array([5e-7]), x, b, free, 'tiny', reference_norm=1.) == 5e-7
try:
    check_residual(Matrix(), np.array([2e-6]), x, b, free, 'reject', reference_norm=1.)
except RuntimeError:
    pass
else:
    raise AssertionError('large true residual accepted')
assert 'radia' not in sys.modules
'''
    env = dict(os.environ, PYTHONPATH=str(source), OPENBLAS_NUM_THREADS='1')
    completed = subprocess.run([sys.executable, '-c', code], env=env,
                               capture_output=True, text=True, timeout=30)
    assert completed.returncode == 0, completed.stdout + completed.stderr
