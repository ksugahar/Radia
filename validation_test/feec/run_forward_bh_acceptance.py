"""Installed-wheel acceptance of forward B(H) plus existing HDiv routes."""
import json
from pathlib import Path
import platform
import time
import pytest

root = Path.cwd()
start = time.perf_counter()
tests = [
    'validation_test/feec/test_hdiv_forward_bh.py',
    'validation_test/feec/test_hdiv_vim_demag_solve.py::test_rt1_linear_solve_is_physically_sane',
    'validation_test/feec/test_hdiv_vim_energy_newton.py::test_default_nonlinear_is_energy_newton_cpp',
]
status = pytest.main(['-q', '--basetemp', str(root/'pytest_tmp'), *tests])
evidence = [json.loads(p.read_text()) for p in (root/'pytest_tmp').rglob('result.json')]
(root/'result.json').write_text(json.dumps(dict(
    passed=status == 0, host=platform.node(), tests=tests,
    wall_s=time.perf_counter()-start, forward_evidence=evidence), indent=2))
raise SystemExit(status)
