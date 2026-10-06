"""Execute the existing HDiv linear/nonlinear acceptance tests as a compute job."""
import json
from pathlib import Path
import pytest

root = Path(__file__).resolve().parents[1]
nodes = [
    "validation_test/feec/test_hdiv_vim_demag_solve.py::test_rt1_linear_solve_is_physically_sane",
    "validation_test/feec/test_hdiv_vim_energy_newton.py::test_default_nonlinear_is_energy_newton_cpp",
]
code = pytest.main(["-q", "-x", "--import-mode=importlib", "-p", "no:cacheprovider",
                    *[str(root / node) for node in nodes]])
(root / "result.json").write_text(json.dumps(dict(
    status="passed" if code == 0 else "failed", pytest_exitcode=int(code),
    purpose="Installed-wheel HDiv linear and nonlinear runtime acceptance; not customer-model accuracy",
    tests=nodes), indent=2))
raise SystemExit(code)
