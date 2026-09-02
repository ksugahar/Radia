"""Fast contracts for the canonical laminated-steel material helper."""

from __future__ import annotations

import cmath
import importlib.util
import math
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "radia" / "lamination.py"
SPEC = importlib.util.spec_from_file_location("radia_lamination_contract", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
LAMINATION = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(LAMINATION)
MU0 = LAMINATION.MU0
laminated_mu_eff = LAMINATION.laminated_mu_eff


def test_laminated_mu_eff_static_limit():
    value = laminated_mu_eff(1000.0, 2.0e6, 0.0, 0.5e-3, fill=0.95)

    assert value == pytest.approx(MU0 * (0.95 * 1000.0 + 0.05))


def test_laminated_mu_eff_matches_closed_form():
    mu_r = 800.0
    sigma = 2.1e6
    omega = 2.0 * math.pi * 1000.0
    thickness = 0.35e-3
    fill = 0.96
    b = 0.5 * thickness * cmath.sqrt(1j * omega * MU0 * mu_r * sigma)
    expected = MU0 * (fill * mu_r * cmath.tanh(b) / b + 1.0 - fill)

    assert laminated_mu_eff(mu_r, sigma, omega, thickness, fill) == pytest.approx(
        expected
    )
    assert expected.imag < 0.0


def test_mcp_laminated_mu_eff_is_a_thin_adapter():
    native_modules = list((ROOT / "src" / "radia").glob("_radia_pybind.*"))
    if not native_modules:
        pytest.skip("MCP adapter parity requires the built Radia extension")
    pytest.importorskip("ngsolve")
    from radia_mcp.radia_ngsolve.solve import laminated_mu_eff as mcp_adapter

    arguments = (500.0, 1.8e6, 2.0 * math.pi * 400.0, 0.5e-3, 0.94)
    assert mcp_adapter(*arguments) == laminated_mu_eff(*arguments)
