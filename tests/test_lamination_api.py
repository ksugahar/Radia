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


@pytest.mark.parametrize("mu_r,fill", [(1.0, 0.6), (5000.0, 0.875), (0.4, 0.35)])
def test_normal_lamination_matches_layer_flux_continuity_and_energy(mu_r, fill):
    # Independently solve the two layer H values from continuity and applied MMF.
    import numpy as np

    thickness = 0.008
    h_average = 250.0
    layer_mu = np.array([MU0 * mu_r, MU0])
    lengths = thickness * np.array([fill, 1.0 - fill])
    h_layers = np.linalg.solve(
        [[mu_r, -1.0], [fill, 1.0 - fill]],
        [0.0, h_average],
    )
    b_layers = layer_mu * h_layers
    effective = LAMINATION.laminated_mu_perpendicular(mu_r, fill)
    np.testing.assert_allclose(b_layers, effective * h_average, rtol=2e-14)
    energy_per_area = 0.5 * np.sum(layer_mu * h_layers**2 * lengths)
    assert energy_per_area == pytest.approx(
        0.5 * effective * h_average**2 * thickness, rel=2e-14
    )


def test_normal_lamination_pure_layers_and_parallel_bound():
    assert LAMINATION.laminated_mu_perpendicular(2000.0, 0.0) == MU0
    assert LAMINATION.laminated_mu_perpendicular(2000.0, 1.0) == MU0 * 2000.0
    normal = LAMINATION.laminated_mu_perpendicular(2000.0, 0.9)
    parallel = laminated_mu_eff(2000.0, 0.0, 0.0, 0.001, 0.9).real
    assert MU0 < normal < parallel < MU0 * 2000.0


@pytest.mark.parametrize("mu_r,fill", [
    (0.0, 0.5), (-1.0, 0.5), (math.nan, 0.5), (math.inf, 0.5),
    (1000.0, -0.1), (1000.0, 1.1), (1000.0, math.nan), (1000.0, math.inf),
])
def test_normal_lamination_rejects_nonphysical_layers(mu_r, fill):
    with pytest.raises(ValueError):
        LAMINATION.laminated_mu_perpendicular(mu_r, fill)


def test_normal_lamination_public_module_import():
    pytest.importorskip('ngsolve')
    from radia.lamination import laminated_mu_perpendicular

    assert laminated_mu_perpendicular(2000.0, 0.9) == LAMINATION.laminated_mu_perpendicular(2000.0, 0.9)
