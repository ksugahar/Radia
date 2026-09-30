"""Explicit accelerator solver requests must not select a different backend."""
import pytest

from radia.em_design import EMDesignSpec
from radia.panels.calc_accel_magnet import _select_accel_solver, build_argparser


@pytest.mark.parametrize("order", [1, 2, 3])
def test_ams_request_is_not_an_alias_for_bddc(order):
    parser = build_argparser()
    with pytest.raises(SystemExit):
        parser.parse_args(["--coil-script", "coil.py",
                           "--solver", "ams", "--formulation", "a",
                           "--fes-order", str(order)])
    with pytest.raises(ValueError, match="Unsupported FE solver"):
        EMDesignSpec(coil_script="coil.py", solver="ams").build_command()
    with pytest.raises(ValueError, match="Unsupported FE solver"):
        _select_accel_solver("ams", "a", order, False, 100)
    assert _select_accel_solver("bddc_ams", "a", order, False, 100) == "bddc_ams"
    assert _select_accel_solver("auto", "a", order, False, 100) == "bddc_ams"


@pytest.mark.parametrize("form,order,periodic", [
    ("a", 1, True), ("a", 4, False), ("omega", 2, False),
])
def test_bddc_ams_rejects_unsupported_spaces(form, order, periodic):
    with pytest.raises(ValueError, match="nonperiodic HCurl"):
        _select_accel_solver("bddc_ams", form, order, periodic, 100)


def test_other_explicit_backends_are_preserved():
    assert _select_accel_solver("sparsecholesky", "omega", 2, False, 900000) == "sparsecholesky"
    assert _select_accel_solver("bddc", "omega", 2, False, 100) == "bddc"
    with pytest.raises(ValueError, match="Unsupported FE solver"):
        EMDesignSpec(coil_script="coil.py", solver="pardiso").build_command()
