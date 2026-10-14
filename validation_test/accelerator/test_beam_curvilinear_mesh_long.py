"""Long order-five curvilinear mesh regression."""

import pytest

from tests.test_beam_curvilinear_mesh import (
    _validate_p5_design_orbit_gauge_zeroes_As_Ay_without_changing_curl,
)

pytestmark = pytest.mark.slow


def test_p5_design_orbit_gauge_zeroes_As_Ay_without_changing_curl():
    _validate_p5_design_orbit_gauge_zeroes_As_Ay_without_changing_curl()
