"""Long EarlyTimes field-map regressions split from the fast test suite."""

import pytest

from tests.test_beam_gridfunction_native import (
    _validate_curvilinear_s_rk_uses_direct_hcurl_a_and_hdiv_b_map,
    _validate_design_orbit_gauge_preserves_the_earlytimes_lie_map,
    _validate_hcurl_volume_recovers_full_p5_xy_jet_and_lie_gradient,
    _validate_tracked_lie_map_self_contained_difference_against_field_rk,
    hcurl_cubic_vector_potential,
    hcurl_p5_xy_vector_potential,
)

pytestmark = pytest.mark.slow


def test_tracked_lie_map_self_contained_difference_against_field_rk(
    hcurl_cubic_vector_potential,
):
    _validate_tracked_lie_map_self_contained_difference_against_field_rk(
        hcurl_cubic_vector_potential
    )


def test_curvilinear_s_rk_uses_direct_hcurl_a_and_hdiv_b_map(
    hcurl_cubic_vector_potential,
):
    _validate_curvilinear_s_rk_uses_direct_hcurl_a_and_hdiv_b_map(
        hcurl_cubic_vector_potential
    )


def test_hcurl_volume_recovers_full_p5_xy_jet_and_lie_gradient(
    hcurl_p5_xy_vector_potential,
):
    _validate_hcurl_volume_recovers_full_p5_xy_jet_and_lie_gradient(
        hcurl_p5_xy_vector_potential
    )


def test_design_orbit_gauge_preserves_the_earlytimes_lie_map():
    _validate_design_orbit_gauge_preserves_the_earlytimes_lie_map()
