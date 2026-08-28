"""Long fourth-order Lie-map optimization regressions."""

import pytest

from tests.test_accelerator_lie_topopt import (
    _validate_direct_As_Ay_polynomials_generate_fourth_order_lie_map_with_ad,
    _validate_fourth_order_lie_objective_and_material_pipeline_control_v_term,
)

pytestmark = pytest.mark.slow


def test_direct_As_Ay_polynomials_generate_fourth_order_lie_map_with_ad():
    _validate_direct_As_Ay_polynomials_generate_fourth_order_lie_map_with_ad()


def test_fourth_order_lie_objective_and_material_pipeline_control_v_term():
    _validate_fourth_order_lie_objective_and_material_pipeline_control_v_term()
