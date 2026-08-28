"""Long native Lie-kernel Jacobian regression."""

import pytest

from tests.test_lie_kernel_parameter_jacobians import (
    _validate_no_jacobian_mode_matches_and_zero_width,
)

pytestmark = pytest.mark.slow


def test_no_jacobian_mode_matches_and_zero_width():
    _validate_no_jacobian_mode_matches_and_zero_width()
