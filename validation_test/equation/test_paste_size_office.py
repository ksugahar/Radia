"""Real-PowerPoint validation for the equation paste-size contract."""

import pytest

from tests.equation.test_paste_size import _validate_powerpoint_pastes_it_at_24_pt

pytestmark = pytest.mark.slow


def test_powerpoint_pastes_it_at_24_pt():
    _validate_powerpoint_pastes_it_at_24_pt()
