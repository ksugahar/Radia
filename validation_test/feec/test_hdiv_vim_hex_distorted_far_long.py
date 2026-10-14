"""Long subprocess A/B regression for distorted HDiv far entries."""

import pytest

from tests.feec.test_hdiv_vim_hex_distorted_far import (
    _validate_far_switch_matches_general_path_entries,
)

pytestmark = pytest.mark.slow


def test_far_switch_matches_general_path_entries():
    _validate_far_switch_matches_general_path_entries()
