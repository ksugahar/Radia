"""Real-Cubit cold-start validation for the toolbar display contract."""

import os

import pytest

from tests.test_cubit_toolbar_smoke import (
    _validate_real_cubit_displays_toolbar_on_two_cold_starts,
)

pytestmark = pytest.mark.slow


@pytest.mark.skipif(
    os.environ.get("RADIA_RUN_CUBIT_GUI_TESTS") != "1",
    reason="set RADIA_RUN_CUBIT_GUI_TESTS=1 to cold-start the real Cubit GUI",
)
def test_real_cubit_displays_toolbar_on_two_cold_starts():
    _validate_real_cubit_displays_toolbar_on_two_cold_starts()
