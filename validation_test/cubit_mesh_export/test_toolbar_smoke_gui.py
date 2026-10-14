"""Real-Cubit cold-start validation for the toolbar display contract."""

import os

import pytest

from cubit_mesh_export.toolbar_smoke import run_smoke_test

pytestmark = pytest.mark.slow


@pytest.mark.skipif(
    os.environ.get("CUBIT_MESH_EXPORT_RUN_GUI_TESTS") != "1",
    reason="set CUBIT_MESH_EXPORT_RUN_GUI_TESTS=1 to cold-start the real Cubit GUI",
)
def test_real_cubit_displays_toolbar_on_two_cold_starts():
    assert run_smoke_test(restarts=2, timeout=45.0, keep=True) == 0
