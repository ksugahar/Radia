"""Long ESIM numerical regressions split from the fast contracts."""

import pytest

from tests.test_esim_envelope import (
    _validate_peak_envelope_matches_bessel,
    cell_solver,
)
from tests.test_esim_integration import _validate_physical_consistency

pytestmark = pytest.mark.slow


def test_peak_envelope_matches_bessel(cell_solver):
    _validate_peak_envelope_matches_bessel(cell_solver)


def test_physical_consistency():
    _validate_physical_consistency()
