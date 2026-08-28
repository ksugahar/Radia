"""Long coil-from-CAD resampling and isolated-import regressions."""

import pytest

from tests.coil_from_cad.test_intra_package_imports import (
    _validate_filaments_from_step_runs_only_src_on_path,
    _validate_fixed_modules_import_only_src_on_path,
)
from tests.coil_from_cad.test_keiko_outsideline_centerline import (
    _validate_keiko_outsideline_centerline_covers_lead,
    _validate_keiko_outsideline_succeeds_with_adaptive_resampling,
)

pytestmark = pytest.mark.slow


def test_fixed_modules_import_only_src_on_path():
    _validate_fixed_modules_import_only_src_on_path()


def test_filaments_from_step_runs_only_src_on_path(monkeypatch):
    _validate_filaments_from_step_runs_only_src_on_path(monkeypatch)


def test_keiko_outsideline_centerline_covers_lead():
    _validate_keiko_outsideline_centerline_covers_lead()


def test_keiko_outsideline_succeeds_with_adaptive_resampling(monkeypatch):
    _validate_keiko_outsideline_succeeds_with_adaptive_resampling(monkeypatch)
