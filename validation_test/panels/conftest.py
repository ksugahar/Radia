"""Shared pytest fixtures for the Radia application validation layer.

Production human interfaces are Simulink blocks. Notebook workbenches and the
old desktop ``radia_*.py`` PySide6 panels remain retired; Cubit's embedded
toolbar is a separate surface.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest


def _have_pyside6() -> bool:
    try:
        from PySide6 import QtCore as _qtcore  # noqa: F401
        return True
    except ImportError:
        return False


_HAVE_PYSIDE6 = _have_pyside6()

_RETIRED_DESKTOP_PANEL_TESTS = [
    "test_build_command_parses.py",   # imports retired panel_qa at module load
    "test_combo_wheel_guard.py",
    "test_ih_panel_qt.py",
    "test_open_gmsh_button.py",
    "test_panel_output_health.py",
    "test_panel_qa.py",               # imports retired panel_qa at module load
    "test_panel_state_restore.py",
    "test_run_button_browse.py",
    "test_streamfunction_panel_qt.py",
    "test_subprocess_failure_ux.py",
]

collect_ignore = list(_RETIRED_DESKTOP_PANEL_TESTS)
if not _HAVE_PYSIDE6:
    # Current Cubit toolbar test.  It is allowed to use PySide6, but normal
    # Radia Python environments do not need that dependency.
    collect_ignore.append("test_radia_export_menu.py")

# Make local panel calc helpers importable for validation tests.
_REPO = Path(__file__).resolve().parents[2]
_RADIA = _REPO / "src" / "radia"
_PANELS = _RADIA / "panels"
_VALIDATION_PANELS = Path(__file__).resolve().parent
for p in (_RADIA, _PANELS, _VALIDATION_PANELS):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))


@pytest.fixture(scope="session")
def qapp():
    """Session-scoped QApplication for the Cubit-toolbar test only."""
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    if not _HAVE_PYSIDE6:
        pytest.skip("PySide6 is not installed; Cubit toolbar validation skipped")
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    yield app
    # Do NOT call app.quit() — pytest may share the QApp across
    # tests and the next file's first instantiation would crash.


@pytest.fixture
def ih_panel(qapp):
    """Retired desktop IH panel fixture."""
    pytest.skip("retired desktop IH panel was removed; use the Simulink block")


@pytest.fixture
def sf_panel(qapp):
    """Retired desktop stream-function panel fixture."""
    pytest.skip("retired desktop stream-function panel was removed; use the Simulink block")


@pytest.fixture
def ih_window(qapp, tmp_path):
    """Retired desktop IH window fixture."""
    pytest.skip("retired desktop IH window was removed; use the Simulink block")
