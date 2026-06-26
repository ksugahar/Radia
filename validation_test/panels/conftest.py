"""
Shared pytest fixtures for the Radia panel validation layer.

The panel tests run real PySide6 widgets headless via the
``offscreen`` Qt platform plugin. A single QApplication is created
session-wide so individual tests can instantiate IHPanel /
EMPanel / etc. without leaking widgets across runs.

Run from the repo root::

    pytest validation_test/panels/

CI sets ``QT_QPA_PLATFORM=offscreen`` automatically; locally the
fixture sets it before importing PySide6 so the tests work on a
machine without an X server.  PySide6 is no longer a required Radia
runtime dependency, so PySide-specific validations are skipped when it
is not installed.
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

if not _HAVE_PYSIDE6:
    collect_ignore = [
        "test_build_command_parses.py",   # imports panel_qa at module load
        "test_combo_wheel_guard.py",
        "test_panel_output_health.py",
        "test_panel_qa.py",               # imports panel_qa at module load
        "test_radia_export_menu.py",
        "test_subprocess_failure_ux.py",
    ]

# Make sure ``import radia_ih`` resolves to ``src/radia/radia_ih.py``
# (the panel modules are not part of the importable ``radia`` package
# — register_toolbar.py adds the panels dir to sys.path at runtime).
_REPO = Path(__file__).resolve().parents[2]
_RADIA = _REPO / "src" / "radia"
_PANELS = _RADIA / "panels"
_VALIDATION_PANELS = Path(__file__).resolve().parent
for p in (_RADIA, _PANELS, _VALIDATION_PANELS):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))


@pytest.fixture(scope="session")
def qapp():
    """Session-scoped QApplication on the offscreen platform.

    Importing PySide6 binds to whichever QPA plugin is in env at
    import time, so we set the env var FIRST and import lazily.

    Applies the lab-standard panel font baseline so panel_qa
    font-size checks see the same QApplication font as the real
    runtime (apply_panel_base_font runs in run_app).
    """
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    if not _HAVE_PYSIDE6:
        pytest.skip("PySide6 is not installed; legacy Qt panel validation skipped")
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    from radia_gui_base import apply_panel_base_font
    apply_panel_base_font(app)
    yield app
    # Do NOT call app.quit() — pytest may share the QApp across
    # tests and the next file's first instantiation would crash.


@pytest.fixture
def ih_panel(qapp):
    """Fresh IHPanel for each test (no leaked widget state)."""
    from radia_ih import IHPanel
    panel = IHPanel()
    yield panel
    panel.deleteLater()


@pytest.fixture
def sf_panel(qapp):
    """Fresh StreamFunctionPanel (Design / Pareto / Manufacture) per test.
    Constructing it does NOT import ngsolve / radia (the calc_streamfunction
    argparser is pure argparse), so this runs inside pytest unlike the
    subprocess calc golden."""
    from radia_streamfunction import StreamFunctionPanel
    panel = StreamFunctionPanel()
    yield panel
    panel.deleteLater()


@pytest.fixture
def ih_window(qapp, tmp_path):
    """Fresh IHWindow with an empty .vol path. Useful when the test
    needs the full Run / Stop / Open GMSH button machinery from
    AnalysisWindow, not just the IHPanel widgets."""
    from radia_ih import IHWindow
    win = IHWindow(vol_path="")
    yield win
    win.deleteLater()
