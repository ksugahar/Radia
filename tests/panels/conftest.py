"""
Shared pytest fixtures for the Radia panel test layer.

The panel tests run real PySide6 widgets headless via the
``offscreen`` Qt platform plugin. A single QApplication is created
session-wide so individual tests can instantiate IHPanel /
AccelPanel / etc. without leaking widgets across runs.

Run from the repo root::

    pytest tests/panels/

CI sets ``QT_QPA_PLATFORM=offscreen`` automatically; locally the
fixture sets it before importing PySide6 so the tests work on a
machine without an X server.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

# Make sure ``import radia_ih`` resolves to ``src/radia/radia_ih.py``
# (the panel modules are not part of the importable ``radia`` package
# — register_toolbar.py adds the panels dir to sys.path at runtime).
_REPO = Path(__file__).resolve().parents[2]
_RADIA = _REPO / "src" / "radia"
_PANELS = _RADIA / "panels"
for p in (_RADIA, _PANELS):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))


@pytest.fixture(scope="session")
def qapp():
    """Session-scoped QApplication on the offscreen platform.

    Importing PySide6 binds to whichever QPA plugin is in env at
    import time, so we set the env var FIRST and import lazily.
    """
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
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
def ih_window(qapp, tmp_path):
    """Fresh IHWindow with an empty .vol path. Useful when the test
    needs the full Run / Stop / Open GMSH button machinery from
    AnalysisWindow, not just the IHPanel widgets."""
    from radia_ih import IHWindow
    win = IHWindow(vol_path="")
    yield win
    win.deleteLater()
