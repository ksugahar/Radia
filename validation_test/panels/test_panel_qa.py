"""Regression tests for panel quality-assurance.

Runs every registered panel x mode combination through the checks in
``panel_qa.py``.  Fails on any single check regression.

Run::

    QT_QPA_PLATFORM=offscreen pytest validation_test/panels/test_panel_qa.py -v

Or via the ``/deploy`` skill (Stage 1 Step 8), which calls the same
``run_all_panel_checks`` entry point.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from panel_qa import (CheckResult, PanelQAError,
                       check_panel_health, get_panel_registry,
                       grab_screenshot)


_REPO = Path(__file__).resolve().parents[2]
_TEMP = _REPO / "temp"


def _force_mode(window, combo_ref, value):
    if combo_ref is None or value is None:
        return
    panel = window._panel
    combo = getattr(panel, combo_ref, None) \
            or (panel._widgets.get(combo_ref)
                if hasattr(panel, "_widgets") else None)
    if combo is not None:
        combo.setCurrentText(value)


@pytest.mark.parametrize("entry", get_panel_registry(),
                          ids=lambda e: e[0])
def test_panel_health(qapp, entry, tmp_path):
    tag, WinCls, combo_ref, value = entry
    window = WinCls("")
    try:
        _force_mode(window, combo_ref, value)
        results = check_panel_health(window)
        # Save screenshot next to the temp/ directory so failures can
        # be inspected after the test.
        _TEMP.mkdir(parents=True, exist_ok=True)
        grab_screenshot(window, _TEMP / f"panel_{tag}.png")
        fails = [r for r in results
                  if not r.ok and r.severity == "error"]
        assert not fails, (
            f"{tag}: failed checks: "
            + "\n".join(f"  - {r.name}: {r.detail}" for r in fails))
    finally:
        window.close()
        window.deleteLater()


def test_all_panels_together(qapp, tmp_path):
    """Sanity: the full runner completes without raising."""
    from panel_qa import run_all_panel_checks
    run_all_panel_checks(screenshot_dir=str(_TEMP), strict=True)
