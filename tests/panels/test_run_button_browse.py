"""Regression tests for the Run button enable / disable rule in
``IHWindow`` after the .vol path changes via Browse... or manual edit.

History
-------
2026-05-08: kubota reported on mdx + 100号機 (PyPI 4.28.0) that the
Run button stayed grayed out for every non-vacuum method (PEEC+BEM,
BEM-A+BEM, PEEC+FEM+Kelvin, Full FEM A-V) when ``radia_ih`` was
launched directly (no ``--vol`` arg) and a .vol picked via the
Browse... dialog.

Two parallel bugs caused the regression:

1. ``AnalysisWindow._browse_vol`` updated the .vol line edit text but
   never re-inspected the file.  ``_vol_mats`` stayed ``None``,
   ``is_runnable()`` returned ``False``, the Run button never enabled.

2. ``IHWindow.__init__`` called ``_reload_vol_info(vol_path)`` using
   the constructor argument, so even when ``_restore_settings()``
   repopulated the line edit from the saved ``radia_ih.json``, the
   QSettings-restored path was never inspected.

Fix shipped in v4.28.1 (commit edf583cb): added overridable
``_on_vol_changed(path)`` hook in ``AnalysisWindow``, invoked from
``_browse_vol`` and ``QLineEdit.editingFinished``.  ``IHWindow``
overrides the hook to re-run ``_reload_vol_info`` +
``_update_run_state``.  ``IHWindow.__init__`` now passes
``self._vol_edit.text()`` (post-restore) to the initial inspect.

The tests in this module pin both halves so a future refactor cannot
silently regress the workflow kubota uses every day.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[2]
SAMPLE_VOL = REPO / "tests" / "panels" / "test_3d_sibc_copper.vol"


@pytest.fixture
def isolated_settings(tmp_path, monkeypatch):
    """Redirect _SETTINGS_DIR to a tmp dir so QSettings round-trip
    tests don't pollute the developer's real ``~/.radia/radia_ih.json``
    and aren't influenced by whatever the last interactive session
    wrote there."""
    import radia_gui_base
    monkeypatch.setattr(radia_gui_base, "_SETTINGS_DIR", str(tmp_path))
    return tmp_path


@pytest.fixture
def ih_window_clean(qapp, isolated_settings):
    """IHWindow with isolated QSettings — fresh state every test."""
    from radia_ih import IHWindow
    win = IHWindow(vol_path="")
    yield win
    win.deleteLater()


# ============================================================
# Run button enable rule — non-vacuum methods need .vol with labels
# ============================================================

class TestRunButtonAfterBrowse:
    """Pins the kubota 2026-05-08 bug: Browse... must re-inspect the
    .vol so the Run button can re-enable when the user changes path."""

    @pytest.mark.parametrize("method_const", [
        "METHOD_PEEC_BEM",
        "METHOD_BEMA_BEM",
        "METHOD_PEEC_FEM_KELVIN",
    ])
    def test_run_disabled_on_fresh_launch_with_empty_vol(
            self, ih_window_clean, method_const):
        """Empty vol_path + non-vacuum method => Run is grayed out.

        This reproduces kubota's symptom (the bug + the expected
        post-fix behaviour at the entry point — both are "Run is
        disabled before any .vol is picked").
        """
        import radia_ih
        method = getattr(radia_ih, method_const)
        ih_window_clean._panel._method_combo.setCurrentText(method)
        assert ih_window_clean._panel._vol_mats is None
        assert ih_window_clean._panel.is_runnable() is False
        assert ih_window_clean._run_btn.isEnabled() is False

    @pytest.mark.parametrize("method_const", [
        "METHOD_PEEC_BEM",
        "METHOD_BEMA_BEM",
        "METHOD_PEEC_FEM_KELVIN",
        "METHOD_FEM_FULL",
    ])
    def test_browse_to_valid_vol_enables_run(
            self, ih_window_clean, method_const):
        """Picking a .vol that has every label required by the method
        (sibc, and source/sink/coil for FEM_FULL) re-enables Run.

        The bug pre-v4.28.1: this transition NEVER happened because
        ``_browse_vol`` did not call ``_reload_vol_info``.  With the
        fix, ``_on_vol_changed`` runs the inspect + ``_update_run_state``
        chain and the button enables on the same click that picks the
        file.
        """
        import radia_ih
        method = getattr(radia_ih, method_const)
        ih_window_clean._panel._method_combo.setCurrentText(method)

        # Simulate the Browse... dialog committing a path: line edit
        # gets the text + the AnalysisWindow hook fires.
        ih_window_clean._vol_edit.setText(str(SAMPLE_VOL))
        ih_window_clean._on_vol_changed(str(SAMPLE_VOL))

        assert ih_window_clean._panel._vol_mats == {"coil", "kelvin", "air"}
        assert ih_window_clean._panel._vol_bnds == {
            "sibc", "source", "sink", "kelvin_int", "kelvin_ext", "default",
        }
        assert ih_window_clean._panel.is_runnable() is True, (
            f"is_runnable() must be True after Browse for method={method!r} "
            "on a .vol with all required labels")
        assert ih_window_clean._run_btn.isEnabled() is True, (
            "Run button must reflect is_runnable() — kubota's symptom is "
            "this assertion failing for every non-vacuum method")

    def test_browse_to_nonexistent_vol_keeps_run_disabled(
            self, ih_window_clean, tmp_path):
        """A path that does not exist on disk must drop ``_vol_mats``
        back to ``None`` and re-disable Run.  Without this, a typo in
        the Browse field could leave a stale OK badge from a previous
        valid path."""
        import radia_ih
        ih_window_clean._panel._method_combo.setCurrentText(
            radia_ih.METHOD_PEEC_BEM)

        # First a valid path (sets _vol_mats).
        ih_window_clean._on_vol_changed(str(SAMPLE_VOL))
        assert ih_window_clean._run_btn.isEnabled() is True

        # Then a nonexistent path — the hook must reset state.
        bogus = tmp_path / "does_not_exist.vol"
        ih_window_clean._on_vol_changed(str(bogus))
        assert ih_window_clean._panel._vol_mats is None
        assert ih_window_clean._run_btn.isEnabled() is False


# ============================================================
# QSettings round-trip — _restore_settings + initial inspect
# ============================================================

class TestQSettingsRoundTrip:
    """Pins the second half of the kubota bug: ``IHWindow.__init__``
    must inspect ``self._vol_edit.text()`` (post-restore) instead of
    the constructor ``vol_path`` argument.

    Pre-fix: even when the user's previous session wrote a valid .vol
    path to ``radia_ih.json``, the next launch ignored it and Run
    stayed disabled until the user clicked Browse...  That made the
    Browse... bug more visible (saved state was never honoured).
    """

    def test_init_uses_restored_vol_path_not_constructor_arg(
            self, qapp, isolated_settings):
        """Write a saved-state JSON pointing at SAMPLE_VOL, then launch
        IHWindow with empty constructor arg.  Run must enable for
        PEEC+BEM because the restored line-edit text is inspected."""
        from radia_ih import IHWindow, METHOD_PEEC_BEM

        # Pre-seed QSettings with the canonical IH save_state schema.
        json_path = Path(isolated_settings) / "radia_ih.json"
        json_path.write_text(json.dumps({
            "vol": str(SAMPLE_VOL),
            "panel": {"method": METHOD_PEEC_BEM},
        }), encoding="utf-8")

        win = IHWindow(vol_path="")  # empty arg, restore-from-disk path
        try:
            assert win._vol_edit.text() == str(SAMPLE_VOL), (
                "_restore_settings must repopulate the line edit")
            assert win._panel._vol_mats == {"coil", "kelvin", "air"}, (
                "_reload_vol_info must inspect the RESTORED path, not "
                "the empty constructor arg (the kubota bug pre-v4.28.1)")
            assert win._run_btn.isEnabled() is True
        finally:
            win.deleteLater()

    def test_save_browse_restore_round_trip(
            self, qapp, isolated_settings):
        """Full lifecycle: launch empty -> Browse to a .vol -> save ->
        relaunch -> Run is enabled without further action."""
        from radia_ih import IHWindow, METHOD_PEEC_BEM

        # Session 1: launch empty, simulate Browse, then save.
        win1 = IHWindow(vol_path="")
        try:
            win1._panel._method_combo.setCurrentText(METHOD_PEEC_BEM)
            win1._vol_edit.setText(str(SAMPLE_VOL))
            win1._on_vol_changed(str(SAMPLE_VOL))
            assert win1._run_btn.isEnabled() is True
            win1._save_settings()
        finally:
            win1.deleteLater()

        # Session 2: relaunch empty, expect saved state honoured.
        win2 = IHWindow(vol_path="")
        try:
            assert win2._vol_edit.text() == str(SAMPLE_VOL)
            assert win2._panel._method_combo.currentText() == METHOD_PEEC_BEM
            assert win2._run_btn.isEnabled() is True
        finally:
            win2.deleteLater()


# ============================================================
# AnalysisWindow base class — _on_vol_changed hook contract
# ============================================================

class TestOnVolChangedHook:
    """The override hook is the entry point for the fix; verify the
    contract so a future refactor cannot remove it without breaking
    callers."""

    def test_analysis_window_default_hook_is_noop(self, qapp):
        """AnalysisWindow.AnalysisWindow's default hook is a no-op so
        non-IH panels (radia_em, radia_pcb, radia_motor) keep their
        existing behaviour."""
        from radia_gui_base import AnalysisWindow
        # The base class declares the hook with default no-op semantics.
        assert hasattr(AnalysisWindow, "_on_vol_changed")
        # Calling it on the base class with any args must not raise.
        win = AnalysisWindow("test", "")
        try:
            win._on_vol_changed("anything")
        finally:
            win.deleteLater()

    def test_ih_window_overrides_hook(self, ih_window_clean):
        """IHWindow must override _on_vol_changed (not inherit the
        no-op).  Without the override the .vol line edit and
        ``_vol_mats`` desynchronise on Browse..."""
        from radia_ih import IHWindow
        from radia_gui_base import AnalysisWindow
        assert (IHWindow._on_vol_changed
                is not AnalysisWindow._on_vol_changed)

    def test_browse_vol_invokes_hook(self, ih_window_clean, monkeypatch):
        """``_browse_vol`` must call the hook with the chosen path so
        subclasses can react.  We monkey-patch QFileDialog to return a
        canned path and assert ``_on_vol_changed`` ran."""
        from PySide6.QtWidgets import QFileDialog
        called_with = []
        original = type(ih_window_clean)._on_vol_changed
        def spy(self, path):
            called_with.append(path)
            return original(self, path)
        monkeypatch.setattr(type(ih_window_clean), "_on_vol_changed", spy)

        canned = str(SAMPLE_VOL)
        monkeypatch.setattr(QFileDialog, "getOpenFileName",
                            staticmethod(lambda *a, **kw: (canned, "")))
        ih_window_clean._browse_vol()

        assert called_with == [canned], (
            "_browse_vol must invoke _on_vol_changed exactly once with "
            "the path the user selected (kubota bug: it didn't)")
        assert ih_window_clean._vol_edit.text() == canned
