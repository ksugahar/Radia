"""
Save / restore tests for ModePanel state.

The 2026-04-12 "empty Method combo" bug came from saving combo
state by index. After dropping "BEM-SIBC (WP)" from the Method
combo, a stale saved index 2 (= the dropped item in the old
3-item combo) became out-of-range on the new 2-item combo and
``setCurrentIndex(2)`` silently set the selection to -1, leaving
the user with a blank widget.

These tests pin the new behaviour:

  1. save_state stores combo selection by TEXT
  2. restore_state with TEXT looks up via findText()
  3. restore_state with legacy INT bounds-checks against
     current item count and falls back to the panel default
"""

from __future__ import annotations

import pytest


class TestComboSaveByText:

    def test_save_returns_text_not_index(self, ih_panel):
        ih_panel._method_combo.setCurrentText("FEM")
        state = ih_panel.save_state()
        assert state["method"] == "FEM"  # text, not 1
        ih_panel._method_combo.setCurrentText("PEEC+FEM")
        state = ih_panel.save_state()
        assert state["method"] == "PEEC+FEM"

    def test_restore_text_lookup(self, ih_panel):
        ih_panel._method_combo.setCurrentText("PEEC+FEM")
        ih_panel.restore_state({"method": "FEM"})
        assert ih_panel._method_combo.currentText() == "FEM"

    def test_restore_unknown_text_keeps_default(self, ih_panel):
        """A saved value that no longer exists in the combo MUST
        leave the widget at its panel-level default rather than
        blanking the selection."""
        before = ih_panel._method_combo.currentText()
        ih_panel.restore_state({"method": "BEM-SIBC (WP)"})
        after = ih_panel._method_combo.currentText()
        assert after == before  # default preserved
        assert after != ""      # never blank


class TestLegacyIndexRestore:
    """The new restore_state must still accept legacy int values
    for backward compatibility, but only when they are in range."""

    def test_legacy_index_in_range(self, ih_panel):
        # New combo: PEEC+FEM=0, FEM=1
        ih_panel.restore_state({"method": 1})
        assert ih_panel._method_combo.currentText() == "FEM"

    def test_legacy_index_out_of_range_keeps_default(self, ih_panel):
        """Index 2 from the old 3-item combo on the new 2-item combo
        must NOT call setCurrentIndex(2). The widget must stay
        non-blank at its panel default."""
        ih_panel._method_combo.setCurrentText("PEEC+FEM")
        ih_panel.restore_state({"method": 2})
        text = ih_panel._method_combo.currentText()
        assert text in ("PEEC+FEM", "FEM")
        assert text != ""

    def test_legacy_negative_index_keeps_default(self, ih_panel):
        ih_panel._method_combo.setCurrentText("PEEC+FEM")
        ih_panel.restore_state({"method": -1})
        assert ih_panel._method_combo.currentText() != ""


class TestRoundTrip:

    def test_save_then_restore_preserves_values(self, ih_panel):
        ih_panel._method_combo.setCurrentText("FEM")
        ih_panel._widgets["workpiece_mode"].setCurrentText("ESIM")
        ih_panel._widgets["wp_sigma"].setText("3.5e7")
        ih_panel._widgets["mu_r"].setText("250")
        state = ih_panel.save_state()

        # Make a fresh panel and restore
        from radia_ih import IHPanel
        p2 = IHPanel()
        p2.restore_state(state)

        assert p2._method_combo.currentText() == "FEM"
        assert p2._widgets["workpiece_mode"].currentText() == "ESIM"
        assert p2._widgets["wp_sigma"].text() == "3.5e7"
        assert p2._widgets["mu_r"].text() == "250"
        p2.deleteLater()

    def test_unknown_key_in_state_is_ignored(self, ih_panel):
        """A saved key that no longer exists in the panel must
        not raise — restore should silently skip unknown widgets."""
        ih_panel.restore_state({
            "method": "PEEC+FEM",
            "removed_widget_xyz": "garbage",
        })
        assert ih_panel._method_combo.currentText() == "PEEC+FEM"
