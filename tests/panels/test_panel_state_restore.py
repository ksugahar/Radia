"""
Save / restore tests for ModePanel state.

The 2026-04-12 "empty Method combo" bug came from saving combo
state by index.  After dropping "BEM-SIBC (WP)" from the Method
combo, a stale saved index 2 (= the dropped item in the old
3-item combo) became out-of-range on the new 2-item combo and
``setCurrentIndex(2)`` silently set the selection to -1, leaving
the user with a blank widget.

These tests pin the new behaviour:

  1. save_state stores combo selection by TEXT
  2. restore_state with TEXT looks up via findText()
  3. restore_state with legacy INT bounds-checks against
     current item count and falls back to the panel default

Refreshed 2026-04-26 to use the canonical method names from
radia_ih.METHOD_* constants (was hardcoded to the obsolete 2-item
combo strings "PEEC+FEM"/"FEM").  Round-trip widgets updated to
the post-2026-04-19 panel structure (impedance_model in place of
the retired workpiece_mode dropdown).
"""

from __future__ import annotations

import pytest


# Pull the canonical method labels from the panel module so the
# tests track the panel's own constants instead of duplicating
# strings that drift.
@pytest.fixture(scope="module")
def methods():
    from radia_ih import (METHOD_PEEC_IND, METHOD_PEEC_BEM,
                          METHOD_PEEC_FEM_KELVIN, METHOD_FEM_FULL)
    return {
        "PEEC_IND":    METHOD_PEEC_IND,
        "PEEC_BEM":    METHOD_PEEC_BEM,
        "PEEC_FEM_K":  METHOD_PEEC_FEM_KELVIN,
        "FEM_FULL":    METHOD_FEM_FULL,
    }


class TestComboSaveByText:

    def test_save_returns_text_not_index(self, ih_panel, methods):
        ih_panel._method_combo.setCurrentText(methods["FEM_FULL"])
        state = ih_panel.save_state()
        assert state["method"] == methods["FEM_FULL"]  # text, not 3
        ih_panel._method_combo.setCurrentText(methods["PEEC_IND"])
        state = ih_panel.save_state()
        assert state["method"] == methods["PEEC_IND"]

    def test_restore_text_lookup(self, ih_panel, methods):
        ih_panel._method_combo.setCurrentText(methods["PEEC_IND"])
        ih_panel.restore_state({"method": methods["FEM_FULL"]})
        assert ih_panel._method_combo.currentText() == methods["FEM_FULL"]

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

    def test_legacy_index_in_range(self, ih_panel, methods):
        target_index = ih_panel._method_combo.findText(methods["FEM_FULL"])
        assert target_index >= 0
        ih_panel.restore_state({"method": target_index})
        assert ih_panel._method_combo.currentText() == methods["FEM_FULL"]

    def test_legacy_index_out_of_range_keeps_default(self, ih_panel):
        """Index 99 (way past the end) must NOT blank the selection.
        The widget stays non-blank at its panel default."""
        before = ih_panel._method_combo.currentText()
        ih_panel.restore_state({"method": 99})
        after = ih_panel._method_combo.currentText()
        assert after == before
        assert after != ""

    def test_legacy_negative_index_keeps_default(self, ih_panel):
        before = ih_panel._method_combo.currentText()
        ih_panel.restore_state({"method": -1})
        after = ih_panel._method_combo.currentText()
        assert after == before
        assert after != ""


class TestRoundTrip:

    def test_save_then_restore_preserves_values(self, ih_panel, methods):
        # PEEC_BEM exposes the workpiece widgets we want to round-trip.
        ih_panel._method_combo.setCurrentText(methods["PEEC_BEM"])
        ih_panel._widgets["impedance_model"].setCurrentIndex(1)  # ESIM
        ih_panel._widgets["wp_sigma"].setText("3.5e7")
        ih_panel._widgets["mu_r"].setText("250")
        state = ih_panel.save_state()

        from radia_ih import IHPanel
        p2 = IHPanel()
        p2.restore_state(state)

        assert p2._method_combo.currentText() == methods["PEEC_BEM"]
        assert p2._widgets["impedance_model"].currentIndex() == 1
        assert p2._widgets["wp_sigma"].text() == "3.5e7"
        assert p2._widgets["mu_r"].text() == "250"
        p2.deleteLater()

    def test_unknown_key_in_state_is_ignored(self, ih_panel, methods):
        """A saved key that no longer exists in the panel must
        not raise -- restore should silently skip unknown widgets."""
        ih_panel.restore_state({
            "method": methods["PEEC_IND"],
            "removed_widget_xyz": "garbage",
        })
        assert ih_panel._method_combo.currentText() == methods["PEEC_IND"]
