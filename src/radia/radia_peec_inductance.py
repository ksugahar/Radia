"""Radia PEEC inductance analysis — coil-only, STEP-only, no mesh.

Surfaces the "PEEC inductance (coil only, STEP)" IH method as a
first-class Analysis choice in the Radia-NGSolve launcher, so the
user does not have to:
  1. pick "Induction Heating" (which forces a .vol export), then
  2. notice the Method combo inside IHWindow, then
  3. ignore the already-exported .vol.

Instead:
  Solve -> Radia-NGSolve -> Analysis = "PEEC inductance (coil only, STEP)"
  -> (no .vol export, no mesh) -> IHWindow opens pre-set to the right Method

Module-level contract with the Cubit C++ launcher (.ccl):
  * TITLE           — shown in the Analysis combo
  * REQUIRED_LABELS — none (PEEC takes STEP, not .vol)
  * OPTIONAL_LABELS — none
  * NEEDS_VOL       — False: the launcher must skip the `.vol` row
                      and the `radia_export netgen` call
  * main()          — entry point, invoked by run_app
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from radia_gui_base import run_app
from radia_ih import IHWindow, METHOD_PEEC_IND


TITLE = "PEEC inductance (coil only, STEP)"
REQUIRED_LABELS = []
OPTIONAL_LABELS = []
OPTIONAL_FILES = {}
# Launcher contract: this analysis takes a STEP file, not a .vol
# mesh.  The C++ launcher (`RadiaComp.cpp`) must not insist on a
# .vol path and must not call `radia_export netgen` before opening
# the window.
NEEDS_VOL = False


class PEECInductanceWindow(IHWindow):
    """IHWindow pre-set to the PEEC-inductance method.

    Inherits all of IHWindow's plumbing (Method combo, widget
    visibility rules, command dispatch).  Just forces the method
    so the user does not have to find and switch it manually.
    """

    def __init__(self, vol_path=""):
        super().__init__(vol_path)
        # Force the PEEC-inductance method and re-run the visibility
        # hooks (workpiece rows hidden, .vol row hidden, n_peri shown).
        self._panel._method_combo.setCurrentText(METHOD_PEEC_IND)
        self._panel._on_method_changed(METHOD_PEEC_IND)
        # Auto-fill the STEP path from the Cubit working directory
        # (the .jou's dir, set by the .ccl launcher via startDetached).
        # If the user just played a .jou that did `export step`, the
        # freshest `.step` there is very likely the coil they want.
        # Falls back silently when nothing is found.
        self._auto_fill_step_from_cwd()

    def _auto_fill_step_from_cwd(self):
        """Populate the STEP Browse field with the newest *.step / *.stp
        in the current working directory, if any.  The user can still
        override via Browse..."""
        import glob
        candidates = sorted(glob.glob("*.step") + glob.glob("*.stp"),
                             key=lambda p: os.path.getmtime(p),
                             reverse=True)
        if not candidates:
            return
        newest = os.path.abspath(candidates[0])
        step_widget = self._panel._widgets.get("peec_step")
        if step_widget is not None:
            step_widget.setText(newest)


def main():
    run_app(PEECInductanceWindow)


if __name__ == "__main__":
    main()
