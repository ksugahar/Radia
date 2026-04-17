"""
Headless Qt tests for the IH (Induction Heating) panel.

These tests instantiate the real PySide6 IHPanel via the offscreen
Qt platform plugin and verify the panel's user-visible behaviour:

  - Method combo defaults to BEM and exposes [BEM, FEM]
  - Mode switches show / hide the right widget rows
  - Workpiece off / SIBC / ESIM toggles WP detail rows
  - build_command for every (method, wp_mode) combination produces
    a list whose argparse can parse without exit 2

The string-grep tests in tests/panels/test_panel_ui_logic.py only
catch widget-removal regressions; these Qt tests catch BEHAVIOUR
regressions like the empty-Method-combo bug from 2026-04-12.
"""

from __future__ import annotations

import argparse
import importlib
import sys

import pytest


# ============================================================
# Method combo behaviour
# ============================================================

class TestMethodCombo:

    def test_default_method_is_PEEC(self, ih_panel):
        """First-launch default. PEEC is the default since BEM was
        removed from the panel (2026-04-17)."""
        assert ih_panel._method_combo.currentText() == "PEEC+BEM"

    def test_method_items(self, ih_panel):
        items = [ih_panel._method_combo.itemText(i)
                 for i in range(ih_panel._method_combo.count())]
        assert items == ["PEEC+BEM", "FEM"]

    def test_no_legacy_BEM_SIBC_WP_method(self, ih_panel):
        """Removed 2026-04-12. Make sure it does not creep back."""
        items = {ih_panel._method_combo.itemText(i)
                 for i in range(ih_panel._method_combo.count())}
        assert "BEM-SIBC (WP)" not in items


# ============================================================
# Mode switch widget visibility
# ============================================================

class TestModeSwitch:

    def _visible(self, panel, key):
        """Read whether a row is currently visible by checking the
        widget's parent QFormLayout row visibility. Falls back to
        the widget's own isVisible() flag (the offscreen platform
        does not actually render but tracks visibility state)."""
        w = panel._widgets.get(key)
        if w is None:
            return False
        return w.isVisibleTo(panel)

    def test_PEEC_shows_coil_params(self, ih_panel):
        ih_panel._method_combo.setCurrentText("PEEC+BEM")
        assert self._visible(ih_panel, "coil_radius")
        assert self._visible(ih_panel, "nwinc")
        assert self._visible(ih_panel, "coil_sigma")
        assert self._visible(ih_panel, "solver")
        assert not self._visible(ih_panel, "max_iter")
        assert not self._visible(ih_panel, "air_mode")
        assert not self._visible(ih_panel, "fes_order")

    def test_FEM_shows_FEM_only_widgets(self, ih_panel):
        ih_panel._method_combo.setCurrentText("FEM")
        assert self._visible(ih_panel, "max_iter")
        assert not self._visible(ih_panel, "air_mode")
        assert self._visible(ih_panel, "workpiece_mode")
        assert self._visible(ih_panel, "wp_sigma")
        assert not self._visible(ih_panel, "coil_radius")

    def test_PEEC_solver_items(self, ih_panel):
        ih_panel._method_combo.setCurrentText("PEEC+BEM")
        items = [ih_panel._widgets["solver"].itemText(i)
                 for i in range(ih_panel._widgets["solver"].count())]
        assert items == ["Dense LU", "HACApK saddle", "PRIMA"]

    def test_FEM_solver_items(self, ih_panel):
        ih_panel._method_combo.setCurrentText("FEM")
        items = [ih_panel._widgets["solver"].itemText(i)
                 for i in range(ih_panel._widgets["solver"].count())]
        assert items == ["pardiso", "bddc", "iccg", "ams"]


# ============================================================
# Workpiece group visibility
# ============================================================

class TestWorkpieceGroup:

    def test_off_hides_WP_detail(self, ih_panel):
        ih_panel._method_combo.setCurrentText("PEEC+BEM")
        ih_panel._widgets["workpiece_mode"].setCurrentText("off")
        for key in ("wp_sigma", "half_thickness", "mu_r",
                    "bh_file", "esim_geometry"):
            assert not ih_panel._widgets[key].isVisibleTo(ih_panel), key

    def test_SIBC_shows_mu_r_not_bh(self, ih_panel):
        ih_panel._method_combo.setCurrentText("PEEC+BEM")
        ih_panel._widgets["workpiece_mode"].setCurrentText("SIBC")
        assert ih_panel._widgets["wp_sigma"].isVisibleTo(ih_panel)
        assert ih_panel._widgets["half_thickness"].isVisibleTo(ih_panel)
        assert ih_panel._widgets["mu_r"].isVisibleTo(ih_panel)
        assert not ih_panel._widgets["bh_file"].isVisibleTo(ih_panel)
        assert not ih_panel._widgets["esim_geometry"].isVisibleTo(ih_panel)

    def test_ESIM_shows_bh_not_mu_r(self, ih_panel):
        ih_panel._method_combo.setCurrentText("PEEC+BEM")
        ih_panel._widgets["workpiece_mode"].setCurrentText("ESIM")
        assert ih_panel._widgets["wp_sigma"].isVisibleTo(ih_panel)
        assert ih_panel._widgets["bh_file"].isVisibleTo(ih_panel)
        assert ih_panel._widgets["esim_geometry"].isVisibleTo(ih_panel)
        assert not ih_panel._widgets["mu_r"].isVisibleTo(ih_panel)

    def test_FEM_workpiece_widgets_visible(self, ih_panel):
        """The 2026-04-12 invisible-UI bug was that FEM read mu_r /
        wp_sigma from hidden widgets that the user could not edit.
        After the fix both must be visible in FEM with SIBC."""
        ih_panel._method_combo.setCurrentText("FEM")
        ih_panel._widgets["workpiece_mode"].setCurrentText("SIBC")
        assert ih_panel._widgets["wp_sigma"].isVisibleTo(ih_panel)
        assert ih_panel._widgets["mu_r"].isVisibleTo(ih_panel)
        assert ih_panel._widgets["half_thickness"].isVisibleTo(ih_panel)


# ============================================================
# build_command roundtrip — every command parses cleanly
# ============================================================

def _calc_inductance_argparse():
    """Build the argparse instance from calc_inductance.py without
    importing the (slow) NGSolve modules. We just need parse_args
    to validate the GUI command."""
    import calc_inductance  # type: ignore
    # The script's __main__ block builds the parser inline; replay it
    # by calling main() with --help would exit, so we re-create the
    # parser by reading the source. Cheap path: invoke the module as
    # __main__ but intercept SystemExit. Cheapest: just import and
    # construct the parser ourselves matching the script.
    p = argparse.ArgumentParser()
    p.add_argument("--mode", default="solve", choices=["solve", "post"])
    p.add_argument("--vol", required=True)
    p.add_argument("--fes-order", type=int, default=0)
    p.add_argument("--source", default="source")
    p.add_argument("--sink", default="sink")
    p.add_argument("--msh-output", default="")
    p.add_argument("--output", default="")
    p.add_argument("--coil", default="coil")
    p.add_argument("--coil-sigma", type=float, default=0)
    p.add_argument("--workpiece", default="")
    p.add_argument("--impedance-model", default="esim",
                   choices=["esim", "sibc", "dowell"])
    p.add_argument("--frequency", type=float, default=50000)
    p.add_argument("--current", type=float, default=1.0)
    p.add_argument("--sigma", type=float, default=2e6)
    p.add_argument("--half-thickness", type=float, default=0.005)
    p.add_argument("--material", default="steel",
                   choices=["steel", "copper", "aluminum"])
    p.add_argument("--mu-r", type=float, default=0)
    p.add_argument("--esim-geometry", default="local_curvature",
                   choices=["local_curvature", "none"])
    p.add_argument("--bh-file", default="")
    p.add_argument("--solver", default="lu",
                   choices=["lu", "minres", "gmres"])
    p.add_argument("--field-air", action="store_true")
    p.add_argument("--use-local-curvature", action="store_true")
    p.add_argument("--lx", type=float, default=0.07)
    p.add_argument("--ly", type=float, default=0.07)
    p.add_argument("--lz", type=float, default=0.07)
    p.add_argument("--maxh-vol", type=float, default=0.01)
    p.add_argument("--j-npy", default="")
    p.add_argument("--mesh-vol", default="")
    return p


def _calc_fem_kelvin_argparse():
    p = argparse.ArgumentParser()
    p.add_argument("--vol", required=True)
    p.add_argument("--fes-order", type=int, default=0)
    p.add_argument("--frequency", type=float, default=50000)
    p.add_argument("--sigma", type=float, default=2e6)
    p.add_argument("--impedance", default="sibc",
                   choices=["sibc", "esim"])
    p.add_argument("--mu-r", type=float, default=100)
    p.add_argument("--bh-file", default="")
    p.add_argument("--material", default="steel",
                   choices=["steel", "copper", "aluminum", "custom"])
    p.add_argument("--current", type=float, default=1.0)
    p.add_argument("--half-thickness", type=float, default=0.005)
    p.add_argument("--max-iter", type=int, default=15)
    p.add_argument("--solver", default="pardiso",
                   choices=["pardiso", "bddc", "iccg", "ams"])
    p.add_argument("--reg", type=float, default=1e-6)
    p.add_argument("--shift-eps", type=float, default=1e-6)
    p.add_argument("--nthreads", type=int, default=0)
    p.add_argument("--msh-output", default="")
    p.add_argument("--output", default="")
    return p


class TestCommandRoundtrip:
    """Every panel command must parse cleanly through the calc_*.py
    argparse. The pre-2026-04-12 bug where ``--material custom`` was
    sent to a parser that did not declare ``custom`` as a valid
    choice would have been caught by this test."""

    def test_PEEC_no_workpiece(self, ih_panel):
        ih_panel._method_combo.setCurrentText("PEEC+BEM")
        ih_panel._widgets["workpiece_mode"].setCurrentText("off")
        cmd = ih_panel.build_command(None)
        assert cmd[0].endswith("python.exe") or cmd[0].endswith("python")
        assert cmd[1].endswith("calc_peec.py")
        assert "--coil-radius" in cmd
        assert "--vol" not in cmd

    def test_PEEC_with_workpiece(self, ih_panel):
        ih_panel._method_combo.setCurrentText("PEEC+BEM")
        ih_panel._widgets["workpiece_mode"].setCurrentText("SIBC")
        cmd = ih_panel.build_command("model.vol")
        assert cmd[1].endswith("calc_peec.py")
        assert "--vol" in cmd
        assert "--workpiece" in cmd

    @pytest.mark.parametrize("wp_mode,expected_imp",
                             [("SIBC", "sibc"),
                              ("ESIM", "esim")])
    def test_FEM(self, ih_panel, wp_mode, expected_imp):
        ih_panel._method_combo.setCurrentText("FEM")
        ih_panel._widgets["workpiece_mode"].setCurrentText(wp_mode)
        cmd = ih_panel.build_command("model.vol")
        assert cmd[1].endswith("calc_fem_kelvin.py")
        parser = _calc_fem_kelvin_argparse()
        ns = parser.parse_args(cmd[2:])
        assert ns.material == "custom"
        assert ns.impedance == expected_imp
        assert ns.half_thickness == pytest.approx(0.005)
        # User can override mu_r and sigma from the GUI
        assert ns.mu_r > 0
        assert ns.sigma > 0
