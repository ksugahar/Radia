"""Headless Qt tests for the IH (Induction Heating) panel.

These tests instantiate the real PySide6 IHPanel via the offscreen
Qt platform plugin and verify the panel's user-visible behaviour:

  - Method combo defaults to PEEC inductance and exposes the 4
    canonical methods (METHOD_PEEC_IND, METHOD_PEEC_BEM,
    METHOD_PEEC_FEM_KELVIN, METHOD_FEM_FULL).
  - Mode switches show / hide the right widget rows.
  - Workpiece impedance toggles (SIBC vs ESIM) drive the right
    widget visibility.
  - build_command for every method produces a list whose argv[1]
    targets the right calc_*.py and whose top-level flags include
    the values the user typed in the panel.

The string-grep tests in tests/panels/test_panel_ui_logic.py only
catch widget-removal regressions; these Qt tests catch BEHAVIOUR
regressions like the empty-Method-combo bug from 2026-04-12.

Refreshed 2026-04-26 to match the post-2026-04-19 IH panel
restructure (4 methods, separate coil_material / wp_material
sections, impedance_model dropdown, n_peri vs nwinc/nhinc split).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src" / "radia"))


# ============================================================
# Method combo behaviour
# ============================================================

class TestMethodCombo:

    def test_default_method_is_PEEC_inductance(self, ih_panel):
        """First-launch default."""
        from radia_ih import METHOD_PEEC_IND
        assert ih_panel._method_combo.currentText() == METHOD_PEEC_IND

    def test_four_methods_present(self, ih_panel):
        from radia_ih import (METHOD_PEEC_IND, METHOD_PEEC_BEM,
                              METHOD_PEEC_FEM_KELVIN, METHOD_FEM_FULL)
        items = [ih_panel._method_combo.itemText(i)
                 for i in range(ih_panel._method_combo.count())]
        assert items == [METHOD_PEEC_IND, METHOD_PEEC_BEM,
                         METHOD_PEEC_FEM_KELVIN, METHOD_FEM_FULL]

    def test_no_legacy_methods(self, ih_panel):
        """Catch creep-back of retired method strings."""
        items = {ih_panel._method_combo.itemText(i)
                 for i in range(ih_panel._method_combo.count())}
        for retired in ("BEM-SIBC (WP)",
                        "PEEC+FEM",   # 2026-04-19 split into 2 methods
                        "FEM"):       # 2026-04-19 became METHOD_FEM_FULL
            assert retired not in items, f"retired method '{retired}' creep"


# ============================================================
# Mode switch widget visibility
# ============================================================

class TestModeSwitch:

    def _visible(self, panel, key):
        w = panel._widgets.get(key)
        return w is not None and w.isVisibleTo(panel)

    def test_PEEC_inductance_minimal_widgets(self, ih_panel):
        from radia_ih import METHOD_PEEC_IND
        ih_panel._method_combo.setCurrentText(METHOD_PEEC_IND)
        # PEEC-IND needs only STEP + freq + current + coil material
        assert self._visible(ih_panel, "peec_step")
        assert self._visible(ih_panel, "peec_n_peri")
        assert self._visible(ih_panel, "freq")
        assert self._visible(ih_panel, "coil_sigma")
        # No workpiece widgets
        assert not self._visible(ih_panel, "wp_sigma")
        assert not self._visible(ih_panel, "mu_r")
        assert not self._visible(ih_panel, "half_thickness")
        assert not self._visible(ih_panel, "impedance_model")

    def test_PEEC_BEM_shows_workpiece(self, ih_panel):
        from radia_ih import METHOD_PEEC_BEM
        ih_panel._method_combo.setCurrentText(METHOD_PEEC_BEM)
        # STEP coil + workpiece widgets
        assert self._visible(ih_panel, "peec_step")
        assert self._visible(ih_panel, "peec_nwinc")
        assert self._visible(ih_panel, "peec_nhinc")
        assert self._visible(ih_panel, "wp_sigma")
        assert self._visible(ih_panel, "mu_r")
        assert self._visible(ih_panel, "half_thickness")
        assert self._visible(ih_panel, "impedance_model")
        # n_peri is PEEC-IND only
        assert not self._visible(ih_panel, "peec_n_peri")

    def test_FEM_full_no_step(self, ih_panel):
        from radia_ih import METHOD_FEM_FULL
        ih_panel._method_combo.setCurrentText(METHOD_FEM_FULL)
        # FEM full uses volumetric coil mesh from .vol -- no STEP
        assert not self._visible(ih_panel, "peec_step")
        # Workpiece widgets visible (coil mesh + WP mesh from .vol)
        assert self._visible(ih_panel, "wp_sigma")
        assert self._visible(ih_panel, "mu_r")
        # coil_mu_r ONLY shown for FEM-full (the only path with a
        # volumetric coil that can carry magnetic material)
        assert self._visible(ih_panel, "coil_mu_r")

    def test_PEEC_FEM_KELVIN_step_AND_workpiece(self, ih_panel):
        from radia_ih import METHOD_PEEC_FEM_KELVIN
        ih_panel._method_combo.setCurrentText(METHOD_PEEC_FEM_KELVIN)
        # PEEC coil filaments (STEP) + FEM workpiece (vol)
        assert self._visible(ih_panel, "peec_step")
        assert self._visible(ih_panel, "peec_nwinc")
        assert self._visible(ih_panel, "wp_sigma")
        assert self._visible(ih_panel, "mu_r")
        # Non-magnetic-coil convention for PEEC-side: hide coil_mu_r
        assert not self._visible(ih_panel, "coil_mu_r")


# ============================================================
# Solver items per method
# ============================================================

class TestSolverItems:

    def test_PEEC_inductance_solver_items(self, ih_panel):
        """PEEC-IND uses dense / hacapk solver list."""
        from radia_ih import METHOD_PEEC_IND
        ih_panel._method_combo.setCurrentText(METHOD_PEEC_IND)
        items = [ih_panel._widgets["solver"].itemText(i)
                 for i in range(ih_panel._widgets["solver"].count())]
        assert items == ["Dense LU (small)", "HACApK (large)"]

    def test_PEEC_BEM_solver_items(self, ih_panel):
        from radia_ih import METHOD_PEEC_BEM
        ih_panel._method_combo.setCurrentText(METHOD_PEEC_BEM)
        items = [ih_panel._widgets["solver"].itemText(i)
                 for i in range(ih_panel._widgets["solver"].count())]
        assert items == ["Dense LU (small)", "HACApK (large)"]

    def test_FEM_solver_items(self, ih_panel):
        """FEM side has 5 solvers (pardiso / AMS / shifted AMS / BDDC / iccg)."""
        from radia_ih import METHOD_FEM_FULL
        ih_panel._method_combo.setCurrentText(METHOD_FEM_FULL)
        items = [ih_panel._widgets["solver"].itemText(i)
                 for i in range(ih_panel._widgets["solver"].count())]
        assert items == [
            "pardiso (direct)",
            "AMS (iterative, p=1)",
            "shifted AMS (iterative, p=1)",
            "BDDC (iterative, p>=2)",
            "iccg (fallback)",
        ]


# ============================================================
# Workpiece impedance model (SIBC vs ESIM)
# ============================================================

class TestImpedanceModel:

    def test_SIBC_hides_ESIM_widgets(self, ih_panel):
        from radia_ih import METHOD_PEEC_BEM
        ih_panel._method_combo.setCurrentText(METHOD_PEEC_BEM)
        ih_panel._widgets["impedance_model"].setCurrentIndex(0)  # SIBC
        # SIBC: no BH file, no ESIM iter controls
        assert not ih_panel._widgets["bh_file"].isVisibleTo(ih_panel)
        assert not ih_panel._widgets["esim_max_iter"].isVisibleTo(ih_panel)
        assert not ih_panel._widgets["esim_tol"].isVisibleTo(ih_panel)
        # SIBC needs mu_r
        assert ih_panel._widgets["mu_r"].isVisibleTo(ih_panel)

    def test_ESIM_shows_BH_and_iter_controls(self, ih_panel):
        from radia_ih import METHOD_PEEC_BEM
        ih_panel._method_combo.setCurrentText(METHOD_PEEC_BEM)
        ih_panel._widgets["impedance_model"].setCurrentIndex(1)  # ESIM
        # ESIM: BH file + iter + tol all visible
        assert ih_panel._widgets["bh_file"].isVisibleTo(ih_panel)
        assert ih_panel._widgets["esim_max_iter"].isVisibleTo(ih_panel)
        assert ih_panel._widgets["esim_tol"].isVisibleTo(ih_panel)


# ============================================================
# build_command roundtrip per method
# ============================================================

class TestBuildCommand:
    """Each method must produce a command targeting the right
    calc_*.py and including the user's panel inputs as flags.
    """

    @pytest.fixture
    def fake_step(self, tmp_path):
        """Materialise a fake STEP file for PEEC paths."""
        f = tmp_path / "fake.step"
        f.write_text("ISO-10303-21;\nENDSEC;\nEND-ISO-10303-21;")
        return str(f)

    def test_PEEC_inductance_command(self, ih_panel, fake_step):
        from radia_ih import METHOD_PEEC_IND
        ih_panel._method_combo.setCurrentText(METHOD_PEEC_IND)
        ih_panel._widgets["peec_step"].setText(fake_step)
        cmd = ih_panel.build_command(None)  # PEEC-IND ignores .vol
        assert cmd[1].endswith("calc_peec_inductance.py"), cmd[1]
        assert "--peec-step" in cmd
        assert fake_step in cmd
        assert "--frequency" in cmd
        assert "--current" in cmd
        assert "--coil-sigma" in cmd
        assert "--peec-n-peri" in cmd

    def test_PEEC_BEM_command(self, ih_panel, fake_step):
        from radia_ih import METHOD_PEEC_BEM
        ih_panel._method_combo.setCurrentText(METHOD_PEEC_BEM)
        ih_panel._widgets["peec_step"].setText(fake_step)
        cmd = ih_panel.build_command("model.vol")
        assert cmd[1].endswith("calc_peec_bem.py"), cmd[1]
        assert "--peec-step" in cmd
        assert "--vol" in cmd
        assert "model.vol" in cmd
        assert "--peec-nwinc" in cmd
        assert "--peec-nhinc" in cmd
        # Workpiece settings present
        assert "--sigma" in cmd
        assert "--mu-r" in cmd
        assert "--half-thickness" in cmd
        # Impedance model
        assert "--impedance-model" in cmd

    def test_PEEC_FEM_KELVIN_command(self, ih_panel, fake_step):
        from radia_ih import METHOD_PEEC_FEM_KELVIN
        ih_panel._method_combo.setCurrentText(METHOD_PEEC_FEM_KELVIN)
        ih_panel._widgets["peec_step"].setText(fake_step)
        cmd = ih_panel.build_command("model.vol")
        assert cmd[1].endswith("calc_fem_kelvin.py"), cmd[1]
        assert "--vol" in cmd
        # Material passed via 'custom' so panel mu_r/sigma decouple
        # from EMMaterial preset names.
        assert "--material" in cmd
        assert cmd[cmd.index("--material") + 1] == "custom"
        assert "--sigma" in cmd
        assert "--mu-r" in cmd

    def test_FEM_FULL_command(self, ih_panel):
        from radia_ih import METHOD_FEM_FULL
        ih_panel._method_combo.setCurrentText(METHOD_FEM_FULL)
        cmd = ih_panel.build_command("model.vol")
        assert cmd[1].endswith("calc_fem_coilmesh.py"), cmd[1]
        assert "--vol" in cmd
        assert "--frequency" in cmd

    def test_PEEC_inductance_requires_STEP(self, ih_panel):
        """build_command without STEP raises ValueError."""
        from radia_ih import METHOD_PEEC_IND
        ih_panel._method_combo.setCurrentText(METHOD_PEEC_IND)
        ih_panel._widgets["peec_step"].setText("")
        with pytest.raises(ValueError, match=r"STEP"):
            ih_panel.build_command(None)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
