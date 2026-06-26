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

The string-grep tests in validation_test/panels/test_panel_ui_logic.py only
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

    def test_nine_methods_present(self, ih_panel):
        """9-method dropdown after thermal split into 3 modes (v4.63.0,
        commit de1d6271).  6 EM-side + 3 thermal-side.  Updated 2026-05-24
        from the original 6-method assertion."""
        from radia_ih import (METHOD_PEEC_IND, METHOD_BEMA_IND,
                              METHOD_PEEC_BEM, METHOD_BEMA_BEM,
                              METHOD_PEEC_FEM_KELVIN, METHOD_FEM_FULL,
                              METHOD_THERMAL_3D_STATIC,
                              METHOD_THERMAL_3D_ROTATING,
                              METHOD_THERMAL_AXISYM)
        items = [ih_panel._method_combo.itemText(i)
                 for i in range(ih_panel._method_combo.count())]
        assert items == [METHOD_PEEC_IND, METHOD_BEMA_IND,
                         METHOD_PEEC_BEM, METHOD_BEMA_BEM,
                         METHOD_PEEC_FEM_KELVIN, METHOD_FEM_FULL,
                         METHOD_THERMAL_3D_STATIC,
                         METHOD_THERMAL_3D_ROTATING,
                         METHOD_THERMAL_AXISYM]

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
        # 4.17.0+: PEEC+BEM uses perimeter-only filaments (n_peri),
        # the volume grid (nwinc/nhinc) was retired -- it is now
        # PEEC+FEM+Kelvin only.
        assert self._visible(ih_panel, "peec_n_peri")
        assert not self._visible(ih_panel, "peec_nwinc")
        assert not self._visible(ih_panel, "peec_nhinc")
        assert self._visible(ih_panel, "wp_sigma")
        assert self._visible(ih_panel, "mu_r")
        assert self._visible(ih_panel, "half_thickness")
        assert self._visible(ih_panel, "impedance_model")

    def test_FEM_full_no_step(self, ih_panel):
        from radia_ih import METHOD_FEM_FULL
        ih_panel._method_combo.setCurrentText(METHOD_FEM_FULL)
        # FEM full uses volumetric coil mesh from .vol -- no STEP
        assert not self._visible(ih_panel, "peec_step")
        # Workpiece widgets visible (coil mesh + WP mesh from .vol)
        assert self._visible(ih_panel, "wp_sigma")
        assert self._visible(ih_panel, "mu_r")

    def test_PEEC_FEM_KELVIN_step_AND_workpiece(self, ih_panel):
        from radia_ih import METHOD_PEEC_FEM_KELVIN
        ih_panel._method_combo.setCurrentText(METHOD_PEEC_FEM_KELVIN)
        # PEEC coil filaments (STEP) + FEM workpiece (vol)
        assert self._visible(ih_panel, "peec_step")
        assert self._visible(ih_panel, "peec_nwinc")
        assert self._visible(ih_panel, "wp_sigma")
        assert self._visible(ih_panel, "mu_r")


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
        """FEM side has 4 solvers (pardiso / AMS / BDDC / iccg).  The
        'shifted AMS' entry was removed as a UI-duplicate map - it
        targeted the same 'ams' backend that AMS itself uses (caught
        by panel-review skill 2026-05-12, see _FEM_SOLVER_MAP).
        Updated 2026-05-24 from the original 5-solver assertion."""
        from radia_ih import METHOD_FEM_FULL
        ih_panel._method_combo.setCurrentText(METHOD_FEM_FULL)
        items = [ih_panel._widgets["solver"].itemText(i)
                 for i in range(ih_panel._widgets["solver"].count())]
        assert items == [
            "pardiso (direct)",
            "AMS (iterative, p=1)",
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

    @pytest.fixture
    def fake_coil_vol(self, tmp_path):
        """Materialise a fake coil .vol for BEM-A paths.

        Only an existence check is exercised by build_command(); the file
        contents are not parsed in these UI tests (NGSolve Mesh() loading
        is exercised in validation_test/cubit/test_inductance_p_convergence.py).
        """
        f = tmp_path / "fake_coil.vol"
        f.write_text("mesh3d 1\n# stub for build_command existence check\n")
        return str(f)

    def test_PEEC_inductance_command(self, ih_panel, fake_step):
        """Vacuum inductance via unified calc_inductance.py (--coil-solver peec)."""
        from radia_ih import METHOD_PEEC_IND
        ih_panel._method_combo.setCurrentText(METHOD_PEEC_IND)
        ih_panel._widgets["peec_step"].setText(fake_step)
        cmd = ih_panel.build_command(None)
        assert cmd[1].endswith("calc_inductance.py"), cmd[1]
        assert "--coil-step" in cmd
        assert fake_step in cmd
        assert "--coil-solver" in cmd
        # PEEC mode: --coil-solver peec; --peec-n-peri present
        i = cmd.index("--coil-solver")
        assert cmd[i + 1] == "peec", f"expected 'peec', got {cmd[i + 1]!r}"
        assert "--peec-n-peri" in cmd
        assert "--coil-maxh" not in cmd
        assert "--frequency" in cmd
        assert "--current" in cmd
        assert "--coil-sigma" in cmd

    def test_BEMA_inductance_command(self, ih_panel, fake_coil_vol):
        """Vacuum inductance via calc_inductance.py (--coil-solver bem-a).

        BEM-A consumes a pre-meshed surface .vol (--coil-vol), NOT a
        CAD .step.  The coil_vol panel row is shown when method=BEM-A
        and replaces peec_step.
        """
        from radia_ih import METHOD_BEMA_IND
        ih_panel._method_combo.setCurrentText(METHOD_BEMA_IND)
        ih_panel._widgets["coil_vol"].setText(fake_coil_vol)
        cmd = ih_panel.build_command(None)
        assert cmd[1].endswith("calc_inductance.py"), cmd[1]
        # BEM-A: --coil-vol present, --coil-step absent.
        assert "--coil-vol" in cmd
        assert "--coil-step" not in cmd
        assert fake_coil_vol in cmd
        assert "--coil-solver" in cmd
        i = cmd.index("--coil-solver")
        assert cmd[i + 1] == "bem-a", f"expected 'bem-a', got {cmd[i + 1]!r}"
        # --coil-maxh retired from the panel CLI (BEM-A reads pre-meshed
        # .vol, no on-the-fly OCC re-mesh).  --peec-n-peri PEEC-only.
        assert "--coil-maxh" not in cmd
        assert "--peec-n-peri" not in cmd

    def test_PEEC_BEM_command(self, ih_panel, fake_step):
        """Weak-coupled PEEC coil + scalar BEM-SIBC."""
        from radia_ih import METHOD_PEEC_BEM
        ih_panel._method_combo.setCurrentText(METHOD_PEEC_BEM)
        ih_panel._widgets["peec_step"].setText(fake_step)
        cmd = ih_panel.build_command("model.vol")
        assert cmd[1].endswith("calc_inductance.py"), cmd[1]
        assert "--coil-step" in cmd
        assert "--coil-solver" in cmd
        i = cmd.index("--coil-solver")
        assert cmd[i + 1] == "peec"
        assert "--vol" in cmd
        assert "model.vol" in cmd
        assert "--peec-n-peri" in cmd
        assert "--coil-maxh" not in cmd
        # Volume-grid nwinc/nhinc retired (4.17.0+)
        assert "--peec-nwinc" not in cmd
        assert "--peec-nhinc" not in cmd
        # Workpiece settings present
        assert "--sigma" in cmd
        assert "--mu-r" in cmd
        assert "--half-thickness" in cmd
        assert "--impedance-model" in cmd

    def test_BEMA_BEM_command(self, ih_panel, fake_coil_vol):
        """Weak-coupled BEM-A coil + scalar BEM-SIBC.

        BEM-A coil consumes --coil-vol (pre-meshed); workpiece consumes
        the AnalysisWindow's main .vol.  Two .vol files (coil + wp) are
        passed independently.
        """
        from radia_ih import METHOD_BEMA_BEM
        ih_panel._method_combo.setCurrentText(METHOD_BEMA_BEM)
        ih_panel._widgets["coil_vol"].setText(fake_coil_vol)
        cmd = ih_panel.build_command("model.vol")
        assert cmd[1].endswith("calc_inductance.py"), cmd[1]
        assert "--coil-vol" in cmd
        assert fake_coil_vol in cmd
        assert "--coil-step" not in cmd
        assert "--coil-solver" in cmd
        i = cmd.index("--coil-solver")
        assert cmd[i + 1] == "bem-a"
        assert "--vol" in cmd
        assert "model.vol" in cmd
        assert "--coil-maxh" not in cmd
        assert "--peec-n-peri" not in cmd
        # Workpiece args still present
        assert "--sigma" in cmd
        assert "--mu-r" in cmd
        assert "--half-thickness" in cmd
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
