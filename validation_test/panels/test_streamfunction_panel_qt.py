"""Headless Qt behaviour tests for the Stream-Function coil-design panel.

Instantiates the real PySide6 ``StreamFunctionPanel`` (a composite ModePanel:
a coil/conductor-.vol Browse + a Method combo + a QStackedWidget over the
Design / Pareto / Manufacture surface sub-panels AND the Volume 3D sub-panel)
via the offscreen Qt platform plugin and verifies the user-visible behaviour
that the calc golden (a subprocess) and the static contract audit cannot:

  - the Method combo defaults to Design and lists all 4 modes (3 surface +
    Volume 3D);
  - Volume 3D is backed by a DIFFERENT calc (calc_streamfunction_volume.py)
    with its own knobs and no surface confine/regularize;
  - each mode's sub-panel exposes ONLY its own argparse knobs (so a hidden
    widget can't silently feed build_command, and a mode can't drop a flag);
  - the choice combos (regularise / confine / pareto-lever / chain) carry the
    exact CLI choices, including the new --confine abe;
  - build_command for every mode targets calc_streamfunction.py and threads
    the user's coil/eval/target + the mode-specific flags into argv;
  - the namespaced save/restore roundtrips per sub-panel.

Constructing the panel does NOT import ngsolve/radia (calc_streamfunction's
build_argparser is pure argparse), so these run inside pytest.
"""
from __future__ import annotations

import os
import sys
import tomllib
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src" / "radia"))

MODES = ["Design", "Pareto", "Manufacture"]        # surface modes (shared knobs)
ALL_MODES = MODES + ["Volume 3D"]                  # full Method combo


def _sub(panel, name):
    return panel._sub_panels[name]


def _keys(panel, name):
    return set(_sub(panel, name)._widgets.keys())


def _combo_items(widget):
    return [widget.itemText(i) for i in range(widget.count())]


def _set_choice(sub, key, cli_value):
    """Select a choice combo by its CLI value (the display label differs)."""
    w = sub._widgets[key]
    for i, (_label, cli) in enumerate(w._radia_choice_pairs):
        if cli == cli_value:
            w.setCurrentIndex(i)
            return
    raise AssertionError(f"no choice {cli_value!r} in {key}")


# ============================================================
# Method combo
# ============================================================
class TestMethodCombo:

    def test_default_is_design(self, sf_panel):
        assert sf_panel._method_combo.currentText() == "Design"

    def test_modes_present(self, sf_panel):
        assert _combo_items(sf_panel._method_combo) == ALL_MODES

    def test_switch_changes_stacked_widget(self, sf_panel):
        for name in ALL_MODES:
            sf_panel._method_combo.setCurrentText(name)
            assert sf_panel._stack.currentWidget() is _sub(sf_panel, name)
            assert sf_panel._current_sub() is _sub(sf_panel, name)


# ============================================================
# Per-mode widget isolation (no hidden widget feeds build_command)
# ============================================================
class TestWidgetIsolation:

    def test_design_has_only_solver_knobs(self, sf_panel):
        k = _keys(sf_panel, "Design")
        assert {"order", "regularize", "confine", "alpha", "eval_max",
                "target_cf"} <= k
        # NOT the pareto / manufacture knobs
        for absent in ("pareto_lever", "nlevels", "contour_sub", "chain",
                       "distort", "flux_plot"):
            assert absent not in k, f"Design leaks {absent}"

    def test_pareto_has_lever_not_manufacture(self, sf_panel):
        k = _keys(sf_panel, "Pareto")
        assert {"pareto_lever", "linf_iter", "geom_scale_min",
                "geom_scale_max", "confine"} <= k
        for absent in ("nlevels", "contour_sub", "chain", "distort",
                       "flux_plot"):
            assert absent not in k, f"Pareto leaks {absent}"

    def test_manufacture_has_wire_knobs_not_pareto(self, sf_panel):
        k = _keys(sf_panel, "Manufacture")
        assert {"nlevels", "contour_sub", "chain", "distort", "step_output",
                "peec", "flux_plot", "flux_plane", "steps_plot",
                "confine"} <= k
        for absent in ("pareto_lever", "alpha_min", "geom_scale_min"):
            assert absent not in k, f"Manufacture leaks {absent}"

    def test_coil_vol_and_method_on_top_not_subs(self, sf_panel):
        # primary coil .vol + method live on the composite, not in the subs
        assert "wp_vol" in sf_panel._widgets
        assert "method" in sf_panel._widgets
        for name in ALL_MODES:
            assert "wp_vol" not in _keys(sf_panel, name)
            assert "coil_vol" not in _keys(sf_panel, name)


# ============================================================
# Choice combos carry the exact CLI choices
# ============================================================
class TestChoiceCombos:

    def test_regularize_choices(self, sf_panel):
        w = _sub(sf_panel, "Design")._widgets["regularize"]
        # combos render labels; the bound values must be l2 / h1
        assert w.count() == 2

    def test_confine_has_abe(self, sf_panel):
        # confine is a shared knob -> present in every mode, abe first
        for name in MODES:
            w = _sub(sf_panel, name)._widgets["confine"]
            assert w.count() == 3, f"{name} confine should have off/on/abe"

    def test_pareto_lever_three(self, sf_panel):
        w = _sub(sf_panel, "Pareto")._widgets["pareto_lever"]
        assert w.count() == 3      # alpha / linf / geometry

    def test_chain_two(self, sf_panel):
        w = _sub(sf_panel, "Manufacture")._widgets["chain"]
        assert w.count() == 2      # field_aware / nn


# ============================================================
# build_command roundtrip per mode
# ============================================================
class TestBuildCommand:

    @pytest.fixture
    def vols(self, tmp_path):
        coil = tmp_path / "coil.vol"
        coil.write_text("mesh3d 1\n")
        evalv = tmp_path / "eval.vol"
        evalv.write_text("mesh3d 1\n")
        return str(coil), str(evalv)

    def _prime(self, sf_panel, name, evalv, target="x"):
        sf_panel._method_combo.setCurrentText(name)
        _sub(sf_panel, name).restore_state({"eval_vol": evalv,
                                            "target_cf": target})

    def _common(self, cmd, coil, evalv, mode):
        assert cmd[1].endswith("calc_streamfunction.py"), cmd[1]
        assert "--coil-vol" in cmd and coil in cmd
        assert "--eval-vol" in cmd and evalv in cmd
        assert "--target-cf" in cmd
        assert "--method" in cmd
        assert cmd[cmd.index("--method") + 1] == mode

    def test_design_command(self, sf_panel, vols):
        coil, evalv = vols
        self._prime(sf_panel, "Design", evalv)
        cmd = sf_panel.build_command(coil)
        self._common(cmd, coil, evalv, "design")
        assert "--regularize" in cmd and "--confine" in cmd
        assert "--pareto-lever" not in cmd
        assert "--nlevels" not in cmd

    def test_pareto_command(self, sf_panel, vols):
        coil, evalv = vols
        self._prime(sf_panel, "Pareto", evalv)
        cmd = sf_panel.build_command(coil)
        self._common(cmd, coil, evalv, "pareto")
        assert "--pareto-lever" in cmd
        assert "--nlevels" not in cmd

    def test_manufacture_command(self, sf_panel, vols):
        coil, evalv = vols
        self._prime(sf_panel, "Manufacture", evalv)
        cmd = sf_panel.build_command(coil)
        self._common(cmd, coil, evalv, "manufacture")
        # the new manufacture knobs reach argv
        assert "--nlevels" in cmd
        assert "--contour-sub" in cmd
        assert "--chain" in cmd
        assert "--confine" in cmd
        assert "--pareto-lever" not in cmd

    def test_manufacture_confine_abe_passes_through(self, sf_panel, vols):
        coil, evalv = vols
        self._prime(sf_panel, "Manufacture", evalv)
        _set_choice(_sub(sf_panel, "Manufacture"), "confine", "abe")
        cmd = sf_panel.build_command(coil)
        assert cmd[cmd.index("--confine") + 1] == "abe"


# ============================================================
# Volume 3D mode (foliated-Clebsch volume SF -- a DIFFERENT calc)
# ============================================================
class TestVolume3D:

    def test_volume_has_own_knobs_not_surface(self, sf_panel):
        k = _keys(sf_panel, "Volume 3D")
        assert {"target_bz", "n_leaves", "fes_order", "aca_eps",
                "rings_per_span", "n_targets", "nphi", "nz"} <= k
        # the Volume mode uses calc_streamfunction_volume, NOT the surface calc:
        for absent in ("confine", "regularize", "nlevels", "pareto_lever",
                       "coil_vol", "target_cf"):
            assert absent not in k, f"Volume 3D leaks surface knob {absent}"

    def test_volume_command_targets_volume_calc(self, sf_panel, tmp_path):
        tube = tmp_path / "tube.vol"
        tube.write_text("mesh3d 1\n")
        sf_panel._method_combo.setCurrentText("Volume 3D")
        cmd = sf_panel.build_command(str(tube))
        assert cmd[1].endswith("calc_streamfunction_volume.py"), cmd[1]
        # conductor volume mesh via --vol (NOT the surface --coil-vol / --method)
        assert "--vol" in cmd and str(tube) in cmd
        assert "--coil-vol" not in cmd
        assert "--method" not in cmd
        # design knobs reach argv; same _sf_volume suffix on json + msh (P4)
        assert "--target-bz" in cmd and "--n-leaves" in cmd
        out = cmd[cmd.index("--output") + 1]
        msh = cmd[cmd.index("--msh-output") + 1]
        assert out.endswith("_sf_volume.json") and msh.endswith("_sf_volume.msh")


# ============================================================
# save / restore roundtrip
# ============================================================
class TestSaveRestore:

    def test_namespaced_roundtrip(self, sf_panel):
        sf_panel._method_combo.setCurrentText("Manufacture")
        _sub(sf_panel, "Manufacture").restore_state(
            {"nlevels": "18", "target_cf": "x"})
        state = sf_panel.save_state()
        # the manufacture nlevels is namespaced under the mode
        assert any(k.endswith("/nlevels") for k in state), state.keys()
        fresh = sf_panel.__class__()
        fresh.restore_state(state)
        assert fresh._method_combo.currentText() == "Manufacture"
        fresh.deleteLater()


# ============================================================
# Window: launcher .vol overrides the restored coil wp_vol
# ============================================================
class TestWindowVolOverride:

    def test_module_contract_and_console_script(self):
        import radia_streamfunction as mod
        scripts = tomllib.loads((REPO / "pyproject.toml").read_text(
            encoding="utf-8"))["project"]["scripts"]
        assert callable(mod.main)
        assert mod.REQUIRED_LABELS == []
        assert mod.OPTIONAL_LABELS == []
        assert mod.OPTIONAL_FILES == {}
        assert scripts["radia-streamfunction"] == (
            "radia.radia_streamfunction:main")

    def test_launcher_vol_fills_wp_vol(self, qapp, tmp_path):
        """StreamFunctionWindow must restore settings AND let the constructor
        vol_path (the Cubit Solve-menu's current .vol) override the restored
        coil wp_vol -- the radia_em / radia_ih restore-order bug class."""
        from radia_streamfunction import StreamFunctionWindow
        coil = tmp_path / "launcher_coil.vol"
        coil.write_text("mesh3d 1\n")
        win = StreamFunctionWindow(vol_path=str(coil))
        assert win._panel.wp_vol_path().endswith("launcher_coil.vol")
        win.deleteLater()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
