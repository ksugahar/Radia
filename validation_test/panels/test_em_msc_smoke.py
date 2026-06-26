"""Static smoke tests for the EM panel's MSC formulation.

These tests catch panel-vs-CLI drift WITHOUT running the heavy
Radia C++ MSC solver.  They verify:

  1. The EM panel's MSC build_command produces flags that
     `calc_accel_msc.py` argparse accepts.
  2. The 4 distinct material modes (linear / BH / hysteresis +
     the implicit `--material custom` from "mu_r (Linear)") all
     parse without error.

The panel was silently broken before 2026-04-26: the MSC mode
sent `--material custom` but `calc_accel_msc.py` registered
`add_material_args(include_custom=False)`, so every panel-built
MSC command failed at argparse.  This test class is the
regression guard.

This is a STATIC test (subprocess --help only); it does not
exercise the Radia C++ MSC solver.  The full physics regression
would require a sample .vol with a 'yoke' block and is covered
by the slow `validation_test/panels/test_em_golden.py` family.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src" / "radia"))


@pytest.fixture
def msc_panel(qapp):
    """The EM panel's MSC sub-panel.

    EMPanel was refactored (committed) from a single panel with a
    'formulation' combo into a Method combo + per-method sub-panels
    (QStackedWidget).  MSC now lives in its own ``_MSCPanel``, reached
    via ``EMPanel._sub_panels[FORM_MSC]``.  Yield that sub-panel so the
    build_command / widget assertions below target the MSC mode (its
    widgets are bound from calc_accel_msc.py's argparser).
    """
    from radia_em import EMPanel, FORM_MSC
    panel = EMPanel()
    panel._method_combo.setCurrentText(FORM_MSC)   # fires _on_method_changed
    msc = panel._sub_panels[FORM_MSC]
    msc._widgets["coil_script"].setText(
        str(REPO / "src" / "radia" / "panels" / "samples"
            / "em_sample_coil.py"))
    yield msc
    panel.deleteLater()


def _calc_msc_path():
    return str(REPO / "src" / "radia" / "panels" / "calc_accel_msc.py")


def _argparse_dry_run(cmd):
    """Run `python calc_accel_msc.py <flags> --help`. Returns (rc, stderr).

    `--help` makes argparse print and exit before any Radia / NGSolve
    import, so this is a sub-second static check.
    """
    # Replace --output / --msh-output values with placeholder paths
    # to keep them inside any valid filesystem path.
    proc = subprocess.run(cmd + ["--help"], capture_output=True, text=True,
                          timeout=30)
    return proc.returncode, proc.stderr


def test_msc_command_uses_calc_accel_msc(msc_panel):
    cmd = msc_panel.build_command("model.vol")
    assert cmd[1].endswith("calc_accel_msc.py"), cmd[1]


def test_msc_linear_material_argparse(msc_panel):
    """Material = 'custom' -> --material custom + --mu-r.

    Select by VALUE (not index): the material combo is bound directly
    from calc_accel_msc.py's --material choices, so the index order is
    an implementation detail; the *value* is the panel<->CLI contract.
    """
    msc_panel._widgets["material"].setCurrentText("custom")
    msc_panel._widgets["mu_r"].setText("1000")
    cmd = msc_panel.build_command("model.vol")
    assert "--material" in cmd
    assert cmd[cmd.index("--material") + 1] == "custom"
    assert "--mu-r" in cmd
    rc, err = _argparse_dry_run(cmd)
    assert rc == 0, f"argparse rejected the panel command:\n{err}"


def test_msc_bh_curve_material_argparse(msc_panel):
    """Material = 'steel' (library BH-curve soft iron) -> --material steel."""
    msc_panel._widgets["material"].setCurrentText("steel")
    msc_panel._widgets["bh_file"].setText(str(
        REPO / "src" / "radia" / "panels" / "samples" / "em_sample_bh.txt"))
    cmd = msc_panel.build_command("model.vol")
    assert "--material" in cmd
    assert cmd[cmd.index("--material") + 1] == "steel"
    rc, err = _argparse_dry_run(cmd)
    assert rc == 0, f"argparse rejected the panel command:\n{err}"


def test_msc_hysteresis_material_argparse(msc_panel):
    """Material = 'hysteresis' -> --material hysteresis."""
    msc_panel._widgets["material"].setCurrentText("hysteresis")
    cmd = msc_panel.build_command("model.vol")
    assert "--material" in cmd
    assert cmd[cmd.index("--material") + 1] == "hysteresis"
    rc, err = _argparse_dry_run(cmd)
    assert rc == 0, f"argparse rejected the panel command:\n{err}"


def test_msc_solver_choice_argparse(msc_panel):
    """All 3 MSC solvers (0=LU / 1=BiCGSTAB / 2=HACApK) parse cleanly."""
    msc_panel._widgets["material"].setCurrentIndex(0)
    msc_panel._widgets["mu_r"].setText("1000")
    for solver_text in ("0 (LU)", "1 (BiCGSTAB)", "2 (HACApK)"):
        msc_panel._widgets["solver"].setCurrentText(solver_text)
        cmd = msc_panel.build_command("model.vol")
        assert "--solver" in cmd
        assert cmd[cmd.index("--solver") + 1] in ("0", "1", "2")
        rc, err = _argparse_dry_run(cmd)
        assert rc == 0, f"argparse rejected solver '{solver_text}':\n{err}"


def test_msc_ima_passthrough(msc_panel):
    """IMA string from the panel survives argparse parsing."""
    msc_panel._widgets["material"].setCurrentIndex(0)
    msc_panel._widgets["mu_r"].setText("1000")
    msc_panel._widgets["ima"].setText("+x-z")
    cmd = msc_panel.build_command("model.vol")
    assert "--ima" in cmd
    assert cmd[cmd.index("--ima") + 1] == "+x-z"
    rc, err = _argparse_dry_run(cmd)
    assert rc == 0, err


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
