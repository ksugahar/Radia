"""Static smoke tests for the EM notebook panel's HDiv-VIM method.

These checks verify the DesignSpec-to-CLI contract without running the
numerical solver.  The command generated for the HDiv-VIM method must call
``calc_accel_hdiv.py`` and must be accepted by argparse.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src" / "radia"))


def _hdiv_spec():
    from radia.em_design import EMDesignSpec, METHOD_HDIV

    return EMDesignSpec(
        method=METHOD_HDIV,
        coil_script=str(REPO / "src" / "radia" / "panels" / "samples" / "em_sample_coil.py"),
        vol="model.vol",
        material="linear",
        mu_r="1000",
    )


def _argparse_dry_run(cmd):
    proc = subprocess.run(cmd + ["--help"], capture_output=True, text=True, timeout=30)
    return proc.returncode, proc.stderr


def test_hdiv_command_uses_calc_accel_hdiv():
    cmd = _hdiv_spec().build_command(python=sys.executable)
    assert cmd[1].endswith("calc_accel_hdiv.py"), cmd[1]


@pytest.mark.parametrize("material,mu_r,bh_file,hys_file", [
    ("custom", "1000", "", ""),
    ("linear", "1000", "", ""),
    ("steel", "", str(REPO / "src" / "radia" / "panels" / "samples" / "em_sample_bh.txt"), ""),
    ("hysteresis", "", "", "dummy.hys"),
])
def test_hdiv_material_modes_parse(material, mu_r, bh_file, hys_file):
    spec = _hdiv_spec()
    spec.material = material
    spec.mu_r = mu_r
    spec.bh_file = bh_file
    spec.hys_file = hys_file
    cmd = spec.build_command(python=sys.executable)
    assert "--material" in cmd
    assert cmd[cmd.index("--material") + 1] == material
    rc, err = _argparse_dry_run(cmd)
    assert rc == 0, f"argparse rejected the HDiv-VIM panel command:\n{err}"


def test_hdiv_solver_choice_argparse():
    spec = _hdiv_spec()
    for solver_id in (0, 1, 2):
        spec.hdiv_solver = solver_id
        cmd = spec.build_command(python=sys.executable)
        assert "--solver" in cmd
        assert cmd[cmd.index("--solver") + 1] == str(solver_id)
        rc, err = _argparse_dry_run(cmd)
        assert rc == 0, f"argparse rejected solver {solver_id}:\n{err}"


def test_hdiv_fes_order_reaches_cli():
    spec = _hdiv_spec()
    spec.fes_order = 2
    cmd = spec.build_command(python=sys.executable)
    assert cmd[cmd.index("--hdiv-order") + 1] == "2"
    rc, err = _argparse_dry_run(cmd)
    assert rc == 0, err


def test_hdiv_order_registry_metadata_is_domain_specific():
    registry_path = REPO / "src" / "radia" / "panels" / "panel_registry.json"
    registry = json.loads(registry_path.read_text(encoding="utf-8"))["panels"]
    for panel_id, cli in (
        ("accel_hdiv", "--hdiv-order"),
        ("motor_hdiv_reduced", "--order"),
    ):
        param = next(p for p in registry[panel_id]["params"] if p["cli"] == cli)
        assert param["ja"] == "HDiv 次数"
        assert param["physics"] == "1=RT1, 2=RT2"


def test_hdiv_ima_passthrough():
    spec = _hdiv_spec()
    spec.ima = "+x-z"
    cmd = spec.build_command(python=sys.executable)
    assert "--ima" in cmd
    assert cmd[cmd.index("--ima") + 1] == "+x-z"
    rc, err = _argparse_dry_run(cmd)
    assert rc == 0, err
