"""Static smoke tests for the Electromagnet block's HDiv-VIM method.

These checks verify the DesignSpec-to-CLI contract without running the
numerical solver.  The command generated for the HDiv-VIM method must call
``calc_accel_hdiv.py`` and must be accepted by argparse.
"""

from __future__ import annotations

import json
import subprocess
import sys
import types
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
RADIA_PACKAGE = REPO / "src" / "radia"

# This is a DesignSpec/CLI contract test; loading the native extension would
# turn it into a build test. Provide only the package namespace needed for
# relative imports inside em_design.py.
if "radia" not in sys.modules:
    radia_package = types.ModuleType("radia")
    radia_package.__path__ = [str(RADIA_PACKAGE)]
    sys.modules["radia"] = radia_package


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


def test_hdiv_command_does_not_offer_fem_iteration_options():
    """HDiv-VIM owns its solver and convergence contract (radia.vim.Solve).

    The DesignSpec must not expose or send options the route cannot apply,
    and the CLI must reject them instead of ignoring them.
    """
    spec = _hdiv_spec()
    assert not {"max_iter", "tol", "relax", "hdiv_solver"} & spec.visible_fields()
    cmd = spec.build_command(python=sys.executable)
    for flag in ("--solver", "--max-iter", "--tol", "--relax"):
        assert flag not in cmd
    rc, err = _argparse_dry_run(cmd)
    assert rc == 0, err
    for extra in (["--solver", "2"], ["--max-iter", "30"], ["--tol", "1e-3"],
                  ["--relax", "0.3"]):
        proc = subprocess.run(cmd + extra, capture_output=True, text=True, timeout=60)
        assert proc.returncode == 2, (extra, proc.stderr)
        assert "unrecognized arguments" in proc.stderr


@pytest.mark.parametrize("name,value", [
    ("tol", "1e-9"), ("max_iter", 200), ("relax", "0.0"), ("newton", True),
    ("solver", "iccg"), ("n_steps", 4),
])
def test_hdiv_spec_rejects_an_explicit_fem_option(name, value):
    """A changed FEM option is an error for HDiv-VIM, never silently dropped."""
    spec = _hdiv_spec()
    setattr(spec, name, value)
    with pytest.raises(ValueError, match=name):
        spec.build_command(python=sys.executable)


def test_hdiv_fes_order_reaches_cli():
    spec = _hdiv_spec()
    spec.fes_order = 2
    cmd = spec.build_command(python=sys.executable)
    assert cmd[cmd.index("--hdiv-order") + 1] == "2"
    rc, err = _argparse_dry_run(cmd)
    assert rc == 0, err


def test_hdiv_headless_entries_are_owned_by_application_manifest():
    manifest_path = (
        REPO / "src" / "radia" / "panels"
        / "application_interface_manifest.json"
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    applications = {item["id"]: item for item in manifest["applications"]}

    assert manifest["schema"] == "radia.application_interface_manifest.v2"
    assert "src/radia/panels/calc_accel_hdiv.py" in applications["radia-em"]["headless"]
    assert (
        "src/radia/panels/calc_motor_hdiv_reduced.py"
        in applications["radia-motor"]["headless"]
    )


def test_hdiv_ima_passthrough():
    spec = _hdiv_spec()
    spec.ima = "+x-z"
    cmd = spec.build_command(python=sys.executable)
    assert "--ima" in cmd
    assert cmd[cmd.index("--ima") + 1] == "+x-z"
    rc, err = _argparse_dry_run(cmd)
    assert rc == 0, err
