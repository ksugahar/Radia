"""The thermal-chain command lines expose the strict options and not the
removed fallback ones."""
from __future__ import annotations

import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PANELS = os.path.join(ROOT, "src", "radia", "panels")

CASES = {
    "calc_heat.py": (["--max-halvings", "--allow-frozen-ht",
                      "--em-reference-temperature", "--power-tolerance",
                      "--em-heat-boundaries", "--rotor-states",
                      "--angle-step-tolerance"],
                     ["--q-phi-average-n", "--surface-label"]),
    "calc_heat_axisym.py": (["--max-halvings", "--allow-frozen-ht",
                             "--n-phi-samples", "--power-tolerance"],
                            ["--q-phi-average-n"]),
    "calc_ih_axisym_coupled.py": (["--max-halvings", "--skin-resolution",
                                   "--newton-max-iter"], []),
    "calc_heat_with_em_table.py": (["--allow-frozen-ht",
                                    "--allow-em-table-extrapolation"],
                                   ["--ht-order",
                                    "--allow-table-extrapolation"]),
}


def _help(args):
    env = dict(os.environ, PYTHONPATH=os.path.join(ROOT, "src"))
    proc = subprocess.run([sys.executable, *args, "--help"],
                          capture_output=True, text=True, timeout=120,
                          env=env)
    assert proc.returncode == 0, proc.stderr
    return proc.stdout


@pytest.mark.parametrize("script", sorted(CASES))
def test_panel_flags(script):
    text = _help([os.path.join(PANELS, script)])
    present, removed = CASES[script]
    for flag in present:
        assert flag in text, f"{script} lacks {flag}"
    for flag in removed:
        assert flag + " " not in text and flag + "\n" not in text, \
            f"{script} still offers {flag}"


def test_post_processing_takes_region_and_geometry_from_the_sidecar():
    text = _help(["-m", "radia.ih_thermal_post", "depth"])
    assert "--span" in text
    for flag in ("--axisymmetric", "--stations-unit", "--region"):
        assert flag not in text


def test_sidecar_command_offers_every_record_field():
    text = _help(["-m", "radia.ih_thermal", "sidecar"])
    for flag in ("--boundaries", "--p-wp", "--frequency", "--definedon",
                 "--geometry"):
        assert flag in text
