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

def test_rotor_states_alone_is_an_accepted_heat_source(tmp_path):
    """--rotor-states replaces --qsurf-sol, so it must pass the source gate
    and reach the next input check."""
    env = dict(os.environ, PYTHONPATH=os.path.join(ROOT, "src"))
    missing = str(tmp_path / "missing_wp.vol")
    base = [sys.executable, os.path.join(PANELS, "calc_heat.py"),
            "--wp-vol", missing]
    with_rotor = subprocess.run(
        base + ["--rotor-states", str(tmp_path / "rotor.json")],
        capture_output=True, text=True, timeout=120, env=env)
    without = subprocess.run(base, capture_output=True, text=True,
                             timeout=120, env=env)
    assert "--rotor-states is required" in without.stdout
    assert "is required" not in with_rotor.stdout
    assert "--wp-vol not found" in with_rotor.stdout


@pytest.mark.parametrize("script,function,kwargs", [
    ("calc_fem_kelvin.py", "solve_fem", {"solver": "pardiso"}),
    ("calc_fem_kelvin.py", "solve_fem", {"solver": "misspelled"}),
    ("calc_fem_coilmesh.py", "solve_fem_coilmesh", dict(vol="missing.vol", frequency=1.,
     I_target=1., coil_sigma=1., wp_sigma=1., wp_mu_r=1., half_thickness=1., solver="bddc")),
])
def test_unsupported_solver_fails_before_loading_a_mesh(script, function, kwargs):
    import ast
    from pathlib import Path
    path = Path(PANELS) / script
    node = next(n for n in ast.parse(path.read_text(encoding="utf-8")).body
                if isinstance(n, ast.FunctionDef) and n.name == function)
    ns = {}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"), ns)
    with pytest.raises(ValueError):
        ns[function](**kwargs)
