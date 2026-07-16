"""Contract tests: the --wp-phi-inc CLI surface of calc_inductance.

The reconstruction math is locked by tests/test_phi_inc_poisson.py
(icosphere goldens); here we pin the CLI wiring: argparse accepts the
mode, the default stays "path", and unsupported combinations fail fast
BEFORE any expensive solve.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
_PANELS = _REPO / "src" / "radia" / "panels"
for p in (_PANELS, _REPO / "validation_test" / "panels"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import calc_inductance as ci
from radia.ih_design import IHDesignSpec, METHOD_PEEC_BEM
from radia.ih_notebook import IH_NOTEBOOK_FIELD_ORDER


def _args(extra):
    p = ci.build_argparser()
    return p.parse_args([
        "--coil-solver", "peec", "--coil-step", "c.step",
        "--vol", "w.vol", "--wp-label", "sibc",
        "--frequency", "7000", "--sigma", "5.8e6",
        "--wp-phi-inc", "poisson",
    ] + extra)


def test_argparse_accepts_mode():
    ns = _args([])
    assert ns.wp_phi_inc == "poisson"
    # default is the validated path-integration route
    p = ci.build_argparser()
    ns0 = p.parse_args(["--coil-solver", "peec", "--coil-step", "c.step",
                        "--frequency", "7000"])
    assert ns0.wp_phi_inc == "path"


@pytest.mark.parametrize("extra,frag", [
    (["--coupling-mode", "strong", "--no-peec-proximity"], "weak-coupling"),
    (["--h1-order", "2"], "P1"),
])
def test_early_guards_fail_fast(extra, frag):
    """Unsupported combinations return an error dict BEFORE the coil
    solve -- the dummy 'c.step' would raise IOError if the guard were
    placed after it."""
    out = ci.run_inductance(_args(extra))
    assert out.get("status") == "error", out
    assert frag in out["error"], out["error"]


def test_composes_with_loop_dof_flags():
    """--wp-phi-inc poisson + --wp-loop-dof parse together."""
    ns = _args(["--wp-loop-dof", "--wp-bem-backend", "intree-dense"])
    assert ns.wp_phi_inc == "poisson" and ns.wp_loop_dof is True


def test_notebook_designspec_exposes_and_emits_both_controls():
    spec = IHDesignSpec(
        method=METHOD_PEEC_BEM,
        peec_step="coil.step",
        wp_vol="workpiece.vol",
        solver="Dense LU (small)",
        fes_order=1,
        wp_phi_inc="poisson",
        wp_loop_dof=True,
    )

    assert {"wp_phi_inc", "wp_loop_dof"} <= spec.visible_fields()
    assert {"wp_phi_inc", "wp_loop_dof"} <= set(IH_NOTEBOOK_FIELD_ORDER)
    cmd = spec.build_command(python="python", panels_dir=_PANELS)
    assert cmd[cmd.index("--wp-phi-inc") + 1] == "poisson"
    assert "--wp-loop-dof" in cmd
