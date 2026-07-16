"""Contract tests: the --wp-phi-inc CLI surface of calc_inductance.

The reconstruction math is locked by tests/test_phi_inc_poisson.py
(icosphere goldens); here we pin the CLI wiring: argparse accepts the
mode, the default stays synchronized with the panel, and unsupported
combinations fail fast BEFORE any expensive solve.
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
    # default is auto (poisson on the weak P1 path, path otherwise)
    p = ci.build_argparser()
    ns0 = p.parse_args(["--coil-solver", "peec", "--coil-step", "c.step",
                        "--frequency", "7000"])
    assert ns0.wp_phi_inc == "auto"


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
    assert ns.wp_phi_inc == "poisson" and ns.wp_loop_dof == "on"


def test_ih_designspec_roundtrips_wp_flags():
    """The notebook-panel IHDesignSpec emits the flags and calc_inductance's
    argparse accepts everything it emits (panel cannot silently drop or
    misspell a flag)."""
    from radia.ih_design import IHDesignSpec, METHOD_BEMA_BEM

    def _calc_argv(cmd):
        """Drop [python, script] and the calc_common-level '--output PATH'
        pair (added by calc_main, not build_argparser); everything else
        must parse STRICTLY."""
        argv, skip = [], False
        for a in cmd[2:]:
            if skip:
                skip = False
                continue
            if a == "--output":
                skip = True
                continue
            argv.append(a)
        return argv

    spec = IHDesignSpec(method=METHOD_BEMA_BEM, coil_vol="c.vol",
                        wp_vol="w.vol", solver="Dense LU (small)",
                        wp_loop_dof="on", wp_phi_inc="poisson")
    cmd = spec.build_command()
    assert cmd[cmd.index("--wp-loop-dof") + 1] == "on"
    assert cmd[cmd.index("--wp-phi-inc") + 1] == "poisson"
    ns = ci.build_argparser().parse_args(_calc_argv(cmd))
    assert ns.wp_loop_dof == "on"
    assert ns.wp_phi_inc == "poisson"
    assert ns.wp_bem_backend == "intree-dense"

    # Panel and CLI defaults agree, with no redundant flag churn.
    base = IHDesignSpec(method=METHOD_BEMA_BEM, coil_vol="c.vol",
                        wp_vol="w.vol", solver="Dense LU (small)")
    assert base.wp_loop_dof == "auto" and base.wp_phi_inc == "auto"
    cmd0 = base.build_command()
    assert "--wp-loop-dof" not in cmd0 and "--wp-phi-inc" not in cmd0
    ns0 = ci.build_argparser().parse_args(_calc_argv(cmd0))
    assert ns0.wp_loop_dof == "auto" and ns0.wp_phi_inc == "auto"

    # forcing the legacy behavior stays expressible and round-trips
    legacy = IHDesignSpec(method=METHOD_BEMA_BEM, coil_vol="c.vol",
                          wp_vol="w.vol", solver="Dense LU (small)",
                          wp_loop_dof="off", wp_phi_inc="path")
    cmdl = legacy.build_command()
    nsl = ci.build_argparser().parse_args(_calc_argv(cmdl))
    assert nsl.wp_loop_dof == "off" and nsl.wp_phi_inc == "path"


def test_ih_designspec_exposes_wp_flags_for_weak_bem_only():
    """visible_fields(): the flags appear on the calc_inductance weak
    methods and stay hidden on strong / FEM-Kelvin (which do not take
    them)."""
    from radia.ih_design import (IHDesignSpec, METHOD_BEMA_BEM,
                                 METHOD_BEMA_BEM_STRONG, METHOD_PEEC_BEM,
                                 METHOD_PEEC_FEM_KELVIN)

    for m in (METHOD_PEEC_BEM, METHOD_BEMA_BEM):
        vis = IHDesignSpec(method=m).visible_fields()
        assert {"wp_loop_dof", "wp_phi_inc"} <= vis, m
    assert {"wp_loop_dof", "wp_phi_inc"} <= set(IH_NOTEBOOK_FIELD_ORDER)
    for m in (METHOD_BEMA_BEM_STRONG, METHOD_PEEC_FEM_KELVIN):
        vis = IHDesignSpec(method=m).visible_fields()
        assert "wp_loop_dof" not in vis and "wp_phi_inc" not in vis, m


@pytest.mark.parametrize("legacy,expected", [(False, "off"), (True, "on")])
def test_ih_designspec_preserves_boolean_loop_dof_constructor(legacy, expected):
    spec = IHDesignSpec(
        method=METHOD_PEEC_BEM,
        peec_step="coil.step",
        wp_vol="workpiece.vol",
        wp_loop_dof=legacy,
    )
    assert spec.wp_loop_dof == expected
    cmd = spec.build_command(python="python", panels_dir=_PANELS)
    assert cmd[cmd.index("--wp-loop-dof") + 1] == expected


def test_ih_designspec_hides_explicit_topology_modes_when_inapplicable():
    esim = IHDesignSpec(method=METHOD_PEEC_BEM, impedance_model="Nonlinear ESIM")
    p2 = IHDesignSpec(method=METHOD_PEEC_BEM, fes_order=2)
    for spec in (esim, p2):
        assert "wp_loop_dof" not in spec.visible_fields()
        assert "wp_phi_inc" not in spec.visible_fields()
