"""Contract tests: the workpiece incident-potential route of calc_inductance.

The reconstruction math is locked by tests/test_phi_inc_poisson.py
(icosphere goldens).  Here we pin the ROUTE contract after the
2026-07-17 removal of the legacy selectable path: the P1 weak route is
surface-Poisson ALWAYS (basis-determined, not a knob), the removed
``--wp-phi-inc`` flag must not silently return, and the notebook-panel
IHDesignSpec emits exactly what calc_inductance's argparse accepts.
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


def test_wp_phi_inc_flag_is_gone():
    """The legacy path-integration route on P1 weak was a known
    branch-cut-wall bias -- it must not be selectable, so the flag
    itself is gone (a superseded twin must not stay reachable)."""
    p = ci.build_argparser()
    with pytest.raises(SystemExit):
        p.parse_args(["--coil-solver", "peec", "--coil-step", "c.step",
                      "--frequency", "7000", "--wp-phi-inc", "path"])


def test_loop_dof_off_is_gone():
    """'off' would re-enable the +25-30% genus-1 legacy solve on a
    mesh where the loop DOF is available -- deliberately removed."""
    p = ci.build_argparser()
    with pytest.raises(SystemExit):
        p.parse_args(["--coil-solver", "peec", "--coil-step", "c.step",
                      "--frequency", "7000", "--wp-loop-dof", "off"])


def test_ih_designspec_roundtrips_wp_flags():
    """The notebook-panel IHDesignSpec emits only flags calc_inductance
    accepts (panel cannot silently drop or misspell a flag), and the
    removed phi_inc knob is gone from the spec too."""
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

    assert not hasattr(IHDesignSpec(), "wp_phi_inc")

    spec = IHDesignSpec(method=METHOD_BEMA_BEM, coil_vol="c.vol",
                        wp_vol="w.vol", solver="Dense LU (small)",
                        wp_loop_dof="on")
    cmd = spec.build_command()
    assert cmd[cmd.index("--wp-loop-dof") + 1] == "on"
    ns = ci.build_argparser().parse_args(_calc_argv(cmd))
    assert ns.wp_loop_dof == "on"
    assert ns.wp_bem_backend == "intree-dense"

    # panel default == CLI default == auto: no flag churn, and the
    # resolved behavior is the FIXED one wherever it applies
    base = IHDesignSpec(method=METHOD_BEMA_BEM, coil_vol="c.vol",
                        wp_vol="w.vol", solver="Dense LU (small)")
    assert base.wp_loop_dof == "auto"
    cmd0 = base.build_command()
    assert "--wp-loop-dof" not in cmd0 and "--wp-phi-inc" not in cmd0
    ns0 = ci.build_argparser().parse_args(_calc_argv(cmd0))
    assert ns0.wp_loop_dof == "auto"


def test_ih_designspec_exposes_loop_dof_for_weak_bem_only():
    """visible_fields(): the loop-DOF mode appears on the
    calc_inductance weak methods and stays hidden on strong /
    FEM-Kelvin (which do not take it)."""
    from radia.ih_design import (IHDesignSpec, METHOD_BEMA_BEM,
                                 METHOD_BEMA_BEM_STRONG, METHOD_PEEC_BEM,
                                 METHOD_PEEC_FEM_KELVIN)

    for m in (METHOD_PEEC_BEM, METHOD_BEMA_BEM):
        vis = IHDesignSpec(method=m).visible_fields()
        assert "wp_loop_dof" in vis, m
        assert "wp_phi_inc" not in vis, m
    for m in (METHOD_BEMA_BEM_STRONG, METHOD_PEEC_FEM_KELVIN):
        vis = IHDesignSpec(method=m).visible_fields()
        assert "wp_loop_dof" not in vis, m


@pytest.mark.parametrize("legacy,expected", [(False, "auto"), (True, "on")])
def test_ih_designspec_preserves_boolean_loop_dof_constructor(legacy, expected):
    spec = IHDesignSpec(
        method=METHOD_PEEC_BEM,
        peec_step="coil.step",
        wp_vol="workpiece.vol",
        wp_loop_dof=legacy,
    )
    assert spec.wp_loop_dof == expected
    cmd = spec.build_command(python="python", panels_dir=_PANELS)
    if expected == "on":
        assert cmd[cmd.index("--wp-loop-dof") + 1] == "on"
    else:
        assert "--wp-loop-dof" not in cmd


def test_ih_designspec_hides_loop_dof_when_inapplicable():
    esim = IHDesignSpec(method=METHOD_PEEC_BEM, impedance_model="Nonlinear ESIM")
    p2 = IHDesignSpec(method=METHOD_PEEC_BEM, fes_order=2)
    for spec in (esim, p2):
        assert "wp_loop_dof" not in spec.visible_fields()
