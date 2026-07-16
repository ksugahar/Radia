"""Contract tests: the --wp-loop-dof CLI surface of calc_inductance.

The numeric physics is locked by the analytic shorted-ring golden
(validation_test/bem/test_loop_extension_ring.py); here we pin the CLI
wiring: argparse accepts the flag, the unsupported combinations fail
fast BEFORE any expensive solve, and a genus-0 workpiece raises with an
actionable message.
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

_SAMPLES = _PANELS / "samples"
_DEMO_COIL_STEP = _SAMPLES / "ih_fem_kelvin_demo_coil.step"
_DEMO_VOL = _SAMPLES / "ih_fem_kelvin_demo.vol"


def _args(extra):
    p = ci.build_argparser()
    return p.parse_args([
        "--coil-solver", "peec", "--coil-step", "c.step",
        "--vol", "w.vol", "--wp-label", "sibc",
        "--frequency", "7000", "--sigma", "5.8e6",
        "--wp-loop-dof",
    ] + extra)


def test_argparse_accepts_flag():
    ns = _args(["--wp-bem-backend", "intree-dense"])
    assert ns.wp_loop_dof is True
    # default off
    p = ci.build_argparser()
    ns0 = p.parse_args(["--coil-solver", "peec", "--coil-step", "c.step",
                        "--frequency", "7000"])
    assert ns0.wp_loop_dof is False


@pytest.mark.parametrize("extra,frag", [
    (["--coupling-mode", "strong", "--no-peec-proximity",
      "--wp-bem-backend", "intree-dense"], "weak-coupling"),
    (["--impedance-model", "esim", "--wp-bem-backend", "intree-dense"],
     "linear SIBC"),
    ([], "intree-dense"),                       # default backend = hacapk
    (["--wp-bem-backend", "intree-dense", "--h1-order", "2"], "P1"),
])
def test_early_guards_fail_fast(extra, frag):
    """Unsupported combinations return an error dict BEFORE the coil
    solve -- the dummy 'c.step' would raise IOError if the guard were
    placed after it."""
    out = ci.run_inductance(_args(extra))
    assert out.get("status") == "error", out
    assert frag in out["error"], out["error"]


@pytest.mark.skipif(
    not (_DEMO_COIL_STEP.is_file() and _DEMO_VOL.is_file()),
    reason="demo fixtures not present")
def test_genus0_workpiece_raises():
    """The demo workpiece is genus 0 (sphere-like): the loop DOF must
    refuse with an actionable message instead of silently no-op'ing."""
    p = ci.build_argparser()
    ns = p.parse_args([
        "--coil-solver", "peec", "--coil-step", str(_DEMO_COIL_STEP),
        "--vol", str(_DEMO_VOL), "--wp-label", "sibc",
        "--frequency", "7000", "--current", "1",
        "--sigma", "5.8e6", "--mu-r", "100", "--half-thickness", "0.005",
        "--wp-bem-backend", "intree-dense", "--wp-loop-dof",
    ])
    with pytest.raises(ValueError, match="genus-1"):
        ci.run_inductance(ns)
