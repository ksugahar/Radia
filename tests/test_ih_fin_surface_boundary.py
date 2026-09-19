"""Boundary contract between the fin-surface PEEC and the production paths.

See docs/induction_heating/FIN_PEEC_PRODUCTION_BOUNDARY.md.  No CAD, no
NGSolve: only the option validation, the CLI handed to calc_inductance and
the argparser acceptance are exercised.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from radia.simulink import ih_operator_assembly as assembly

ROOT = Path(__file__).resolve().parents[1]


def test_options_accept_fin_surface_for_weak_coupling_only():
    opts = assembly.IHOperatorAssemblyOptions(coil_step_solver="fin-surface")
    assert opts.checked().coil_step_solver == "fin-surface"
    with pytest.raises(ValueError):
        assembly.IHOperatorAssemblyOptions(
            coil_step_solver="fin-surface", coupling_mode="strong").checked()
    with pytest.raises(ValueError):
        assembly.IHOperatorAssemblyOptions(coil_step_solver="lumped").checked()
    with pytest.raises(ValueError):
        assembly.IHOperatorAssemblyOptions(fin_n_lanes=4).checked()
    with pytest.raises(ValueError):
        assembly.IHOperatorAssemblyOptions(fin_lane_grading="tip").checked()


def test_coil_solver_name_follows_geometry_and_option():
    peec = assembly.IHOperatorAssemblyOptions()
    fin = assembly.IHOperatorAssemblyOptions(coil_step_solver="fin-surface")
    assert assembly._coil_solver_name("peec", peec) == "peec"
    assert assembly._coil_solver_name("peec", fin) == "fin-surface"
    # a .vol coil is BEM-A whatever the STEP option says
    assert assembly._coil_solver_name("bem-a", fin) == "bem-a"


def _argv(backend, **kw):
    opts = assembly.IHOperatorAssemblyOptions(**kw)
    return assembly._unit_current_argv(
        Path("wp.vol"), Path("coil.step" if backend == "peec" else "coil.vol"),
        backend, opts, Path("out.msh"))


def test_unit_current_argv_for_fin_surface():
    argv = _argv("peec", coil_step_solver="fin-surface", fin_n_lanes=96,
                 fin_n_stations=12, fin_lane_grading="uniform", fin_tip_lanes=6)
    joined = " ".join(argv)
    assert "--coil-solver fin-surface" in joined
    assert "--coil-step coil.step" in joined
    assert "--fin-n-lanes 96 --fin-n-stations 12 --fin-lane-grading uniform --fin-tip-lanes 6" in joined
    assert "--peec-n-peri" not in joined and "--peec-proximity" not in joined
    assert "--coupling-mode weak" in joined
    assert "--impedance-model sibc" in joined


def test_unit_current_argv_for_series_peec_and_bema_are_unchanged():
    peec = " ".join(_argv("peec", peec_perimeter_filaments=24, peec_proximity=False))
    assert "--coil-solver peec" in peec and "--peec-n-peri 24 --no-peec-proximity" in peec
    assert "--fin-n-lanes" not in peec
    bema = " ".join(_argv("bem-a"))
    assert "--coil-solver bem-a" in bema and "--coil-vol coil.vol" in bema
    assert "--coil-source-name source --coil-sink-name sink" in bema


def test_calc_inductance_argparser_accepts_fin_surface():
    if importlib.util.find_spec("radia.panels.calc_inductance") is None:
        pytest.skip("panels package not importable here")
    from radia.panels.calc_inductance import build_argparser

    args = build_argparser().parse_args(_argv("peec", coil_step_solver="fin-surface"))
    assert args.coil_solver == "fin-surface"
    assert args.fin_route == "auto"
    # Nothing pinned by the caller stays None so the solver measures it
    # from the STEP; a pinned value still reaches the parser.
    assert (args.fin_n_lanes, args.fin_n_stations, args.fin_n_outline,
            args.fin_lane_grading, args.fin_tip_lanes) == (
        None, None, None, None, None)
    pinned = build_argparser().parse_args(
        _argv("peec", coil_step_solver="fin-surface", fin_n_lanes=96))
    assert pinned.fin_n_lanes == 96 and pinned.fin_n_stations is None
    with pytest.raises(SystemExit):
        build_argparser().parse_args(["--coil-solver", "fin", "--coil-step", "x.step"])


def test_boundary_document_exists_and_names_the_contract():
    doc = ROOT / "docs" / "induction_heating" / "FIN_PEEC_PRODUCTION_BOUNDARY.md"
    text = doc.read_text(encoding="utf-8")
    for needle in ("fin-surface", "source_type", "eddy_solver", "weak"):
        assert needle in text
