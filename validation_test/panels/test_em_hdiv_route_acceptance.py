"""Electromagnet block HDiv-VIM route: real solve against the production solve.

``calc_accel_hdiv.solve_hdiv`` is what the Simulink Electromagnet block runs
for the HDiv-VIM method.  It must report the solver's own convergence data
(not placeholders) and say which panel options the HDiv route does not use.
The independent reference is ``radia.vim.Solve`` called directly on the same
mesh, B-H table and coil field; the soft iron must also strengthen the coil
field at a nearby point.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
PANELS = REPO / "src" / "radia" / "panels"
SAMPLE_BH = PANELS / "samples" / "em_sample_bh.txt"

rad = pytest.importorskip("radia")
ngsolve = pytest.importorskip("ngsolve")
if str(PANELS) not in sys.path:
    sys.path.insert(0, str(PANELS))

COIL = """
import os, sys
sys.path.insert(0, {radia_src!r})
from coil_builder import CoilBuilder


def build_coil():
    return (CoilBuilder(current={current})
            .set_start([0.04, 0.0, 0.0])
            .set_cross_section(0.004, 0.004)
            .add_arc(0.04, 360))
"""


def _iron_block_vol(path):
    from netgen.occ import Box, OCCGeometry, Pnt

    block = Box(Pnt(-0.01, -0.01, 0.015), Pnt(0.01, 0.01, 0.035))
    block.mat("yoke")
    with ngsolve.TaskManager():
        mesh = OCCGeometry(block).GenerateMesh(maxh=0.006)
    mesh.Save(str(path))


def _coil_script(path, current):
    path.write_text(
        COIL.format(radia_src=str(REPO / "src" / "radia"), current=current), encoding="utf-8"
    )


@pytest.mark.parametrize("current", [2000.0, 200000.0])
def test_hdiv_panel_reports_the_production_solve(tmp_path, current):
    from calc_accel_hdiv import solve_hdiv
    from calc_common import EMMaterial

    vol = tmp_path / "yoke.vol"
    coil = tmp_path / "coil.py"
    _iron_block_vol(vol)
    _coil_script(coil, current)
    material = EMMaterial.from_name("steel", bh_file=str(SAMPLE_BH))
    bh_table, _ = material.get_bh_curve()

    result = solve_hdiv(
        coil_script=str(coil),
        vol_file=str(vol),
        mat=material,
        solver=2,
        max_iter=30,
        tol=1e-3,
        relax=0.3,
    )
    assert "error" not in result, result.get("error")
    assert result["converged"] is True
    assert result["residual_kind"] == "nonlinear_final_relative_residual"
    assert result["residual"] is not None and result["residual"] > 0.0
    assert result["residual"] <= result["residual_tolerance"]
    assert result["panel_options_not_applied"] == {
        "solver": "HACApK",
        "max_iter": 30,
        "tol": 1e-3,
        "relax": 0.3,
    }

    # Independent reference: the production solve on the same inputs.
    import radia.vim as vim
    import importlib.util

    rad.UtiDelAll()
    spec = importlib.util.spec_from_file_location("acceptance_coil", coil)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    source = rad.ObjCnt(module.build_coil().to_radia())
    mesh = ngsolve.Mesh(str(vol))
    with ngsolve.TaskManager():
        direct = vim.Solve(mesh, H_ext=rad.RadiaField(source, "h"), bh_table=bh_table)
    assert np.allclose(result["M_avg"], direct["M_avg"], rtol=1e-9, atol=1e-6)
    assert result["iterations"] == int(direct["iters"])
    # Both final residuals are the solver's own and meet its tolerance; at
    # round-off level they differ between threaded runs, so compare the gate.
    assert float(direct["nonlinear_final_relative_residual"]) <= result["residual_tolerance"]
    assert result["residual_tolerance"] == direct["nonlinear_residual_tolerance"]

    # Physical check: the magnetized block raises the coil's field at the
    # loop centre, which lies on the block axis just below it.
    coil_only = np.linalg.norm(rad.Fld(source, "b", [0.0, 0.0, 0.0]))
    assert result["B_center_mag"] > coil_only
    rad.UtiDelAll()
