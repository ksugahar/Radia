"""Contracts of the coupled axisymmetric eddy-current + heat solver."""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "src", "radia"))
sys.path.insert(0, os.path.join(ROOT, "src", "radia", "panels"))

from radia import ih_axisym_coupled as C  # noqa: E402
from radia.analytical_formulas.induction_heating import (  # noqa: E402
    cylinder_axial_eddy_loss,
)

A, L, RC, TC, ROUT = 0.02, 0.01, 0.024, 5e-4, 0.03


def _slice(delta):
    """Axial slice of an infinite cylinder inside a current sheet.

    With only the axis Dirichlet, the natural condition H_t = 0 on the ends
    and on the outer radius reproduces the infinite solenoid exactly.
    """
    from netgen.occ import Glue, MoveTo, OCCGeometry, X
    from ngsolve import Mesh

    sk = min(4 * delta, A / 2)
    core = MoveTo(0, 0).Rectangle(A - sk, L).Face()
    core.faces.name = "wp"
    core.edges.Min(X).name = "axis"
    # the skin depth opens up to millimetres through the Curie band
    core.maxh = 6e-4
    shell = MoveTo(A - sk, 0).Rectangle(sk, L).Face()
    shell.faces.name = "wp"
    shell.maxh = min(delta / 3, 1e-3)
    coil = MoveTo(RC, 0).Rectangle(TC, L).Face()
    coil.faces.name = "coil"
    coil.maxh = TC / 2
    air = MoveTo(0, 0).Rectangle(ROUT, L).Face()
    air.faces.name = "air"
    air.edges.Min(X).name = "axis"
    air = air - core - shell - coil
    return Mesh(OCCGeometry(Glue([air, core, shell, coil]), dim=2)
                .GenerateMesh(maxh=2e-3))


@pytest.mark.parametrize("sigma, mu_r, rtol", [(5.8e6, 100.0, 5e-4),
                                               (1.0e6, 1.0, 1e-5)])
def test_eddy_loss_matches_the_bessel_solution(sigma, mu_r, rtol):
    f = 7000.0
    delta = math.sqrt(2 / (2 * math.pi * f * C.MU0 * mu_r * sigma))
    mesh = _slice(delta)
    em = C.AxisymEddyCurrent(mesh, frequency=f, workpiece="wp",
                             coils={"coil": 1000.0}, dirichlet="axis")
    n = len(em.workpiece_elements)
    rec = em.solve(np.full(n, sigma), np.full(n, mu_r))
    exact = cylinder_axial_eddy_loss(1000.0 / L, A, 2 * math.pi * f, sigma,
                                     mu_r) * L
    assert rec["P_joule_W"] == pytest.approx(exact, rel=rtol)
    assert abs(rec["power_balance_relative_error"]) < 1e-7


def test_dirichlet_must_include_the_axis():
    mesh = _slice(1e-3)
    with pytest.raises(ValueError, match="axis"):
        C.AxisymEddyCurrent(mesh, frequency=1e3, workpiece="wp",
                            coils={"coil": 1.0}, dirichlet="outer")


def _materials(curie=True):
    from radia import ih_thermal_material as itm
    em = C.EMMaterialTable(
        T=[0.0, 400.0, 700.0, 770.0, 1500.0],
        sigma=[6.0e6, 2.4e6, 1.3e6, 1.15e6, 0.85e6],
        mu_r=[100.0, 90.0, 80.0, 1.0, 1.0] if curie else [1.0] * 5,
        source="illustrative")
    th = itm.ThermalMaterial.constant(7800.0, 600.0, 35.0)
    return em, th


def test_staggered_run_conserves_energy_through_the_curie_band():
    from radia import ih_heat_transient as iht
    em_mat, th = _materials()
    mesh = _slice(2.5e-4)
    result, gfT, em = C.run_coupled(
        mesh, frequency=7000.0, workpiece="wp", coils={"coil": 3000.0},
        dirichlet="axis", em_material=em_mat, thermal_material=th,
        boundaries=iht.HeatBoundaryTerms(), dt=0.05, t_end=1.5)
    nl = result["nonlinear_transient"]
    assert abs(nl["energy_balance_relative_error"]) < 1e-6
    assert nl["halvings"] == 0
    # the stored energy equals the time integral of the EM Joule power
    joule = sum(h["P_W"] for h in result["history"]) * 0.05
    assert nl["energy_in_J"] == pytest.approx(joule, rel=1e-6)
    hist = result["em_history"]
    assert all(abs(h["power_balance_relative_error"]) < 1e-7 for h in hist)
    assert max(h["T_elem_max_C"] for h in hist) > 770.0
    # the skin depth opens up once the surface is above the Curie band
    assert hist[-1]["skin_depth_max_m"] > 20 * hist[0]["skin_depth_min_m"]


def test_command_line_writes_a_reloadable_temperature(tmp_path):
    import calc_ih_axisym_coupled as cli
    from radia import ih_thermal, ih_thermal_post

    mesh = _slice(2.45e-3)                 # skin depth of 6e6 S/m at 7 kHz
    vol = tmp_path / "slice.vol"
    mesh.ngmesh.Save(str(vol))
    table = tmp_path / "em.csv"
    table.write_text("T_C,sigma_S_m,mu_r\n0,6e6,1\n2000,6e6,1\n",
                     encoding="utf-8")
    out = tmp_path / "T.sol"
    result = cli.solve_ih_axisym_coupled(
        str(vol), workpiece="wp", coils={"coil": 2000.0}, frequency=7000.0,
        em_dirichlet="axis", em_material_table=str(table), material="custom",
        rho=7800.0, cp=600.0, k=35.0, dt=0.1, t_end=0.3,
        exposure_thresholds=[100.0], temperature_output=str(out))
    assert "error" not in result, result
    ex = result["thermal_exposure"]
    assert ex["region"] == "wp"
    assert ex["T_min_C"] >= 20.0 - 1e-6        # air is not counted
    mesh2, gf, audit = ih_thermal.load_field(
        str(out), quantity=ih_thermal.TEMPERATURE_QUANTITY)
    assert audit["provenance"] == "sidecar-verified"
    again = ih_thermal_post.thermal_exposure(mesh2, gf, [100.0],
                                             axisymmetric=True, region="wp")
    assert again["T_max_C"] == pytest.approx(ex["T_max_C"], rel=1e-12)
    # the post-processing command reads region and geometry from the sidecar
    report = tmp_path / "exposure.json"
    assert ih_thermal_post.main(["exposure", "--temperature", str(out),
                                 "--thresholds", "100",
                                 "--output", str(report)]) == 0
    import json
    cli_ex = json.loads(report.read_text(encoding="utf-8"))
    assert cli_ex["region"] == "wp"
    assert cli_ex["T_max_C"] == pytest.approx(ex["T_max_C"], rel=1e-12)


def test_unnamed_axis_segments_are_refused():
    """A_phi = 0 must hold on the whole axis, not only where it is named."""
    from netgen.occ import Glue, MoveTo, OCCGeometry, X
    from ngsolve import Mesh

    core = MoveTo(0, 0).Rectangle(A, L).Face()
    core.faces.name = "wp"                   # its axis edge stays unnamed
    air = MoveTo(0, 0).Rectangle(ROUT, L).Face()
    air.faces.name = "air"
    air.edges.Min(X).name = "axis"
    coil = MoveTo(RC, 0).Rectangle(TC, L).Face()
    coil.faces.name = "coil"
    air = air - core - coil
    mesh = Mesh(OCCGeometry(Glue([air, core, coil]), dim=2)
                .GenerateMesh(maxh=3e-3))
    with pytest.raises(ValueError, match="r = 0 axis"):
        C.AxisymEddyCurrent(mesh, frequency=1e3, workpiece="wp",
                            coils={"coil": 1.0}, dirichlet="axis")


def test_material_table_starting_above_zero_is_not_tripped_by_air(tmp_path):
    """The temperature lives on the workpiece only; air vertices (value 0 in
    the definedon space) must not count as leaving a 20 C.. table."""
    from radia import ih_heat_transient as iht, ih_thermal_material as itm
    em_mat, _ = _materials(curie=False)
    th = itm.ThermalMaterial(rho=7800.0, T=[20.0, 1500.0], k=[40.0, 30.0],
                             cp=[450.0, 700.0], source="20C table")
    result, gfT, em = C.run_coupled(
        _slice(2.45e-3), frequency=7000.0, workpiece="wp",
        coils={"coil": 1500.0}, dirichlet="axis", em_material=em_mat,
        thermal_material=th, boundaries=iht.HeatBoundaryTerms(), dt=0.1,
        t_end=0.3)
    assert result["history"][-1]["T_max_nodal_C"] > 20.0
    assert result["nonlinear_transient"]["table_extrapolation_C"] == 0.0