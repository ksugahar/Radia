"""Contracts of the temperature-dependent surface source (--em-table).

q_i(T) = q_EM,i * q_tab(s |H_t|, T) / q_tab(|H_t|, T_ref): at the reference
temperature the EM solution is reproduced exactly, the table supplies the
local temperature (and current) dependence, and a table that does not
describe the EM run is refused.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src", "radia"))
sys.path.insert(0, os.path.join(ROOT, "src", "radia", "panels"))

import ih_thermal  # noqa: E402

RE_Z0 = 3.0e-4          # Ohm, below the Curie band


def _re_z(T):
    """Illustrative Re Z_s(T): constant, then falling to 30 % over 700-800 C."""
    T = np.asarray(T, float)
    return RE_Z0 * np.interp(T, [700.0, 800.0], [1.0, 0.3])


def _table(path, T_max=1500.0, scale=1.0):
    H = np.logspace(2, 6, 25)
    T = np.linspace(0.0, T_max, 61)
    q = 0.5 * scale * _re_z(T)[None, :] * H[:, None] ** 2
    np.savez(path, H_grid=H, T_grid=T, Zs_re=_re_z(T)[None, :].repeat(25, 0),
             Zs_im=_re_z(T)[None, :].repeat(25, 0), q_surf=q,
             meta=np.asarray(json.dumps({"frequency": 7000.0,
                                         "material": "test"})))
    return str(path)


@pytest.fixture(scope="module")
def case(tmp_path_factory):
    """Plate heated on one face; EM mesh = thermal mesh; q and |H_t| saved."""
    from netgen.occ import Box, OCCGeometry, Pnt, Z
    from ngsolve import GridFunction, H1, Mesh, x

    d = tmp_path_factory.mktemp("coupled")
    solid = Box(Pnt(0, 0, 0), Pnt(0.02, 0.02, 0.004))
    solid.faces.name = "edge"
    solid.faces.Max(Z).name = "heated"
    vol = str(d / "plate.vol")
    OCCGeometry(solid).GenerateMesh(maxh=0.004).Save(vol)
    mesh = Mesh(vol)
    H = GridFunction(H1(mesh, order=1))
    H.Set(1.0e5 * (1.0 + 20.0 * x))               # |H_t| varies over the face
    Hv = H.vec.FV().NumPy()
    q = GridFunction(H1(mesh, order=1))
    q.vec.FV().NumPy()[:] = 0.5 * RE_Z0 * Hv ** 2
    qs, hs = str(d / "q.sol"), str(d / "Ht.sol")
    q.Save(qs)
    H.Save(hs)
    for sol, quantity, unit in ((qs, ih_thermal.QSURF_QUANTITY,
                                 ih_thermal.QSURF_UNIT),
                                (hs, ih_thermal.HT_QUANTITY,
                                 ih_thermal.HT_UNIT)):
        ih_thermal.write_field_sidecar(
            sol, mesh_path=vol, mesh=mesh, fes_order=1, quantity=quantity,
            unit=unit, boundaries=["heated"],
            extra={"frequency_Hz": 7000.0})
    return {"dir": d, "vol": vol, "q": qs, "ht": hs, "mesh": mesh}


def _solve(case, table, **kw):
    import calc_heat
    from ngsolve import Mesh
    args = dict(material="steel", h_conv=0.0, heat_flux_boundaries="heated",
                qsurf_sol=case["q"], em_vol=case["vol"], dt=0.25, t_end=3.0,
                t_initial=20.0, fes_order=1, _wp_mesh=Mesh(case["vol"]),
                _write_solution=False)
    args.update(kw)
    if table:
        args.update(em_table=table, ht_sol=case["ht"], allow_frozen_ht=True)
    return calc_heat.solve_heat("<plate>", **args)


def test_below_the_curie_band_the_source_equals_the_em_solution(case):
    table = _table(case["dir"] / "t.npz")
    fixed = _solve(case, None, t_end=0.5)
    coupled = _solve(case, table, t_end=0.5)
    assert "error" not in coupled, coupled
    audit = coupled["qsurf_projection"]["temperature_dependent_source"]
    assert abs(audit["consistency"]["relative_error"]) < 1e-3
    assert coupled["T_max_C"] < 700.0
    assert coupled["Q_input_J"] == pytest.approx(fixed["Q_input_J"], rel=1e-6)
    assert coupled["T_max_C"] == pytest.approx(fixed["T_max_C"], abs=0.05)


def test_a_falling_surface_resistance_limits_the_heating(case):
    table = _table(case["dir"] / "t.npz")
    fixed = _solve(case, None, t_end=8.0)
    coupled = _solve(case, table, t_end=8.0)
    assert "error" not in coupled, coupled
    assert fixed["T_max_C"] > 900.0                   # would run away
    assert coupled["T_max_C"] < fixed["T_max_C"] - 100.0
    assert coupled["Q_input_J"] < fixed["Q_input_J"]
    nl = coupled["nonlinear_transient"]
    assert abs(nl["energy_balance_relative_error"]) < 1e-5
    assert nl["halvings"] == 0


def test_frozen_ht_needs_an_explicit_acknowledgement(case):
    import calc_heat
    from ngsolve import Mesh
    table = _table(case["dir"] / "t.npz")
    result = calc_heat.solve_heat(
        "<plate>", material="steel", h_conv=0.0,
        heat_flux_boundaries="heated", qsurf_sol=case["q"],
        em_vol=case["vol"], dt=0.25, t_end=0.25, fes_order=1,
        _wp_mesh=Mesh(case["vol"]), _write_solution=False, em_table=table,
        ht_sol=case["ht"])
    assert "--allow-frozen-ht" in result["error"]
    assert "coupled_curie_cylinder_frozen_ht" in result["error"]


def test_a_table_for_another_material_is_refused(case):
    wrong = _table(case["dir"] / "wrong.npz", scale=2.0)
    result = _solve(case, wrong, t_end=0.25)
    assert "impedance table at the reference temperature" in result["error"]


def test_leaving_the_table_is_an_error(case):
    short = _table(case["dir"] / "short.npz", T_max=300.0)
    result = _solve(case, short)
    assert "leaves the impedance table" in result["error"]


def test_current_scale_uses_the_table_not_q_scale(case):
    table = _table(case["dir"] / "t.npz")
    result = _solve(case, table, q_scale=0.5, t_end=0.25)
    assert "--ht-scale" in result["error"]
    half = _solve(case, table, ht_scale=0.5, t_end=0.25)
    full = _solve(case, table, t_end=0.25)
    assert half["Q_input_J"] == pytest.approx(0.25 * full["Q_input_J"],
                                              rel=1e-6)


def test_axisymmetric_solver_uses_ring_samples(tmp_path):
    """A 3D EM cylinder source drives a 2D (r, z) model through the table."""
    import calc_heat_axisym
    from netgen.geom2d import SplineGeometry
    from netgen.occ import Axes, Cylinder, OCCGeometry, Pnt, Z
    from ngsolve import GridFunction, H1, Mesh, x

    solid = Cylinder(Axes(Pnt(0, 0, 0), Z), r=0.02, h=0.01)
    solid.faces.name = "side"
    solid.faces.Max(Z).name = "top"
    solid.faces.Min(Z).name = "bottom"
    vol = str(tmp_path / "em.vol")
    OCCGeometry(solid).GenerateMesh(maxh=0.003).Save(vol)
    mesh = Mesh(vol)
    H = GridFunction(H1(mesh, order=1))
    H.Set(2.0e5 * (1.0 + 10.0 * x))              # azimuthally varying |H_t|
    q = GridFunction(H1(mesh, order=1))
    q.vec.FV().NumPy()[:] = 0.5 * RE_Z0 * H.vec.FV().NumPy() ** 2
    qs, hs = str(tmp_path / "q.sol"), str(tmp_path / "Ht.sol")
    q.Save(qs)
    H.Save(hs)
    for sol, quantity, unit in ((qs, ih_thermal.QSURF_QUANTITY,
                                 ih_thermal.QSURF_UNIT),
                                (hs, ih_thermal.HT_QUANTITY,
                                 ih_thermal.HT_UNIT)):
        ih_thermal.write_field_sidecar(
            sol, mesh_path=vol, mesh=mesh, fes_order=1, quantity=quantity,
            unit=unit, boundaries=["side"], extra={"frequency_Hz": 7000.0})
    geo = SplineGeometry()
    geo.AddRectangle((0, 0), (0.02, 0.01), bcs=("bottom", "side", "top",
                                                  "axis"))
    mesh2d = Mesh(geo.GenerateMesh(maxh=0.001))
    table = _table(tmp_path / "t.npz")
    common = dict(material="steel", h_conv=0.0, heat_flux_boundaries="side",
                  qsurf_sol=qs, em_vol=vol, n_phi_samples=32, dt=0.25,
                  t_end=6.0, t_initial=20.0, fes_order=2, _wp_mesh=mesh2d,
                  _write_solution=False)
    fixed = calc_heat_axisym.solve_heat_axisym("<m>", **common)
    coupled = calc_heat_axisym.solve_heat_axisym(
        "<m>", em_table=table, ht_sol=hs, allow_frozen_ht=True, **common)
    assert "error" not in coupled, coupled
    tds = coupled["qsurf_projection"]["temperature_dependent_source"]
    assert tds["azimuth_samples"] == 64
    assert abs(tds["consistency"]["relative_error"]) < 0.02
    assert fixed["T_max_C"] > 900.0
    assert coupled["T_max_C"] < fixed["T_max_C"] - 100.0
    assert coupled["nonlinear_transient"]["halvings"] == 0