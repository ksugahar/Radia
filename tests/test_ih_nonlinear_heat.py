"""Contracts of the temperature-dependent heat integrator.

radia.ih_thermal_material.ThermalMaterial + radia.ih_heat_transient
.NonlinearHeatStepper, and their use by calc_heat --material-table.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src", "radia"))
sys.path.insert(0, os.path.join(ROOT, "src", "radia", "panels"))

from radia import ih_heat_transient as iht
from radia import ih_thermal_material as itm


def _plate(maxh=0.004, t=0.004):
    """0.02 x 0.02 x t plate; both large faces are 'heated'."""
    from netgen.occ import Box, OCCGeometry, Pnt, Z
    from ngsolve import Mesh

    solid = Box(Pnt(0, 0, 0), Pnt(0.02, 0.02, t))
    solid.faces.name = "edge"
    solid.faces.Max(Z).name = "heated"
    solid.faces.Min(Z).name = "heated"
    return Mesh(OCCGeometry(solid).GenerateMesh(maxh=maxh))


def _stepper(mesh, material, q, order=2, **kw):
    from ngsolve import CF, GridFunction, H1
    gf = GridFunction(H1(mesh, order=order))
    gf.Set(CF(20.0))
    st = iht.NonlinearHeatStepper(
        gf, material, iht.HeatBoundaryTerms(heat_flux="heated", **kw),
        q_source=CF(q))
    return gf, st


def _steel_like():
    # an illustrative table (not alloy data): k falls, cp peaks, latent heat
    return itm.ThermalMaterial(
        rho=7800.0, T=[0.0, 700.0, 760.0, 1400.0],
        k=[50.0, 30.0, 27.0, 30.0], cp=[450.0, 750.0, 900.0, 650.0],
        latent_heat=6.0e4, latent_range=(720.0, 770.0), source="test")


def test_constant_material_reproduces_linear_backward_euler():
    from ngsolve import (BilinearForm, CF, GridFunction, H1, LinearForm, ds,
                         dx, grad, InnerProduct)

    mesh = _plate()
    mat = itm.ThermalMaterial(rho=7800.0, T=[0.0, 2000.0], k=[40.0, 40.0001],
                              cp=[460.0, 460.0], source="test")
    gf, st = _stepper(mesh, mat, 2.0e6)
    for _ in range(5):
        st.advance(0.1)

    fes = H1(mesh, order=2)
    u, v = fes.TnT()
    ref = GridFunction(fes)
    ref.Set(CF(20.0))
    a = BilinearForm(fes, symmetric=True)
    a += 7800 * 460 / 0.1 * u * v * dx + 40.0 * InnerProduct(grad(u), grad(v)) * dx
    a.Assemble()
    inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    for _ in range(5):
        f = LinearForm(fes)
        f += 7800 * 460 / 0.1 * ref * v * dx + 2.0e6 * v * ds("heated")
        f.Assemble()
        ref.vec.data = inv * f.vec
    diff = np.max(np.abs(gf.vec.FV().NumPy() - ref.vec.FV().NumPy()))
    assert diff < 0.01                    # k differs by 2.5e-6 relative
    assert st.audit.max_newton_iterations <= 3


def test_energy_balance_is_exact_with_latent_heat():
    mesh = _plate()
    gf, st = _stepper(mesh, _steel_like(), 3.0e6, convection="edge",
                      h_conv=50.0, t_ext=20.0)
    for _ in range(20):
        st.advance(0.25)
    audit = st.audit.as_dict()
    assert np.max(gf.vec.FV().NumPy()) > 770.0      # crossed the latent band
    assert abs(audit["energy_balance_relative_error"]) < 1e-5
    assert audit["energy_loss_J"] > 0.0
    assert audit["halvings"] == 0


def test_uniform_heating_follows_the_enthalpy_curve():
    """A thin, highly conductive plate heated on both faces is isothermal,
    so T(t) must follow H(T) = H(T0) + q A t / V through the latent band."""
    mat = _steel_like()
    fast = itm.ThermalMaterial(rho=mat.rho, T=mat.T, k=mat.k * 1.0e4,
                               cp=mat.cp, latent_heat=mat.latent_heat,
                               latent_range=mat.latent_range, source="fast")
    mesh = _plate(maxh=0.005, t=0.002)
    q = 1.0e6                    # 1e9 W/m^3: about 900 C after 5 s
    gf, st = _stepper(mesh, fast, q, order=1)
    area, vol = 2 * 0.02 * 0.02, 0.02 * 0.02 * 0.002
    t = 0.0
    for _ in range(20):
        st.advance(0.25)
        t += 0.25
        H_expected = float(fast.enthalpy(20.0)) + q * area * t / vol
        grid = np.linspace(0.0, 1400.0, 140001)
        T_expected = float(np.interp(H_expected, fast.enthalpy(grid), grid))
        T = gf.vec.FV().NumPy()
        assert np.ptp(T) < 0.5
        assert T.mean() == pytest.approx(T_expected, abs=0.5)
    assert T.mean() > 770.0


def test_leaving_the_table_is_an_error_unless_allowed():
    mesh = _plate()
    short = itm.ThermalMaterial(rho=7800.0, T=[0.0, 100.0], k=[40.0, 41.0],
                                cp=[460.0, 470.0], source="short")
    gf, st = _stepper(mesh, short, 5.0e6)
    with pytest.raises(ValueError, match="material table range"):
        for _ in range(20):
            st.advance(0.5)
    ok = itm.ThermalMaterial(rho=7800.0, T=[0.0, 100.0], k=[40.0, 41.0],
                             cp=[460.0, 470.0], source="short",
                             allow_extrapolation=True)
    gf, st = _stepper(mesh, ok, 5.0e6)
    for _ in range(4):
        st.advance(0.5)
    assert st.audit.table_extrapolation_C > 0.0


def test_calc_heat_material_table_end_to_end(tmp_path):
    import calc_heat

    table = tmp_path / "steel.csv"
    table.write_text("T_C,k_W_mK,cp_J_kgK\n0,50,450\n700,30,750\n"
                     "760,27,900\n1400,30,650\n", encoding="utf-8")
    result = calc_heat.solve_heat(
        "<in-memory-plate>", material="steel", h_conv=0.0,
        heat_flux_boundaries="heated", q_uniform=3.0e6, dt=0.25, t_end=2.5,
        t_initial=20.0, fes_order=2, _wp_mesh=_plate(), _write_solution=False,
        material_table=str(table), latent_heat=6.0e4,
        latent_range=(720.0, 770.0))
    assert "error" not in result, result
    audit = result["nonlinear_transient"]
    assert audit["steps"] == 10 and audit["halvings"] == 0
    assert abs(audit["energy_balance_relative_error"]) < 1e-5
    assert result["thermal_material"]["temperature_dependent"]
    # the stored energy equals the heat input (no losses here)
    assert audit["energy_stored_J"] == pytest.approx(
        result["Q_input_J"], rel=1e-5)


def test_axisymmetric_material_table_conserves_revolved_energy(tmp_path):
    import calc_heat_axisym
    from netgen.geom2d import SplineGeometry
    from ngsolve import Mesh

    geo = SplineGeometry()
    geo.AddRectangle((0, 0), (0.02, 0.01), bcs=("bottom", "outer", "top",
                                                  "axis"))
    mesh = Mesh(geo.GenerateMesh(maxh=0.002))
    table = tmp_path / "steel.csv"
    table.write_text("T_C,k_W_mK,cp_J_kgK\n0,50,450\n700,30,750\n"
                     "760,27,900\n1400,30,650\n", encoding="utf-8")
    result = calc_heat_axisym.solve_heat_axisym(
        "<in-memory-disc>", material="steel", h_conv=20.0, t_ext=20.0,
        heat_flux_boundaries="outer", convection_boundaries="top|bottom",
        q_uniform=4.0e6, dt=0.25, t_end=2.5, t_initial=20.0, fes_order=2,
        _wp_mesh=mesh, _write_solution=False, material_table=str(table),
        latent_heat=6.0e4, latent_range=(720.0, 770.0))
    assert "error" not in result, result
    audit = result["nonlinear_transient"]
    assert abs(audit["energy_balance_relative_error"]) < 1e-5
    # the residual and the matrix share one quadrature; with the 2 pi r
    # weight a mismatch made Newton diverge and every step halve
    assert audit["halvings"] == 0 and audit["max_newton_iterations"] <= 8
    # revolved heat input: q * 2 pi R H * t
    assert audit["energy_in_J"] == pytest.approx(
        4.0e6 * 2 * np.pi * 0.02 * 0.01 * 2.5, rel=1e-9)