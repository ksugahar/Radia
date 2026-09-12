"""Regression tests for independent IH thermal boundary roles.

The mesh is built in memory.  In particular, this test never calls the
problematic Netgen 2D ``NgMesh.Save()`` path.
"""
from __future__ import annotations

import math
from pathlib import Path
import sys

import pytest


PANELS = Path(__file__).resolve().parents[2] / "src" / "radia" / "panels"
if str(PANELS) not in sys.path:
    sys.path.insert(0, str(PANELS))

import calc_heat  # noqa: E402
import calc_heat_axisym  # noqa: E402
import calc_fem_kelvin  # noqa: E402


def _annulus_section_mesh(maxh=0.16):
    from netgen.geom2d import SplineGeometry
    from ngsolve import Mesh

    geometry = SplineGeometry()
    geometry.AddRectangle(
        (1.0, 0.0),
        (2.0, 1.0),
        bcs=("bottom", "heated", "top", "cooled"),
    )
    return Mesh(geometry.GenerateMesh(maxh=maxh))


def _axis_touching_section_mesh(
    maxh=0.16, *, outer_label="outer", axis_label="axis"
):
    from netgen.geom2d import SplineGeometry
    from ngsolve import Mesh

    geometry = SplineGeometry()
    geometry.AddRectangle(
        (0.0, 0.0),
        (1.0, 1.0),
        bcs=("bottom", outer_label, "top", axis_label),
    )
    return Mesh(geometry.GenerateMesh(maxh=maxh))


def _slab_mesh(maxh=0.35):
    from netgen.occ import Box, OCCGeometry, Pnt
    from ngsolve import Mesh

    solid = Box(Pnt(0.0, 0.0, 0.0), Pnt(1.0, 1.0, 1.0))
    solid.name = "workpiece"
    for face in solid.faces:
        if abs(float(face.center.x) - 1.0) < 1.0e-12:
            face.name = "heated"
        elif abs(float(face.center.x)) < 1.0e-12:
            face.name = "cooled"
        else:
            face.name = "insulated"
    return Mesh(OCCGeometry(solid).GenerateMesh(maxh=maxh))


def test_boundary_role_selector_resolves_expression_and_rejects_ambiguity():
    mesh = _annulus_section_mesh()

    selector, matched = calc_heat._resolve_boundary_role(
        mesh, "heated|top", "--heat-flux-boundaries", required=True
    )
    assert selector == "heated|top"
    assert matched == ["heated", "top"]

    with pytest.raises(ValueError, match="is required"):
        calc_heat._resolve_boundary_role(
            mesh, "", "--heat-flux-boundaries", required=True
        )
    with pytest.raises(ValueError, match="Available boundaries"):
        calc_heat._resolve_boundary_role(
            mesh, "missing", "--heat-flux-boundaries", required=True
        )


def test_axisym_heat_separates_heating_and_cooling_boundaries():
    """A radial annulus has a closed-form steady solution.

    Heat flux q enters at r=b and convection removes it at r=a.  With
    insulated z ends,

        T(r) = T_inf + q*b/(a*h) + q*b/k*log(r/a).

    The actual transient solver is run long enough to reach this state.  This
    catches accidental reuse of the heating selector for convection because
    that changes both the power balance and the temperature profile.
    """
    mesh = _annulus_section_mesh()
    q_flux = 10.0
    h_conv = 5.0
    ambient = 20.0
    conductivity = 1.0
    probe_r = 1.5

    result = calc_heat_axisym.solve_heat_axisym(
        "<in-memory-annulus>",
        material="custom",
        rho=1.0,
        cp=1.0,
        k=conductivity,
        h_conv=h_conv,
        t_ext=ambient,
        t_initial=ambient,
        emissivity=0.0,
        heat_flux_boundaries="heated",
        convection_boundaries="cooled",
        radiation_boundaries="",
        q_uniform=q_flux,
        dt=0.05,
        t_end=8.0,
        fes_order=2,
        probe_point=(probe_r, 0.5),
        _wp_mesh=mesh,
        _write_solution=False,
    )
    assert "error" not in result, result

    expected_heat_area = 2.0 * math.pi * 2.0 * 1.0
    expected_cooling_area = 2.0 * math.pi * 1.0 * 1.0
    expected_power = q_flux * expected_heat_area
    assert result["surface_area_m2"] == pytest.approx(
        expected_heat_area, rel=1.0e-10
    )
    assert result["q_surf_int_W"] == pytest.approx(
        expected_power, rel=1.0e-10
    )

    audit = result["boundary_audit"]
    assert audit["heat_flux"]["matched_boundaries"] == ["heated"]
    assert audit["convection"]["matched_boundaries"] == ["cooled"]
    assert audit["radiation"]["matched_boundaries"] == []
    assert audit["convection"]["area_m2"] == pytest.approx(
        expected_cooling_area, rel=1.0e-10
    )

    expected_probe = (
        ambient
        + q_flux * 2.0 / (1.0 * h_conv)
        + q_flux * 2.0 / conductivity * math.log(probe_r / 1.0)
    )
    assert result["T_probe_history_C"][-1] == pytest.approx(
        expected_probe, abs=0.02
    )


def test_axisym_heat_rejects_center_axis_as_heat_inflow_boundary():
    mesh = _axis_touching_section_mesh(
        outer_label="surface", axis_label="surface"
    )

    result = calc_heat_axisym.solve_heat_axisym(
        "<in-memory-axis-touching-section>",
        material="custom",
        rho=1.0,
        cp=1.0,
        k=1.0,
        h_conv=0.0,
        emissivity=0.0,
        heat_flux_boundaries="surface",
        q_uniform=10.0,
        dt=0.1,
        t_end=0.1,
        fes_order=2,
        _wp_mesh=mesh,
        _write_solution=False,
    )

    assert "error" in result
    assert "--heat-flux-boundaries" in result["error"]
    assert "natural symmetry boundary" in result["error"]
    assert "zero-revolved-area" in result["error"]


@pytest.mark.parametrize(
    ("option_name", "role_args"),
    [
        (
            "--convection-boundaries",
            {"h_conv": 5.0, "convection_boundaries": "axis"},
        ),
        (
            "--radiation-boundaries",
            {"emissivity": 0.5, "radiation_boundaries": "axis"},
        ),
    ],
)
def test_axisym_heat_rejects_center_axis_for_active_exchange_roles(
    option_name, role_args
):
    mesh = _axis_touching_section_mesh()
    kwargs = {
        "h_conv": 0.0,
        "emissivity": 0.0,
        **role_args,
    }

    result = calc_heat_axisym.solve_heat_axisym(
        "<in-memory-axis-touching-section>",
        material="custom",
        rho=1.0,
        cp=1.0,
        k=1.0,
        heat_flux_boundaries="outer",
        q_uniform=10.0,
        dt=0.1,
        t_end=0.1,
        fes_order=2,
        _wp_mesh=mesh,
        _write_solution=False,
        **kwargs,
    )

    assert "error" in result
    assert option_name in result["error"]
    assert "natural symmetry boundary" in result["error"]


def test_axisym_heat_outer_surface_excludes_axis_and_has_expected_power():
    mesh = _axis_touching_section_mesh()

    result = calc_heat_axisym.solve_heat_axisym(
        "<in-memory-axis-touching-section>",
        material="custom",
        rho=1.0,
        cp=1.0,
        k=1.0,
        h_conv=0.0,
        emissivity=0.0,
        heat_flux_boundaries="outer",
        q_uniform=10.0,
        dt=0.1,
        t_end=0.1,
        fes_order=2,
        _wp_mesh=mesh,
        _write_solution=False,
    )

    assert "error" not in result, result
    assert result["surface_area_m2"] == pytest.approx(
        2.0 * math.pi, rel=1.0e-10
    )
    assert result["q_surf_int_W"] == pytest.approx(
        20.0 * math.pi, rel=1.0e-10
    )


def test_3d_heat_separates_heating_and_cooling_boundaries():
    """The 3D solver reproduces the linear insulated-slab solution."""
    mesh = _slab_mesh()
    result = calc_heat.solve_heat(
        "<in-memory-slab>",
        material="custom",
        rho=1.0,
        cp=1.0,
        k=1.0,
        h_conv=5.0,
        t_ext=20.0,
        t_initial=20.0,
        emissivity=0.0,
        heat_flux_boundaries="heated",
        convection_boundaries="cooled",
        q_uniform=10.0,
        dt=0.05,
        t_end=8.0,
        fes_order=1,
        probe_point=(0.5, 0.5, 0.5),
        _wp_mesh=mesh,
        _write_solution=False,
    )
    assert "error" not in result, result
    assert result["q_surf_int_W"] == pytest.approx(10.0, rel=1.0e-10)
    assert result["boundary_audit"]["heat_flux"][
        "matched_boundaries"
    ] == ["heated"]
    assert result["boundary_audit"]["convection"][
        "matched_boundaries"
    ] == ["cooled"]
    assert result["T_probe_history_C"][-1] == pytest.approx(27.0, abs=0.03)


def test_spatial_qsurf_high_order_transfer_fails_before_loading_files():
    assert calc_fem_kelvin.QSURF_HANDOFF_ORDER == 1
    assert calc_heat.QSURF_HANDOFF_ORDER == 1
    with pytest.raises(ValueError, match="--qsurf-order 1 only"):
        calc_heat._validate_qsurf_transfer_order(2)
