"""Axisymmetric current-drive normalization and symmetric-solver regression.

The imposed current is the amperes crossing a meridional conductor section,
``I = integral J_phi dr dz``.  It is not the r-weighted azimuthal volume
integral ``2*pi*integral J_phi*r dr dz``, whose dimensions are A*m.
"""
import math
from pathlib import Path
import sys

import pytest


ng = pytest.importorskip("ngsolve")
from netgen.meshing import (
    Element1D,
    Element2D,
    FaceDescriptor,
    Mesh as NetgenMesh,
    MeshPoint,
    Pnt,
)

_SRC = str(Path(__file__).resolve().parents[2] / "packages/radia-mcp/src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import solve_axi_eddy


def _annular_section_strip(r_inner, r_outer, height, radial_cells=16):
    """One-cell-high Q1 strip, with every magnetic-potential DOF fixed."""
    ngmesh = NetgenMesh(dim=2)
    ngmesh.SetMaterial(1, "conductor")
    names = ("bottom", "outer", "top", "inner")
    for boundary, name in enumerate(names, 1):
        ngmesh.Add(FaceDescriptor(surfnr=boundary, domin=0, bc=boundary))
        ngmesh.SetBCName(boundary - 1, name)
    radii = [
        r_inner + i * (r_outer - r_inner) / radial_cells
        for i in range(radial_cells + 1)
    ]
    bottom = [ngmesh.Add(MeshPoint(Pnt(r, 0.0, 0.0))) for r in radii]
    top = [ngmesh.Add(MeshPoint(Pnt(r, height, 0.0))) for r in radii]
    for i in range(radial_cells):
        ngmesh.Add(Element2D(1, [bottom[i], bottom[i + 1], top[i + 1], top[i]]))
        ngmesh.Add(Element1D([bottom[i], bottom[i + 1]], index=1))
        ngmesh.Add(Element1D([top[i], top[i + 1]], index=3))
    ngmesh.Add(Element1D([bottom[-1], top[-1]], index=2))
    ngmesh.Add(Element1D([top[0], bottom[0]], index=4))
    return ng.Mesh(ngmesh)


def test_axi_current_drive_matches_annular_cross_section_reference():
    """With A_phi fixed to zero, J_phi=sigma*Vc/r has a closed-form port law."""
    r_inner, r_outer, height = 0.01, 0.02, 0.006
    conductivity = 5.8e7
    omega = 2.0 * math.pi * 1000.0
    imposed_current = 3.25 - 0.4j
    mesh = _annular_section_strip(r_inner, r_outer, height)

    with ng.TaskManager():
        solution = solve_axi_eddy(
            mesh,
            ng.CoefficientFunction(1.0 / (4.0e-7 * math.pi)),
            ng.CoefficientFunction(conductivity),
            omega,
            driven_region="conductor",
            total_current=imposed_current,
            order=1,
            dirichlet="bottom|outer|top|inner",
        )

    vector_potential, voltage_coordinate = solution.components
    vc = complex(voltage_coordinate(mesh(0.015, height / 2.0)))
    # Independent calculus reference:
    # I = int_0^h int_ri^ro sigma*Vc/r dr dz
    #   = sigma*Vc*h*log(ro/ri).
    vc_reference = imposed_current / (
        conductivity * height * math.log(r_outer / r_inner)
    )
    assert abs(vc / vc_reference - 1.0) < 1.0e-8

    current_density = conductivity * (
        -1j * omega * vector_potential + voltage_coordinate / ng.x
    )
    recovered_current = ng.Integrate(current_density, mesh)
    assert abs(complex(recovered_current) / imposed_current - 1.0) < 1.0e-8


def test_axi_current_drive_rejects_zero_frequency_scaled_form():
    mesh = _annular_section_strip(0.01, 0.02, 0.006)
    with pytest.raises(ValueError, match="omega > 0"):
        solve_axi_eddy(
            mesh,
            ng.CoefficientFunction(1.0),
            ng.CoefficientFunction(1.0),
            0.0,
            driven_region="conductor",
            total_current=1.0,
            order=1,
            dirichlet="bottom|outer|top|inner",
        )


def test_axi_current_drive_rejects_conductor_touching_axis():
    mesh = _annular_section_strip(0.0, 0.02, 0.006)
    with pytest.raises(ValueError, match="touching r=0"):
        solve_axi_eddy(mesh, ng.CoefficientFunction(1.0),
                       ng.CoefficientFunction(1.0), 1000.,
                       driven_region="conductor", total_current=1., order=2)


@pytest.mark.parametrize("drive, message", [
    ({"total_current": 1.}, "requires driven_region"),
    ({"driven_region": "conductor"}, "requires total_current or applied_Vc"),
    ({"driven_region": "conductor", "total_current": 1., "applied_Vc": 2.},
     "mutually exclusive"),
])
def test_axi_rejects_incomplete_or_conflicting_drive(drive, message):
    # Argument errors must be rejected before touching any mesh or assembling.
    with pytest.raises(ValueError, match=message):
        solve_axi_eddy(None, None, None, 1000., **drive)


def test_axi_rejects_impressed_current_overlapping_current_constraint():
    mesh = _annular_section_strip(.01, .02, .006)
    with ng.TaskManager(), pytest.raises(ValueError, match="Jr must vanish"):
        solve_axi_eddy(mesh, ng.CoefficientFunction(1.), ng.CoefficientFunction(1.),
                       1000., driven_region="conductor", total_current=1.,
                       Jr=ng.CoefficientFunction(1j))
