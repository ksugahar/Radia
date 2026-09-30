"""Current drive stays on ``driven_region`` when other conductors are present.

The impressed field Vc and the net-current constraint act on the driven
conductor only.  A second conductor in ``sigma`` carries induced eddy current
(Vc = 0 there) and must not share the imposed current.
"""
import math

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
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue, X, Y

from radia_mcp.radia_ngsolve.solve import solve_axi_eddy, solve_planar_eddy


MU0 = 4.0e-7 * math.pi


def _two_wire_mesh(radius=1.0e-3, spacing=3.0e-3, r_far=15.0e-3):
    wires = []
    for name, xc in (("driven", -spacing / 2.0), ("passive", spacing / 2.0)):
        wire = WorkPlane().Circle(xc, 0.0, radius).Face()
        wire.faces.name = name
        wire.maxh = radius / 4.0
        wires.append(wire)
    box = MoveTo(-r_far, -r_far).Rectangle(2.0 * r_far, 2.0 * r_far).Face()
    air = box - wires[0] - wires[1]
    air.faces.name = "air"
    for edge in (air.edges.Max(X), air.edges.Min(X), air.edges.Max(Y), air.edges.Min(Y)):
        edge.name = "outer"
    geometry = OCCGeometry(Glue([air] + wires), dim=2)
    return ng.Mesh(geometry.GenerateMesh(maxh=r_far / 6.0))


def test_planar_current_drive_is_confined_to_driven_region():
    conductivity, omega, imposed_current = 5.8e7, 2.0 * math.pi * 2000.0, 2.0 - 0.5j
    mesh = _two_wire_mesh()
    sigma = mesh.MaterialCF({"driven": conductivity, "passive": conductivity}, default=0.0)
    nu = ng.CoefficientFunction(1.0 / MU0)

    with ng.TaskManager():
        solution = solve_planar_eddy(mesh, nu, sigma, omega, driven_region="driven",
                                     total_current=imposed_current, order=2)
        az, vc = solution.components
        vc_value = complex(vc(mesh(-1.5e-3, 0.0)))
        driven_current = complex(ng.Integrate(
            sigma * (-1j * omega * az + vc), mesh, definedon=mesh.Materials("driven")))
        passive_current = complex(ng.Integrate(
            sigma * (-1j * omega * az), mesh, definedon=mesh.Materials("passive")))
        # Independent route: the prescribed-field solve with the recovered Vc
        # on the driven wire alone must give the same potential.
        az_voltage = solve_planar_eddy(mesh, nu, sigma, omega, driven_region="driven",
                                       applied_Ez=vc_value, order=2)
        mismatch = math.sqrt(ng.Integrate(ng.Norm(az - az_voltage) ** 2, mesh))
        scale = math.sqrt(ng.Integrate(ng.Norm(az) ** 2, mesh))

    assert abs(driven_current / imposed_current - 1.0) < 1.0e-8
    # The induced return current is real physics but must stay a fraction of I.
    assert 0.0 < abs(passive_current) < 0.5 * abs(imposed_current)
    assert mismatch / scale < 1.0e-8


@pytest.mark.parametrize("region, conductivities, message", [
    ("missing", {"driven": 1.0e6}, "matches no mesh material"),
    ("passive", {"driven": 1.0e6}, "zero conductivity"),
])
def test_planar_current_drive_rejects_bad_region(region, conductivities, message):
    mesh = _two_wire_mesh()
    sigma = mesh.MaterialCF(conductivities, default=0.0)
    with pytest.raises(ValueError, match=message):
        solve_planar_eddy(mesh, ng.CoefficientFunction(1.0 / MU0), sigma, 1.0e3,
                          driven_region=region, total_current=1.0, order=1)


def _split_annular_strip(r_inner, r_split, r_outer, height, cells=8):
    """One-cell-high Q1 strip, driven | passive, with every A DOF fixed."""
    ngmesh = NetgenMesh(dim=2)
    ngmesh.SetMaterial(1, "driven")
    ngmesh.SetMaterial(2, "passive")
    for boundary, name in enumerate(("bottom", "outer", "top", "inner"), 1):
        ngmesh.Add(FaceDescriptor(surfnr=boundary, domin=0, bc=boundary))
        ngmesh.SetBCName(boundary - 1, name)
    radii = [r_inner + i * (r_split - r_inner) / cells for i in range(cells)]
    radii += [r_split + i * (r_outer - r_split) / cells for i in range(cells + 1)]
    bottom = [ngmesh.Add(MeshPoint(Pnt(r, 0.0, 0.0))) for r in radii]
    top = [ngmesh.Add(MeshPoint(Pnt(r, height, 0.0))) for r in radii]
    for i in range(len(radii) - 1):
        ngmesh.Add(Element2D(1 if i < cells else 2,
                             [bottom[i], bottom[i + 1], top[i + 1], top[i]]))
        ngmesh.Add(Element1D([bottom[i], bottom[i + 1]], index=1))
        ngmesh.Add(Element1D([top[i], top[i + 1]], index=3))
    ngmesh.Add(Element1D([bottom[-1], top[-1]], index=2))
    ngmesh.Add(Element1D([top[0], bottom[0]], index=4))
    return ng.Mesh(ngmesh)


def test_axi_current_drive_is_confined_to_driven_region():
    """With A fixed to zero, J = sigma*Vc/r on the driven ring and 0 elsewhere."""
    r_inner, r_split, r_outer, height = 0.01, 0.015, 0.02, 0.006
    conductivity, omega, imposed_current = 5.8e7, 2.0 * math.pi * 1000.0, 3.25 - 0.4j
    mesh = _split_annular_strip(r_inner, r_split, r_outer, height)
    sigma = ng.CoefficientFunction(conductivity)   # both rings conduct

    with ng.TaskManager():
        solution = solve_axi_eddy(mesh, ng.CoefficientFunction(1.0 / MU0), sigma, omega,
                                  driven_region="driven", total_current=imposed_current,
                                  order=1, dirichlet="bottom|outer|top|inner")
        _, vc = solution.components
        vc_value = complex(vc(mesh(0.0125, height / 2.0)))

    # I = int_0^h int_ri^rs sigma*Vc/r dr dz = sigma*Vc*h*log(rs/ri).
    vc_reference = imposed_current / (conductivity * height * math.log(r_split / r_inner))
    assert abs(vc_value / vc_reference - 1.0) < 1.0e-8
