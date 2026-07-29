from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
for source_root in (ROOT / "src", ROOT / "packages" / "radia-mcp" / "src"):
    source_path = str(source_root)
    if source_path not in sys.path:
        sys.path.insert(0, source_path)

ng = pytest.importorskip("ngsolve")
pytest.importorskip("radia.axifem")

from netgen.meshing import (  # noqa: E402
    Element1D,
    Element2D,
    FaceDescriptor,
    Mesh as NetgenMesh,
    MeshPoint,
    Pnt,
)
from radia.axifem import H1Henrotte  # noqa: E402
from radia_mcp.radia_ngsolve.solve import (  # noqa: E402
    NU0,
    axisymmetric_point_potential_constraint_contract,
    axisymmetric_ring_current_load_contract,
    solve_axi_magnetostatic,
    solve_axi_magnetostatic_nonlinear,
)


def _structured_mesh():
    mesh = NetgenMesh(dim=2)
    mesh.SetMaterial(1, "domain")
    for bc in range(1, 5):
        mesh.Add(FaceDescriptor(surfnr=bc, domin=0, bc=bc))
    mesh.SetBCName(0, "axis")
    mesh.SetBCName(1, "outer")
    mesh.SetBCName(2, "outer")
    mesh.SetBCName(3, "outer")

    radii = (0.0, 1.0, 2.0, 3.0)
    axial = (-1.0, 0.0, 1.0)
    points = {}
    for j, z_value in enumerate(axial):
        for i, radius in enumerate(radii):
            points[i, j] = mesh.Add(MeshPoint(Pnt(radius, z_value, 0.0)))

    for j in range(len(axial) - 1):
        for i in range(len(radii) - 1):
            mesh.Add(
                Element2D(
                    1,
                    [
                        points[i, j],
                        points[i + 1, j],
                        points[i + 1, j + 1],
                        points[i, j + 1],
                    ],
                )
            )
    for j in range(len(axial) - 1):
        mesh.Add(Element1D([points[0, j], points[0, j + 1]], index=1))
        mesh.Add(Element1D([points[3, j], points[3, j + 1]], index=2))
    for i in range(len(radii) - 1):
        mesh.Add(Element1D([points[i, 2], points[i + 1, 2]], index=3))
        mesh.Add(Element1D([points[i, 0], points[i + 1, 0]], index=4))
    return ng.Mesh(mesh)


def test_ring_contract_accepts_numpy_rows_and_rejects_malformed_rows():
    mesh = _structured_mesh()
    rows = np.array([[1.0, 0.0, 2.5], [0.0, 0.0, 7.0]])

    evidence = axisymmetric_ring_current_load_contract(mesh, rows)

    assert evidence[0]["weak_load_amplitude_a_m"] == pytest.approx(2.5)
    assert evidence[1]["axis_annihilated"] is True
    with pytest.raises(ValueError, match="exactly three numeric values"):
        axisymmetric_ring_current_load_contract(mesh, [[1.0, 0.0]])


def test_point_constraint_rejects_nonzero_existing_dirichlet_dof():
    mesh = _structured_mesh()
    fes = H1Henrotte(mesh, order=2, dirichlet="axis|outer")

    interior = axisymmetric_point_potential_constraint_contract(
        fes, mesh, np.array([[2.0, 0.0, 3.0e-7]])
    )
    assert interior[0]["already_dirichlet"] is False

    compatible = axisymmetric_point_potential_constraint_contract(
        fes, mesh, [[3.0, 0.0, 0.0]]
    )
    assert compatible[0]["already_dirichlet"] is True
    with pytest.raises(ValueError, match="existing Dirichlet DOF"):
        axisymmetric_point_potential_constraint_contract(
            fes, mesh, [[3.0, 0.0, 1.0e-7]]
        )


def test_linear_and_nonlinear_exact_sources_agree_for_constant_reluctivity():
    mesh = _structured_mesh()
    rings = [[1.0, 0.0, 4.0]]
    potentials = [[2.0, 0.0, 2.0e-7]]
    magnets = {"domain": (2.0, 90.0)}

    linear = solve_axi_magnetostatic(
        mesh,
        ng.CoefficientFunction(NU0),
        magnets=magnets,
        order=2,
        ring_currents=np.asarray(rings),
        point_potentials=np.asarray(potentials),
    )
    nonlinear = solve_axi_magnetostatic_nonlinear(
        mesh,
        lambda _field: ng.CoefficientFunction(NU0),
        None,
        2,
        "axis|outer",
        1.0,
        3,
        1.0e-12,
        2,
        magnets=magnets,
        ring_currents=(row for row in rings),
        point_potentials=(row for row in potentials),
    )

    point = axisymmetric_point_potential_constraint_contract(
        nonlinear.space, mesh, potentials
    )[0]
    assert nonlinear.vec[point["dof_index"]] == pytest.approx(2.0e-7, abs=1.0e-15)
    np.testing.assert_allclose(
        nonlinear.vec.FV().NumPy(),
        linear.vec.FV().NumPy(),
        rtol=2.0e-12,
        atol=2.0e-15,
    )


def test_order_one_point_constraints_fail_before_space_construction():
    mesh = _structured_mesh()
    with pytest.raises(NotImplementedError, match="order-1"):
        solve_axi_magnetostatic(
            mesh,
            ng.CoefficientFunction(NU0),
            order=1,
            point_potentials=np.array([[1.0, 0.0, 0.0]]),
        )
