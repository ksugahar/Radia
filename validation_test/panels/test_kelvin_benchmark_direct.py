"""Analytical magnetic sphere through the production Kelvin benchmark entry."""
import pytest


def test_kelvin_sphere_sparsecholesky(tmp_path):
    import ngsolve as ng
    from netgen.occ import Sphere, Pnt, Vertex, Glue, OCCGeometry, IdentificationType
    from radia.panels.calc_kelvin_benchmark import solve_kelvin_benchmark

    iron = Sphere(Pnt(0, 0, 0), .05)
    iron.mat("magnetic")
    iron.maxh = .015
    iron.faces.name = "sphere"
    inner = Sphere(Pnt(0, 0, 0), .2)
    inner.faces.name = "kelvin_int"
    air = inner - iron
    air.mat("air")
    outer = Sphere(Pnt(0, 0, .6), .2)
    outer.mat("kelvin")
    outer.faces.name = "kelvin_ext"
    ground = Vertex(Pnt(0, 0, .6))
    ground.name = "GND"
    shape = Glue([air, iron, outer, ground])
    inside = [f for f in shape.faces if f.name == "kelvin_int"][0]
    outside = [f for f in shape.faces if f.name == "kelvin_ext"][0]
    inside.Identify(outside, "kelvin_periodic", IdentificationType.PERIODIC)
    with ng.TaskManager():
        mesh = OCCGeometry(shape).GenerateMesh(maxh=.035)
    path = tmp_path / "sphere.vol"
    mesh.Save(str(path))
    result = solve_kelvin_benchmark(str(path), fes_order=2)
    assert result["linear_solver"] == "sparsecholesky"
    assert result["linear_relative_residual"] < 1e-10
    assert result["slaved_dofs"] > 0
    assert result["converged"]
    assert result["Hi_origin"] == pytest.approx(3 / 102, rel=.01)
