"""Production SIBC HCurl solver parity; this small lane is not a speed benchmark."""
import json
import sys
from pathlib import Path

import ngsolve as ng
import pytest
from radia._residual_gate import RELATIVE_LIMIT

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "radia" / "panels"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import calc_fem_kelvin as solver
import fem_sibc_geometry as geometry
from radia.em_material import EMMaterial


@pytest.mark.parametrize("order", [1, 2, 3])
def test_ams_matches_direct_sibc_loss_and_inductance(tmp_path, order):
    with ng.TaskManager():
        mesh = geometry.air_mesh(
            lambda: geometry.cylinder_part(.025, .025), r_air=.09,
            ring_r=.045, ring_h=.05, maxh_part=.018, maxh_ring=.025,
            maxh_far=.04, order=3)
    vol = tmp_path / "air.vol"
    mesh.ngmesh.Save(str(vol))
    arguments = dict(
        vol_file=str(vol), fes_order=order, frequency=1000,
        mat=EMMaterial(name="copper", sigma=5.8e7, mu_r=1),
        half_thickness=.025, nthreads=1,
        filaments=([geometry.loop_paths((0, 0, 0), (0, 0, 1), .035,
                                        n_seg=48)], [1.]))
    direct = solver.solve_fem(**arguments, solver="sparsecholesky")
    iterative = solver.solve_fem(**arguments, solver="auto")
    record = dict(ngsolve_version=ng.__version__, order=order,
                  direct=direct, iterative=iterative)
    (tmp_path / "result.json").write_text(
        json.dumps(record, default=str, indent=2), encoding="utf-8")
    assert "error" not in direct, direct
    assert "error" not in iterative, iterative
    assert iterative["linear_solver"] == ("ams" if order == 1 else "bddc")
    for result in (direct, iterative):
        assert result["linear_true_residual_limit"] == RELATIVE_LIMIT
        assert result["linear_true_relative_residual"] <= RELATIVE_LIMIT
    assert iterative["P_total"] == pytest.approx(direct["P_total"], rel=1e-6)
    assert iterative["L"] == pytest.approx(direct["L"], rel=1e-6)

    if order == 1:
        scattered = solver.solve_fem(**arguments, solver="sparsecholesky", formulation="scattered")
        assert "error" not in scattered, scattered
        assert scattered["linear_true_relative_residual"] <= RELATIVE_LIMIT
        assert scattered["P_total"] == pytest.approx(direct["P_total"], rel=1e-6)

@pytest.mark.parametrize("order", [1, 2, 3])
def test_periodic_kelvin_auto_uses_the_supported_direct_route(tmp_path, order):
    from netgen.occ import Sphere, Pnt, Vertex, Glue, OCCGeometry, IdentificationType
    workpiece = geometry.cylinder_part(.025, .025)
    workpiece.faces.name = 'sibc'
    interior = Sphere(Pnt(0, 0, 0), .09)
    interior.faces.name = 'kelvin_int'
    air = interior - workpiece
    air.mat('air')
    exterior = Sphere(Pnt(.3, 0, 0), .09)
    exterior.mat('kelvin')
    exterior.faces.name = 'kelvin_ext'
    ground = Vertex(Pnt(.3, 0, 0))
    ground.name = 'GND'
    shape = Glue([air, exterior, ground])
    inside = next(f for f in shape.faces if f.name == 'kelvin_int')
    outside = next(f for f in shape.faces if f.name == 'kelvin_ext')
    inside.Identify(outside, 'kelvin_periodic', IdentificationType.PERIODIC)
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(shape).GenerateMesh(maxh=.02))
        mesh.Curve(max(2, order))
    vol = tmp_path / 'periodic.vol'
    mesh.ngmesh.Save(str(vol))
    args = dict(vol_file=str(vol), fes_order=order, frequency=1000,
                mat=EMMaterial(name='copper', sigma=5.8e7, mu_r=1),
                half_thickness=.025, nthreads=1,
                filaments=([geometry.loop_paths((0,0,0), (0,0,1), .035, n_seg=48)], [1.]))
    result = solver.solve_fem(**args, solver='auto')
    assert 'error' not in result, result
    (tmp_path / 'result.json').write_text(json.dumps(result, default=str, indent=2), encoding='utf-8')
    assert result['kelvin_gradient_gauge_dofs'] > 0
    assert result['linear_solver_requested'] == 'auto'
    assert result['linear_solver'] == 'sparsecholesky'
    assert result['linear_true_relative_residual'] <= RELATIVE_LIMIT
    if order == 1:
        # A closed current must do no work on a discrete scalar gradient.
        # A rule spanning element faces violates this even though curl=0.
        import numpy as np
        with ng.TaskManager():
            h1 = ng.H1(mesh, order=1, definedon=mesh.Materials('air'),
                       dirichlet='kelvin_int|sibc')
            potential = ng.GridFunction(h1)
            free = np.array(list(h1.FreeDofs()), dtype=bool)
            potential.vec.FV().NumPy()[free] = np.random.default_rng(17).normal(size=int(free.sum()))
            gradient = ng.grad(potential)
            nodes, weights = np.polynomial.legendre.leggauss(4)
            work = 0.0
            for start, end in args['filaments'][0][0]:
                start = np.asarray(start)
                direction = np.asarray(end)-start
                for point, weight in solver._element_line_quadrature(
                        mesh, start, direction, .5*(nodes+1), .5*weights):
                    work += weight * np.dot(gradient(point), direction)
            assert abs(work) < 1e-8
    rejected = solver.solve_fem(**args, solver='bddc')
    assert 'not validated' in rejected['error']
