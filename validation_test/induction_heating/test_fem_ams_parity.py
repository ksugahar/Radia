"""Production SIBC HCurl solver parity; this small lane is not a speed benchmark."""
import json
import sys
from pathlib import Path

import ngsolve as ng
import pytest

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
    for run in (direct, iterative):
        assert run["linear_true_relative_residual"] <= run["linear_true_residual_limit"]
    assert direct["linear_krylov_iterations"] == []
    assert iterative["P_total"] == pytest.approx(direct["P_total"], rel=1e-6)
    assert iterative["L"] == pytest.approx(direct["L"], rel=1e-6)
