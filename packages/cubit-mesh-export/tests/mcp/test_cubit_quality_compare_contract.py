"""Fast contracts for cubit_netgen_quality_compare (no Cubit, no
netgen needed -- error paths and input validation only; the real
three-way comparison lives in validation_test/radia_mcp)."""

import json

import pytest

from cubit_mesh_export.mcp import server


cubit_netgen_quality_compare = server.cubit_netgen_quality_compare


def test_missing_step_is_input_error(tmp_path):
    out = json.loads(cubit_netgen_quality_compare(
        str(tmp_path / "nope.step")))
    assert out["status"] == "error"
    assert out["kind"] == "input"


def test_unknown_scheme_is_input_error(tmp_path):
    step = tmp_path / "x.step"
    step.write_text("dummy", encoding="utf-8")
    out = json.loads(cubit_netgen_quality_compare(
        str(step), schemes=["netgen", "voxels"]))
    assert out["status"] == "error"
    assert out["kind"] == "input"
    assert "voxels" in out["error"]


@pytest.mark.parametrize("kwargs", [
    {"order": "bad"},
    {"order": 4},
    {"netgen_maxh": -1.0, "cubit_size": 1.0},
    {"netgen_maxh": 1.0, "cubit_size": float("nan")},
    {"threshold": float("inf")},
    {"schemes": []},
    {"schemes": ["netgen", "netgen"]},
])
def test_invalid_options_fail_before_meshing(tmp_path, kwargs):
    step = tmp_path / "x.step"
    step.write_text("dummy", encoding="ascii")
    out = json.loads(cubit_netgen_quality_compare(str(step), **kwargs))
    assert out["status"] == "error"
    assert out["kind"] == "input"


def test_cubit_referee_exception_stays_in_route_row(tmp_path, monkeypatch):
    step = tmp_path / "x.step"
    step.write_text("dummy", encoding="ascii")
    monkeypatch.setattr(
        server, "_run_batch",
        lambda *_args, **_kwargs: {"status": "ok", "summary": {}})

    from cubit_mesh_export.mcp._support import mesh_quality as msh_inspect

    def fail_referee(*_args, **_kwargs):
        raise RuntimeError("broken quality file")

    monkeypatch.setattr(msh_inspect, "mesh_quality", fail_referee)
    out = json.loads(cubit_netgen_quality_compare(
        str(step), netgen_maxh=1.0, cubit_size=1.0,
        schemes=["cubit_tet"]))
    assert out["status"] == "ok"
    assert out["all_routes_completed"] is False
    assert out["rows"][0]["status"] == "error"
    assert out["rows"][0]["kind"] == "referee"
    assert "broken quality file" in out["rows"][0]["error"]


@pytest.mark.parametrize('order', [1, 2, 3])
def test_netgen_reference_preserves_curved_volume_without_radia(tmp_path, order):
    pytest.importorskip('gmsh')
    from netgen.occ import OCCGeometry, Sphere, Pnt
    from ngsolve import Mesh, Integrate, TaskManager
    from cubit_mesh_export.mcp._support.netgen_compare import _write_tet_msh
    from cubit_mesh_export.mcp._support.mesh_quality import mesh_total_volume, mesh_quality
    from cubit_mesh_export.mcp.gmsh_v41 import summarize_gmsh_v41_ascii
    target = tmp_path / 'sphere.msh'
    # Netgen's sphere tessellation is not deterministic across hosts. A coarse
    # third-order curved mesh can genuinely invert, which the exporter must
    # reject. Retry with finer source meshes rather than accepting a bad Jacobian.
    for maxh in (0.4, 0.3, 0.2):
        with TaskManager():
            mesh = Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1)).GenerateMesh(maxh=maxh))
        try:
            report = _write_tet_msh(mesh, target, order)
            break
        except ValueError as exc:
            if 'Invalid Netgen reference Jacobian' not in str(exc) or maxh == 0.2:
                raise
    assert report['min_jacobian_det'] > 0
    assert '$Entities' in target.read_text()
    inventory = summarize_gmsh_v41_ascii(target.read_text())
    assert inventory['status'] == 'ok'
    volume = mesh_total_volume(target, quadrature='Gauss10')
    assert volume['ok'], volume
    assert volume['min_jacobian_det'] > 0
    with TaskManager():
        expected = Integrate(1, mesh, order=10)
    assert volume['total_volume'] == pytest.approx(expected, rel=1e-8)
    quality = mesh_quality(target, threshold=0.0)
    assert quality['ran'] and quality['total_negative'] == 0, quality
    assert {row['order'] for row in quality['by_type']} == {order}


@pytest.mark.parametrize('minimum', [-1.0, 0.0, float('nan')])
def test_invalid_reference_jacobians_fail_loudly(tmp_path, monkeypatch, minimum):
    from netgen.occ import OCCGeometry, Sphere, Pnt
    from ngsolve import Mesh
    from cubit_mesh_export.mcp._support import netgen_compare, mesh_quality
    mesh = Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1)).GenerateMesh(maxh=0.8))
    monkeypatch.setattr(mesh_quality, 'mesh_total_volume',
                        lambda *_a, **_k: {'ok': True, 'min_jacobian_det': minimum})
    with pytest.raises(ValueError, match='Jacobian'):
        netgen_compare._write_tet_msh(mesh, tmp_path / 'bad.msh')


def test_reference_nodes_are_cached_by_supported_order(monkeypatch):
    from cubit_mesh_export.mcp._support import netgen_compare
    calls = []
    def reference(_script, args, **kwargs):
        calls.append(args)
        return {'ok': True, 'points': []}
    netgen_compare._reference_nodes.cache_clear()
    monkeypatch.setattr(netgen_compare, 'run_gmsh_json_subprocess', reference)
    try:
        netgen_compare._reference_nodes(1)
        netgen_compare._reference_nodes(1)
        netgen_compare._reference_nodes(2)
        assert len(calls) == 2
        with pytest.raises(ValueError):
            netgen_compare._reference_nodes(4)
    finally:
        netgen_compare._reference_nodes.cache_clear()
