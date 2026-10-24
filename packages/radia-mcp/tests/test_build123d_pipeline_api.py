"""Guards: Netgen .vol is the analysis boundary; MSH 4.2 is post output.
Assumes: the stored result fixture has independently recorded provenance.
"""
from pathlib import Path
import sys

import pytest


def test_build123d_pipeline_public_api_imports_without_geometry_runtime():
    from radia_mcp.build123d import pipeline

    assert pipeline.__all__ == ["run_pipeline", "run_pipeline_multi", "save_record"]
    for name in pipeline.__all__:
        assert callable(getattr(pipeline, name))


@pytest.mark.parametrize('entry', ['_load_checked_vol', '_stage_post', '_stage_post_multi'])
def test_pipeline_rejects_msh_analysis_before_loading_native_dependencies(entry, tmp_path, monkeypatch):
    from radia_mcp.build123d import pipeline

    monkeypatch.setitem(sys.modules, 'ngsolve', None)
    monkeypatch.setitem(sys.modules, 'cubit_mesh_export', None)
    args = [tmp_path / 'display.msh']
    if entry == '_stage_post':
        args += [tmp_path, 'display']
    elif entry == '_stage_post_multi':
        args += [['core'], tmp_path, 'display']
    with pytest.raises(ValueError, match='Netgen .vol'):
        getattr(pipeline, entry)(*args)


def test_pipeline_rejects_non_si_vol_before_loading_dependencies(tmp_path, monkeypatch):
    import json
    from radia_mcp.build123d import pipeline

    path = tmp_path / 'solver.vol'
    path.with_suffix('.labels.json').write_text(json.dumps({'units': 'mm'}), encoding='utf-8')
    monkeypatch.setitem(sys.modules, 'ngsolve', None)
    with pytest.raises(ValueError, match='SI units'):
        pipeline._load_checked_vol(path)


def test_pipeline_inspects_stored_msh_only_as_post_output():
    from radia_mcp.build123d import pipeline

    path = Path(__file__).parent / 'fixtures/gmsh/unit_cube_tet.msh'
    result = pipeline._validate_post_output(path)
    assert result['version'] == '4.2'
    assert result['nodes']['count'] == 8
    assert result['elements']['count'] == 6


def test_pipeline_requires_current_exporter_before_creating_msh(monkeypatch):
    from types import ModuleType
    from radia_mcp.build123d import pipeline

    package = ModuleType('radia')
    package.__path__ = []
    backend = ModuleType('radia.gmsh_post_export')

    def forbidden(*args):
        raise AssertionError('unsupported writer was called')

    backend.GmshPostExport = forbidden
    monkeypatch.setitem(sys.modules, 'radia', package)
    monkeypatch.setitem(sys.modules, 'radia.gmsh_post_export', backend)
    with pytest.raises(RuntimeError, match='current Radia MSH 4.2 exporter'):
        pipeline._post_exporter(object())


@pytest.mark.parametrize('entry', ['run_pipeline', 'run_pipeline_multi'])
def test_pipeline_units_are_required_before_creating_artifacts(entry, tmp_path):
    from radia_mcp.build123d import pipeline

    destination = tmp_path / 'result'
    with pytest.raises(TypeError, match='units'):
        getattr(pipeline, entry)(None, destination, 'case')
    with pytest.raises(ValueError, match="explicitly declared"):
        getattr(pipeline, entry)(None, destination, 'case', units='unspecified')
    assert not destination.exists()


def test_native_solver_mesh_save_requires_explicit_si_before_dependencies(tmp_path, monkeypatch):
    from radia_mcp.build123d import pipeline

    monkeypatch.setitem(sys.modules, 'ngsolve', None)
    with pytest.raises(ValueError, match='already use SI'):
        pipeline._save_solver_mesh(None, tmp_path, 'case', units='mm')


@pytest.mark.parametrize('units,side', [('m', 1.0), ('mm', 1000.0)])
@pytest.mark.parametrize('route', ['brep', 'step'])
def test_explicit_source_units_produce_the_same_si_cube(units, side, route, tmp_path):
    """A prescribed one-meter cube survives both CAD routes and source units."""
    import json
    import math

    build123d = pytest.importorskip('build123d')
    ng = pytest.importorskip('ngsolve')
    pytest.importorskip('cubit_mesh_export.check')
    from radia_mcp.build123d import pipeline

    part = build123d.Box(side, side, side).solid()
    if route == 'step':
        cad = pipeline._stage_cad_multi([(part, 'domain')], tmp_path, 'unit_cube', units=units)
        result = pipeline._stage_mesh_multi(
            Path(cad['step']), ['domain'], tmp_path, 'unit_cube',
            maxh=side / 2.0, units=units,
        )
    else:
        cad = pipeline._stage_cad(part, tmp_path, 'unit_cube')
        result = pipeline._stage_mesh(
            Path(cad['brep']), tmp_path, 'unit_cube',
            maxh=side / 2.0, units=units,
        )
    mesh = ng.Mesh(result['vol'])
    coordinates = [vertex.point for vertex in mesh.vertices]
    bound = 64.0 * sys.float_info.epsilon
    for axis in range(3):
        assert math.isclose(min(p[axis] for p in coordinates), -0.5,
                            rel_tol=0.0, abs_tol=bound)
        assert math.isclose(max(p[axis] for p in coordinates), 0.5,
                            rel_tol=0.0, abs_tol=bound)
    with ng.TaskManager():
        volume = float(ng.Integrate(ng.CoefficientFunction(1.0), mesh))
    assert math.isclose(volume, 1.0, rel_tol=1e-9, abs_tol=bound), volume
    assert result['input_units'] == units
    assert result['maxh_m'] == 0.5
    assert json.loads(Path(result['label_contract']).read_text())['units'] == 'm'
    assert b'\r' not in Path(result['vol']).read_bytes()
