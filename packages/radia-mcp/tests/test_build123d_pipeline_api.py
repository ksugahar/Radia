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
