"""Validate explicit source-unit conversion against an analytic SI cube.

The one-meter cube has coordinates +/-0.5 m and volume 1 m^3.
"""
from pathlib import Path
import sys
import pytest


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
