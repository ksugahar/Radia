"""Independent Gmsh-generated field blocks, mesh embeddings and reader boundary."""
from pathlib import Path

import gmsh
import pytest

from radia_mcp.gmsh.msh_inspect import inspect_msh, read_msh_data, validate_msh
from radia_mcp.gmsh.runtime import qualify_runtime, require_msh42


def test_native_api_generated_42_fields_and_embeddings(tmp_path):
    identity = qualify_runtime()
    assert identity['wheel_sha256'] and identity['msh_version'] == 4.2
    path = tmp_path / 'independent.msh'
    expected = {
        'B': ('NodeData', {17: [1., 2., -3.], 19: [-4., 5., 6.]}),
        'I': ('ElementData', {37: [-7.], 41: [9.]}),
        'T': ('ElementNodeData', {37: [1., 2.], 41: [3., 4.]}),
    }
    gmsh.initialize(['-noconfig'])
    try:
        gmsh.model.add('oracle')
        gmsh.model.addDiscreteEntity(0, 3)
        gmsh.model.addDiscreteEntity(1, 1)
        gmsh.model.mesh.embed(0, [3], 1, 1)
        gmsh.model.mesh.addNodes(0, 3, [23], [.5, 0, 0])
        gmsh.model.mesh.addNodes(1, 1, [17, 19, 21], [0, 0, 0, 1, 0, 0, 2, 0, 0])
        gmsh.model.mesh.addElementsByType(1, 1, [37, 41], [17, 19, 19, 21])
        gmsh.write(str(path))
        gmsh.option.setNumber('PostProcessing.SaveMesh', 0)
        for name, (kind, values) in expected.items():
            view = gmsh.view.add(name)
            gmsh.view.addModelData(view, 0, 'oracle', kind,
                                   list(values), list(values.values()),
                                   time=.125, numComponents=3 if name == 'B' else 1)
            gmsh.view.write(view, str(path), append=True)
    finally:
        gmsh.finalize()
    actual = read_msh_data(path, include_elements=True)
    assert actual['elements'] == {37: {'type': 1, 'nodes': [17, 19]},
                                  41: {'type': 1, 'nodes': [19, 21]}}
    for view in actual['views']:
        kind, values = expected.pop(view['name'])
        assert view['section'] == kind and view['rows'] == values
        assert view['time'] == .125
    assert not expected
    assert inspect_msh(path)['ok']
    result = validate_msh(path)
    assert result['ok'], result['errors']


@pytest.mark.parametrize('version', ['4.1', '4.0', '2.2', '4.9'])
def test_other_versions_fail_both_readers_and_native_preflight(tmp_path, version):
    path = tmp_path / 'other.msh'
    path.write_text(f'$MeshFormat\n{version} 0 8\n$EndMeshFormat\n', encoding='utf-8')
    assert not inspect_msh(path)['ok']
    assert not validate_msh(path)['ok']
    with pytest.raises(ValueError):
        read_msh_data(path)
    with pytest.raises(ValueError):
        require_msh42(path)


def test_time_series_writer_is_readable_by_native_reference(tmp_path):
    from radia_mcp.gmsh.post_process import time_series

    qualify_runtime()
    inputs = []
    for step, values in enumerate(([1., 6.], [5., 2.])):
        path = tmp_path / f'step{step}.msh'
        gmsh.initialize(['-noconfig'])
        try:
            gmsh.model.add('series')
            gmsh.model.addDiscreteEntity(1, 1)
            gmsh.model.mesh.addNodes(1, 1, [17, 19], [0, 0, 0, 1, 0, 0])
            gmsh.model.mesh.addElementsByType(1, 1, [37], [17, 19])
            gmsh.write(str(path))
            view = gmsh.view.add('f')
            gmsh.view.addModelData(view, 0, 'series', 'NodeData',
                                   [17, 19], [[v] for v in values],
                                   numComponents=1)
            gmsh.option.setNumber('PostProcessing.SaveMesh', 0)
            gmsh.view.write(view, str(path), append=True)
        finally:
            gmsh.finalize()
        inputs.append(path)
    output = tmp_path / 'stats.msh'
    result = time_series(inputs, view='f', times=[0., 1.],
                         stats=('min', 'max', 'argmax_time'), out_file=output)
    assert result['ok'], result
    expected = {'f_min': {17: [1.], 19: [2.]},
                'f_max': {17: [5.], 19: [6.]},
                'f_argmax_time': {17: [1.], 19: [0.]}}
    gmsh.initialize(['-noconfig'])
    try:
        gmsh.open(str(output))
        for view in gmsh.view.getTags():
            name = gmsh.view.option.getString(view, 'Name')
            kind, tags, values, time, ncomp = gmsh.view.getModelData(view, 0)
            assert kind == 'NodeData' and ncomp == 1 and time == 0.
            assert dict(zip(map(int, tags), map(list, values))) == expected.pop(name)
        assert not expected
    finally:
        gmsh.finalize()
