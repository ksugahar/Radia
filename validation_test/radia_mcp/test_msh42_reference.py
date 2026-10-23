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
