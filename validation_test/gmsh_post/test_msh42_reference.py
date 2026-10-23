"""MSH 4.2 serialization against Gmsh itself, independently prescribed data.

Load the Python-only writer directly so this format validation requires no
Radia native solver build. Gmsh is a validation dependency, not a core runtime.
"""
import importlib.util
import ast
import math
from pathlib import Path
import sys
import types

import gmsh
import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[2]


def _source_module(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


package = types.ModuleType('_msh42_core')
package.__path__ = [str(ROOT / 'src/radia')]
sys.modules[package.__name__] = package
writer = _source_module('_msh42_core.gmsh_post_export', 'src/radia/gmsh_post_export.py')
reader = _source_module('_msh42_core.peec_mesh_import', 'src/radia/peec_mesh_import.py')


def _views(path):
    reader.validate_native_msh_fields(path)
    gmsh.initialize(['-noconfig'])
    try:
        gmsh.option.setNumber('General.Terminal', 0)
        gmsh.open(str(path))
        result = {}
        for tag in gmsh.view.getTags():
            name = gmsh.option.getString(f'View[{gmsh.view.getIndex(tag)}].Name')
            count = int(gmsh.option.getNumber(f'View[{gmsh.view.getIndex(tag)}].NbTimeStep'))
            for step in range(count):
                kind, tags, data, time, components = gmsh.view.getModelData(tag, step)
                result[name, step] = (kind, dict(zip(map(int, tags), map(list, data))), time, components)
        gmsh.option.setNumber('Mesh.MshFileVersion', 4.2)
        gmsh.option.setNumber('Mesh.Binary', 0)
        target = path.with_name(path.stem + '_roundtrip.msh')
        gmsh.write(str(target))
        gmsh.option.setNumber('PostProcessing.SaveMesh', 0)
        for tag in gmsh.view.getTags():
            gmsh.view.write(tag, str(target), append=True)
        return result, target
    finally:
        gmsh.finalize()


def test_payload_fields_roundtrip(tmp_path):
    start = 17
    expected = {start: [1., -2., 3.], start + 1: [4., 5., -6.]}
    payload = writer.PayloadV42(
        groups=[writer.PhysicalGroup(1, 1, 'conductor')],
        entities=[writer.Entity(1, 1, (0, 0, 0, 1, 0, 0), [1])],
        node_blocks=[writer.NodeBlock(1, 1, start, [(0, 0, 0), (1, 0, 0)])],
        element_blocks=[writer.ElementBlock(1, 1, 1, [[start, start + 1]])],
        data_fields=[writer.DataField('$NodeData', 'B', 3, list(expected.items()), time=.25),
                     writer.DataField('$NodeData', 'B', 3, [(t, [2*x for x in v]) for t, v in expected.items()], time=.5, timestep=1),
                     writer.DataField('$ElementData', 'I', 1, [(1, -7.)])])
    path = tmp_path / 'payload.msh'
    with path.open('w', encoding='utf-8') as stream:
        writer._emit_v42_msh(stream, payload)
    actual, target = _views(path)
    assert actual['B', 0] == ('NodeData', expected, .25, 3)
    assert actual['B', 1][1] == {t: [2*x for x in v] for t, v in expected.items()}
    assert actual['I', 0] == ('ElementData', {1: [-7.]}, 0., 1)
    reopened, _ = _views(target)
    assert reopened == actual
    nodes, edges = reader.read_gmsh_centerline(path)
    assert nodes == {start: [0., 0., 0.], start+1: [1., 0., 0.]}
    assert edges == [[start, start+1]]


def test_large_tags_remain_supported_for_mesh_only_and_pure_python(tmp_path):
    start = 2147483665
    payload = writer.PayloadV42(
        entities=[writer.Entity(1, 1, (0, 0, 0, 1, 0, 0))],
        node_blocks=[writer.NodeBlock(1, 1, start, [(0, 0, 0), (1, 0, 0)])],
        element_blocks=[writer.ElementBlock(1, 1, 1, [[start, start + 1]])])
    path = tmp_path / 'large_mesh.msh'
    with path.open('w', encoding='utf-8') as stream:
        writer._emit_v42_msh(stream, payload)
    reader.validate_native_msh_fields(path)
    gmsh.initialize(['-noconfig'])
    try:
        gmsh.open(str(path))
        assert list(map(int, gmsh.model.mesh.getNodes()[0])) == [start, start + 1]
        assert list(map(int, gmsh.model.mesh.getElementsByType(1)[1])) == [start, start + 1]
    finally:
        gmsh.finalize()
    # Pure serialization/parsing retains size_t indices even with field data.
    # Do not pass this field-bearing file to the native reference.
    with path.open('a', encoding='utf-8') as stream:
        writer._emit_data_section(stream, writer.DataField(
            '$NodeData', 'B', 3, [(start, [1., -2., 3.]), (start + 1, [4., 5., -6.])]))
    nodes, edges = reader.read_gmsh_centerline(path)
    assert nodes == {start: [0., 0., 0.], start + 1: [1., 0., 0.]}
    assert edges == [[start, start + 1]]


def test_separate_line_writer_preserves_fields(tmp_path):
    path = tmp_path / 'lines.msh'
    stats = writer._write_line_element_msh(
        str(path), [((0, 0, 0), (1, 0, 0)), ((1, 0, 0), (2, 0, 0))],
        [0, 0], ['coil'], [{'name': 'I', 'ncomp': 1, 'values': [2., -3.]}])
    actual, target = _views(path)
    assert actual['I', 0][1] == dict(zip(stats['element_tags'], [[2.], [-3.]]))
    assert _views(target)[0] == actual


def test_streamfunction_wire_writer_opens_in_gmsh(tmp_path):
    # Execute the production function without importing unrelated solver kernels.
    tree = ast.parse((ROOT / 'src/radia/streamfunction_volume.py').read_text(encoding='utf-8'))
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'write_wires_msh')
    namespace = {'math': math}
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(ROOT / 'src/radia/streamfunction_volume.py'), 'exec'), namespace)
    path = tmp_path / 'wires.msh'
    assert namespace['write_wires_msh']([[(0, 0, 0), (1, 0, 0), (1, 1, 0)]], path) == (3, 2)
    gmsh.initialize(['-noconfig'])
    try:
        gmsh.open(str(path))
        assert list(map(int, gmsh.model.mesh.getNodes()[0])) == [1, 2, 3]
        assert list(map(int, gmsh.model.mesh.getElementsByType(1)[1])) == [1, 2, 2, 3]
    finally:
        gmsh.finalize()


def test_combined_surface_volume_writer_roundtrips_fields(tmp_path):
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    mesh = MakeStructured3DMesh(hexes=False, nx=1, ny=1, nz=1)
    vol = tmp_path / 'combined_source.vol'
    mesh.ngmesh.Save(str(vol))
    checker = _source_module('msh42_combined_vol_checker', 'packages/cubit-mesh-export/src/cubit_mesh_export/check.py')
    contract = {'schema': checker.LABEL_CONTRACT_SCHEMA,
                'materials': {'required': list(mesh.GetMaterials())},
                'boundaries': {'required': list(mesh.GetBoundaries())}}
    assert checker.check_label_contract(vol, contract=contract)['passed']
    mesh = ng.Mesh(str(vol))
    coordinates = [list(v.point) for v in mesh.vertices]
    elements = [[v.nr for v in e.vertices] for e in mesh.Elements(ng.VOL)]
    path = tmp_path / 'combined.msh'
    with ng.TaskManager():
        writer.write_combined_msh_v42(str(path), mesh, [[1, 2, 3]]*mesh.nv,
                                     coordinates, elements, [[4, 5, 6]]*mesh.nv)
    actual, target = _views(path)
    assert set(actual['B', 0][1]) == set(range(1, mesh.nv+1))
    assert all(v == [4., 5., 6.] for v in actual['B', 0][1].values())
    assert all(v == [1., 2., 3.] for v in actual['J', 0][1].values())
    assert _views(target)[0] == actual


def test_checked_vol_mesh_field_export_and_animation(tmp_path):
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    # Persist/reload the solver mesh; MSH is a derived post-processing artifact.
    mesh = MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1)
    vol = tmp_path / 'cube.vol'
    mesh.ngmesh.Save(str(vol))
    checker = _source_module('msh42_vol_checker', 'packages/cubit-mesh-export/src/cubit_mesh_export/check.py')
    contract = {'schema': checker.LABEL_CONTRACT_SCHEMA,
                'materials': {'required': list(mesh.GetMaterials())},
                'boundaries': {'required': list(mesh.GetBoundaries())}}
    assert checker.check_label_contract(vol, contract=contract)['passed']
    mesh = ng.Mesh(str(vol))
    assert mesh.GetMaterials() and mesh.ne == 1
    with ng.TaskManager():
        post = writer.GmshPostExport(mesh)
        post.add_vector_field('B', ng.CF((ng.x, 2*ng.y, -ng.z)))
        path = tmp_path / 'volume.msh'
        post.write(str(path))
        actual, target = _views(path)
        coords = {v.nr+1: list(v.point) for v in mesh.vertices}
        expected = {t: [x, 2*y, -z] for t, (x, y, z) in coords.items()}
        assert actual['B', 0][1] == expected
        assert _views(target)[0] == actual
        active = tmp_path / 'activation.msh'
        writer.export_element_activation_animation(mesh, [[True]], str(active))
        assert _views(active)[0]
        moving = tmp_path / 'moving.msh'
        displacement = np.zeros((mesh.nv, 3))
        writer.export_nodal_deformation_animation(mesh, [displacement], str(moving))
        assert all(v == [0., 0., 0.] for v in next(iter(_views(moving)[0].values()))[1].values())


@pytest.mark.parametrize('header', ['4.1 0 8', '4.0 0 8', '2.2 0 8', '4.9 0 8', '4.2 1 8', '4.2'])
def test_core_readers_reject_other_formats(tmp_path, header):
    path = tmp_path / 'unsupported.msh'
    path.write_text(f'$MeshFormat\n{header}\n$EndMeshFormat\n', encoding='utf-8')
    with pytest.raises(ValueError):
        reader.read_gmsh_centerline(path)
    with pytest.raises(ValueError):
        reader._load_gmsh(path)
