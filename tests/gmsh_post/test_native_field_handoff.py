"""Fast native-boundary checks; no Gmsh import or large native allocation."""
import importlib.util
import hashlib
import json
from pathlib import Path
import sys
import types

import pytest


ROOT = Path(__file__).resolve().parents[2]
package = types.ModuleType('_native_field_source')
package.__path__ = [str(ROOT / 'src/radia')]
sys.modules[package.__name__] = package
spec = importlib.util.spec_from_file_location(
    package.__name__ + '._gmsh_native_fields', ROOT / 'src/radia/_gmsh_native_fields.py')
reader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reader)


def _field_file(path, kind, tag, declared_max=None):
    block = '1 2 1' if kind == 'ElementNodeData' else '1'
    values = '1.25 -2.5' if kind == 'ElementNodeData' else '1.25'
    maximum = tag if declared_max is None else declared_max
    path.write_text(
        '$MeshFormat\n4.2 0 8\n$EndMeshFormat\n'
        f'${kind}\n1\n"field"\n1\n0\n3\n0\n1\n1\n'
        f'1 1 {tag} {maximum}\n{block}\n{tag}\n{values}\n$End{kind}\n',
        encoding='utf-8')


@pytest.mark.parametrize('kind', ['NodeData', 'ElementData', 'ElementNodeData'])
def test_large_field_tag_rejected_before_native_access(tmp_path, monkeypatch, kind):
    path = tmp_path / 'field.msh'
    _field_file(path, kind, 2147483665)
    native_calls = []

    class NativeMustNotRun:
        def initialize(self):
            native_calls.append('initialize')
            raise AssertionError('native Gmsh was accessed before tag rejection')

    monkeypatch.setitem(sys.modules, 'gmsh', NativeMustNotRun())
    with pytest.raises(ValueError, match='2147483665 exceeds 2147483647'):
        reader.validate_native_msh_fields(path)
    assert native_calls == []


@pytest.mark.parametrize('tag', [268435456, 2147483647])
@pytest.mark.parametrize('kind', ['NodeData', 'ElementData', 'ElementNodeData'])
def test_supported_field_tag_preflight_is_small_and_allocation_free(tmp_path, kind, tag):
    path = tmp_path / 'field.msh'
    _field_file(path, kind, tag)
    reader.validate_native_msh_fields(path)


def test_large_actual_tag_cannot_hide_behind_small_header(tmp_path):
    path = tmp_path / 'field.msh'
    _field_file(path, 'NodeData', 2147483665, declared_max=17)
    with pytest.raises(ValueError, match='2147483665 exceeds 2147483647'):
        reader.validate_native_msh_fields(path)


def test_geo_native_handoff_checks_merged_fields(tmp_path):
    _field_file(tmp_path / 'field.msh', 'NodeData', 2147483665)
    geo = tmp_path / 'display.geo'
    geo.write_text('Merge "field.msh";\n', encoding='utf-8')
    with pytest.raises(ValueError, match='2147483665 exceeds 2147483647'):
        reader.validate_native_msh_fields(geo)


@pytest.mark.parametrize('name', ['point_vector.msh', 'line_vector.msh'])
def test_native_vector_element_node_fixture_passes_core_preflight(name):
    """Guards: component-aware field consumption; assumes native oracle fixtures."""
    fixtures = Path(__file__).parent / 'fixtures/native_element_node_data/inputs'
    oracle = json.loads((fixtures / 'native-oracle.json').read_text(encoding='utf-8'))
    expected = next(f for f in oracle['fixtures'] if f['file'] == name)
    path = fixtures / name
    assert hashlib.sha256(path.read_bytes()).hexdigest() == expected['sha256']
    assert expected['components'] == 3
    assert all(len(row) == 3 * expected['samples_per_element']
               for row in expected['values'])
    reader.validate_native_msh_fields(path)


@pytest.mark.parametrize('scheme,samples', [('', 0), ('', 1), ('', 3), ('', 5), ('named', 0)])
def test_invalid_element_node_samples_rejected_before_native(tmp_path, monkeypatch, scheme, samples):
    """Guards: sample shape before native access; assumes a two-node line."""
    path = tmp_path / 'field.msh'
    path.write_text('$MeshFormat\n4.2 0 8\n$EndMeshFormat\n$ElementNodeData\n'
                    + f'2\n"field"\n"{scheme}"\n1\n0\n3\n0\n1\n1\n'
                    + f'1 1 17 17\n1 {samples} 1\n17\n' + '1.25 ' * samples
                    + '\n$EndElementNodeData\n', encoding='utf-8')
    native_calls = []
    monkeypatch.setitem(sys.modules, 'gmsh', types.SimpleNamespace(
        initialize=lambda: native_calls.append('initialize')))
    with pytest.raises(ValueError, match='sample count/type mismatch'):
        reader.validate_native_msh_fields(path)
    assert native_calls == []


def test_named_interpolation_allows_positive_sample_count(tmp_path):
    """Guards: named interpolation differs from geometry; no evaluation implied."""
    path = tmp_path / 'field.msh'
    path.write_text('$MeshFormat\n4.2 0 8\n$EndMeshFormat\n$ElementNodeData\n'
                    '2\n"field"\n"named"\n1\n0\n3\n0\n3\n1\n'
                    '1 1 17 17\n1 1 1\n17\n1.25 -2.5 3.75\n$EndElementNodeData\n',
                    encoding='utf-8')
    reader.validate_native_msh_fields(path)


@pytest.mark.parametrize('etype', [34, 999])
@pytest.mark.parametrize('scheme', ['', 'named'])
def test_unknown_element_type_rejected_before_native(tmp_path, etype, scheme):
    """Guards: explicit unknown-type rejection; assumes no native type lookup."""
    path = tmp_path / 'unknown.msh'
    path.write_text('$MeshFormat\n4.2 0 8\n$EndMeshFormat\n$ElementNodeData\n'
                    + f'2\n"field"\n"{scheme}"\n1\n0\n3\n0\n1\n1\n'
                    + f'1 1 17 17\n{etype} 2 1\n17\n1.25 -2.5\n'
                    + '$EndElementNodeData\n', encoding='utf-8')
    with pytest.raises(ValueError, match='sample count/type mismatch'):
        reader.validate_native_msh_fields(path)
