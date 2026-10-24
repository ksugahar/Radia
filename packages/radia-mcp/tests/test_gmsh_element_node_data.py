"""Guards: native MSH 4.2 sample counts, component ordering and preflight.

Assumes: committed post-only fixtures and getModelData oracle records generated
by the qualified Gmsh API; pure parsing needs no native Gmsh installation.
"""
import hashlib
import io
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from radia_mcp.gmsh import msh_inspect as msh, runtime

ROOT = Path(__file__).parent / 'fixtures/gmsh/element_node_data'
PROVENANCE = json.loads((ROOT / 'provenance.json').read_text(encoding='utf-8'))


@pytest.mark.parametrize('fixture', PROVENANCE['fixtures'], ids=lambda f: f['file'])
def test_native_element_node_data_roundtrip_matches_oracle(fixture):
    path = ROOT / fixture['file']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == fixture['sha256']
    runtime.require_msh42(path)
    checked = msh.validate_msh(path)
    assert checked['ok'], checked['errors']
    parsed = msh.read_msh_data(path)
    view, = parsed['views']
    assert (view['section'], view['components'], view['time']) == (
        fixture['data_type'], fixture['components'], fixture['time'])
    expected = dict(zip(fixture['tags'], fixture['values']))
    assert view['rows'] == expected, (view['rows'], expected)
    assert all(len(row) == fixture['samples_per_element'] * fixture['components']
               for row in view['rows'].values())


def test_native_preflight_discards_samples_but_validates_every_value():
    # Interpolation sample counts need not match geometric connectivity.
    # A large record must not be retained by the native eligibility iterator.
    samples = 10000
    prefix = f'1 1 17 17\n1 {samples} 1\n17\n'
    text = prefix + '1.25 -2.5 3.75 ' * samples
    rows = list(msh._iter_data_blocks('ElementNodeData', io.StringIO(text),
                lambda: io.StringIO(text), 3, 1, native_limit=2**31-1,
                collect_values=False, scheme_present=True))
    assert rows == [(17, samples, [])]
    # Replace the last numeric value inside the declared record.
    broken = text.rsplit('3.75 ', 1)[0] + 'not-a-number'
    with pytest.raises(ValueError, match='could not convert'):
        list(msh._iter_data_blocks('ElementNodeData', io.StringIO(broken),
             lambda: io.StringIO(broken), 3, 1, native_limit=2**31-1,
             collect_values=False, scheme_present=True))


@pytest.mark.parametrize('tag', [2**28, 2147483665])
def test_vector_element_node_tags_checked_before_native(tmp_path, tag):
    path = tmp_path / 'field.msh'
    path.write_text('$MeshFormat\n4.2 0 8\n$EndMeshFormat\n$ElementNodeData\n'
                    f'1\n"v"\n1\n0\n3\n0\n3\n1\n1 1 {tag} {tag}\n'
                    f'15 1 1\n{tag}\n1.25 -2.5 3.75\n$EndElementNodeData\n', encoding='utf-8')
    calls = []
    native = SimpleNamespace(open=lambda p: calls.append(p), merge=lambda p: None,
                             write=lambda p: None, view=SimpleNamespace(write=lambda *a: None))
    runtime.install_io_guards(native)
    if tag > 2**31 - 1:
        with pytest.raises(ValueError, match='exceeds 2147483647'):
            native.open(path)
        assert calls == []
    else:
        native.open(path)
        assert calls == [path]


def test_pure_reader_keeps_size_t_element_field_tag():
    tag = 2**63 + 17
    rows = msh._decode_data_blocks('ElementNodeData',
             [f'1 1 {tag} {tag}', '15 1 1', str(tag), '1.25 -2.5 3.75'], 3, 1)
    assert rows == [f'{tag} 1 1.25 -2.5 3.75']


@pytest.mark.parametrize('count,values,message', [
    (-1, '1 2 3', 'sample count/type mismatch'),
    (2, '1 2 3', 'sample count/type mismatch'),
    (1, '1 2 3 4', 'trailing data tokens'),
])
def test_bad_sample_counts_and_value_widths_fail(count, values, message):
    with pytest.raises(ValueError, match=message):
        msh._decode_data_blocks('ElementNodeData',
             ['1 1 17 17', f'15 {count} 1', '17', values], 3, 1)


@pytest.mark.parametrize('scheme,count', [('', 0), ('', 1), ('', 3), ('', 5), ('named', 0)])
def test_invalid_element_samples_rejected_by_all_entry_points(tmp_path, scheme, count):
    strings = '2\n"v"\n"' + scheme + '"\n'
    path = tmp_path / 'invalid.msh'
    path.write_text('$MeshFormat\n4.2 0 8\n$EndMeshFormat\n$ElementNodeData\n'
                    + strings + '1\n0\n3\n0\n1\n1\n1 1 17 17\n'
                    + f'1 {count} 1\n17\n' + '1.25 ' * count
                    + '\n$EndElementNodeData\n', encoding='utf-8')
    with pytest.raises(ValueError, match='sample count/type mismatch'):
        msh._validate_native_fields(path)
    checked = msh.validate_msh(path)
    assert not checked['ok'], checked
    assert any('sample count/type mismatch' in error for error in checked['errors']), checked
    with pytest.raises(ValueError, match='sample count/type mismatch'):
        msh.read_msh_data(path)
