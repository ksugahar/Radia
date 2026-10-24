"""Guards: bounded field reads and text-cookie tag replay across chunks.

Assumes: ASCII MSH 4.2 data with UTF-8 headers and optional CRLF newlines.
This is a parser contract, not numerical validation or a timing gate.
"""
import io

import pytest

from radia_mcp.gmsh import msh_inspect as msh


@pytest.mark.parametrize('chunk', [1, 2, 7, 65536])
def test_tag_cursor_replays_split_tokens_and_delimiters(chunk, monkeypatch):
    monkeypatch.setattr(msh._DataTokens, '_CHUNK', chunk)
    stream = io.StringIO('日本語\r\n17 268435456\t29\n1.25 -2.5 3.75')
    stream.readline()
    tokens = msh._DataTokens(stream)
    cursor = tokens.position()
    assert [next(tokens) for _ in range(3)] == ['17', '268435456', '29']
    replay = msh._DataTokens(io.StringIO(stream.getvalue()))
    replay.seek(cursor)
    assert [next(replay) for _ in range(3)] == ['17', '268435456', '29']
    assert list(tokens) == ['1.25', '-2.5', '3.75']


@pytest.mark.parametrize('newline', ['\n', '\r\n'])
@pytest.mark.parametrize('kind', ['NodeData', 'ElementData', 'ElementNodeData'])
def test_bounded_single_line_field_and_utf8_cookie(tmp_path, monkeypatch, newline, kind):
    monkeypatch.setattr(msh._DataTokens, '_CHUNK', 31)
    n = 2000
    block = f'15 1 {n}' if kind == 'ElementNodeData' else str(n)
    header = f'${kind}\n1\n"日本語"\n1\n0\n3\n0\n1\n{n}\n1 {n} 1 {n}\n{block}\n'
    text = header + ' '.join(map(str, range(1, n + 1))) + '\n'
    text += ' '.join(str(-i) for i in range(1, n + 1)) + f'\n$End{kind}\n'
    path = tmp_path / 'field.msh'
    path.write_bytes(text.replace('\n', newline).encode('utf-8'))

    class BoundedBody(msh._SectionBody):
        def readline(self, size=-1):
            assert 0 < size <= 31, f'unbounded token read: {size}'
            return super().readline(size)

    def reopen(): return path.open(encoding='utf-8')
    with reopen() as source:
        for _ in header.splitlines()[:-2]: source.readline()
        records = list(msh._iter_data_blocks(
            kind, BoundedBody(source, kind), reopen, 1, n, native_limit=2**31 - 1))
    nodes = 1 if kind == 'ElementNodeData' else None
    assert records == [(i, nodes, [str(-i)]) for i in range(1, n + 1)]


def test_closure_after_chunked_whitespace_still_required(tmp_path, monkeypatch):
    monkeypatch.setattr(msh._DataTokens, '_CHUNK', 7)
    prefix = '$NodeData\n1\n"v"\n1\n0\n3\n0\n1\n1\n1 1 1 1\n1\n1\n7\n'
    path = tmp_path / 'field.msh'
    path.write_text(prefix + ' ' * 50 + '$EndNodeData\n', encoding='utf-8')
    msh._validate_native_fields(path)
    path.write_text(prefix + ' ' * 50 + '$EndOther\n', encoding='utf-8')
    with pytest.raises(ValueError, match='not closed'):
        msh._validate_native_fields(path)


def test_one_value_may_span_multiple_chunks(monkeypatch):
    monkeypatch.setattr(msh._DataTokens, '_CHUNK', 17)
    value = '0' * 1000 + '1.25'
    tokens = msh._DataTokens(io.StringIO(value + '\t2'))
    assert list(tokens) == [value, '2']


@pytest.mark.parametrize('chunk', [7, 65536])
@pytest.mark.parametrize('suffix', [' x', 'X'])
def test_closure_marker_must_end_after_its_name(tmp_path, monkeypatch, chunk, suffix):
    monkeypatch.setattr(msh._DataTokens, '_CHUNK', chunk)
    path = tmp_path / 'field.msh'
    path.write_text('$NodeData\n1\n"v"\n1\n0\n3\n0\n1\n1\n'
                    '1 1 1 1\n1\n1\n7\n$EndNodeData' + suffix + '\n', encoding='utf-8')
    with pytest.raises(ValueError):
        msh._validate_native_fields(path)


@pytest.mark.parametrize('chunk', [1, 7])
def test_closure_marker_cannot_end_a_data_line(tmp_path, monkeypatch, chunk):
    monkeypatch.setattr(msh._DataTokens, '_CHUNK', chunk)
    path = tmp_path / 'field.msh'
    path.write_text('$NodeData\n1\n"v"\n1\n0\n3\n0\n1\n1\n'
                    '1 1 1 1\n1\n1\n7 $EndNodeData\n', encoding='utf-8')
    with pytest.raises(ValueError):
        msh._validate_native_fields(path)
