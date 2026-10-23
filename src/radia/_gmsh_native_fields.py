"""Bound native Gmsh field ingestion without limiting mesh-only size_t tags."""
from pathlib import Path
import json
import re


NATIVE_FIELD_TAG_MAX = 2**31 - 1


def validate_native_field_tags(tags, association):
    """Reject unsupported native field indices; never renumber or substitute."""
    for tag in tags:
        if int(tag) > NATIVE_FIELD_TAG_MAX:
            raise ValueError(
                f'Gmsh native {association} field tag {tag} exceeds '
                f'{NATIVE_FIELD_TAG_MAX}; large-tag native field data is '
                'unsupported. Re-export with compact producer numbering; '
                'tags are not remapped and no fallback is used. '
                'Pure-Python readers and mesh-only tags are unaffected.')


def validate_native_msh_fields(path, _seen=None):
    """Preflight MSH 4.2 field tags before importing/initializing native Gmsh.

    Scan only field associations, not mesh tags: a mesh-only large index does
    not allocate dense field data. Static GEO Merge/Include targets are checked
    recursively; dynamic targets cannot be preflighted and fail explicitly.
    """
    path = Path(path).resolve()
    seen = set() if _seen is None else _seen
    if path in seen:
        return
    seen.add(path)
    if path.suffix.lower() not in ('.msh', '.geo'):
        return  # Other explicit CAD/list-data formats have no mesh field tags.
    text = path.read_text(encoding='utf-8')
    if path.suffix.lower() == '.geo':
        text = re.sub(r'/\*.*?\*/', '', text, flags=re.S)
        text = re.sub(r'(?m)^\s*//.*$', '', text)
        for match in re.finditer(r'\b(?:Merge|Include)\s+([^;]+);', text):
            try:
                target = json.loads(match[1].strip())
            except (ValueError, TypeError) as exc:
                raise ValueError(f'{path}: dynamic Merge/Include cannot be checked '
                                 'before native field access') from exc
            if not isinstance(target, str):
                raise ValueError(f'{path}: Merge/Include target must be a static string')
            validate_native_msh_fields(path.parent / target, seen)
        return
    rows = text.splitlines()
    if len(rows) < 3 or rows[0].strip() != '$MeshFormat' or rows[1].split() != ['4.2', '0', '8']:
        raise ValueError(f'{path}: native field hand-off requires ASCII MSH 4.2 (size 8)')
    lines = iter(rows)
    for line in lines:
        kind = line.strip()
        if kind not in ('$NodeData', '$ElementData', '$ElementNodeData'):
            continue
        association = 'node' if kind == '$NodeData' else 'element'
        try:
            strings = int(next(lines))
            if strings < 0:
                raise ValueError('negative string tag count')
            for _ in range(strings):
                next(lines)
            body = []
            end = kind.replace('$', '$End')
            for row in lines:
                if row.strip() == end:
                    break
                body.append(row)
            else:
                raise ValueError(f'missing {end}')
            tokens = iter(' '.join(body).split())
            nreal = int(next(tokens))
            if nreal < 0:
                raise ValueError('negative real tag count')
            for _ in range(nreal):
                next(tokens)
            ninteger = int(next(tokens))
            if ninteger < 3:
                raise ValueError('field requires three integer tags')
            integers = [int(next(tokens)) for _ in range(ninteger)]
            ncomp, declared = integers[1:3]
            blocks, count, minimum, maximum = [int(next(tokens)) for _ in range(4)]
            if ncomp <= 0 or min(blocks, count, declared) < 0 or count != declared:
                raise ValueError('invalid field count/component header')
            if count:
                validate_native_field_tags((maximum,), association)
            consumed = 0
            for _ in range(blocks):
                width = ncomp
                if kind == '$ElementNodeData':
                    int(next(tokens))  # element type
                    width = int(next(tokens))
                n = int(next(tokens))
                if n < 0 or n > count - consumed or width <= 0:
                    raise ValueError('invalid field block width/count')
                for _ in range(n):
                    validate_native_field_tags((int(next(tokens)),), association)
                for _ in range(n * width):
                    next(tokens)
                consumed += n
            if consumed != count or next(tokens, None) is not None:
                raise ValueError('field block count/trailing token mismatch')
        except StopIteration as exc:
            raise ValueError(f'{path}: truncated {kind} before native hand-off') from exc
