"""Small test-only MSH 4.2 field fixture builder, independent of production."""


def field_section(kind, name, ncomp, rows, time=0., step=0):
    rows = list(rows)
    rows = [(r[0], *r[1]) if len(r) == 2 and isinstance(r[1], (list, tuple)) else r for r in rows]
    tags = [r[0] for r in rows]
    lines = ['$'+kind, '1', f'"{name}"', '1', str(time), '3',
             str(step), str(ncomp), str(len(rows))]
    if rows:
        lines += [f'1 {len(tags)} {min(tags)} {max(tags)}', str(len(tags))]
        lines += [str(t) for t in tags]
        lines += [' '.join(str(v) for v in r[1:]) for r in rows]
    else:
        lines += ['0 0 0 0']
    return '\n'.join(lines + ['$End'+kind, ''])
