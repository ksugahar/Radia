"""Use the existing BibTeX/TeX4ht contract; never execute notebook code."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'packages/radia-mcp/src'))
from radia_mcp.bibliography.plans import T14_canonical as canonical


def sha_lf(value: str) -> str:
    return hashlib.sha256(value.replace('\r\n', '\n').encode('utf-8')).hexdigest()


def render(path: Path) -> dict:
    path = path.resolve()
    path.relative_to(ROOT / 'docs')
    raw = path.read_bytes()
    notebook = json.loads(raw)
    before = [c for c in notebook['cells'] if c['cell_type'] == 'code']
    meta = notebook.setdefault('metadata', {}).setdefault('radia', {})
    declared = meta.get('citation_keys', [])
    if not declared or len(set(declared)) != len(declared):
        raise ValueError(f'Explicit unique citation_keys required: {path}')
    keys = []
    for cell in notebook['cells']:
        if cell['cell_type'] != 'markdown' or cell.get('id') == 'radia-canonical-references':
            continue
        for key in re.findall(r'@([A-Za-z0-9_:./+-]+)', ''.join(cell['source'])):
            if key in declared and key not in keys:
                keys.append(key)
    keys.extend(k for k in declared if k not in keys)
    meta['bibliography'] = {'keys': keys, 'style': 'IEEEtran'}
    make4ht = shutil.which('make4ht')
    if not make4ht:
        raise RuntimeError('TeX4ht is required; no hand-authored display fallback')
    version = subprocess.run([make4ht, '--version'], capture_output=True, timeout=20)
    renderer = version.stdout.decode('utf-8', errors='replace').strip()
    if version.returncode:
        raise RuntimeError('Cannot identify TeX4ht renderer')
    temp_root = Path('C:/temp')
    temp_root.mkdir(exist_ok=True)
    tempfile.tempdir = str(temp_root)
    # Publish explicit metadata first: the existing generator snapshots it.
    path.write_text(json.dumps(notebook, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    receipt = canonical.bibliography_make_bbl(str(path))
    if not receipt.startswith('bibliography_make_bbl:'):
        raise RuntimeError(receipt)
    bbl = path.with_suffix('.bbl')
    selected = canonical._citation_source_sha256(keys, canonical.CANONICAL.read_bytes())
    if f'selected_source_sha256: {selected}' not in receipt:
        raise RuntimeError('Selected bibliography changed after generation')
    with tempfile.TemporaryDirectory(prefix='radia-docs-bbl-', dir=temp_root) as work_name:
        work = Path(work_name)
        (work / 'references.bbl').write_bytes(bbl.read_bytes())
        cjk = bool(re.search(r'[\u3000-\u9fff]', bbl.read_text(encoding='utf-8')))
        (work / 'references.tex').write_text(
            '\\documentclass{article}\n\\usepackage[T1]{fontenc}\n'
            '\\usepackage[utf8]{inputenc}\n\\usepackage{url}\n\\usepackage{amsmath,amssymb}\n'
            + ('\\usepackage{CJKutf8}\n' if cjk else '')
            + '\\begin{document}\n'
            + ('\\begin{CJK}{UTF8}{min}\n' if cjk else '')
            + '\\input{references.bbl}\n'
            + ('\\end{CJK}\n' if cjk else '')
            + '\\end{document}\n', encoding='utf-8')
        result = subprocess.run([make4ht, '-u', 'references.tex'], cwd=work,
                                capture_output=True, timeout=90)
        html_path = work / 'references.html'
        if result.returncode or not html_path.exists():
            raise RuntimeError('TeX4ht failed: ' + result.stdout.decode('utf-8', errors='replace')[-4000:])
        html = html_path.read_text(encoding='utf-8')
        match = re.search(r'<body[^>]*>(.*?)</body>', html, re.S | re.I)
        if not match:
            raise RuntimeError('TeX4ht output has no complete body')
        body = match[1].strip()
    source = f'## References\n\nGenerated from the canonical bibliography: [{bbl.name}]({bbl.name}).\n\n' + body + '\n'
    display = {
        'cell_type': 'markdown', 'id': 'radia-canonical-references',
        'metadata': {'radia_bibliography': {
            'bbl': bbl.name, 'bbl_sha256_lf': sha_lf(bbl.read_text(encoding='utf-8')),
            'selected_source_sha256': selected, 'style': 'IEEEtran',
            'renderer': renderer, 'display_sha256_lf': sha_lf(source),
        }}, 'source': source.splitlines(keepends=True),
    }
    notebook['cells'] = [c for c in notebook['cells'] if c.get('id') != display['id']
                         and 'radia_bibliography' not in c.get('metadata', {})] + [display]
    assert before == [c for c in notebook['cells'] if c['cell_type'] == 'code']
    path.write_text(json.dumps(notebook, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    return {'notebook': str(path.relative_to(ROOT)), 'keys': len(keys), 'receipt': receipt}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('notebooks', nargs='+', type=Path)
    args = parser.parse_args()
    for notebook in args.notebooks:
        print(json.dumps(render(notebook), ensure_ascii=False), flush=True)
