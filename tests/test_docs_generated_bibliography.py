"""Canonical generated-reference contract, integrated from mcp-development-main.

This gate checks tracked public notebooks only. It does not certify the
scientific adequacy of a reference or execute a numerical cell.
"""
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def test_tracked_notebook_bibliographies_keep_generated_bbl_and_display():
    from radia_mcp.bibliography.plans.T14_canonical import _citation_source_sha256, _generated_keys

    canonical = ROOT / 'packages/radia-mcp/src/radia_mcp/bibliography/data/references.bib'
    source = canonical.read_bytes()
    paths = subprocess.check_output(['git', 'ls-files', 'docs'], cwd=ROOT, text=True).splitlines()
    notebooks = [ROOT / p for p in paths if p.endswith('.ipynb') and (ROOT / p).exists()]
    assert notebooks
    for path in notebooks:
        notebook = json.loads(path.read_text(encoding='utf-8'))
        radia = notebook.get('metadata', {}).get('radia', {})
        declaration = radia.get('bibliography')
        assert declaration is not None, path
        keys = declaration['keys']
        assert keys and len(keys) == len(set(keys)), path
        assert set(keys) == set(radia['citation_keys']), path
        bbl = path.with_suffix('.bbl')
        raw = bbl.read_bytes()
        generated = _generated_keys(raw)
        assert len(generated) == len(keys) and set(generated) == set(keys), path
        cells = [c for c in notebook['cells'] if 'radia_bibliography' in c.get('metadata', {})]
        assert len(cells) == 1, path
        display = cells[0]
        provenance = display['metadata']['radia_bibliography']
        assert provenance['style'] == declaration.get('style', 'IEEEtran'), path
        assert provenance['selected_source_sha256'] == _citation_source_sha256(keys, source), path
        assert provenance['bbl'] == bbl.name, path
        assert provenance['bbl_sha256_lf'] == hashlib.sha256(raw.replace(b'\r\n', b'\n')).hexdigest(), path
        rendered = ''.join(display['source'])
        assert provenance['display_sha256_lf'] == hashlib.sha256(rendered.replace('\r\n', '\n').encode('utf-8')).hexdigest(), path
        assert 'doc-bibliography' in rendered, path
        assert all(f"id='X{key}'" in rendered for key in keys), path
