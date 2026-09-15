"""Saved ESIM evidence checks; no numerical solver is launched on LAB."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / 'docs/ih_esim_benchmark'


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def test_esim_spatial_scenes_and_provenance():
    path = DOC / 'esim_spatial_demo.ipynb'
    nb = json.loads(path.read_bytes())
    meta = nb['metadata']['radia']
    assert meta['notebook_role'] == 'example'
    assert meta['webgui_required'] and meta['webgui_field_required']
    assert set(meta['citation_keys']) == {'hollaus2026nonlinear','yuferev2009surface'}
    from radia_mcp.bibliography.plans.T14_canonical import _generated_keys, _citation_source_sha256
    bib = ROOT/'packages/radia-mcp/src/radia_mcp/bibliography/data/references.bib'
    assert set(_generated_keys(path.with_suffix('.bbl').read_bytes())) == set(meta['citation_keys'])
    refs=next(c for c in nb['cells'] if c.get('id')=='radia-canonical-references')
    assert refs['metadata']['radia_bibliography']['selected_source_sha256'] == _citation_source_sha256(meta['citation_keys'],bib.read_bytes())
    state = nb['metadata']['widgets']['application/vnd.jupyter.widget-state+json']['state']
    views = []
    for cell in nb['cells']:
        if cell['cell_type'] != 'code':
            continue
        assert cell['execution_count'] is not None
        for output in cell['outputs']:
            assert output['output_type'] != 'error'
            view = output.get('data', {}).get('application/vnd.jupyter.widget-view+json')
            if view:
                payload = json.dumps(state[view['model_id']]['state'])
                assert 'webgui' in payload.lower() and len(payload) > 1000
                views.append(view)
    assert len(views) == 10  # assembly, mesh, 3 Ht, 3 qlocal, ReZs, qtransfer
    side = json.loads((DOC/'esim_spatial_demo_result.json').read_bytes())
    assert side['notebook_sha256'] == digest(path)
    assert side['execution']['hostname'].lower() in {'mdx1','mdx2'}
    assert 'pending' in side['execution']['visual_qa'].lower()


def test_esim_fresh_cases_and_mesh_scope():
    folder = DOC / 'spatial_demo_data'
    record = json.loads((folder/'compute_record.json').read_bytes())
    assert record['hostname'].lower() in {'mdx1','mdx2'}
    for name, expected in record['inputs'].items():
        assert digest(ROOT/'src/radia/panels/samples'/name) == expected
    assert [r['frequency_hz'] for r in record['cases']] == [10000,50000,100000]
    for item in record['cases']:
        assert item['returncode'] == 0
        path = folder / f"f{item['frequency_hz']}.json"
        assert digest(path) == item['result_sha256']
        d = json.loads(path.read_bytes())
        assert d['esim_converged'] and d['current_A'] == 100
        assert d['coupling_mode'] == 'weak' and d['wp_basis_order'] == 1
        assert len(d['esim_per_panel_H_t']) == len(d['esim_per_panel_Z_s_real']) == d['wp_ndof']
        assert min(d['esim_per_panel_Z_s_real']) > 0
    check = json.loads((folder/'vol_check.json').read_bytes())
    assert check['passed'] and check['threshold_pct'] == 10
    for group,name in [('materials','workpiece'),('boundaries','sibc')]:
        assert abs(next(x for x in check[group] if x['name'] == name)['error_pct']) < 1
    metrics = json.loads((folder/'spatial_metrics.json').read_bytes())
    assert len(metrics) == 3
    for r in metrics:
        assert r['P_BIE_W'] > 0 and r['P_local_W'] > 0 and r['P_transfer_W'] > 0
        # Conservation check only: this does NOT validate the spatial pattern.
        assert r['transfer_relative_gap'] < 1e-5


def test_esim_manual_distinguishes_loss_definitions():
    from radia_mcp.ih.esim_knowledge import ESIM_USAGE_OVERVIEW
    assert 'esim_spatial_demo.ipynb' in ESIM_USAGE_OVERVIEW
    assert 'qsurf_sol' in ESIM_USAGE_OVERVIEW
    assert 'hollaus2026nonlinear' in ESIM_USAGE_OVERVIEW
    assert 'not a relaxed FEM acceptance gate' in ESIM_USAGE_OVERVIEW
