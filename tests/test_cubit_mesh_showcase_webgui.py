"""Saved-scene acceptance, without importing a solver or executing a notebook."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / 'docs/cubit_mesh_export'


def test_curved_hex_showcase_has_four_recoverable_scenes():
    nb = json.loads((DOC / 'cubit_mesh_export_showcase.ipynb').read_bytes())
    assert nb['metadata']['radia']['webgui_required'] is True
    assert nb['metadata']['radia']['webgui_field_required'] is True
    states = nb['metadata']['widgets']['application/vnd.jupyter.widget-state+json']['state']
    scenes = []
    for cell in nb['cells']:
        if cell['cell_type'] != 'code':
            continue
        assert cell['execution_count'] is not None
        for output in cell['outputs']:
            assert output['output_type'] != 'error'
            widget = output.get('data', {}).get('application/vnd.jupyter.widget-view+json')
            if widget:
                state = states[widget['model_id']]['state']
                # No orphan model IDs: the actual WebGUI scene must survive kernel exit.
                assert 'webgui' in json.dumps(state).lower()
                assert len(json.dumps(state)) > 1000
                scenes.append(cell)
    assert len(scenes) == 4
    assert any('Draw(solution, mesh,' in ''.join(c['source']) for c in scenes)


def test_curved_hex_evidence_matches_committed_inputs():
    result = json.loads((DOC / 'cubit_mesh_webgui_results.json').read_bytes())
    assert result['host'].lower() in {'mdx1', 'mdx2'}
    assert [row['order'] for row in result['meshes']] == [1, 2, 3]
    for row in result['meshes']:
        assert row['check_vol']['passed'] is True
        assert row['elements'] == {'ET.HEX': 56}
        assert hashlib.sha256((DOC / row['mesh_path']).read_bytes()).hexdigest() == row['sha256']
    assert result['poisson']['l2_error'] > 0
    assert result['poisson']['ndof'] > 0


def test_showcase_sidecar_tracks_notebook_and_compute_host():
    path = DOC / 'cubit_mesh_export_showcase.ipynb'
    result = json.loads(path.with_name(path.stem + '_result.json').read_bytes())
    assert result['notebook_sha256'] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert result['execution_provenance']['host'].lower() in {'mdx1', 'mdx2'}
