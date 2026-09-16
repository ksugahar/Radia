"""Protect independent operating guidance and executable recipe artifacts."""
import json
import runpy

from cubit_mesh_export.mcp import server
from cubit_mesh_export.mcp.knowledge import scripting


def test_curated_python_roundtrips_json_values(monkeypatch, tmp_path):
    pool = tmp_path / 'learned.jsonl'
    record = {'signature': {'volumes': 1, 'surf_per_vol': 6, 'closed': True},
              'recipe': ['volume all scheme tetmesh'],
              'quality': {'min': 0.7, 'optional': None, 'failed': False},
              'source': 'local-test'}
    pool.write_text(json.dumps(record), encoding='utf-8')
    monkeypatch.setattr(server, '_learned_recipes_path', lambda: pool)
    output = tmp_path / 'curated.py'
    response = json.loads(server.cubit_curate_learned_recipes(str(output)))
    assert response['status'] == 'ok'
    restored = runpy.run_path(str(output))['CURATED'][0]
    for key in record:
        assert restored[key] == record[key]
    assert 'consented via' not in output.read_text()
    assert 'radia-mcp' not in response['next_step']


def test_manual_does_not_reintroduce_radia_install_or_compute_host_deploy():
    text = '\n'.join(value for value in vars(scripting).values()
                     if isinstance(value, str))
    assert 'pip install radia[cubit]' not in text
    assert 'Radia adds optional toolbar integration' not in text
    assert '`release_quad.py phase8e` installs PyPI wheels for `radia` and' not in text
    assert 'pip install cubit-mesh-export' in text
