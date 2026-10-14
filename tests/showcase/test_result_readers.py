"""Portable numerical-record and result-reader contracts; no solver imports."""
import ast
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
LANE=ROOT/"validation_test/showcase"
CASES={
    'complex_coil':'complex_coil_geometry/complex_coil.ipynb',
    'motor':'electric_machine/cogging_skew_demo.ipynb',
    'winding':'stream_function/theory.ipynb',
    'lift':'maglev/maglev_showcase.ipynb',
    'particles':'gmsh_post/em_particle_orbits.ipynb',
    'heating':'esim_spatial/esim_spatial_demo.ipynb',
}


def test_readers_do_not_import_solvers_and_keep_saved_plots():
    for case,path in CASES.items():
        nb=json.loads((ROOT/'docs'/path).read_text(encoding='utf-8'))
        cells=[c for c in nb['cells'] if c['cell_type']=='code']
        code='\n'.join(''.join(c['source']) for c in cells)
        assert f'showcase/results/{case}.json' in code
        for node in ast.walk(ast.parse(code)):
            if isinstance(node,ast.Import):
                assert all(a.name.split('.')[0] in {'json','numpy','matplotlib'} for a in node.names)
            if isinstance(node,ast.ImportFrom):
                assert node.module=='pathlib'
        outputs=[o for c in cells for o in c.get('outputs',[])]
        assert not any(o['output_type']=='error' for o in outputs)
        assert any('image/png' in o.get('data',{}) for o in outputs)


def test_numerical_records_retain_acceptance_quantities():
    records={k:json.loads((LANE/'results'/f'{k}.json').read_text()) for k in CASES}
    for key,record in records.items():
        assert record['case']==key and record['generated_at_utc'] and record['versions']['radia']
        assert record['result']['checks'] and all(record['result']['checks'].values())
    assert records['complex_coil']['result']['current_reversal_relative_error']<1e-10
    assert records['winding']['result']['continuous_relative_residual']<1e-6
    assert records['winding']['result']['wire_relative_rms']<.06
    motor=records['motor']['result']
    assert motor['dominant_order']==2 and motor['mean_fraction']<.15
    for row in motor['skew_rows']:
        assert math.isclose(row['amplitude_ratio'],row['predicted_ratio'],rel_tol=.02,abs_tol=1e-12)
    lift=records['lift']['result']
    assert 0.99<lift['lift_N'][-1]/lift['perfect_conductor_lift_N']<1.01
    for item in records['particles']['result']['tracks']:
        assert item['track']['success']
        assert item['track']['maximum_relative_kinetic_energy_drift']<1e-6
    for row in records['heating']['result']['cases']:
        assert row['converged'] and math.isfinite(row['power_W']) and row['power_W']>0
        old=json.loads((ROOT/'validation_test/esim_spatial'/f"f{row['frequency_hz']}.json").read_text())
        assert math.isclose(row['power_W'],old['P_wp_W'],rel_tol=1e-6)


def test_recorded_result_hashes_are_unchanged():
    import hashlib
    campaign = json.loads((LANE/'campaign.json').read_text())
    for name, expected in campaign['result_sha256'].items():
        assert hashlib.sha256((LANE/'results'/name).read_bytes()).hexdigest() == expected


def test_archived_study_narrative_links_exist():
    import re
    from urllib.parse import unquote
    for path in (LANE/'studies').glob('*/*.ipynb'):
        notebook = json.loads(path.read_text(encoding='utf-8'))
        markdown = '\n'.join(''.join(c['source']) for c in notebook['cells'] if c['cell_type']=='markdown')
        for match in re.finditer(r'\]\(([^)\s]+)\)|(?:src|href)="([^"]+)"', markdown):
            link = next(g for g in match.groups() if g)
            if link.startswith(('https:', 'http:', '#', 'mailto:', 'data:')):
                continue
            target = path.parent / unquote(link.split('#')[0])
            assert target.exists(), (path, link)


def test_scratch_directory_never_targets_repository(tmp_path):
    import importlib.util
    import pytest
    spec = importlib.util.spec_from_file_location('showcase_recording', LANE/'recording.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with pytest.raises(ValueError, match='outside the repository'):
        module.scratch_directory(ROOT/'validation_test/showcase/results/intermediate')
    assert module.scratch_directory(tmp_path/'scratch') == tmp_path/'scratch'


def test_heating_driver_scratch_and_failed_run_preserves_record(tmp_path, monkeypatch):
    """Exercise orchestration with a fake child solver, without numerical imports."""
    import importlib.util
    import runpy
    import sys
    import types
    import subprocess
    import pytest
    lane = LANE
    monkeypatch.syspath_prepend(str(lane))
    fake_package = tmp_path/'fake_package'
    (fake_package/'panels').mkdir(parents=True)
    (fake_package/'panels/calc_inductance.py').write_text('# fake child solver')
    monkeypatch.setitem(sys.modules, 'radia', types.SimpleNamespace(__file__=str(fake_package/'__init__.py')))
    output = tmp_path/'accepted.json'
    output.write_text('accepted record')
    scratch = tmp_path/'scratch'
    monkeypatch.setattr(sys, 'argv', ['run_heating.py', '--workdir', str(scratch), '--output', str(output)])
    destinations = []
    def failed_child(command, **kwargs):
        target = Path(command[command.index('--output')+1])
        destinations.append(target)
        target.write_text(json.dumps({'status': 'ok', 'esim_converged': True, 'P_wp_W': float('nan')}))
        return types.SimpleNamespace(returncode=0)
    monkeypatch.setattr(subprocess, 'run', failed_child)
    with pytest.raises(AssertionError, match='convergence/power'):
        runpy.run_path(str(lane/'run_heating.py'), run_name='__main__')
    assert output.read_text() == 'accepted record'
    assert destinations and all(p.parent == scratch for p in destinations)
