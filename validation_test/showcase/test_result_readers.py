"""Portable numerical-record and result-reader contracts; no solver imports."""
import ast
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
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
    records={k:json.loads((Path(__file__).parent/'results'/f'{k}.json').read_text()) for k in CASES}
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
