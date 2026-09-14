"""Cross-host scientific receipts must not mix designs or reuse stale success."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    'shape_reacceptance', ROOT / 'validation_test/isochronous_topopt/run_shape_reacceptance.py')
DRIVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DRIVER)


def receipts(tmp_path):
    preparation = {'inputs': {'design.stl': 'design-a'}}
    (tmp_path / 'prepare.json').write_text(json.dumps(preparation), encoding='utf-8')
    meshes = {'prepare_sha256': DRIVER.digest(tmp_path / 'prepare.json'),
              'input_sha256': preparation['inputs'].copy()}
    for key in ('mesh_results', 'checks', 'vol_sha256'):
        meshes[key] = dict.fromkeys(DRIVER.MESH_NAMES, {})
    return preparation, meshes


def test_mesh_parent_binds_exact_preparation_and_inputs(tmp_path):
    preparation, meshes = receipts(tmp_path)
    DRIVER.verify_mesh_parent(tmp_path, preparation, meshes)
    meshes['input_sha256'] = {'design.stl': 'design-b'}
    with pytest.raises(ValueError, match='inputs differ'):
        DRIVER.verify_mesh_parent(tmp_path, preparation, meshes)
    meshes['input_sha256'] = preparation['inputs']
    (tmp_path / 'prepare.json').write_text('other design')
    with pytest.raises(ValueError, match='another preparation'):
        DRIVER.verify_mesh_parent(tmp_path, preparation, meshes)


@pytest.mark.parametrize('field', ['mesh_results', 'checks', 'vol_sha256'])
def test_empty_mesh_set_is_not_vacuous_success(tmp_path, field):
    preparation, meshes = receipts(tmp_path)
    meshes[field] = {}
    with pytest.raises(ValueError, match='mesh set'):
        DRIVER.verify_mesh_parent(tmp_path, preparation, meshes)


def test_failed_phase_cannot_reuse_directory(tmp_path):
    with pytest.raises(RuntimeError, match='failure'):
        with DRIVER.phase_receipt(tmp_path, 'evaluate'):
            raise RuntimeError('injected early failure')
    state = json.loads((tmp_path / 'evaluate.state.json').read_text())
    assert state['status'] == 'failed'
    assert state['run_id']
    with pytest.raises(FileExistsError):
        with DRIVER.phase_receipt(tmp_path, 'evaluate'):
            pytest.fail('a failed run directory must not be reused')


def test_old_result_without_state_is_not_reused(tmp_path):
    (tmp_path / 'field.json').write_text('{"completed":true}')
    with pytest.raises(FileExistsError):
        with DRIVER.phase_receipt(tmp_path, 'evaluate'):
            pytest.fail('old success must not be overwritten')


def test_success_receipt_finishes_only_after_body(tmp_path):
    with DRIVER.phase_receipt(tmp_path, 'prepare') as receipt:
        assert json.loads((tmp_path / 'prepare.state.json').read_text())['status'] == 'in_progress'
        run_id = receipt['run_id']
    assert json.loads((tmp_path / 'prepare.state.json').read_text()) == {
        'run_id': run_id, 'phase': 'prepare', 'status': 'completed'}


@pytest.mark.parametrize('mutation', [dict(run_id='other'), dict(phase='mesh'), dict(status='failed')])
def test_parent_receipt_rejects_other_run_or_phase(tmp_path, mutation):
    result = {'run_id': 'expected'}
    state = {'run_id': 'expected', 'phase': 'prepare', 'status': 'completed'}
    (tmp_path / 'prepare.json').write_text(json.dumps(result))
    (tmp_path / 'prepare.state.json').write_text(json.dumps(state))
    DRIVER.verify_parent_receipt(tmp_path, 'prepare')
    state.update(mutation)
    (tmp_path / 'prepare.state.json').write_text(json.dumps(state))
    with pytest.raises(ValueError, match='receipt mismatch'):
        DRIVER.verify_parent_receipt(tmp_path, 'prepare')


def completed_chain(out):
    (out / 'design.stl').write_bytes(b'first-party synthetic input')
    preparation = {'run_id': 'prepare-id', 'staircase': {'J': 0.75},
                   'inputs': {'design.stl': DRIVER.digest(out / 'design.stl')}}
    (out / 'prepare.json').write_text(json.dumps(preparation))
    meshes = {'run_id': 'mesh-id', 'prepare_sha256': DRIVER.digest(out / 'prepare.json'),
              'input_sha256': preparation['inputs'], 'vol_sha256': {},
              'mesh_results': dict.fromkeys(DRIVER.MESH_NAMES, {}),
              'checks': dict.fromkeys(DRIVER.MESH_NAMES, {'passed': True})}
    for name in DRIVER.MESH_NAMES:
        (out / f'{name}.vol').write_bytes(name.encode())
        meshes['vol_sha256'][name] = DRIVER.digest(out / f'{name}.vol')
    (out / 'mesh.json').write_text(json.dumps(meshes))
    field = {'run_id': 'evaluate-id', 'completed': True, 'preparation': preparation,
             'mesh_identity': meshes, 'mesh_receipt_sha256': DRIVER.digest(out / 'mesh.json'),
             'J_staircase': 0.75, 'provenance': {'driver_sha256': 'original-driver'}}
    (out / 'field.json').write_text(json.dumps(field))
    for phase in ('prepare', 'mesh', 'evaluate'):
        (out / f'{phase}.state.json').write_text(json.dumps({
            'phase': phase, 'run_id': f'{phase}-id', 'status': 'completed'}))
    return field


def test_post_execution_audit_preserves_numerical_provenance(tmp_path):
    completed_chain(tmp_path)
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    result = DRIVER.audit_completed_receipts(tmp_path)
    assert result['passed'] is True
    assert result['numerical_recomputed'] is False
    assert result['original_driver_sha256'] == 'original-driver'
    assert result['validator_sha256'] == DRIVER.digest(DRIVER.__file__)
    assert before == {p.name: p.read_bytes() for p in tmp_path.iterdir()}


@pytest.mark.parametrize('mutation', ['uuid', 'mesh_receipt', 'input', 'vol', 'staircase'])
def test_post_execution_audit_rejects_changed_evidence(tmp_path, mutation):
    field = completed_chain(tmp_path)
    if mutation == 'uuid':
        field['run_id'] = 'other'
    elif mutation == 'mesh_receipt':
        field['mesh_receipt_sha256'] = 'other'
    elif mutation == 'staircase':
        field['J_staircase'] = 0.5
    elif mutation == 'input':
        (tmp_path / 'design.stl').write_bytes(b'changed')
    else:
        (tmp_path / 'hex_fine.vol').write_bytes(b'changed')
    (tmp_path / 'field.json').write_text(json.dumps(field))
    with pytest.raises(ValueError):
        DRIVER.audit_completed_receipts(tmp_path)
