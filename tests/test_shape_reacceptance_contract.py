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
