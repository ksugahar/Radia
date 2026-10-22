"""Tree reuse and impact contracts need no solver, credentials or network."""
import importlib.util
import io
import json
from pathlib import Path
import zipfile

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('tree_receipt', ROOT / 'tools/ci_tree_receipt.py')
receipt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(receipt)
spec = importlib.util.spec_from_file_location('impact', ROOT / 'tools/sparsesolv_ci_impact.py')
impact = importlib.util.module_from_spec(spec)
spec.loader.exec_module(impact)


def fixture():
    current = dict(schema='radia.ci-tree.v1', lane='sparsesolv', repository='owner/repo', commit='main', tree='tree', run_id=22)
    current['test_set'] = receipt.numerical_test_set(ROOT)
    prior = dict(current, commit='merge-ref', run_id=11, impact_required=True)
    data = {
        'commits/main/pulls': [dict(number=7, merged_at='today', merge_commit_sha='main',
            base=dict(ref='main', repo=dict(full_name='owner/repo')),
            head=dict(sha='head', repo=dict(full_name='owner/repo')))],
        'branches/main/protection/required_status_checks': dict(contexts=['ams-regression'], checks=[]),
        'rules/branches/main': [],
        'actions/workflows/sparsesolv.yml/runs?event=pull_request&head_sha=head': [dict(id=11, status='completed', conclusion='success', event='pull_request', head_sha='head', path='.github/workflows/sparsesolv.yml', head_repository=dict(full_name='owner/repo'), pull_requests=[dict(number=7)])],
        'actions/runs/11/artifacts': [dict(id=33, name='sparsesolv-tree-identity', expired=False)],
        'git/commits/merge-ref': dict(tree=dict(sha='tree'), parents=[dict(sha='base'), dict(sha='head')]),
        'commits/head/check-runs?filter=latest': [dict(name='ams-regression', status='completed', conclusion='success', details_url='https://github.com/owner/repo/actions/runs/11/job/12')],
        'commits/head/statuses': [],
        'actions/runs/11/jobs?filter=latest': [dict(name='ams-regression', conclusion='success')],
    }
    return current, prior, data


def evaluate(current, prior, data):
    def get(path, binary=False):
        if binary:
            stream = io.BytesIO()
            with zipfile.ZipFile(stream, 'w') as archive:
                archive.writestr('ci-tree-receipt.json', json.dumps(prior))
            return stream.getvalue()
        return data[path]
    return receipt.reuse('sparsesolv', {}, current, get, lambda path, key=None: data[path])


def test_exact_tree_receipt_records_pr_and_run_provenance():
    current, prior, data = fixture()
    result = evaluate(current, prior, data)
    assert result['mode'] == 'tree already validated'
    assert (result['tree'], result['pr'], result['run_ids'], result['artifact_id']) == ('tree', 7, [11], 33)


def test_rulesets_and_pinned_check_applications_are_required():
    current, prior, data = fixture()
    data['rules/branches/main'] = [dict(type='required_status_checks', parameters=dict(
        required_status_checks=[dict(context='ams-regression', integration_id=123)]))]
    assert evaluate(current, prior, data) is None
    data['commits/head/check-runs?filter=latest'][0]['app'] = dict(id=123)
    assert evaluate(current, prior, data)


def test_unprotected_branch_still_requires_successful_observed_checks_and_lane():
    current, prior, data = fixture()
    data['branches/main/protection/required_status_checks'] = dict(contexts=[], checks=[])
    assert evaluate(current, prior, data)
    data['commits/head/check-runs?filter=latest'][0]['conclusion'] = 'failure'
    assert evaluate(current, prior, data) is None


def test_skipped_pr_cannot_certify_full_postmerge_lane():
    current, prior, data = fixture()
    prior['impact_required'] = False
    assert evaluate(current, prior, data) is None


def test_solver_pr_only_receipt_cannot_certify_main():
    current, prior, data = fixture()
    prior['test_set'] = {'solver-pr': ['tests/test_ih_thermal.py']}
    assert evaluate(current, prior, data) is None
    prior['test_set'] = {**current['test_set'], 'solver-numerics': []}
    assert evaluate(current, prior, data) is None
    del prior['test_set']
    assert evaluate(current, prior, data) is None


@pytest.mark.parametrize('change', ['tree', 'direct', 'head', 'failed', 'expired', 'required', 'job', 'foreign', 'parent', 'run', 'pending'])
def test_unproven_tree_never_reuses(change):
    current, prior, data = fixture()
    if change == 'tree': prior['tree'] = 'other'
    if change == 'direct': data['commits/main/pulls'] = []
    if change == 'head': data['commits/main/pulls'][0]['merge_commit_sha'] = 'other'
    if change == 'failed': data['actions/workflows/sparsesolv.yml/runs?event=pull_request&head_sha=head'][0]['conclusion'] = 'failure'
    if change == 'expired': data['actions/runs/11/artifacts'][0]['expired'] = True
    if change == 'required': data['branches/main/protection/required_status_checks']['contexts'].append('missing')
    if change == 'job': data['actions/runs/11/jobs?filter=latest'][0]['conclusion'] = 'skipped'
    if change == 'foreign': prior['repository'] = 'foreign/repo'
    if change == 'parent': data['git/commits/merge-ref']['parents'][1]['sha'] = 'other'
    if change == 'run': prior['run_id'] = 99
    if change == 'pending': data['commits/head/check-runs?filter=latest'][0]['status'] = 'in_progress'
    assert evaluate(current, prior, data) is None


@pytest.mark.parametrize('lane', receipt.LANES)
def test_unrelated_package_never_selects_heavy_lane(lane):
    for path in ('packages/eqnedit64/src/change.py', 'docs/guide.md'):
        assert not impact.lane_required(lane, {path}, None, impact.manifest())[0]


@pytest.mark.parametrize('lane', receipt.LANES)
def test_every_declared_existing_input_selects_its_lane(lane):
    data = impact.manifest()
    projection = data['lanes'][lane]['events']['pull_request']
    for path in impact.CONTRACT.tracked(ROOT):
        if path not in impact.CONTRACT.ROOT_CONFIGS and impact.CONTRACT.selected(path, projection):
            assert impact.lane_required(lane, {path}, None, data)[0], path
    for path in impact.CONTRACT.ROOT_CONFIGS:
        assert impact.lane_required(lane, {path}, lambda _: (None, None), data)[0]


def test_numerical_tests_are_retained_and_full_release_is_mandatory():
    tiers = json.loads((ROOT / 'tests/test_tier_manifest.json').read_text())['profiles']
    # Preserve the original 91 routes and explicitly register the added validations.
    esim_routes = {
        'tests/test_esim_adaptive_depth.py',
        'validation_test/induction_heating/test_esim_adaptive_depth_validation.py',
    }
    direct_routes = {
        'tests/test_planar_eddy_direct_inverse.py',
        'tests/test_vector_potential_direct_inverse.py',
    }
    topology_routes = {
        'validation_test/induction_heating/test_ngsbem_loop_star_topology_validation.py',
    }
    mmm_routes = {
        'tests/test_mmm_topology.py',
    }
    added_routes = esim_routes | direct_routes | topology_routes | mmm_routes
    assert added_routes <= set(tiers['solver-numerics']['paths'])
    assert len(tiers['solver-numerics']['paths']) == 91 + len(added_routes)
    assert 'solver-pr' not in tiers
    workflow = yaml.load((ROOT / '.github/workflows/sparsesolv.yml').read_text(), Loader=yaml.BaseLoader)
    assert workflow['on']['schedule'] and 'workflow_dispatch' in workflow['on']
    full = next(s for s in workflow['jobs']['ams-regression']['steps'] if s.get('name') == 'Run required NGSolve numerical contracts')
    assert "github.event_name" not in full['if']
    native = (ROOT / '.github/workflows/build-test.yml').read_text()
    assert '--profile solver-numerics' in native
    assert 'Run full release numerical contracts\n        if:' not in native


def test_evidence_changes_do_not_build_gpu_but_executable_helpers_do():
    data = impact.manifest()
    for path in ('validation_test/ams_cuda/run.py', 'validation_test/ams_cuda/setup_benchmark.py', 'validation_test/radia_optuna/validate_cuda.m'):
        lane = 'radia-optuna' if 'radia_optuna' in path else 'radia-simulink'
        assert impact.lane_required(lane, {path}, None, data)[0]
    assert impact.lane_required('radia-simulink', {'validation_test/ams_cuda/input.json'}, None, data)[0]


@pytest.mark.parametrize('path,expected', [
    ('src/radia/ih_axisym_coupled.py', True),
    ('src/radia/ih_thermal.py', True),
    ('src/radia/vim/_nonlinear.py', True),
    ('src/core/rad_parallel.h', True),
    ('Build.ps1', True),
    ('packages/cubit-mesh-export/src/cubit_mesh_export/toolbar_install.py', True),
    ('packages/radia-simulink/matlab/change.m', False),
    ('validation_test/maglev/test_demo_evidence.py', True),
])
def test_full_pr_dependency_closure(path, expected):
    def git(*args):
        if args[0] == 'merge-base': return 'base'
        if args[0] == 'diff': return 'M\0' + path + '\0'
        raise ValueError('missing blob selects root files conservatively')
    assert impact.required_for_event('sparsesolv', 'pull_request',
        {'pull_request': {'base': {'sha': 'base'}}}, git, impact.manifest())[0] is expected


def test_manual_and_api_failure_execute_instead_of_reusing(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv('GITHUB_REPOSITORY', 'owner/repo')
    monkeypatch.setenv('GITHUB_RUN_ID', '22')
    def unavailable(*args): raise OSError('offline')
    monkeypatch.setattr(receipt, 'reuse', unavailable)
    for event in ('workflow_dispatch', 'schedule', 'push'):
        assert not receipt.route('sparsesolv', event, {'ref': 'refs/heads/main'}, lambda *args: 'a' * 40)
        assert json.loads(Path('ci-tree-receipt.json').read_text())['mode'] == 'execute'


def test_python_closure_and_all_native_inputs_have_workflow_admission():
    spec = importlib.util.spec_from_file_location('ams_inputs', ROOT / 'tools/ci_ams_inputs.py')
    inputs = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(inputs)
    data = impact.manifest()
    lane = data['lanes']['sparsesolv']
    paths = inputs.python_inputs(ROOT)
    for pattern in lane['pr_heavy_inputs']:
        paths.update(p for p in impact.CONTRACT.tracked(ROOT) if not p.endswith('.md') and impact.CONTRACT.matches(p, pattern))
    for path in paths:
        assert impact.CONTRACT.selected(path, lane['events']['pull_request']), path


def test_full_profile_changes_cannot_be_exempted_as_unrelated():
    before = (ROOT / 'tests/test_tier_manifest.json').read_text()
    after = json.loads(before)
    after['profiles']['solver-numerics']['paths'].append('tests/new_numerical_contract.py')
    assert not impact.ams_tier_projection_unchanged(before, json.dumps(after))


@pytest.mark.parametrize('path', [
    'src/radia/p1_newton.py', 'src/radia/p1_linear.py',
    'src/radia/topology_optimization.py', 'src/radia/kelvin_solver.py',
    'pyproject.toml', 'packages/radia-mcp/pyproject.toml',
    'requirements-ci.txt', '.github/workflows/radia.yml',
    'tests/fixtures/input.h5',
])
def test_full_numerics_pr_admission(path):
    def git(*args):
        if args[0] == 'merge-base': return 'base'
        if args[0] == 'diff': return 'M\0' + path + '\0'
        raise ValueError('missing blob')
    required, _ = impact.required_for_event('sparsesolv', 'pull_request',
        {'pull_request': {'base': {'sha': 'base'}}}, git, impact.manifest())
    assert required
    workflow = yaml.load((ROOT / '.github/workflows/sparsesolv.yml').read_text(), Loader=yaml.BaseLoader)
    steps = workflow['jobs']['ams-regression']['steps']
    full = next(s for s in steps if s.get('name') == 'Run required NGSolve numerical contracts')
    assert 'github.event_name' not in full['if']
    assert '--profile solver-numerics' in full['run']
    assert all('--profile solver-pr' not in s.get('run', '') for s in steps)


def test_nightly_is_independent_and_reports_failures():
    workflow = yaml.load((ROOT / '.github/workflows/sparsesolv.yml').read_text(), Loader=yaml.BaseLoader)
    assert 'github.event_name' in workflow['concurrency']['group']
    job = workflow['jobs']['nightly-failure']
    assert "github.event_name == 'schedule'" in job['if']
    assert 'issues.create' in job['steps'][0]['with']['script']


def test_closure_follows_cli_loaders_and_non_python_fixture(tmp_path):
    spec = importlib.util.spec_from_file_location('fixture_inputs', ROOT / 'tools/ci_ams_inputs.py')
    inputs = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(inputs)
    files = {
        'tests/test_tier_manifest.json': json.dumps({'profiles': {
            'sparsesolv': {'paths': []}, 'solver-heavy': {'paths': []},
            'native-build-contracts': {'paths': []},
            'solver-numerics': {'extends': 'base', 'paths': []},
            'base': {'paths': ['tests/test_entry.py']}}}),
        'tests/test_entry.py': "import tests.helper\nimport subprocess\nsubprocess.run(['python', 'tools/child.py'])\n",
        'tests/helper.py': 'VALUE = 1\n',
        'tools/child.py': "from pathlib import Path\nPath('tests/fixtures/model.custom').read_bytes()\n",
        'tests/fixtures/model.custom': 'fixture',
    }
    for path, content in files.items():
        target = tmp_path / path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(content)
    closure = inputs.python_inputs(tmp_path)
    assert {'tests/test_entry.py', 'tests/helper.py', 'tools/child.py', 'tests/fixtures/model.custom'} <= closure
