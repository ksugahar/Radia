"""Guards: profile, ownership or environment edits cannot silently bypass admission.

Assumes: source-contract payloads; parent distribution admission is conservative,
while the separately verified immutable plan selects actual numerical consumers.
"""
import copy
import importlib.util
import json
from pathlib import Path
import shutil

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('ams_projection', ROOT / 'tools/sparsesolv_ci_impact.py')
impact = importlib.util.module_from_spec(spec)
spec.loader.exec_module(impact)
TIER = 'tests/test_tier_manifest.json'
OWNERS = 'tools/ci_distribution_manifest.json'


@pytest.fixture(scope='module')
def real_inputs():
    spec = importlib.util.spec_from_file_location('real_ams_inputs', ROOT / 'tools/ci_ams_inputs.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(autouse=True)
def reuse_graph_cache_for_immutable_checkout(monkeypatch, real_inputs):
    # Reuse the real helper's existing ROOT-keyed cache across cases. Neither
    # traversal nor its results are mocked; temporary consumer trees use other keys.
    monkeypatch.setattr(impact, '_ams_inputs', lambda: real_inputs)


def documents():
    data = json.loads((ROOT / OWNERS).read_text(encoding='utf-8'))
    before = copy.deepcopy(data)
    before['rules'] = [rule for rule in before['rules']
                       if 'validation_records/' not in rule['pattern']]
    tier = json.loads((ROOT / TIER).read_text(encoding='utf-8'))
    unrelated = copy.deepcopy(tier)
    unrelated['profiles']['repo-contracts']['description'] += ' records-only contract'
    return {OWNERS: (json.dumps(before), json.dumps(data)),
            TIER: (json.dumps(tier), json.dumps(unrelated))}, data


def decide(blobs, data, event='pull_request', extra=None):
    paths = list(blobs) + ([extra] if extra else [])
    def git(*args):
        if args[0] == 'merge-base':
            return 'base'
        if args[0] == 'diff':
            return ''.join('M\0' + path + '\0' for path in paths)
        if args[0] == 'show':
            ref, path = args[1].split(':', 1)
            return blobs[path][0 if ref == 'base' else 1]
        raise AssertionError(args)
    return impact.required_for_event('sparsesolv', event,
        {'before': 'base', 'pull_request': {'base': {'sha': 'base'}}}, git, data)[0]


def test_records_rules_and_unrelated_profile_data_do_not_start_pr_ams():
    blobs, data = documents()
    assert not decide(blobs, data)
    assert decide(blobs, data, event='push')
    assert decide(blobs, data, event='workflow_dispatch')


@pytest.mark.parametrize('extra', [
    'src/radia/new_solver.py', 'src/core/new_header.h', 'Build.ps1',
    'CMakeLists.txt', 'cmake/toolchain.cmake', 'tests/fixtures/new.vol',
    'tools/ci_ams_inputs.py', '.github/workflows/sparsesolv.yml',
])
def test_data_exemptions_never_hide_mixed_numeric_build_or_gate_changes(extra):
    blobs, data = documents()
    assert decide(blobs, data, extra=extra)


@pytest.mark.parametrize('profile', ['sparsesolv', 'solver-numerics', 'solver-heavy', 'native-build-contracts'])
@pytest.mark.parametrize('mutation', ['path', 'budget', 'ancestor'])
def test_every_numerical_profile_and_ancestor_stays_an_execution_input(profile, mutation):
    blobs, data = documents()
    before = json.loads(blobs[TIER][0])
    after = copy.deepcopy(before)
    if mutation == 'path':
        after['profiles'][profile]['paths'].append('tests/new_numerical_contract.py')
    elif mutation == 'budget':
        after['profiles'][profile]['max_elapsed_seconds'] = 12345
    else:
        for value in (before, after):
            value['profiles']['new-parent'] = {'paths': []}
            value['profiles'][profile]['extends'] = 'new-parent'
        after['profiles']['new-parent']['paths'].append('tests/new_numerical_contract.py')
    blobs[TIER] = json.dumps(before), json.dumps(after)
    assert decide(blobs, data)
    assert decide(blobs, data, event='push')


@pytest.mark.parametrize('mutation', ['unknown', 'cycle', 'missing-parent', 'bool-budget',
                                    'bad-impact', 'bad-paths', 'unknown-profile-field'])
def test_invalid_or_future_tier_semantics_fail_closed(mutation):
    blobs, data = documents()
    after = json.loads(blobs[TIER][1])
    profile = after['profiles']['repo-contracts']
    if mutation == 'unknown': after['new_execution_option'] = True
    elif mutation == 'cycle': profile['extends'] = 'repo-contracts'
    elif mutation == 'missing-parent': profile['extends'] = 'missing'
    elif mutation == 'bool-budget': profile['max_elapsed_seconds'] = True
    elif mutation == 'bad-impact': after['impact_rules'] = []
    elif mutation == 'bad-paths': profile['paths'] = 'tests/new.py'
    else: profile['new_execution_option'] = True
    blobs[TIER] = blobs[TIER][0], json.dumps(after)
    assert decide(blobs, data)


@pytest.mark.parametrize('mutation', ['native-rule', 'lane', 'other-lane', 'attribution',
                                    'include-closure', 'unknown', 'records-field'])
def test_records_exemption_keeps_all_executable_projections_and_unknown_data_closed(mutation):
    blobs, data = documents()
    after = json.loads(blobs[OWNERS][1])
    if mutation == 'native-rule':
        after['rules'].append({'pattern': 'src/new_solver/**', 'owners': ['radia']})
    elif mutation in ('lane', 'other-lane'):
        name = 'sparsesolv' if mutation == 'lane' else 'radia-optuna'
        after['lanes'][name]['events']['pull_request']['paths'].append('src/new_solver/**')
    elif mutation == 'attribution':
        after['root_content_diff']['unknown'] = 'ignore'
    elif mutation == 'include-closure':
        after['shared_include_closures'].append({'new': ['src/new_solver.h']})
    elif mutation == 'unknown':
        after['future_execution_option'] = True
    else:
        next(rule for rule in after['rules'] if rule['pattern'] == 'validation_records/**')['future'] = True
    blobs[OWNERS] = blobs[OWNERS][0], json.dumps(after)
    assert decide(blobs, data)


@pytest.mark.parametrize('path', [TIER, OWNERS])
@pytest.mark.parametrize('value', ['{invalid', '{"schema": 1, "schema": 2}', '{"x": NaN}'])
def test_bad_json_cannot_exempt_a_pr(path, value):
    blobs, data = documents()
    blobs[path] = blobs[path][0], value
    assert decide(blobs, data)


def test_parent_admission_and_internal_case_plan_remain_separate():
    workflow = yaml.load((ROOT / '.github/workflows/sparsesolv.yml').read_text(encoding='utf-8'), Loader=yaml.BaseLoader)
    parent = [step['run'] for step in workflow['jobs']['ci-impact']['steps']
              if 'tools/sparsesolv_ci_impact.py' in step.get('run', '')]
    assert len(parent) == 1 and '--lane sparsesolv' in parent[0]
    internal = [step['run'] for step in workflow['jobs']['ams-regression']['steps']
                if step.get('id') == 'native-impact']
    assert len(internal) == 1 and 'tools/ci_numerical_plan.py --github-event' in internal[0]
    assert '--github-output' in internal[0]


@pytest.mark.parametrize('consumed', [TIER, OWNERS])
def test_a_numerical_consumer_of_metadata_prevents_its_exemption(tmp_path, monkeypatch, consumed):
    blobs, data = documents()
    (tmp_path / 'tools').mkdir()
    (tmp_path / 'tests').mkdir()
    shutil.copy2(ROOT / 'tools/ci_ams_inputs.py', tmp_path / 'tools/ci_ams_inputs.py')
    (tmp_path / OWNERS).write_text(json.dumps(data), encoding='utf-8')
    (tmp_path / TIER).write_text(json.dumps({'profiles': {
        name: {'paths': ['tests/test_consumer.py']}
        for name in ('sparsesolv', 'solver-numerics', 'solver-heavy', 'native-build-contracts')}}), encoding='utf-8')
    (tmp_path / 'tests/test_consumer.py').write_text(
        "from pathlib import Path\nPath(" + repr(consumed) + ").read_text()\n", encoding='utf-8')
    monkeypatch.setattr(impact, 'ROOT', tmp_path)
    assert decide(blobs, data)


@pytest.mark.parametrize('mutation', ['beam-impact', 'empty-native', 'native-budget'])
def test_review_native_contract_counterexamples_remain_required(mutation):
    blobs, data = documents()
    before = json.loads(blobs[TIER][0])
    after = copy.deepcopy(before)
    if mutation == 'beam-impact':
        after['impact_rules']['validation_test/accelerator/benchmark_beam_transfer_mex_python.py'].remove(
            'tests/test_beam_transfer_mex_pybind_validation.py')
    elif mutation == 'empty-native':
        after['profiles']['native-build-contracts']['paths'] = []
    else:
        after['profiles']['native-build-contracts']['max_elapsed_seconds'] = 0.001
    blobs[TIER] = json.dumps(before), json.dumps(after)
    assert decide(blobs, data)
    assert decide(blobs, data, event='push')
    assert not impact.ams_tier_projection_unchanged(*blobs[TIER])


@pytest.mark.parametrize('package', ['radia_simulink', 'radia_mcp', 'radia_mcp/radia_ngsolve'])
@pytest.mark.parametrize('version_change', [False, True])
def test_review_parent_initializers_are_real_execution_inputs(package, version_change):
    blobs, data = documents()
    initializer = ('packages/radia-simulink/src/' if package == 'radia_simulink'
                   else 'packages/radia-mcp/src/') + package + '/__init__.py'
    if version_change:
        blobs = {'pyproject.toml': ('[project]\nversion="1"\n', '[project]\nversion="2"\n')}
    else:
        blobs = {}
    # Supply actual initializer blobs: a missing blob must not be the reason
    # this regression selects AMS.
    blobs[initializer] = ('__version__ = "1"\n',
                          '__version__ = "1"\nimport numpy\nnumpy.seterr(all="raise")\n')
    assert decide(blobs, data)



def test_version_exemption_does_not_hide_a_numerical_pyproject_reader(tmp_path, monkeypatch):
    _, data = documents()
    (tmp_path / 'tools').mkdir()
    (tmp_path / 'tests').mkdir()
    shutil.copy2(ROOT / 'tools/ci_ams_inputs.py', tmp_path / 'tools/ci_ams_inputs.py')
    (tmp_path / OWNERS).write_text(json.dumps(data), encoding='utf-8')
    (tmp_path / TIER).write_text(json.dumps({'profiles': {
        name: {'paths': ['tests/test_consumer.py']} for name in
        ('sparsesolv', 'solver-numerics', 'solver-heavy', 'native-build-contracts')}}), encoding='utf-8')
    (tmp_path / 'tests/test_consumer.py').write_text(
        "from pathlib import Path\nPath('pyproject.toml').read_text()\n", encoding='utf-8')
    before, after = '[project]\nversion="1"\n', '[project]\nversion="2"\n'
    (tmp_path / 'pyproject.toml').write_text(after, encoding='utf-8')
    monkeypatch.setattr(impact, 'ROOT', tmp_path)
    assert decide({'pyproject.toml': (before, after)}, data)
