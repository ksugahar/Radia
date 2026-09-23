"""Exercise the checked test-tier selector without running native solvers."""

from __future__ import annotations

import copy
import importlib.util
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "tools" / "run_test_tier.py"


def runner_module():
    spec = importlib.util.spec_from_file_location("run_test_tier", SOURCE)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_native_smoke_extends_fast_contracts_without_duplicates():
    runner = runner_module()
    fast_paths, fast_budget = runner.load_profile("fast-contracts")
    native_paths, native_budget = runner.load_profile("native-smoke")

    assert set(fast_paths) <= set(native_paths)
    assert len(native_paths) == len(set(native_paths))
    assert fast_budget == 60
    assert native_budget == 120


def test_unknown_test_tier_fails_before_pytest_is_started():
    runner = runner_module()
    with pytest.raises(ValueError, match="unknown test tier"):
        runner.load_profile("does-not-exist")


def test_changed_release_tool_selects_its_regressions():
    runner = runner_module()
    paths, _ = runner.load_profile('fast-contracts')
    selected = runner.select_impact_tests(paths, ['tools/release_quad.py'])
    assert 'tests/test_release_quad_state.py' in selected
    assert 'tests/test_release_quad_editable_source.py' in selected
    assert len(selected) == len(set(selected))


def test_changed_test_selects_itself_but_not_unrelated_regressions():
    runner = runner_module()
    selected = runner.select_impact_tests([], ['tests/test_ci_preflight_mdx.py'])
    assert selected == ['tests/test_ci_preflight_mdx.py']
    assert runner.select_impact_tests([], ['docs/intro.md']) == []


def test_unknown_base_selects_all_registered_impacts():
    runner = runner_module()
    selected = runner.select_impact_tests([], None)
    assert 'tests/test_release_quad_state.py' in selected
    assert 'tests/test_ci_preflight_mdx.py' in selected


def test_manifest_change_without_comparison_checks_all_registered_impacts():
    runner = runner_module()
    assert runner.select_impact_tests([], ['tests/test_tier_manifest.json']) == runner.select_impact_tests([], None)


def test_additive_manifest_selects_new_rules_and_ordinary_impacts(monkeypatch, tmp_path):
    runner = runner_module()
    previous = json.loads(runner.MANIFEST.read_text(encoding='utf-8'))
    current = copy.deepcopy(previous)
    # Source need not change: adding a mapping must exercise its tests anyway.
    current['impact_rules']['src/new_gate.py'] = ['tests/test_ci_monitor.py']
    manifest = tmp_path / 'manifest.json'
    manifest.write_text(json.dumps(current), encoding='utf-8')
    monkeypatch.setattr(runner, 'MANIFEST', manifest)
    selected = runner.select_impact_tests(
        ['tests/test_test_tier_policy.py'],
        ['tests/test_tier_manifest.json', 'tools/ci_preflight_mdx.py'],
        previous_manifest=previous,
    )
    assert set(selected) == {
        'tests/test_test_tier_policy.py', 'tests/test_ci_monitor.py',
        'tests/test_ci_preflight_mdx.py',
    }


@pytest.mark.parametrize('change', ['remove_rule', 'replace_tests', 'append_test'])
def test_changed_manifest_rule_selects_only_old_and_new_tests(change, monkeypatch, tmp_path):
    runner = runner_module()
    current = json.loads(runner.MANIFEST.read_text(encoding='utf-8'))
    previous = copy.deepcopy(current)
    source = next(iter(current['impact_rules']))
    old_tests = list(previous['impact_rules'][source])
    if change == 'remove_rule':
        del current['impact_rules'][source]
        expected = set(old_tests)
    elif change == 'replace_tests':
        current['impact_rules'][source] = ['tests/test_ci_monitor.py']
        expected = {*old_tests, 'tests/test_ci_monitor.py'}
    else:
        current['impact_rules'][source].append('tests/test_ci_monitor.py')
        expected = {*old_tests, 'tests/test_ci_monitor.py'}
    manifest = tmp_path / 'manifest.json'
    manifest.write_text(json.dumps(current), encoding='utf-8')
    monkeypatch.setattr(runner, 'MANIFEST', manifest)
    assert set(runner.select_impact_tests(
        [], ['tests/test_tier_manifest.json'], previous_manifest=previous,
    )) == expected


@pytest.mark.parametrize('change', ['profile', 'schema', 'missing_rules', 'inheritance',
                                  'description', 'new_profile_key', 'profile_added',
                                  'invalid_paths', 'duplicate_paths'])
def test_structural_manifest_change_stays_broad(change, monkeypatch, tmp_path):
    runner = runner_module()
    current = json.loads(runner.MANIFEST.read_text(encoding='utf-8'))
    previous = copy.deepcopy(current)
    if change == 'profile':
        previous['profiles']['fast-contracts']['max_elapsed_seconds'] = 59
    elif change == 'schema':
        previous['schema'] = 'other-schema'
    elif change == 'missing_rules':
        del previous['impact_rules']
    elif change == 'inheritance':
        previous['profiles']['native-smoke']['extends'] = 'other-profile'
    elif change == 'description':
        previous['profiles']['fast-contracts']['description'] = 'other-description'
    elif change == 'new_profile_key':
        previous['profiles']['native-smoke']['unknown_policy'] = True
    elif change == 'profile_added':
        previous['profiles']['another-profile'] = {'paths': []}
    elif change == 'invalid_paths':
        previous['profiles']['fast-contracts']['paths'] = 'not-a-list'
    else:
        previous['profiles']['fast-contracts']['paths'] *= 2
    manifest = tmp_path / 'manifest.json'
    manifest.write_text(json.dumps(current), encoding='utf-8')
    monkeypatch.setattr(runner, 'MANIFEST', manifest)
    assert runner.select_impact_tests(
        [], ['tests/test_tier_manifest.json'], previous_manifest=previous,
    ) == runner.select_impact_tests([], None)


@pytest.mark.parametrize('change', ['add', 'remove', 'replace', 'reorder'])
def test_profile_membership_only_selects_changed_tests(change, monkeypatch, tmp_path):
    runner = runner_module()
    current = json.loads(runner.MANIFEST.read_text(encoding='utf-8'))
    previous = copy.deepcopy(current)
    paths = current['profiles']['fast-contracts']['paths']
    before = set(paths)
    if change == 'add':
        paths.append('tests/test_ci_preflight_mdx.py')
    elif change == 'remove':
        paths.remove('tests/test_ci_monitor.py')
    elif change == 'replace':
        paths.remove('tests/test_ci_monitor.py')
        paths.append('tests/test_ci_preflight_mdx.py')
    else:
        paths.reverse()
    manifest = tmp_path / 'manifest.json'
    manifest.write_text(json.dumps(current), encoding='utf-8')
    monkeypatch.setattr(runner, 'MANIFEST', manifest)
    selected = runner.select_impact_tests([], ['tests/test_tier_manifest.json'], previous_manifest=previous)
    assert set(selected) == before ^ set(paths)


@pytest.mark.parametrize('renamed', [False, True])
def test_removed_test_is_not_passed_to_pytest(renamed, monkeypatch, tmp_path):
    runner = runner_module()
    previous = {'schema': 'v1', 'impact_rules': {}, 'profiles': {'fast': {'paths': ['old.py']}}}
    current = copy.deepcopy(previous)
    current['profiles']['fast']['paths'] = ['new.py'] if renamed else []
    if renamed:
        (tmp_path / 'new.py').write_text('')
    manifest = tmp_path / 'manifest.json'
    manifest.write_text(json.dumps(current), encoding='utf-8')
    monkeypatch.setattr(runner, 'MANIFEST', manifest)
    monkeypatch.setattr(runner, 'ROOT', tmp_path)
    selected = runner.select_impact_tests([], ['tests/test_tier_manifest.json'], previous_manifest=previous)
    assert selected == (['new.py'] if renamed else [])


def test_missing_test_still_referenced_by_current_manifest_fails(monkeypatch, tmp_path):
    runner = runner_module()
    previous = {'schema': 'v1', 'impact_rules': {}, 'profiles': {'fast': {'paths': ['old.py']}}}
    current = copy.deepcopy(previous)
    current['profiles']['fast']['paths'] = []
    current['impact_rules']['source.py'] = ['old.py']
    manifest = tmp_path / 'manifest.json'
    manifest.write_text(json.dumps(current), encoding='utf-8')
    monkeypatch.setattr(runner, 'MANIFEST', manifest)
    monkeypatch.setattr(runner, 'ROOT', tmp_path)
    with pytest.raises(ValueError, match='missing test: old.py'):
        runner.select_impact_tests([], ['tests/test_tier_manifest.json'], previous_manifest=previous)


def test_new_missing_profile_member_fails_before_pytest(monkeypatch, tmp_path):
    runner = runner_module()
    previous = {'schema': 'v1', 'impact_rules': {}, 'profiles': {'fast': {'paths': []}}}
    current = copy.deepcopy(previous)
    current['profiles']['fast']['paths'] = ['typo.py']
    manifest = tmp_path / 'manifest.json'
    manifest.write_text(json.dumps(current), encoding='utf-8')
    monkeypatch.setattr(runner, 'MANIFEST', manifest)
    monkeypatch.setattr(runner, 'ROOT', tmp_path)
    with pytest.raises(ValueError, match='missing test: typo.py'):
        runner.select_impact_tests([], ['tests/test_tier_manifest.json'], previous_manifest=previous)


@pytest.mark.parametrize('changed', [
    'src/radia/vim/_nonlinear.py', 'src/radia/panels/samples/em_sample_bh.txt',
    'validation_test/feec/bh_saturation_audit.py',
    'validation_test/feec/test_hdiv_vim_energy_newton.py',
    'validation_test/feec/run_energy_newton_audit.py',
    'validation_test/feec/results/bh_saturation_20260914/qualified/result.json',
    '.github/workflows/radia-fast.yml', 'tests/test_hdiv_bh_table_saturation_contract.py',
])
def test_bh_dependency_changes_select_bh_contract(changed):
    runner = runner_module()
    base, _ = runner.load_profile('fast-contracts')
    contract = 'tests/test_hdiv_bh_table_saturation_contract.py'
    assert contract not in base
    assert contract in runner.select_impact_tests(base, [changed])
    native, _ = runner.load_profile('native-smoke')
    assert contract in native


def test_fast_workflow_checks_dependencies_without_unrelated_impacts():
    runner = runner_module()
    base, _ = runner.load_profile('fast-contracts')
    selected = runner.select_impact_tests(base, ['.github/workflows/radia-fast.yml'])
    assert set(base) <= set(selected)
    assert 'tests/test_hdiv_bh_table_saturation_contract.py' in selected
    assert 'tests/test_simulink_release_package.py' not in selected
    assert 'tests/test_build_failure_propagation.py' not in selected
    assert 'tests/test_hdiv_bh_table_saturation_contract.py' not in runner.select_impact_tests(base, ['docs/intro.md'])


@pytest.mark.parametrize('payload,code', [(b'not json', 0), (b'[]', 0), (b'', 1),
                                         (b'\xff', 0)])
def test_unavailable_previous_manifest_remains_unknown(monkeypatch, payload, code):
    runner = runner_module()
    monkeypatch.setattr(runner.subprocess, 'run',
                        lambda *args, **kwargs: SimpleNamespace(stdout=payload, returncode=code))
    assert runner.read_previous_manifest('base') is None


def test_cli_passes_exact_base_manifest_to_selection(monkeypatch):
    runner = runner_module()
    previous = json.loads(runner.MANIFEST.read_text(encoding='utf-8'))
    calls = []

    def fake_run(command, **kwargs):
        calls.append(command)
        if 'diff' in command:
            return SimpleNamespace(returncode=0, stdout=b'tests/test_tier_manifest.json\0')
        if 'show' in command:
            assert command[-1] == 'exact-base:tests/test_tier_manifest.json'
            return SimpleNamespace(returncode=0, stdout=json.dumps(previous).encode())
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(runner.subprocess, 'run', fake_run)
    assert runner.main(['--since', 'exact-base', '--collect-only']) == 0
    selected = calls[-1]
    assert 'tests/test_ci_preflight_mdx.py' not in selected
    assert 'tests/test_run_test_tier.py' in selected


def test_previous_manifest_timeout_remains_unknown(monkeypatch):
    runner = runner_module()

    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired('git', 30)

    monkeypatch.setattr(runner.subprocess, 'run', timeout)
    assert runner.read_previous_manifest('base') is None


@pytest.mark.parametrize('changed,expected', [
    ('matlab/+radia/+sparsesolv/AMS.m', 'tests/test_sparsesolv_matlab_contract.py'),
    ('matlab/+radia/+sparsesolv/private/setup.m', 'tests/test_sparsesolv_matlab_contract.py'),
    ('matlab/+radia/+python/sparsesolv.m', 'tests/test_matlab_python_parity_manifest.py'),
])
def test_directory_rules_select_descendant_changes(changed, expected):
    assert runner_module().select_impact_tests([], [changed]) == [expected]


@pytest.mark.parametrize('changed', [
    'matlab/+radia/+sparsesolv_extra/AMS.m',
    'matlab/+radia/+python_extra/sparsesolv.m',
    'tools/release_quad.py.bak',
])
def test_rules_do_not_select_similarly_named_siblings(changed):
    assert runner_module().select_impact_tests([], [changed]) == []


@pytest.mark.parametrize('changed', [
    '.github/workflows/ltspice-data-safety.yml',
    '.github/workflows/radia-optuna.yml',
    'tests/test_ltspice_ci_scope.py',
])
def test_ltspice_scope_regression_is_reachable_from_ci(changed):
    assert 'tests/test_ltspice_ci_scope.py' in runner_module().select_impact_tests([], [changed])


def test_adding_numerical_profile_does_not_expand_fast_lane(monkeypatch, tmp_path):
    runner = runner_module()
    current = json.loads(runner.MANIFEST.read_text(encoding='utf-8'))
    previous = copy.deepcopy(current)
    del previous['profiles']['solver-numerics']
    manifest = tmp_path / 'manifest.json'
    manifest.write_text(json.dumps(current), encoding='utf-8')
    monkeypatch.setattr(runner, 'MANIFEST', manifest)
    paths, _ = runner.load_profile('fast-contracts')
    assert runner.select_impact_tests(
        paths, ['tests/test_tier_manifest.json'], previous_manifest=previous,
        profile_name='fast-contracts') == paths


def test_scoped_native_lane_still_checks_changed_parent_profile():
    runner = runner_module()
    current = json.loads(runner.MANIFEST.read_text(encoding='utf-8'))
    previous = copy.deepcopy(current)
    previous['profiles']['fast-contracts']['paths'].remove('tests/test_ci_monitor.py')
    assert runner.changed_impact_tests(
        current, previous, profile_name='native-smoke') == {'tests/test_ci_monitor.py'}
