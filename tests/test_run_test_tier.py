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


@pytest.mark.parametrize('change', ['profile', 'schema', 'missing_rules'])
def test_structural_manifest_change_stays_broad(change, monkeypatch, tmp_path):
    runner = runner_module()
    current = json.loads(runner.MANIFEST.read_text(encoding='utf-8'))
    previous = copy.deepcopy(current)
    if change == 'profile':
        previous['profiles']['fast-contracts']['max_elapsed_seconds'] = 59
    elif change == 'schema':
        previous['schema'] = 'other-schema'
    else:
        del previous['impact_rules']
    manifest = tmp_path / 'manifest.json'
    manifest.write_text(json.dumps(current), encoding='utf-8')
    monkeypatch.setattr(runner, 'MANIFEST', manifest)
    assert runner.select_impact_tests(
        [], ['tests/test_tier_manifest.json'], previous_manifest=previous,
    ) == runner.select_impact_tests([], None)


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
