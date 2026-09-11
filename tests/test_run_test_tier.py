"""Exercise the checked test-tier selector without running native solvers."""

from __future__ import annotations

import importlib.util
from pathlib import Path

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


def test_manifest_change_checks_all_registered_impacts():
    runner = runner_module()
    assert runner.select_impact_tests([], ['tests/test_tier_manifest.json']) == runner.select_impact_tests([], None)


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
