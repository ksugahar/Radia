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
