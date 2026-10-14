"""Heavy validation uses explicit isolated job roles rather than hostnames."""
import pytest
from validation_test import conftest as validation_conftest

@pytest.mark.parametrize("role", ["cpu_validation", "gpu_evaluation"])
def test_isolated_validation_roles_are_recognized(monkeypatch, role):
    monkeypatch.setenv("RADIA_VALIDATION_HOST_ROLE", role)
    assert validation_conftest._is_compute_host()

@pytest.mark.parametrize("role", ["", "ci", "developer", "retired_worker"])
def test_other_roles_are_rejected(monkeypatch, role):
    monkeypatch.setenv("RADIA_VALIDATION_HOST_ROLE", role)
    assert not validation_conftest._is_compute_host()

def test_validation_requires_explicit_opt_in(monkeypatch):
    monkeypatch.delenv("RADIA_VALIDATION_HOST_ROLE", raising=False)
    assert not validation_conftest._is_compute_host()
