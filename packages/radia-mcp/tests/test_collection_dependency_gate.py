"""Regression checks for the lightweight-matrix collection gate."""

import conftest
import pytest


@pytest.mark.parametrize("filename", [
    "test_urn_fit_contract.py", "test_paper_writing_conclusion_first_use.py",
    "test_paper_writing_sentence_endings.py",
])
def test_lightweight_contract_is_collectable_in_minimal_ci(monkeypatch, filename):
    """Pure input/text checks must not disappear behind optional dependencies."""
    monkeypatch.setattr(conftest, "_FORCE_MINIMAL", True)
    conftest._PROJECT_IMPORT_CACHE.clear()
    source = (conftest._TEST_ROOT / filename).read_text(
        encoding="utf-8"
    )
    assert "numpy" in conftest._MINIMAL_BASELINE
    for module in conftest._imported_modules(source):
        assert not conftest._module_absent(module), module
        assert not conftest._project_import_has_absent_dependency(module), module


def test_project_dependency_probe_does_not_import_server(monkeypatch):
    """Dependency discovery must not instantiate FastMCP during collection."""

    conftest._PROJECT_IMPORT_CACHE.clear()

    def fail_runtime_import(module_name):
        raise AssertionError(f"unexpected runtime import: {module_name}")

    monkeypatch.setattr(conftest.importlib, "import_module", fail_runtime_import)

    assert not conftest._project_import_has_absent_dependency(
        "radia_mcp.accelerator.server"
    )
