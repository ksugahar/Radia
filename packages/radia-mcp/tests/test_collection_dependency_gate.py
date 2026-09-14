"""Regression checks for the lightweight-matrix collection gate."""

import conftest


def test_solver_free_urn_contract_is_collectable_in_minimal_ci(monkeypatch):
    """CSV boundary coverage must not disappear behind optional solver imports."""
    monkeypatch.setattr(conftest, "_FORCE_MINIMAL", True)
    conftest._PROJECT_IMPORT_CACHE.clear()
    source = (conftest._TEST_ROOT / "test_urn_fit_contract.py").read_text(
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
