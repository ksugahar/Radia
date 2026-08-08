"""Regression checks for the lightweight-matrix collection gate."""

import conftest


def test_project_dependency_probe_does_not_import_server(monkeypatch):
    """Dependency discovery must not instantiate FastMCP during collection."""

    conftest._PROJECT_IMPORT_CACHE.clear()

    def fail_runtime_import(module_name):
        raise AssertionError(f"unexpected runtime import: {module_name}")

    monkeypatch.setattr(conftest.importlib, "import_module", fail_runtime_import)

    assert not conftest._project_import_has_absent_dependency(
        "radia_mcp.accelerator.server"
    )
