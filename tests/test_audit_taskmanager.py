from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


@pytest.fixture()
def audit_module(monkeypatch, tmp_path):
    script = Path(__file__).parents[1] / "tools" / "audit_taskmanager.py"
    spec = importlib.util.spec_from_file_location("audit_taskmanager_test", script)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "ROOT", tmp_path)
    return module


@pytest.mark.parametrize(
    ("relative_path", "expected"),
    [
        ("src/radia/solver.py", "helper"),
        ("src/radia/panels/calc_heat.py", "caller"),
        ("validation_test/heat/check_energy.py", "caller"),
        ("tests/test_heat.py", "caller"),
        ("docs/heat/notebook_helper.py", "caller"),
        ("examples/retired.py", None),
    ],
)
def test_classification_matches_repository_lanes(
    audit_module, tmp_path, relative_path, expected
):
    path = tmp_path / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("", encoding="utf-8")
    assert audit_module._classify(path) == expected


def test_caller_without_region_is_reported(audit_module, tmp_path):
    path = tmp_path / "tests" / "test_solver.py"
    path.parent.mkdir(parents=True)
    path.write_text(
        "def test_solve():\n"
        "    form = BilinearForm(space)\n"
        "    form.Assemble()\n",
        encoding="utf-8",
    )

    findings = audit_module._audit_caller(path)

    assert len(findings) == 1
    assert findings[0].kind == "caller-missing-wrap"


def test_mesh_constructor_candidate_is_not_skipped(audit_module, tmp_path):
    path = tmp_path / "tests" / "test_mesh.py"
    path.parent.mkdir(parents=True)
    path.write_text(
        "def test_mesh():\n"
        "    mesh = Mesh(make_mesh())\n",
        encoding="utf-8",
    )

    findings = audit_module._audit_caller(path)

    assert len(findings) == 1
    assert findings[0].snippet.startswith("[Mesh(...)]")


def test_caller_region_satisfies_minimum_gate(audit_module, tmp_path):
    path = tmp_path / "validation_test" / "solver.py"
    path.parent.mkdir(parents=True)
    path.write_text(
        "def main():\n"
        "    with TaskManager():\n"
        "        form = BilinearForm(space)\n"
        "        form.Assemble()\n",
        encoding="utf-8",
    )

    assert audit_module._audit_caller(path) == []


def test_helper_owned_region_is_reported(audit_module, tmp_path):
    path = tmp_path / "src" / "radia" / "solver.py"
    path.parent.mkdir(parents=True)
    path.write_text(
        "def assemble():\n"
        "    with TaskManager():\n"
        "        return operator.Assemble()\n",
        encoding="utf-8",
    )

    findings = audit_module._audit_helper(path)

    assert len(findings) == 1
    assert findings[0].kind == "helper-wraps"
