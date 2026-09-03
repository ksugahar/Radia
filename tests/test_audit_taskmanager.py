from __future__ import annotations

import importlib.util
import warnings
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
        ("src/radia/simulink/ih_operator_assembly.py", "caller"),
        ("validation_test/heat/check_energy.py", "caller"),
        ("tests/test_heat.py", "caller"),
        ("docs/heat/notebook_helper.py", "caller"),
        ("tests/_ngsolve_2606.py", "helper"),
        ("tests/axifem/_vol_mesh.py", "helper"),
        ("validation_test/cubit/cubit_202512_helpers.py", "helper"),
        ("validation_test/feec/conftest.py", "helper"),
        ("validation_test/ngsolve_matlab_parity/extended_catalog.py", "helper"),
        ("validation_test/stream_function/regcoil_fusion_helpers.py", "helper"),
        ("docs/electric_machine/planar_vim_motor_helpers.py", "helper"),
        ("validation_test/cubit/_ho_volume_worker.py", "caller"),
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


def test_pytest_module_may_request_shared_taskmanager_fixture(
    audit_module, tmp_path
):
    path = tmp_path / "validation_test" / "solver" / "test_parallel.py"
    path.parent.mkdir(parents=True)
    (path.parent / "conftest.py").write_text(
        "def ngsolve_taskmanager():\n"
        "    with TaskManager():\n"
        "        yield\n",
        encoding="utf-8",
    )
    path.write_text(
        "import pytest\n"
        "pytestmark = pytest.mark.usefixtures('ngsolve_taskmanager')\n"
        "def test_solve():\n"
        "    form = BilinearForm(space)\n"
        "    form.Assemble()\n",
        encoding="utf-8",
    )

    assert audit_module._audit_caller(path) == []


def test_fixture_marker_without_fixture_is_reported(audit_module, tmp_path):
    path = tmp_path / "tests" / "test_parallel.py"
    path.parent.mkdir(parents=True)
    path.write_text(
        "import pytest\n"
        "pytestmark = pytest.mark.usefixtures('ngsolve_taskmanager')\n"
        "form = BilinearForm(space)\n"
        "form.Assemble()\n",
        encoding="utf-8",
    )

    findings = audit_module._audit_caller(path)

    assert len(findings) == 1
    assert findings[0].kind == "caller-missing-wrap"


def test_shared_fixture_marker_does_not_exempt_standalone_generator(
    audit_module, tmp_path
):
    path = tmp_path / "validation_test" / "solver" / "generate_result.py"
    path.parent.mkdir(parents=True)
    path.write_text(
        "import pytest\n"
        "pytestmark = pytest.mark.usefixtures('ngsolve_taskmanager')\n"
        "form = BilinearForm(space)\n"
        "form.Assemble()\n",
        encoding="utf-8",
    )

    findings = audit_module._audit_caller(path)

    assert len(findings) == 1
    assert findings[0].kind == "caller-missing-wrap"


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


def test_syntax_warning_identifies_source_file(audit_module, tmp_path):
    path = tmp_path / "tests" / "test_invalid_escape.py"
    path.parent.mkdir(parents=True)
    path.write_text(
        '"""invalid \\i escape"""\nform = BilinearForm(space)\n',
        encoding="utf-8",
    )

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", SyntaxWarning)
        audit_module._audit_caller(path)

    assert caught
    assert caught[0].filename == str(path)
