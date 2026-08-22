"""Regression tests for Cubit's startup hook and GUI ownership boundary."""

from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTER_TOOLBAR = ROOT / "src" / "radia" / "panels" / "register_toolbar.py"
AUDIT_SCRIPT = ROOT / "tools" / "audit_pyside6_only.py"

_AUDIT_SPEC = importlib.util.spec_from_file_location(
    "radia_audit_pyside6_only", AUDIT_SCRIPT
)
assert _AUDIT_SPEC is not None and _AUDIT_SPEC.loader is not None
_AUDIT_MODULE = importlib.util.module_from_spec(_AUDIT_SPEC)
_AUDIT_SPEC.loader.exec_module(_AUDIT_MODULE)
check_deployed_panel_source = _AUDIT_MODULE.check_deployed_panel_source


def test_deployment_audit_rejects_startup_from_another_checkout(tmp_path):
    expected = tmp_path / "current" / "src" / "radia" / "panels" \
        / "register_toolbar.py"
    foreign = tmp_path / "old-release" / "src" / "radia" / "panels" \
        / "register_toolbar.py"
    expected.parent.mkdir(parents=True)
    foreign.parent.mkdir(parents=True)
    expected.write_text("# current\n", encoding="utf-8")
    foreign.write_text("# stale\n", encoding="utf-8")

    startup = tmp_path / "radia_startup.py"
    startup.write_text(
        f"exec(open(r'{foreign.as_posix()}').read())\n",
        encoding="utf-8",
    )
    cubit_file = tmp_path / ".cubit"
    cubit_file.write_text(
        "## BEGIN radia toolbar\n"
        f'play "{startup.as_posix()}"\n'
        "## END radia toolbar\n",
        encoding="utf-8",
    )

    status, issues = check_deployed_panel_source(cubit_file, expected)
    assert status == "checked"
    assert len(issues) == 1
    assert "different checkout" in issues[0]
    assert expected.as_posix().lower() in issues[0]

    startup.write_text(
        f"exec(open(r'{expected.as_posix()}').read())\n",
        encoding="utf-8",
    )
    status, issues = check_deployed_panel_source(cubit_file, expected)
    assert status == "checked"
    assert issues == []


def test_cubit_startup_does_not_inject_a_qmenu():
    """``~/.cubit`` must leave persistent GUI ownership to Coreform."""
    source = REGISTER_TOOLBAR.read_text(encoding="utf-8")

    assert "def _install_radia_export_menu" not in source
    assert "radia_export_menu.install_menu()" not in source
    assert '"Export Mesh", "Radia Export"' in source
    assert "official WorkflowToolbar package" in source
