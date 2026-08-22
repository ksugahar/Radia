"""Regression tests for Cubit's startup hook and GUI ownership boundary."""

from __future__ import annotations

from pathlib import Path

from tools.audit_pyside6_only import check_deployed_panel_source


ROOT = Path(__file__).resolve().parents[1]
REGISTER_TOOLBAR = ROOT / "src" / "radia" / "panels" / "register_toolbar.py"


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
