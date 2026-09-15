"""Regression tests for the real-Cubit Radia Export display contract."""

from __future__ import annotations

import copy
import importlib.util
import tomllib
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
PACKAGE = ROOT / "packages/cubit-mesh-export"
MODULE_PATH = PACKAGE / "src/cubit_mesh_export/toolbar_smoke.py"
SPEC = importlib.util.spec_from_file_location("cubit_toolbar_smoke_test", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
SMOKE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SMOKE)


def _healthy_payload():
    actions = list(SMOKE.EXPECTED_ACTIONS)
    return {
        "schema": SMOKE.PROBE_SCHEMA,
        "ok": True,
        "main_window_visible": True,
        "toolbar_count": 1,
        "toolbar_visible": True,
        "toolbar_visible_region_nonempty": True,
        "toolbar_size": [480, 32],
        "toolbar_actions": actions,
        "action_visible": {name: True for name in actions},
        "action_enabled": {name: True for name in actions},
        "toolbar_menu_has_radia_export": True,
        "toolbar_owner": "WorkflowToolbar",
    }


def test_display_contract_accepts_one_visible_complete_toolbar():
    assert SMOKE.validate_probe_result(_healthy_payload()) == []


@pytest.mark.parametrize(
    ("mutation", "expected_issue"),
    [
        (lambda data: data.update(toolbar_count=0), "exactly one"),
        (lambda data: data.update(toolbar_count=2), "exactly one"),
        (lambda data: data.update(toolbar_visible=False), "not visible"),
        (
            lambda data: data.update(toolbar_visible_region_nonempty=False),
            "no visible screen region",
        ),
        (
            lambda data: data["action_visible"].update({"GMSH (.msh)": False}),
            "action is not visible: GMSH",
        ),
        (
            lambda data: data["action_enabled"].update({"VTK (.vtk)": False}),
            "action is not enabled: VTK",
        ),
        (
            lambda data: data.update(toolbar_menu_has_radia_export=False),
            "absent from Cubit's toolbar menu",
        ),
        (
            lambda data: data.update(toolbar_owner="QToolBar"),
            "not owned by Coreform",
        ),
    ],
)
def test_display_contract_rejects_missing_hidden_or_duplicate_ui(
    mutation, expected_issue,
):
    payload = copy.deepcopy(_healthy_payload())
    mutation(payload)

    issues = SMOKE.validate_probe_result(payload)

    assert any(expected_issue in issue for issue in issues), issues


def test_display_contract_rejects_missing_or_reordered_actions():
    payload = _healthy_payload()
    payload["toolbar_actions"] = list(reversed(payload["toolbar_actions"]))

    issues = SMOKE.validate_probe_result(payload)

    assert any("toolbar actions differ" in issue for issue in issues)


def test_cubit_probe_checks_runtime_visibility_and_enabled_state():
    source = (
        PACKAGE / "src/cubit_mesh_export/cubit_gui/toolbar_probe.py"
    ).read_text(encoding="utf-8")

    assert "main.isVisible()" in source
    assert "toolbar.isVisible()" in source
    assert "toolbar.visibleRegion().isEmpty()" in source
    assert "action.isVisible()" in source
    assert "action.isEnabled()" in source
    assert "toolbar.toggleViewAction()" in source
    assert "main.menuBar()" not in source
    assert "main.createPopupMenu()" not in source
    assert 'cubit.cmd("exit 0")' in source
    assert 'temporary_result = RESULT_PATH + ".tmp"' in source
    assert "os.replace(temporary_result, RESULT_PATH)" in source


def test_gui_smoke_is_an_installed_release_gate():
    project = tomllib.loads((PACKAGE / "pyproject.toml").read_text(encoding="utf-8"))
    assert project["project"]["scripts"]["cubit-toolbar-smoke-test"] == (
        "cubit_mesh_export.toolbar_smoke:main"
    )
    radia = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert "cubit-toolbar-smoke-test" not in radia["project"]["scripts"]
    assert SMOKE._probe_path().is_file()

    release_source = (ROOT / "tools" / "release_quad.py").read_text(
        encoding="utf-8"
    )
    lab_start = release_source.index("def _deploy_lab():")
    lab_end = release_source.index("def _deploy_editable_remote", lab_start)
    lab_deploy = release_source[lab_start:lab_end]
    assert 'run(["cubit-toolbar-smoke-test", "--restarts", "2"])' in lab_deploy


@pytest.mark.parametrize('returncode', [0, 1, 2, -1073740791])
def test_good_probe_requires_clean_process_exit(monkeypatch, tmp_path, returncode):
    class Process:
        def poll(self):
            return returncode
        def wait(self, timeout):
            return returncode
    def launch(*args, **kwargs):
        assert args[0] == ['cubit.exe', '-nojournal', '-commandplugindir',
                           'plugins', str(tmp_path / 'probe_bootstrap.py')]
        assert kwargs['stdout'].name == str(tmp_path / 'launcher.log')
        assert kwargs['stderr'] == SMOKE.subprocess.STDOUT
        (tmp_path / 'result.json').write_text(json.dumps(_healthy_payload()))
        return Process()
    monkeypatch.setattr(SMOKE, '_cubit_pids', lambda: set())
    monkeypatch.setattr(SMOKE.subprocess, 'Popen', launch)
    result = SMOKE._run_one(Path('cubit.exe'), tmp_path, 1)
    assert result['ok'] is (returncode == 0)


def test_timeout_terminates_only_owned_process(monkeypatch, tmp_path):
    calls = []
    class Process:
        closed = False
        def poll(self):
            return 1 if self.closed else None
        def terminate(self):
            calls.append('terminate-owned')
            self.closed = True
        def wait(self, timeout):
            return 1
    process = Process()
    monkeypatch.setattr(SMOKE, '_cubit_pids', lambda: set())
    monkeypatch.setattr(SMOKE.subprocess, 'Popen', lambda *a, **k: process)
    result = SMOKE._run_one(Path('cubit.exe'), tmp_path, 0)
    assert not result['ok']
    assert calls == ['terminate-owned']
