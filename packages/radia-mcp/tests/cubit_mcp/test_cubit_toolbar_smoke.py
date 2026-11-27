"""Guards: the toolbar identifiers, shipped definitions and probe payloads agree,
and the GUI smoke judges only what the Cubit-side probe reports.

Assumes: offline source contracts; ``subprocess.Popen`` is replaced, no Cubit
GUI starts.  Cubit owns Qt and the optional live GUI probe.
"""
from __future__ import annotations

import ast
import copy
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from radia_mcp.cubit import toolbar_smoke as SMOKE

GUI = Path(SMOKE.__file__).resolve().parent / "gui"


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
        "toolbar_menu_has_entry": True,
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
        (lambda data: data.update(toolbar_visible_region_nonempty=False),
         "no visible screen region"),
        (lambda data: data["action_visible"].update({"Export Netgen (.vol)...": False}),
         "action is not visible: Export Netgen"),
        (lambda data: data["action_enabled"].update({"Export Netgen (.vol)...": False}),
         "action is not enabled: Export Netgen"),
        (lambda data: data.update(toolbar_menu_has_entry=False),
         "absent from Cubit's toolbar menu"),
        (lambda data: data.update(toolbar_owner="QToolBar"), "not owned by Coreform"),
        (lambda data: data.update(schema="other.toolbar-probe.v2"),
         "unexpected probe schema"),
    ],
)
def test_display_contract_rejects_missing_hidden_or_duplicate_ui(mutation, expected_issue):
    payload = copy.deepcopy(_healthy_payload())
    mutation(payload)
    issues = SMOKE.validate_probe_result(payload)
    assert any(expected_issue in issue for issue in issues), issues


def test_display_contract_rejects_missing_or_extra_actions():
    payload = _healthy_payload()
    payload["toolbar_actions"] = payload["toolbar_actions"] + ["Gmsh (.msh)"]
    assert any("toolbar actions differ" in i for i in SMOKE.validate_probe_result(payload))
    payload = _healthy_payload()
    payload["toolbar_actions"] = []
    assert any("toolbar actions differ" in i for i in SMOKE.validate_probe_result(payload))


def test_cubit_probe_checks_runtime_visibility_and_enabled_state():
    source = (GUI / "toolbar_probe.py").read_text(encoding="utf-8")
    for required in ("main.isVisible()", "toolbar.isVisible()", "toolbar.visibleRegion().isEmpty()",
                     "action.isVisible()", "action.isEnabled()", "toolbar.toggleViewAction()",
                     'cubit.cmd("exit 0")', 'temporary_result = RESULT_PATH + ".tmp"',
                     "os.replace(temporary_result, RESULT_PATH)"):
        assert required in source, required
    for forbidden in ("main.menuBar()", "main.createPopupMenu()", "import radia"):
        assert forbidden not in source, forbidden


def test_probe_and_runner_share_schema_title_and_environment():
    probe = ast.parse((GUI / "toolbar_probe.py").read_text(encoding="utf-8"))
    constants = {t.id: ast.literal_eval(n.value) for n in probe.body if isinstance(n, ast.Assign)
                 for t in n.targets if isinstance(t, ast.Name) and t.id in ("SCHEMA", "EXPECTED_ACTIONS")}
    assert constants == {"SCHEMA": SMOKE.PROBE_SCHEMA, "EXPECTED_ACTIONS": SMOKE.EXPECTED_ACTIONS}
    source = (GUI / "toolbar_probe.py").read_text(encoding="utf-8")
    assert f'== "{SMOKE.TOOLBAR_TITLE}"' in source
    smoke_source = Path(SMOKE.__file__).read_text(encoding="utf-8")
    for variable in ("RADIA_MCP_CUBIT_TOOLBAR_PROBE_RESULT", "RADIA_MCP_CUBIT_TOOLBAR_PROBE_TIMEOUT"):
        assert variable in source and variable in smoke_source, variable
    assert SMOKE._probe_path().is_file()


@pytest.mark.parametrize('returncode', [0, 1, 2, -1073740791])
@pytest.mark.parametrize('executable', ['cubit.exe', 'clarox.exe'])
def test_good_probe_requires_clean_process_exit(monkeypatch, tmp_path, returncode, executable):
    class Process:
        def poll(self):
            return returncode

        def wait(self, timeout):
            return returncode

    def launch(*args, **kwargs):
        expected = [executable]
        if executable == 'clarox.exe':
            expected += ['-prompt', 'Cubit']
        expected += ['-nojournal']
        if executable != 'clarox.exe':
            expected += ['-commandplugindir', 'plugins']
        expected += [str(tmp_path / 'probe_bootstrap.py')]
        assert args[0] == expected
        assert kwargs['stdout'].name == str(tmp_path / 'launcher.log')
        assert kwargs['stderr'] == SMOKE.subprocess.STDOUT
        assert kwargs['env']['RADIA_MCP_CUBIT_TOOLBAR_PROBE_RESULT'] == str(tmp_path / 'result.json')
        (tmp_path / 'result.json').write_text(json.dumps(_healthy_payload()))
        return Process()

    monkeypatch.setattr(SMOKE, '_cubit_pids', lambda: set())
    monkeypatch.setattr(SMOKE.subprocess, 'Popen', launch)
    result = SMOKE._run_one(Path(executable), tmp_path, 1)
    assert result['ok'] is (returncode == 0), result


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


def test_shipped_toolbar_identifiers_match_both_probe_consumers():
    template = GUI / "cubit_toolbar/toolbars/netgen_export_toolbar.ttb.tmpl"
    root = ET.parse(template).getroot()
    assert root.attrib["name"] == SMOKE.TOOLBAR_TITLE, root.attrib
    actions = [node.attrib["name"] for node in root.iter("WAction")]
    assert actions == SMOKE.EXPECTED_ACTIONS, actions
    menu = ast.parse((GUI / "cubit_export_menu.py").read_text(encoding="utf-8"))
    names = {t.id: n.value.value for n in menu.body if isinstance(n, ast.Assign)
             for t in n.targets if isinstance(t, ast.Name) and isinstance(n.value, ast.Constant)}
    assert [names["_MENU_LABEL"]] == actions, (names.get("_MENU_LABEL"), actions)
    scripts = [Path(node.text.split("/")[-1]) for node in root.iter("filename")]
    assert scripts == [Path("export_netgen.py")], scripts
    tree = ast.parse((GUI / "cubit_toolbar/scripts" / scripts[0]).read_text(encoding="utf-8"))
    calls = [getattr(c.func, "id", "") for c in ast.walk(tree) if isinstance(c, ast.Call)]
    assert "launch_netgen_export" in calls, calls
    installed = sorted(p.name for p in (GUI / "cubit_toolbar/scripts").glob("*.py"))
    assert installed == sorted(s.name for s in scripts), (installed, scripts)
    icons = [Path(node.text).name for node in root.iter("icon")]
    assert all((GUI / "cubit_toolbar/icons" / icon).is_file() for icon in icons), icons
