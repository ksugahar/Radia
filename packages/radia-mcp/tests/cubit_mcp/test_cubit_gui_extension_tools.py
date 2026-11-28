"""In-Cubit GUI extension knowledge and generators (toolbar, dialog, docs).

Guards: the toolbar/C++ SDK guides and cubit_docs panel_*/cpp_sdk_* routes
stay reachable, generated scripts follow the in-Cubit conventions, and the
scaffold refuses to overwrite.
Assumes: no Cubit; pure text and file checks.
"""
import ast
import json

import pytest

from radia_mcp.cubit import server
from radia_mcp.cubit.knowledge.cpp_sdk import get_cpp_sdk_documentation
from radia_mcp.cubit.knowledge.custom_toolbar import (
    TOOLBAR_EXTERNAL_JOBS,
    TOOLBAR_MENU_BAR_RECIPE,
)
from radia_mcp.cubit.knowledge.panel_conventions import LABEL_GUIDE, PANEL_CONVENTIONS


@pytest.mark.parametrize("topic,body", [
    ("menu_bar", TOOLBAR_MENU_BAR_RECIPE),
    ("external_jobs", TOOLBAR_EXTERNAL_JOBS),
])
def test_toolbar_guide_serves_menu_and_job_recipes(topic, body):
    assert server.cubit_toolbar_guide(topic) == body
    assert body in server.cubit_toolbar_guide("all")


def test_menu_recipe_compiles_and_keeps_its_lifecycle():
    code = TOOLBAR_MENU_BAR_RECIPE.split("```python\n", 1)[1].split("```", 1)[0]
    tree = ast.parse(code)
    funcs = {n.name for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}
    assert {"find_claro", "install_menu", "remove_menu", "_deferred_attach"} <= funcs, funcs
    assert "aboutToQuit.connect(remove_menu)" in code
    assert "QTimer.singleShot(0" in code
    assert ".lower() == \"claro\"" in code


def test_external_job_recipe_is_async_and_isolated():
    assert '"-I", "-B"' in TOOLBAR_EXTERNAL_JOBS
    launch = TOOLBAR_EXTERNAL_JOBS.split("## Launch", 1)[1].split("## Poll", 1)[0]
    code = launch.split("```python\n", 1)[1].split("```", 1)[0]
    ast.parse(code)
    assert "Popen(" in code and ".wait(" not in code and ".communicate(" not in code


def test_generated_dialog_follows_in_cubit_conventions():
    script = server.cubit_generate_dialog(
        "Brick",
        json.dumps([{"name": "size", "label": "Size", "default": "1.0", "kind": "float"}]),
        json.dumps(["create brick x {size}"]),
    )
    ast.parse(script)
    assert script.startswith("#!python")
    assert "\nimport cubit" not in script
    assert "_coreform_cubit" in script
    assert "objectName().lower() == 'claro'" in script


def test_scaffold_toolbar_writes_once(tmp_path):
    buttons = json.dumps([{"label": "Mesh", "script_basename": "mesh",
                           "body": "cubit.cmd('mesh volume all')"}])
    out = server.cubit_scaffold_toolbar("tb", buttons, str(tmp_path))
    assert out.startswith("Toolbar scaffold written to"), out
    assert (tmp_path / "tb" / "toolbar.xml").is_file()
    again = server.cubit_scaffold_toolbar("tb", buttons, str(tmp_path))
    assert again.startswith("ERROR"), again


@pytest.mark.parametrize("topic,expected", [
    ("panel_conventions", PANEL_CONVENTIONS),
    ("panel_labels", LABEL_GUIDE),
    ("cpp_sdk", get_cpp_sdk_documentation("overview")),
    ("cpp_sdk_custom_commands", get_cpp_sdk_documentation("custom_commands")),
    ("cpp_sdk_command_plugin", get_cpp_sdk_documentation("command_plugin")),
])
def test_cubit_docs_routes_gui_extension_topics(topic, expected):
    assert server.cubit_docs(topic) == expected


def test_gui_extension_tools_are_classified():
    tools = server.mcp._tool_manager
    for name in ("cubit_toolbar_guide", "cubit_cpp_sdk_guide", "cubit_generate_dialog"):
        assert tools.get_tool(name).annotations.readOnlyHint is True, name
    assert "cubit_scaffold_toolbar" in server._WRITING_TOOLS
