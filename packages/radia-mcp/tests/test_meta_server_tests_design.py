"""Read-only tests-design knowledge contracts.

Guards: Missing authoring sections, mutable shared responses, silent selector
fallback and a public MCP schema/dispatch that differs from the guide.
Assumes: Packaged pure-Python knowledge and the real MCP SDK. These are source
and interface contracts; prose review and numerical lock certification are
independent, and no solver/native execution is claimed.
"""

import asyncio
import json
import subprocess
import sys

import pytest

from radia_mcp.meta.knowledge.tests_design import get_tests_design_guide


def test_sections_are_complete_and_responses_do_not_share_mutable_state():
    overview = get_tests_design_guide()
    sections = [row["id"] for row in overview["sections"]]
    assert set(sections) == {
        "overview", "feature", "quantity", "lock", "declarations", "mutation",
        "coverage", "contracts", "review", "regeneration", "checklist", "adoption", "inference", "layers",
    }, sections
    for section in sections:
        response = get_tests_design_guide(section)
        assert response["section"] == section, response
        assert response["content"] and response["title"], response
        assert response["read_only"] is True, response
    overview["sections"].clear()
    overview["related_tools"].append("caller-owned-mutation")
    fresh = get_tests_design_guide()
    assert len(fresh["sections"]) == len(sections), fresh
    assert "caller-owned-mutation" not in fresh["related_tools"], fresh


def test_unknown_section_fails_loudly_without_overview_fallback():
    with pytest.raises(ValueError, match="expected one of"):
        get_tests_design_guide("unrecognized-section")


def test_real_mcp_schema_and_dispatch_expose_the_read_only_guide():
    from radia_mcp.meta import catalog
    from radia_mcp.meta.server import mcp

    async def exercise():
        tools = await mcp.list_tools()
        tool = next(item for item in tools if item.name == "tests_design_guide")
        assert tool.inputSchema["properties"]["section"]["default"] == "overview"
        for arguments, section in (({}, "overview"), ({"section": "mutation"}, "mutation")):
            result = await mcp.call_tool("tests_design_guide", arguments)
            content = result[0] if isinstance(result, tuple) else result
            actual = json.loads(content[0].text)
            assert actual == get_tests_design_guide(section), actual

    asyncio.run(exercise())
    assert "tests_design_guide" in catalog.get("meta")["primary_tools"]


def test_packaged_guide_imports_without_numerical_runtime():
    # The real child interpreter and installed package are used, not fake SDKs
    # or handwritten numerical-library objects.
    script = """
import json, sys
from radia_mcp.meta.knowledge.tests_design import get_tests_design_guide
assert not ({'radia', 'ngsolve', 'torch', 'numpy'} & set(sys.modules))
print(json.dumps(get_tests_design_guide('checklist')))
"""
    process = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True,
        encoding="utf-8", check=True,
    )
    assert json.loads(process.stdout)["section"] == "checklist"
