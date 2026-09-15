from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_PATH = ROOT / ".github/workflows/radia-mcp-matrix.yml"


def test_minimum_sdk_workflow_derives_floor_from_package_metadata():
    workflow = WORKFLOW_PATH.read_text(encoding="utf-8")
    assert 'pip install -e "packages/radia-mcp[maintenance]"' in workflow
    assert 'minimum-sdk:' in workflow
    assert '"mcp==$SDK_MIN"' in workflow
    assert 'RADIA_MCP_EXPECTED_SDK' in workflow
