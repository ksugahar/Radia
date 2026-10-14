from pathlib import Path
import shlex

ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_PATH = ROOT / ".github/workflows/radia-mcp-matrix.yml"


def test_minimum_sdk_workflow_derives_floor_from_package_metadata():
    workflow = WORKFLOW_PATH.read_text(encoding="utf-8")
    install_commands = [shlex.split(line.strip()) for line in workflow.splitlines()
                        if line.strip().startswith("python -m pip install ")]
    assert any(
        ["-e", "packages/radia-mcp[maintenance]"] == tokens[i:i + 2]
        for tokens in install_commands for i in range(len(tokens) - 1)
    )
    assert 'minimum-sdk:' in workflow
    assert '"mcp==$SDK_MIN"' in workflow
    assert 'RADIA_MCP_EXPECTED_SDK' in workflow
