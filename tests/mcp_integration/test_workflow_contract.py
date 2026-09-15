from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_PATH = ROOT / ".github/workflows/radia-mcp-matrix.yml"


def test_minimum_sdk_workflow_derives_floor_from_package_metadata():
    workflow = WORKFLOW_PATH.read_text(encoding="utf-8")
    assert 'pip install -e "packages/radia-mcp[maintenance]"' in workflow
    assert 'minimum-sdk:' in workflow
    assert '"mcp==$SDK_MIN"' in workflow
    assert 'RADIA_MCP_EXPECTED_SDK' in workflow

def test_workflow_discovers_selected_tests_from_the_test_directory():
    workflow = WORKFLOW_PATH.read_text(encoding="utf-8")

    assert 'env["RADIA_MCP_CI_SELECTION_JSON"]' in workflow
    assert '"-m", "not xval and not slow", "tests"' in workflow
    package_step = workflow.split('package_root = Path("packages/radia-mcp")', 1)[1].split("\n          PY", 1)[0]
    assert '*targets' not in package_step
    assert '"--confcutdir=tests/mcp_integration", *targets' in workflow
