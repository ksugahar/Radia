"""The two distributions ship one canonical numerical acceptance policy."""
from pathlib import Path


def test_vendored_gate_matches_single_canonical_source():
    root = Path(__file__).resolve().parents[1]
    canonical = root / "src/radia/_residual_gate.py"
    vendored = root / "packages/radia-mcp/src/radia_mcp/radia_ngsolve/_vendor/residual_gate.py"
    assert canonical.read_bytes() == vendored.read_bytes(), (
        "Do not edit the MCP distribution copy; run tools/sync_mcp_residual_gate.py")
