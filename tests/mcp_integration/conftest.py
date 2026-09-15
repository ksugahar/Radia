"""Repository evidence contracts: no solver/native imports or shared installs."""

from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "packages/radia-mcp/src"))


@pytest.fixture(autouse=True)
def selected_repository(monkeypatch):
    """Never inspect a different editable checkout through a user's environment."""
    monkeypatch.setenv("RADIA_REPO_ROOT", str(ROOT))
