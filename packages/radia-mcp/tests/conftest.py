"""pytest fixtures + path setup for radia-mcp tests."""

import sys
from pathlib import Path

# Ensure tests resolve `radia_mcp` to THIS checkout's src/, not whatever
# `pip install -e` happens to point at on the editable-install machine.
_SRC = Path(__file__).resolve().parent.parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))
