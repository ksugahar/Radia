# Cubit-owned implementation; intentionally maintained independently of Radia MCP.
# Derived support retains BSD-3-Clause terms: see LICENSE-BSD-3-Clause.txt.
"""Private Cubit MCP runtime support; not a separate distribution."""
from importlib import import_module


def __getattr__(name):
    if name in {"register_status_tool", "build_status_payload"}:
        return getattr(import_module(".status", __name__), name)
    raise AttributeError(name)
