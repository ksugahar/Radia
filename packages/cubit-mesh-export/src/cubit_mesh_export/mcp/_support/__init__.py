"""Private Cubit MCP runtime support; not a separate distribution."""
from importlib import import_module


def __getattr__(name):
    if name in {"register_status_tool", "build_status_payload"}:
        return getattr(import_module(".status", __name__), name)
    raise AttributeError(name)
