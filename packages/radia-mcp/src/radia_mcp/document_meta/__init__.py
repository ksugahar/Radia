"""document_meta sub-module — cross-cutting authoring helpers."""
from . import tools as _tools


def register(mcp) -> int:
    """Register document_meta_* tools with the given FastMCP instance."""
    count = 0
    for name in dir(_tools):
        if name.startswith("document_meta_") and callable(getattr(_tools, name)):
            mcp.tool()(getattr(_tools, name))
            count += 1
    return count
