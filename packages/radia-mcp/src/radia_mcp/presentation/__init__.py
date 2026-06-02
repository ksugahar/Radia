"""presentation sub-module — 学会発表スライドの作文技術."""
from . import tools as _tools


def register(mcp) -> int:
    """Register presentation_* tools with the given FastMCP instance."""
    count = 0
    for name in dir(_tools):
        if name.startswith("presentation_") and callable(getattr(_tools, name)):
            mcp.tool()(getattr(_tools, name))
            count += 1
    return count
