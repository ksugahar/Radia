"""document_meta sub-module — cross-cutting authoring helpers."""
from . import tools as _tools


PUBLIC_TOOLS = (
    "document_meta_deadline_countdown",
    "document_meta_diff_versions",
    "document_meta_template_loader",
    "document_meta_lint_all",
    "document_meta_notebook_result_audit",
)


def register(mcp) -> int:
    """Register the intentional public document-meta surface."""
    for name in PUBLIC_TOOLS:
        mcp.tool()(getattr(_tools, name))
    return len(PUBLIC_TOOLS)
