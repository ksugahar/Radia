"""radia_mcp.literature_index: PDF literature search across W:/03_文献・論文/.

Meta-MCP that indexes all PDFs in the lab literature folder and
provides search capabilities:
- `lit_search(query)` — find papers by filename keywords
- `lit_by_folder(folder)` — list papers in a specific topic folder
- `lit_stats` — coverage statistics

Index is built on-demand by walking W:/.../00_電磁界解析/ recursively.
Cached as JSON for fast subsequent access.

Covers ~3,889 files / 51 GB across 25 numbered top-level folders.
"""

from .index_tools import (
    lit_search,
    lit_by_folder,
    lit_stats,
    lit_folder_tree,
)

__all__ = [
    "lit_search",
    "lit_by_folder",
    "lit_stats",
    "lit_folder_tree",
]
