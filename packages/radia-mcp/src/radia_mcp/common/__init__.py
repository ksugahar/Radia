"""Shared utilities across radia-mcp servers.

Existing utilities:
  - `examples`, `failure_log`, `utf8_stdout`, `web_docs`

New (2026-05-24, learned from COMSOL MCP analysis):
  - `prompts_loader` — load knowledge as .md files (avoids docstring
    inside r-string parser errors hit in sibc/hoibc)
  - `async_runner` — threading-based long-running command wrapper
    with progress + cancellation
  - `chroma_retriever` — ChromaDB + sentence-transformers RAG
    infrastructure (optional dependency)
"""

from .prompts_loader import load_prompt, list_prompts
from .async_runner import AsyncRunner, RunStatus
from .status import register_status_tool, build_status_payload
from .topics import register_topics_tool
from .chroma_retriever import (
    ChromaRetriever,
    extract_pdf_chunks,
    detect_filename_language,
    find_chapters,
    CHAPTER_PATTERNS,
    DEFAULT_EMBEDDING_MODEL,
)

__all__ = [
    "load_prompt",
    "list_prompts",
    "AsyncRunner",
    "RunStatus",
    "register_status_tool",
    "build_status_payload",
    "register_topics_tool",
    # ChromaDB RAG helpers (multilingual support added 2026-05-25)
    "ChromaRetriever",
    "extract_pdf_chunks",
    "detect_filename_language",
    "find_chapters",
    "CHAPTER_PATTERNS",
    "DEFAULT_EMBEDDING_MODEL",
]
