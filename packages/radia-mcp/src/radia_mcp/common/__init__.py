"""Shared utilities across radia-mcp servers.

The package exports a convenient flat API, but resolves each symbol lazily.
Almost every MCP server imports one status helper; eagerly importing the RAG,
PDF, and asynchronous execution helpers at that point made every server pay
for unrelated infrastructure during cold start.
"""

from __future__ import annotations

from importlib import import_module
from typing import Any


_LAZY_ATTRS = {
    "load_prompt": (".prompts_loader", "load_prompt"),
    "list_prompts": (".prompts_loader", "list_prompts"),
    "AsyncRunner": (".async_runner", "AsyncRunner"),
    "RunStatus": (".async_runner", "RunStatus"),
    "register_status_tool": (".status", "register_status_tool"),
    "build_status_payload": (".status", "build_status_payload"),
    "register_topics_tool": (".topics", "register_topics_tool"),
    "CoarseToolRegistry": (".tool_group", "CoarseToolRegistry"),
    "selected_tool_profile": (".tool_group", "selected_tool_profile"),
    "lazy_callable": (".lazy_call", "lazy_callable"),
    "ChromaRetriever": (".chroma_retriever", "ChromaRetriever"),
    "extract_pdf_chunks": (".chroma_retriever", "extract_pdf_chunks"),
    "detect_filename_language": (
        ".chroma_retriever",
        "detect_filename_language",
    ),
    "find_chapters": (".chroma_retriever", "find_chapters"),
    "CHAPTER_PATTERNS": (".chroma_retriever", "CHAPTER_PATTERNS"),
    "DEFAULT_EMBEDDING_MODEL": (
        ".chroma_retriever",
        "DEFAULT_EMBEDDING_MODEL",
    ),
}

_LAZY_MODULES = {
    "async_runner",
    "chroma_retriever",
    "examples",
    "failure_log",
    "learning_quality",
    "mcp_contract",
    "prompts_loader",
    "server_hardening",
    "status",
    "tool_group",
    "topics",
    "utf8_stdout",
    "web_docs",
}


def __getattr__(name: str) -> Any:
    target = _LAZY_ATTRS.get(name)
    if target is not None:
        module_name, attribute = target
        value = getattr(import_module(module_name, __name__), attribute)
        globals()[name] = value
        return value
    if name in _LAZY_MODULES:
        value = import_module(f".{name}", __name__)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(_LAZY_ATTRS) | _LAZY_MODULES)

__all__ = [
    "load_prompt",
    "list_prompts",
    "AsyncRunner",
    "RunStatus",
    "register_status_tool",
    "build_status_payload",
    "register_topics_tool",
    "CoarseToolRegistry",
    "selected_tool_profile",
    "lazy_callable",
    # ChromaDB RAG helpers (multilingual support added 2026-05-25)
    "ChromaRetriever",
    "extract_pdf_chunks",
    "detect_filename_language",
    "find_chapters",
    "CHAPTER_PATTERNS",
    "DEFAULT_EMBEDDING_MODEL",
]
