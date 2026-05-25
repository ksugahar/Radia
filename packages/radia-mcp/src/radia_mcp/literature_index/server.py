"""MCP Server: radia_mcp.literature_index

Meta-MCP for searching the lab literature corpus.

Indexes W:/03_文献・論文/00_電磁界解析/ (3,889 files / 51 GB).

Two search modes (2026-05-24 update):
  - `literature_search` — keyword/filename match (fast, no embedding)
  - `literature_semantic_search` — ChromaDB + sentence-transformers
    RAG over full PDF text (requires `pip install chromadb
    sentence-transformers pymupdf`); first-time index build runs
    asynchronously via radia_mcp.common.AsyncRunner.

Cache: %LOCALAPPDATA%/radia_mcp_literature_index/lit_index.json
Vector store: %LOCALAPPDATA%/radia_mcp_literature_index/chroma/

Usage:
    mcp-server-literature-index              # stdio
    mcp-server-literature-index --selftest   # self-test
"""
import os
import sys
from pathlib import Path

from mcp.server.fastmcp import FastMCP

from .index_tools import lit_search, lit_by_folder, lit_folder_tree, lit_stats
from ..common import AsyncRunner, register_status_tool
from ..common.chroma_retriever import ChromaRetriever, extract_pdf_chunks

mcp = FastMCP("mcp-server-literature-index")


# Persistent ChromaDB location (per-user, off the NAS)
_CHROMA_DIR = Path(
    os.environ.get("LOCALAPPDATA",
                    os.path.expanduser("~/.cache"))
) / "radia_mcp_literature_index" / "chroma"

_retriever: ChromaRetriever | None = None
_index_runner = AsyncRunner()


def _get_retriever() -> ChromaRetriever:
    """Lazy-init the ChromaDB retriever (only when first asked)."""
    global _retriever
    if _retriever is None:
        _retriever = ChromaRetriever(
            db_dir=_CHROMA_DIR,
            collection="lit_em",
        )
    return _retriever


@mcp.tool()
def literature_search(query: str, limit: int = 30, folder_filter: str = "") -> str:
    """
    Search the lab literature corpus by filename keywords.

    Args:
        query: Case-insensitive substring(s) — all must match (AND).
               Special: query='rebuild' forces index rebuild.
        limit: Max number of matches to return (default 30)
        folder_filter: Optional top-folder name substring filter
                       (e.g. '30_磁気特性', '99_アプリ')

    Examples:
        literature_search("RWG basis")
        literature_search("hysteresis vector", folder_filter="30_磁気特性")
        literature_search("rebuild")  # force re-scan W: drive
    """
    return lit_search(query, limit=limit,
                      folder_filter=folder_filter or None)


@mcp.tool()
def literature_by_folder(folder: str, limit: int = 50) -> str:
    """
    List papers in a specific top-level / nested folder.

    Args:
        folder: Folder path substring (e.g. '30_磁気特性',
                '99_アプリケーション/01_誘導加熱', '11_BEM_モーメント法')
        limit: Max number of files (default 50)

    Returns:
        Markdown file list grouped by sub-folder.
    """
    return lit_by_folder(folder, limit=limit)


@mcp.tool()
def literature_folder_tree() -> str:
    """List all top-level folders with file count + size."""
    return lit_folder_tree()


@mcp.tool()
def literature_stats() -> str:
    """Index statistics + cache info."""
    return lit_stats()


# ============================================================
# Semantic search (RAG) — added 2026-05-24
# Pattern adopted from wjc9011/COMSOL_Multiphysics_MCP
# ============================================================

@mcp.tool()
def literature_semantic_search(
    query: str,
    n_results: int = 5,
    source_filter: str = "",
    language_filter: str = "",
) -> str:
    """Semantic search over indexed PDF text via ChromaDB + sentence-
    transformers (all-MiniLM-L6-v2).

    Returns markdown-formatted hits with source, page, similarity score,
    language tag (if any), and the matching text chunk.

    First-time use: call `literature_build_vector_index` to populate
    the vector store from W:/.../00_電磁界解析/ PDFs. Without an index,
    this returns an empty result + hint.

    Args:
        query: free-text query (Japanese or English; the embedding model
               is multilingual-friendly).
        n_results: top N hits to return (default 5).
        source_filter: optional filename substring filter (exact match
                       on source filename).
        language_filter: optional ISO 639-1 tag (e.g. "ja", "en", "zh")
                       to restrict hits to one language. Only effective
                       when the index was built with
                       ``auto_detect_language=True`` (the default for
                       ``literature_build_vector_index``) or with an
                       explicit ``default_language=...``.
                       Use "ja" to filter the lab's Japanese textbooks,
                       "en" to filter the IEEE / Elsevier paper subset.

    Examples:
        literature_semantic_search("rotational core loss measurement")
        literature_semantic_search("Hayaho MCTS PM motor", n_results=10)
        literature_semantic_search("ヒステリシス", language_filter="ja")
    """
    rag = _get_retriever()
    stats = rag.stats()
    if stats.get("n_docs", 0) == 0:
        return (
            "**No documents indexed yet.**\n\n"
            "Run `literature_build_vector_index` to index PDFs from "
            "W:/03_文献・論文/00_電磁界解析/ into the vector store. "
            "First-time build can take 30-90 min depending on # PDFs."
        )

    where = {"source": {"$eq": source_filter}} if source_filter else None
    hits = rag.search(
        query,
        n_results=n_results,
        where=where,
        language_filter=language_filter or None,
    )
    if not hits:
        suffix = ""
        if language_filter:
            suffix = (f" (language='{language_filter}' -- did you index "
                      f"with `auto_detect_language=True`?)")
        return f"No semantic hits for `{query}`{suffix}."

    title = f"Semantic search: `{query}`"
    if language_filter:
        title += f" [lang={language_filter}]"
    out = [f"# {title} ({len(hits)} hits)\n"]
    for i, hit in enumerate(hits, 1):
        meta = hit.get("metadata", {})
        lang_tag = f" [{meta['language']}]" if meta.get("language") else ""
        out.append(f"## {i}. score={hit['score']:.3f}  "
                   f"`{meta.get('source', '?')}` p.{meta.get('page', '?')}{lang_tag}")
        text = hit["text"].strip()
        snippet = (text[:600] + "...") if len(text) > 600 else text
        out.append(f"\n{snippet}\n")
    return "\n".join(out)


@mcp.tool()
def literature_build_vector_index(
    folder: str = "30_磁気特性",
    max_pdfs: int = 50,
    chunk_size: int = 1000,
    overlap: int = 100,
    default_language: str = "",
    auto_detect_language: bool = True,
) -> str:
    """Build / extend the ChromaDB vector index from PDFs.

    Runs ASYNCHRONOUSLY (via radia_mcp.common.AsyncRunner). Returns
    immediately with a job ID; poll with `literature_index_job_status`.

    Args:
        folder: subfolder of W:/.../00_電磁界解析/ to ingest
                (e.g. "30_磁気特性", "11_BEM_モーメント法",
                or "" for ENTIRE corpus -- multi-hour build).
        max_pdfs: cap for this run (default 50; safety against
                  accidental full-corpus runs).
        chunk_size: characters per chunk (default 1000).
        overlap: chunk overlap (default 100).
        default_language: ISO 639-1 tag attached to EVERY chunk in this
                  run (e.g. "ja" for an entirely-Japanese folder,
                  "en" for an English papers folder). Overrides
                  auto-detect. "" = unset (use auto-detect if enabled).
        auto_detect_language: when True (default), guess each PDF's
                  language from its filename via CJK Unicode heuristic
                  (radia_mcp.common.detect_filename_language). Enables
                  ``literature_semantic_search(..., language_filter=)``
                  filtering on the bilingual lab corpus.

    Returns:
        Job status message.
    """
    if _index_runner.is_running:
        s = _index_runner.get_status()
        return (f"Already running: {s.get('label')}, "
                f"progress {s.get('progress', 0):.1%}, "
                f"elapsed {s.get('elapsed_seconds', 0):.0f}s.")

    base = Path(r"W:\03_文献・論文\00_電磁界解析")
    target = base / folder if folder else base
    if not target.exists():
        return f"Target folder does not exist: {target}"

    lang_default = default_language or None

    def _build(_progress_cb=None, _cancel_check=None):
        pdfs = list(target.rglob("*.pdf"))
        pdfs = pdfs[:max_pdfs]
        if not pdfs:
            return {"added": 0, "msg": "No PDFs found."}
        rag = _get_retriever()
        rag.initialize()
        total_added = 0
        failed = []
        for i, p in enumerate(pdfs):
            if _cancel_check and _cancel_check():
                break
            if _progress_cb:
                _progress_cb(i / len(pdfs), f"Processing {p.name}")
            try:
                chunks = extract_pdf_chunks(
                    p,
                    chunk_size=chunk_size,
                    overlap=overlap,
                    default_language=lang_default,
                    auto_detect_language=auto_detect_language,
                )
                if chunks:
                    total_added += rag.add_chunks(chunks)
            except Exception as e:
                failed.append((p.name, str(e)))
        return {
            "added": total_added,
            "n_pdfs": len(pdfs),
            "n_failed": len(failed),
            "failed_samples": failed[:5],
        }

    started = _index_runner.start(
        target=_build,
        label=f"build_index({folder or 'ALL'},max={max_pdfs},"
              f"lang={lang_default or ('auto' if auto_detect_language else 'none')})",
    )
    if started:
        return (f"Started indexing: folder={folder or 'ALL'}, "
                f"max_pdfs={max_pdfs}, language="
                f"{lang_default or ('auto-detect' if auto_detect_language else 'none')}. "
                f"Poll with `literature_index_job_status`.")
    return "Failed to start (another job running?)."


@mcp.tool()
def literature_index_job_status() -> str:
    """Status of the most recent indexing job (idle / running /
    completed / failed)."""
    s = _index_runner.get_status()
    lines = [
        f"**Status**: {s['status']}",
        f"**Label**: {s.get('label') or '(none)'}",
        f"**Progress**: {s['progress']:.1%}",
        f"**Message**: {s.get('message', '')}",
        f"**Elapsed**: {s['elapsed_seconds']:.1f} s",
    ]
    if s.get("error"):
        lines.append(f"**Error**: {s['error']}")
    # If completed, also show the result payload
    if s["status"] == "completed":
        final = _index_runner.wait(timeout=0)   # already done
        if "result" in final and final["result"]:
            r = final["result"]
            lines.append(f"**Added chunks**: {r.get('added', 0)} "
                         f"from {r.get('n_pdfs', 0)} PDFs "
                         f"({r.get('n_failed', 0)} failed)")
    # Vector store stats
    try:
        stats = _get_retriever().stats()
        lines.append(f"**Vector store**: {stats.get('n_docs', 0)} chunks "
                     f"in `{stats.get('collection')}`")
    except Exception as e:
        lines.append(f"**Vector store**: not available ({e})")
    return "\n".join(lines)


@mcp.tool()
def literature_index_cancel() -> str:
    """Cancel the running indexing job (cooperative; checks every PDF)."""
    if _index_runner.cancel():
        return "Cancellation signalled (will stop after current PDF)."
    return "No running job to cancel."


register_status_tool(
    mcp,
    server_name="mcp-server-literature-index",
    description=("Meta-MCP: keyword + semantic search across 3,889 lab "
                  "literature files (W:/03_文献・論文/00_電磁界解析)."),
    subpackage="radia_mcp.literature_index",
    related_servers=["radia-meta", "radia-mcmc", "radia-team-benchmark"],
    optional_deps=["chromadb", "sentence_transformers", "fitz"],
)


def main():
    if "--selftest" in sys.argv:
        print("literature_index MCP server self-test:")
        s = lit_stats()
        # Print first 5 lines only
        for L in s.split("\n")[:8]:
            print(f"  {L}")
        print("OK")
        return
    mcp.run()


if __name__ == "__main__":
    main()
