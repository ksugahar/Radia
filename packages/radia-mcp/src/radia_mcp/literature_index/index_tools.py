"""Index + search tools for a user-provided literature corpus."""
import json
import os
from pathlib import Path
from typing import Optional

# Source root. Public wheels do not embed lab-local NAS paths.  Set
# RADIA_LIT_ROOT explicitly on machines that own a literature corpus; otherwise
# the index degrades gracefully to empty.
_LIT_ROOT_RAW = os.environ.get("RADIA_LIT_ROOT", "")
LIT_ROOT = Path(_LIT_ROOT_RAW) if _LIT_ROOT_RAW else Path("__RADIA_LIT_ROOT_NOT_SET__")

# Cache location
CACHE_DIR = Path(os.environ.get("LOCALAPPDATA", str(Path.home() / ".cache"))) / "radia_mcp_literature_index"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
INDEX_PATH = CACHE_DIR / "lit_index.json"


def _build_index() -> list[dict]:
    """Walk LIT_ROOT and build a flat index of every PDF / DOC / etc."""
    if not LIT_ROOT.exists():
        return []
    out = []
    for p in LIT_ROOT.rglob("*"):
        if not p.is_file():
            continue
        if p.suffix.lower() not in (".pdf", ".docx", ".doc", ".pptx", ".ppt",
                                     ".xlsx", ".bib", ".md"):
            continue
        try:
            size_kb = p.stat().st_size // 1024
        except OSError:
            continue
        rel = str(p.relative_to(LIT_ROOT))
        # Top-level folder for sorting
        parts = rel.split(os.sep)
        top = parts[0] if parts else "."
        out.append({
            "name": p.name,
            "rel_path": rel,
            "top_folder": top,
            "size_kb": size_kb,
            "suffix": p.suffix.lower(),
        })
    return out


def _load_or_build_index(force_rebuild: bool = False) -> list[dict]:
    """Load cached index or rebuild from disk."""
    if INDEX_PATH.exists() and not force_rebuild:
        try:
            return json.loads(INDEX_PATH.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            pass
    idx = _build_index()
    INDEX_PATH.write_text(json.dumps(idx, ensure_ascii=False), encoding="utf-8")
    return idx


def lit_search(query: str, limit: int = 30, folder_filter: Optional[str] = None) -> str:
    """Search literature index by filename keywords.

    Args:
        query: case-insensitive substring(s) to match in filename.
               Multiple keywords separated by space, all must match (AND).
               Use 'rebuild' to force index rebuild.
        limit: max number of matches to return
        folder_filter: optional substring to filter by top folder name

    Returns:
        Markdown-formatted match list with rel_path + size.
    """
    if query.strip().lower() == "rebuild":
        idx = _load_or_build_index(force_rebuild=True)
        return f"Index rebuilt: {len(idx)} files cached at {INDEX_PATH}"

    idx = _load_or_build_index()
    if not idx:
        return ("No index available. Run with query='rebuild' to scan "
                f"{LIT_ROOT}. Note: the configured literature root must be accessible.")

    kws = [k.lower() for k in query.split() if k.strip()]
    if not kws:
        return "No search keywords provided."

    matches = []
    for entry in idx:
        name_lower = entry["name"].lower()
        rel_lower = entry["rel_path"].lower()
        if all(kw in name_lower or kw in rel_lower for kw in kws):
            if folder_filter and folder_filter.lower() not in entry["top_folder"].lower():
                continue
            matches.append(entry)

    if not matches:
        return f"No matches for query '{query}' in {len(idx)} indexed files."

    lines = [f"# {len(matches)} matches for '{query}'"]
    if folder_filter:
        lines[0] += f" (filtered by folder: {folder_filter})"
    if len(matches) > limit:
        lines[0] += f" (showing first {limit})"
    lines.append("")
    for m in matches[:limit]:
        size_mb = m["size_kb"] / 1024
        lines.append(f"- `{m['top_folder']}/.../`**{m['name']}**  ({size_mb:.1f} MB)")
        lines.append(f"   {m['rel_path']}")
    return "\n".join(lines)


def lit_by_folder(folder: str, limit: int = 50) -> str:
    """List papers in a specific top-level folder.

    Args:
        folder: substring matching top-level folder (e.g. '30_磁気特性',
                '11_BEM', '99_アプリケーション/02_加速器')
        limit: max number of files to list

    Returns:
        Markdown-formatted file list grouped by sub-folder.
    """
    idx = _load_or_build_index()
    if not idx:
        return ("No index available. Run lit_search('rebuild') first.")

    matches = [e for e in idx if folder.lower() in e["rel_path"].lower()]
    if not matches:
        return f"No papers in folders matching '{folder}'."

    # Group by parent folder
    from collections import defaultdict
    grouped = defaultdict(list)
    for m in matches:
        parent = os.path.dirname(m["rel_path"])
        grouped[parent].append(m)

    lines = [f"# Papers in folders matching '{folder}': {len(matches)} total"]
    if len(matches) > limit:
        lines[0] += f" (showing first {limit})"
    lines.append("")
    shown = 0
    for parent in sorted(grouped.keys()):
        if shown >= limit:
            break
        files = grouped[parent]
        lines.append(f"\n## {parent or '(root)'}")
        for m in files:
            if shown >= limit:
                break
            size_mb = m["size_kb"] / 1024
            lines.append(f"- {m['name']}  ({size_mb:.1f} MB)")
            shown += 1
    return "\n".join(lines)


def lit_folder_tree() -> str:
    """List all top-level folders + file counts.

    Returns:
        Markdown table of top folder, file count, total size.
    """
    idx = _load_or_build_index()
    if not idx:
        return "No index available."

    from collections import defaultdict
    totals = defaultdict(lambda: {"count": 0, "size_kb": 0})
    for e in idx:
        totals[e["top_folder"]]["count"] += 1
        totals[e["top_folder"]]["size_kb"] += e["size_kb"]

    lines = ["# 00_電磁界解析 top-level folders\n",
             "| Folder | Files | Size (MB) |",
             "|--------|-------|-----------|"]
    for top, info in sorted(totals.items()):
        lines.append(f"| {top} | {info['count']} | {info['size_kb']/1024:.0f} |")
    total_files = sum(t["count"] for t in totals.values())
    total_size = sum(t["size_kb"] for t in totals.values()) / 1024 / 1024
    lines.append(f"| **TOTAL** | **{total_files}** | **{total_size:.1f} GB** |")
    return "\n".join(lines)


def lit_stats() -> str:
    """Index statistics + cache info."""
    idx = _load_or_build_index()
    if not idx:
        return f"No index. LIT_ROOT={LIT_ROOT} (exists={LIT_ROOT.exists()})"

    from collections import Counter
    suffix_count = Counter(e["suffix"] for e in idx)
    total_size_gb = sum(e["size_kb"] for e in idx) / 1024 / 1024

    lines = [
        f"# Literature Index Statistics",
        "",
        f"- Source: `{LIT_ROOT}`",
        f"- Index cache: `{INDEX_PATH}`",
        f"- Total files: {len(idx)}",
        f"- Total size: {total_size_gb:.1f} GB",
        "",
        "## File type breakdown:",
    ]
    for suf, n in suffix_count.most_common():
        lines.append(f"- `{suf}`: {n}")
    return "\n".join(lines)
