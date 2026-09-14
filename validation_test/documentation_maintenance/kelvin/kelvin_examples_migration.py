"""Inventory helpers for promoting Kelvin examples into docs or validation_test."""

from __future__ import annotations

import datetime as _dt
import hashlib
import importlib.metadata as _metadata
import json
import platform
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path


TEXT_SUFFIXES = {
    ".py", ".pyi", ".md", ".rst", ".txt", ".toml", ".json", ".ps1", ".yml", ".yaml"
}
SCAN_ROOTS = ("validation_test", "tests", "src", "packages", "docs")
SKIP_TEXT_PARTS = {
    ".git", ".pytest_cache", ".ipynb_checkpoints", "__pycache__",
    "build", "dist", "radia.egg-info",
}
SKIP_TEXT_NAMES = {
    "kelvin_examples_migration_results.json",
    "kelvin_examples_migration_result.json",
    "kelvin_classic_demos_results.json",
    "kelvin_classic_demos_result.json",
}

_LOCAL_PATH_PATTERNS = (
    re.compile(r"\\\\192\.168\.11\.100\\work\\[^\s'\"<>()\],}]+", re.IGNORECASE),
    re.compile(r"//192\.168\.11\.100/work/[^\s'\"<>()\],}]+", re.IGNORECASE),
    re.compile(r"\b[A-Za-z]:\\[^\s'\"<>()\],}]+"),
    re.compile(r"\b[A-Za-z]:/[^\s'\"<>()\],}]+"),
)


def public_text(text: str) -> str:
    """Redact LAB-local absolute paths from public docs archives."""
    out = text
    for pattern in _LOCAL_PATH_PATTERNS:
        out = pattern.sub("<LOCAL_PATH>", out)
    return out


def find_repo_root(start: str | Path | None = None) -> Path:
    here = Path(start) if start else Path(__file__)
    here = here if here.is_dir() else here.parent
    for cand in (here, *here.parents):
        if (cand / ".git").exists() and (cand / "docs" / "kelvin").is_dir():
            return cand
    raise FileNotFoundError("Could not locate Radia repository root")


def rel(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.as_posix()


def package_versions() -> dict:
    versions = {
        "python_version": sys.version.split()[0],
        "python_executable": sys.executable,
        "platform": platform.platform(),
    }
    for package in ("radia", "radia-mcp", "ngsolve"):
        try:
            versions[package.replace("-", "_") + "_version"] = _metadata.version(package)
        except Exception:
            versions[package.replace("-", "_") + "_version"] = None
    return versions


def git_tracked_files(root: Path) -> set[str]:
    try:
        proc = subprocess.run(
            ["git", "ls-files", "--", "examples/kelvin_transformation"],
            cwd=str(root),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    except Exception:
        return set()
    if proc.returncode != 0:
        return set()
    return {line.strip().replace("\\", "/") for line in proc.stdout.splitlines() if line.strip()}


def top_group(py_file: Path, examples_root: Path) -> str:
    relative = py_file.relative_to(examples_root)
    if len(relative.parts) == 1:
        return relative.name
    return relative.parts[0]


def is_validation_named(path: Path) -> bool:
    name = path.name.lower()
    stem = path.stem.lower()
    return (
        name.startswith(("validation_", "validate_"))
        or stem.endswith(("_validation", "_cross_validation"))
        or "validation" in {part.lower() for part in path.parts}
    )


def migration_lane(py_file: Path, examples_root: Path, validation_ref_count: int) -> str:
    parts = {part.lower() for part in py_file.relative_to(examples_root).parts}
    group = top_group(py_file, examples_root)
    if group == "Cubit_1_4_p_convergence":
        return "validation_test_locked_sample"
    if validation_ref_count:
        return "validation_test_or_src_api_locked"
    if is_validation_named(py_file):
        return "validation_test_candidate"
    if "dtn_spectrum" in parts:
        return "open_boundary_src_api_or_validation"
    if "adaptivemesh" in parts:
        return "collapse_to_notebook_or_memory"
    if group in {"A-formulation", "H-formulation", "Omega_ReducedOmega"}:
        return "docs_classic_demo_candidate"
    return "manual_review"


def iter_text_files(root: Path):
    for base_name in SCAN_ROOTS:
        base = root / base_name
        if not base.exists():
            continue
        paths = [base] if base.is_file() else base.rglob("*")
        for path in paths:
            if path.name in SKIP_TEXT_NAMES:
                continue
            if path.is_file() and path.suffix.lower() in TEXT_SUFFIXES:
                if {part.lower() for part in path.parts} & SKIP_TEXT_PARTS:
                    continue
                yield path


def reference_hits(root: Path, target_rel: str, max_hits: int = 8) -> list[dict]:
    needles = {
        target_rel,
        target_rel.replace("/", "\\"),
    }
    hits: list[dict] = []
    for path in iter_text_files(root):
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            try:
                text = path.read_text(encoding="cp932")
            except Exception:
                continue
        except Exception:
            continue
        for lineno, line in enumerate(text.splitlines(), start=1):
            if any(needle in line for needle in needles):
                hits.append({
                    "path": rel(path, root),
                    "line": lineno,
                    "text": line.strip()[:220],
                })
                if len(hits) >= max_hits:
                    return hits
    return hits


def reference_hit_map(root: Path,
                      target_rels: list[str],
                      max_hits_per_target: int = 8) -> dict[str, list[dict]]:
    """Scan text files once and return reference hits for many target paths."""
    needles_by_target = {
        target: (target, target.replace("/", "\\"))
        for target in target_rels
    }
    hits = {target: [] for target in target_rels}
    for path in iter_text_files(root):
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            try:
                text = path.read_text(encoding="cp932")
            except Exception:
                continue
        except Exception:
            continue
        for lineno, line in enumerate(text.splitlines(), start=1):
            if (
                "examples/kelvin_transformation" not in line
                and "examples\\kelvin_transformation" not in line
            ):
                continue
            for target, needles in needles_by_target.items():
                if len(hits[target]) >= max_hits_per_target:
                    continue
                if any(needle in line for needle in needles):
                    hits[target].append({
                        "path": rel(path, root),
                        "line": lineno,
                        "text": line.strip()[:220],
                    })
    return hits


def sibling_artifacts(py_file: Path) -> list[str]:
    names = []
    for suffix in (".json", ".png", ".pdf", ".html", ".vtu", ".vol", ".geo", ".mat"):
        peer = py_file.with_suffix(suffix)
        if peer.exists():
            names.append(peer.name)
    return names


def source_preview(py_file: Path, max_chars: int = 500) -> str:
    try:
        text = py_file.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        text = py_file.read_text(encoding="cp932")
    except Exception as exc:
        return f"<read failed: {exc}>"
    return public_text(text)[:max_chars]


def source_text(py_file: Path) -> str:
    try:
        return py_file.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return py_file.read_text(encoding="cp932")


def source_record(py_file: Path, root: Path, include_source: bool = True) -> dict:
    text = source_text(py_file)
    safe_text = public_text(text)
    record = {
        "path": rel(py_file, root),
        "bytes": py_file.stat().st_size,
        "sha256": hashlib.sha256(text.encode("utf-8", errors="replace")).hexdigest(),
        "line_count": len(text.splitlines()),
    }
    if include_source:
        record["source_text"] = safe_text
    return record


def _git_show_text(root: Path, rel_path: str) -> str:
    """Read a tracked file from the worktree, falling back to HEAD."""
    path = root / rel_path
    if path.exists():
        return source_text(path)
    proc = subprocess.run(
        ["git", "show", f"HEAD:{rel_path}"],
        cwd=str(root),
        text=True,
        encoding="utf-8",
        errors="replace",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if proc.returncode != 0:
        raise FileNotFoundError(f"cannot read tracked source {rel_path}: {proc.stderr}")
    return proc.stdout


def _git_ls_files(root: Path, prefix: str) -> list[str]:
    proc = subprocess.run(
        ["git", "ls-files", "--", prefix],
        cwd=str(root),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    if proc.returncode != 0:
        return []
    return [line.strip().replace("\\", "/") for line in proc.stdout.splitlines() if line.strip()]


def build_git_source_archive(repo_root: str | Path | None = None,
                             prefix: str = "examples/kelvin_transformation/AdaptiveMesh",
                             lane: str = "collapse_to_notebook_or_memory",
                             include_source: bool = True) -> dict:
    """Archive tracked source under a prefix, including files deleted in worktree.

    This is used before pruning repetitive example runners: the source is read
    from the current worktree when present and from ``HEAD:<path>`` after a file
    has already been deleted but not yet committed.
    """
    root = find_repo_root(repo_root)
    files = [p for p in _git_ls_files(root, prefix) if p.endswith(".py")]
    selected = []
    for rel_path in files:
        text = _git_show_text(root, rel_path)
        safe_text = public_text(text)
        path = root / rel_path
        try:
            top = Path(rel_path).relative_to("examples/kelvin_transformation").parts[0]
        except ValueError:
            top = Path(rel_path).parts[0]
        record = {
            "path": rel_path,
            "top_group": top,
            "migration_lane": lane,
            "validation_named": is_validation_named(Path(rel_path)),
            "reference_hit_count_reported": 0,
            "sibling_artifacts": [],
            "source": {
                "path": rel_path,
                "bytes": len(text.encode("utf-8", errors="replace")),
                "sha256": hashlib.sha256(text.encode("utf-8", errors="replace")).hexdigest(),
                "line_count": len(text.splitlines()),
                "exists_in_worktree": path.exists(),
            },
        }
        if include_source:
            record["source"]["source_text"] = safe_text
        selected.append(record)

    return {
        "schema": "radia.docs.kelvin_source_archive.v1",
        "generated_at_utc": _dt.datetime.now(
            _dt.timezone.utc
        ).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "versions": package_versions(),
        "selection": {
            "source": "git ls-files + worktree/HEAD fallback",
            "prefix": prefix,
            "lane": lane,
            "include_source": include_source,
        },
        "summary": {
            "archived_files": len(selected),
            "archived_lines": sum(row["source"]["line_count"] for row in selected),
            "archived_bytes": sum(row["source"]["bytes"] for row in selected),
            "by_lane": dict(Counter(row["migration_lane"] for row in selected).most_common()),
            "by_group": dict(Counter(row["top_group"] for row in selected).most_common()),
            "missing_from_worktree": sum(
                1 for row in selected if not row["source"]["exists_in_worktree"]
            ),
        },
        "files": selected,
    }


def build_source_archive(report: dict,
                         lanes: list[str],
                         path_contains: str = "",
                         include_source: bool = True) -> dict:
    """Build a full-source archive from selected migration-report rows."""
    root = Path(report["repo_root"])
    selected = []
    lane_set = set(lanes)
    for row in report["files"]:
        if row["migration_lane"] not in lane_set:
            continue
        if path_contains and path_contains not in row["path"]:
            continue
        path = root / row["path"]
        if not path.exists():
            continue
        selected.append({
            **{k: row[k] for k in (
                "path", "top_group", "migration_lane", "validation_named",
                "reference_hit_count_reported", "sibling_artifacts",
            )},
            "source": source_record(path, root, include_source=include_source),
        })
    return {
        "schema": "radia.docs.kelvin_source_archive.v1",
        "generated_at_utc": _dt.datetime.now(
            _dt.timezone.utc
        ).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "versions": package_versions(),
        "source_report_schema": report.get("schema"),
        "selection": {
            "lanes": lanes,
            "path_contains": path_contains,
            "include_source": include_source,
        },
        "summary": {
            "archived_files": len(selected),
            "archived_lines": sum(row["source"]["line_count"] for row in selected),
            "archived_bytes": sum(row["source"]["bytes"] for row in selected),
            "by_lane": dict(Counter(row["migration_lane"] for row in selected).most_common()),
            "by_group": dict(Counter(row["top_group"] for row in selected).most_common()),
        },
        "files": selected,
    }


def build_migration_report(repo_root: str | Path | None = None) -> dict:
    root = find_repo_root(repo_root)
    examples_root = root / "examples" / "kelvin_transformation"
    tracked = git_tracked_files(root)
    py_files = sorted(examples_root.rglob("*.py"))
    target_rels = [rel(py_file, root) for py_file in py_files]
    hits_by_target = reference_hit_map(root, target_rels)
    rows = []
    for py_file in py_files:
        target_rel = rel(py_file, root)
        hits = hits_by_target.get(target_rel, [])
        lane = migration_lane(py_file, examples_root, len(hits))
        rows.append({
            "path": target_rel,
            "top_group": top_group(py_file, examples_root),
            "name": py_file.name,
            "bytes": py_file.stat().st_size,
            "tracked": target_rel in tracked if tracked else None,
            "validation_named": is_validation_named(py_file),
            "reference_hit_count_reported": len(hits),
            "sample_reference_hits": hits[:3],
            "sibling_artifacts": sibling_artifacts(py_file),
            "migration_lane": lane,
            "source_preview": source_preview(py_file),
        })

    by_group = Counter(row["top_group"] for row in rows)
    by_lane = Counter(row["migration_lane"] for row in rows)
    artifact_counts = Counter()
    for path in examples_root.rglob("*"):
        if path.is_file():
            artifact_counts[path.suffix.lower() or "<none>"] += 1

    return {
        "schema": "radia.docs.kelvin_examples_migration.v1",
        "generated_at_utc": _dt.datetime.now(
            _dt.timezone.utc
        ).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "versions": package_versions(),
        "repo_root": str(root),
        "examples_root": rel(examples_root, root),
        "summary": {
            "python_files": len(rows),
            "tracked_python_files": sum(1 for row in rows if row["tracked"] is True),
            "validation_named_python_files": sum(1 for row in rows if row["validation_named"]),
            "referenced_python_files_reported": sum(
                1 for row in rows if row["reference_hit_count_reported"] > 0
            ),
            "top_groups": dict(by_group.most_common()),
            "migration_lanes": dict(by_lane.most_common()),
            "artifact_suffix_counts": dict(artifact_counts.most_common()),
        },
        "policy": {
            "docs": "classic explanatory demos promote to result-bearing docs/kelvin notebooks with synchronized JSON",
            "validation_test": "referenced or validation-named solver checks promote to validation_test or remain protected validation corpus",
            "src": "DtN/open-boundary behavior already used by tests/docs should promote to src API before pruning prototypes",
            "memory": "superseded adaptive/scratch studies should be distilled to memory before deletion",
        },
        "recommended_batches": [
            {
                "name": "classic_kelvin_demo_notebook",
                "lane": "docs_classic_demo_candidate",
                "description": "A/H/Omega formulation demos rendered as one or two docs notebooks.",
            },
            {
                "name": "cubit_p_convergence",
                "lane": "validation_test_locked_sample",
                "description": "Keep as validation-backed sample; docs may point at the result but should not replace it.",
            },
            {
                "name": "dtn_spectrum",
                "lane": "open_boundary_src_api_or_validation",
                "description": "Port stable act scripts into radia.open_boundary/src APIs or validation_test before pruning prototypes.",
            },
            {
                "name": "adaptive_mesh",
                "lane": "collapse_to_notebook_or_memory",
                "description": "Summarize convergence families into docs notebook/JSON, then prune superseded per-order scripts only after distilling lessons.",
            },
        ],
        "files": rows,
    }


def write_report_json(report: dict, output_path: str | Path) -> Path:
    def _public_copy(obj):
        if isinstance(obj, dict):
            return {key: _public_copy(value) for key, value in obj.items()
                    if key != "repo_root"}
        if isinstance(obj, list):
            return [_public_copy(value) for value in obj]
        if isinstance(obj, str):
            return public_text(obj)
        return obj

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(_public_copy(report), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return out


def markdown_table(rows: list[dict], columns: list[str], max_rows: int = 30) -> str:
    shown = rows[:max_rows]
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join("---" for _ in columns) + " |"
    body = []
    for row in shown:
        values = [str(row.get(col, "")).replace("|", "\\|") for col in columns]
        body.append("| " + " | ".join(values) + " |")
    if len(rows) > max_rows:
        body.append(f"| ... | {len(rows) - max_rows} more rows | | | |")
    return "\n".join([header, sep, *body])


def short_summary(report: dict) -> str:
    summary = report["summary"]
    lanes = ", ".join(
        f"{name}={count}" for name, count in summary["migration_lanes"].items()
    )
    groups = ", ".join(
        f"{name}={count}" for name, count in summary["top_groups"].items()
    )
    return "\n".join([
        f"python files: {summary['python_files']}",
        f"tracked python files: {summary['tracked_python_files']}",
        f"validation-named python files: {summary['validation_named_python_files']}",
        f"referenced python files (reported): {summary['referenced_python_files_reported']}",
        f"migration lanes: {lanes}",
        f"top groups: {groups}",
    ])
