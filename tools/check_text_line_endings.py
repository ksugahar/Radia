"""Check tracked text for CRLF/mixed endings, honoring attributes.

Git's index EOL classifier detects binary data; explicit -text/binary evidence
is preserved. The default checks the index and modifies nothing. --worktree
checks the checked-out bytes of LF index text: a long-lived checkout (a CI
runner workspace) can keep CRLF files from an older checkout because Git never
rewrites stat-clean files, which breaks hash-bound fixtures. --repair rewrites
only those files from the index, then checks again. UTF-8/LF is the authoring
policy; this guard checks line endings only.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess


def check_index(root: Path) -> dict:
    result = subprocess.run(
        ["git", "-C", str(root), "ls-files", "--eol", "-z"],
        check=True, capture_output=True)
    violations = []
    checked = preserved = binary = 0
    for record in result.stdout.split(b"\0"):
        if not record:
            continue
        metadata, path = record.split(b"\t", 1)
        fields = metadata.decode("ascii").split()
        indexed = fields[0].removeprefix("i/")
        # Attribute override wins over automatic text classification. This
        # includes historical hash-bound JSON and source-form evidence.
        if "attr/-text" in fields:
            preserved += 1
            continue
        if indexed == "-text":
            binary += 1
            continue
        checked += 1
        if indexed in {"crlf", "mixed"}:
            violations.append({"path": path.decode("utf-8"), "index_eol": indexed})
    return {"schema": "radia.tracked-text-eol.v1", "passed": not violations,
            "checked_text_files": checked, "preserved_attribute_files": preserved,
            "auto_binary_files": binary, "violations": violations}


def check_worktree(root: Path) -> dict:
    result = subprocess.run(
        ["git", "-C", str(root), "ls-files", "--eol", "-z"],
        check=True, capture_output=True)
    violations = []
    checked = 0
    for record in result.stdout.split(b"\0"):
        if not record:
            continue
        metadata, path = record.split(b"\t", 1)
        fields = metadata.decode("ascii").split()
        if "attr/-text" in fields or fields[0] != "i/lf":
            continue
        checked += 1
        worktree = fields[1].removeprefix("w/")
        if worktree in {"crlf", "mixed"}:
            violations.append({"path": path.decode("utf-8"), "worktree_eol": worktree})
    return {"schema": "radia.worktree-text-eol.v1", "passed": not violations,
            "checked_lf_index_files": checked, "violations": violations}


def edited_beyond_line_endings(root: Path, violations: list[dict]) -> list[str]:
    """Paths whose bytes differ from the index blob by more than CRLF -> LF."""
    paths = [v["path"] for v in violations]
    specs = b"".join(b":" + p.encode("utf-8") + b"\n" for p in paths)
    blobs = subprocess.run(["git", "-C", str(root), "cat-file", "--batch"],
                           input=specs, check=True, capture_output=True).stdout
    edited, offset = [], 0
    for path in paths:
        header_end = blobs.index(b"\n", offset)
        header = blobs[offset:header_end].split()
        if len(header) != 3:
            raise RuntimeError(f"index blob unavailable: {path}")
        size = int(header[2])
        indexed = blobs[header_end + 1:header_end + 1 + size]
        offset = header_end + 1 + size + 1
        if (root / path).read_bytes().replace(b"\r\n", b"\n") != indexed:
            edited.append(path)
    return edited


def repair_worktree(root: Path, violations: list[dict]) -> None:
    if not violations:
        return
    # Rewrite only line-ending drift; a real edit must never be discarded.
    edited = edited_beyond_line_endings(root, violations)
    if edited:
        raise SystemExit("refusing --repair: content differs from the index beyond "
                         "CRLF -> LF: " + ", ".join(edited))
    # checkout-index skips a stat-clean file even with --force, so remove the
    # stale copy first; the index still holds every byte to write back.
    for v in violations:
        (root / v["path"]).unlink()
    paths = b"".join(v["path"].encode("utf-8") + b"\0" for v in violations)
    subprocess.run(["git", "-C", str(root), "checkout-index", "--force", "--index", "-z", "--stdin"],
                   input=paths, check=True)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--worktree", action="store_true",
                        help="check checked-out bytes of LF index text instead of the index")
    parser.add_argument("--repair", action="store_true",
                        help="with --worktree, rewrite violating files from the index")
    args = parser.parse_args(argv)
    if args.repair and not args.worktree:
        parser.error("--repair requires --worktree")
    root = args.root.resolve()
    if args.worktree:
        result = check_worktree(root)
        if args.repair and result["violations"]:
            repaired = len(result["violations"])
            repair_worktree(root, result["violations"])
            result = check_worktree(root)
            result["repaired_files"] = repaired
    else:
        result = check_index(root)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
