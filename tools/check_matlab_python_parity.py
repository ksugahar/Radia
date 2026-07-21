"""Validate MATLAB ownership for every tracked Radia Python implementation file."""

from __future__ import annotations

import argparse
import fnmatch
import json
from collections import Counter
from pathlib import Path


def audit(repo_root: Path) -> dict:
    repo_root = repo_root.resolve()
    manifest_path = repo_root / "matlab" / "python_api_parity_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    source_root = repo_root / manifest["source_root"]
    allowed = set(manifest["classifications"])
    rules = manifest["rules"]
    backlog = manifest.get("native_promotion_backlog", [])
    errors: list[str] = []
    assignments: list[dict] = []

    rule_ids = [rule["id"] for rule in rules]
    if len(rule_ids) != len(set(rule_ids)):
        errors.append("rule ids must be unique")

    priorities = [item.get("priority") for item in backlog]
    if priorities != list(range(1, len(backlog) + 1)):
        errors.append("native promotion priorities must be contiguous and ordered")
    families = [item.get("family") for item in backlog]
    if len(families) != len(set(families)):
        errors.append("native promotion families must be unique")
    for item in backlog:
        family = item.get("family", "<missing>")
        for field in ("current", "next_native_boundary", "retained_python_boundary"):
            if not item.get(field):
                errors.append(f"{family}: missing native promotion field {field}")
        gates = item.get("gate", [])
        if not gates:
            errors.append(f"{family}: native promotion entry has no gate")
        for relative in gates:
            if not (repo_root / relative).is_file():
                errors.append(f"{family}: missing native promotion gate {relative}")

    for rule in rules:
        if rule.get("classification") not in allowed:
            errors.append(f"{rule['id']}: unknown classification")
        matlab_files = rule.get("matlab", [])
        if rule.get("classification") != "private/not-applicable" and not matlab_files:
            errors.append(f"{rule['id']}: non-private rule has no MATLAB owner")
        for relative in matlab_files:
            owner = repo_root / relative
            if not owner.is_file():
                errors.append(f"{rule['id']}: missing MATLAB owner {relative}")
            elif rule.get("backend_contract") == "in-process-python":
                text = owner.read_text(encoding="utf-8", errors="ignore")
                if "radia.internal.callPython" not in text:
                    errors.append(
                        f"{rule['id']}: {relative} does not use the checked Python fallback"
                    )
            elif rule.get("backend_contract") == "process-python":
                text = owner.read_text(encoding="utf-8", errors="ignore")
                if "PythonExecutable" not in text:
                    errors.append(
                        f"{rule['id']}: {relative} lacks the explicit process-Python contract"
                    )

    for path in sorted(source_root.rglob("*.py")):
        relative = path.relative_to(source_root).as_posix()
        matches = [
            rule for rule in rules
            if any(fnmatch.fnmatchcase(relative, pattern) for pattern in rule["patterns"])
        ]
        if not matches:
            errors.append(f"unclassified Python module: {relative}")
            continue
        rule = matches[0]
        assignments.append(
            {
                "python": relative,
                "rule": rule["id"],
                "classification": rule["classification"],
                "matlab": rule.get("matlab", []),
            }
        )

    counts = Counter(item["classification"] for item in assignments)
    fallback_rules = sorted(
        {
            item["rule"]
            for item in assignments
            if item["classification"] == "python-fallback"
        }
    )
    return {
        "schema": manifest["schema"],
        "ok": not errors,
        "manifest": str(manifest_path),
        "python_file_count": len(list(source_root.rglob("*.py"))),
        "classified_file_count": len(assignments),
        "counts": dict(sorted(counts.items())),
        "python_fallback_families": fallback_rules,
        "native_promotion_backlog": backlog,
        "binary_extensions": manifest.get("binary_extensions", []),
        "assignments": assignments,
        "errors": errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--compact", action="store_true")
    args = parser.parse_args()
    result = audit(args.repo_root)
    print(json.dumps(result, indent=None if args.compact else 2, ensure_ascii=True))
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
