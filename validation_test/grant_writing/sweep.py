"""Run every grant-writing detector over a corpus of real proposals.

Fixing this suite one document at a time found one false positive per
session. Running every detector over every document at once and adjudicating
the result found eight in a single pass, so the sweep itself is the tool and
this module is it.

The corpus is **not** in the repository. Real proposals belong to their
authors, and several are colleagues' work; the repository is public. Point
``GRANT_WRITING_CORPUS`` at a manifest kept outside the tree:

```json
{
  "documents": [
    {"label": "adopted-kiban", "path": "texts/adopted_kiban.txt",
     "outcome": "adopted", "program": "kaken_oss"}
  ]
}
```

Relative paths resolve against the manifest's directory. ``outcome`` is
recorded but never scored: four measurements have found no relationship
between these checks and adoption.

Usage::

    python validation_test/grant_writing/sweep.py
    python validation_test/grant_writing/sweep.py --write-baseline
"""

from __future__ import annotations

import argparse
import collections
import json
import os
import pathlib

from radia_mcp.grant_writing import tools as gw

MANIFEST_ENV = "GRANT_WRITING_CORPUS"
BASELINE_NAME = "baseline.json"


def manifest_path() -> pathlib.Path | None:
    raw = os.environ.get(MANIFEST_ENV, "").strip()
    if not raw:
        return None
    path = pathlib.Path(raw)
    return path if path.is_file() else None


def load_corpus(manifest: pathlib.Path) -> list[dict]:
    data = json.loads(manifest.read_text(encoding="utf-8"))
    documents = []
    for entry in data.get("documents", []):
        path = pathlib.Path(entry["path"])
        if not path.is_absolute():
            path = manifest.parent / path
        if not path.is_file():
            raise FileNotFoundError(f"corpus document not found: {path}")
        documents.append({
            "label": entry["label"],
            "path": path,
            "outcome": entry.get("outcome", "unknown"),
            "program": entry.get("program", "generic"),
        })
    if not documents:
        raise ValueError(f"corpus manifest lists no documents: {manifest}")
    return documents


def measure(document: dict) -> dict:
    text = document["path"].read_text(encoding="utf-8")
    report = gw.grant_writing_health_report(text, program=document["program"])
    patterns: collections.Counter = collections.Counter()
    for key, result in report["detailed_results"].items():
        if key not in gw._DETECTOR_RESULT_KEYS:
            continue
        for risk in (result.get("risks") or []):
            patterns[f"{key}/{risk.get('type') or '?'}"] += 1
        for issue in (result.get("issues") or []):
            patterns[f"{key}/{issue.get('rule') or '?'}"] += 1
    prose = gw._prose_for_lint(text)
    return {
        "finding_count": len(report["findings"]),
        "prose_chars": len(prose),
        "patterns": dict(sorted(patterns.items())),
    }


def sweep(documents: list[dict]) -> dict:
    return {d["label"]: measure(d) for d in documents}


def render(results: dict) -> str:
    lines = ["%-28s %8s %7s" % ("document", "prose", "findings")]
    for label, row in results.items():
        lines.append("%-28s %8d %7d"
                     % (label, row["prose_chars"], row["finding_count"]))
    totals: collections.Counter = collections.Counter()
    for row in results.values():
        totals.update(row["patterns"])
    lines.append("")
    lines.append("patterns: %d   total findings: %d"
                 % (len(totals), sum(totals.values())))
    for pattern, count in totals.most_common():
        lines.append("  %-52s %3d" % (pattern, count))
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-baseline", action="store_true",
                        help="record the current counts as the expected ones")
    args = parser.parse_args()

    manifest = manifest_path()
    if manifest is None:
        print(f"set {MANIFEST_ENV} to a corpus manifest; nothing to sweep")
        return 2

    results = sweep(load_corpus(manifest))
    print(render(results))

    if args.write_baseline:
        target = manifest.parent / BASELINE_NAME
        target.write_text(
            json.dumps(results, ensure_ascii=False, indent=1) + "\n",
            encoding="utf-8",
        )
        print(f"\nbaseline written: {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
