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
        pdf = entry.get("pdf")
        if pdf:
            pdf = pathlib.Path(pdf)
            if not pdf.is_absolute():
                pdf = manifest.parent / pdf
            if not pdf.is_file():
                raise FileNotFoundError(f"compiled proposal not found: {pdf}")
        documents.append({
            "label": entry["label"],
            "path": path,
            "pdf": pdf,
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

    # A page limit is a property of the rendered document, so it is checked
    # against the compiled PDF when the manifest names one. It is the only
    # defect class that gets a proposal returned before anyone reads it,
    # which makes it the one worth locking even while it reports nothing.
    pages = None
    if document.get("pdf"):
        limits = gw.grant_writing_page_limit_check(str(document["pdf"]))
        for risk in limits["risks"]:
            patterns[f"page_limit/{risk['type'] if 'type' in risk else risk['severity']}"] += 1
        pages = {
            field["field"]: [field["used_pages"], field["declared_max_pages"]]
            for field in limits["fields"]
        }

    prose = gw._prose_for_lint(text)
    measured = {
        "finding_count": len(report["findings"]),
        "prose_chars": len(prose),
        "patterns": dict(sorted(patterns.items())),
    }
    if pages is not None:
        measured["pages"] = pages
    return measured


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


# Result shapes differ per check: some return risks, some issues, some a bare
# count. A check that never reports anything on any real proposal is either
# correctly quiet or quietly broken, and only the inventory tells them apart.
_FINDING_LISTS = (
    "risks", "issues", "findings", "over_threshold_examples", "variants",
    "misuses", "undefined_acronyms", "weak_expressions", "top_fixes",
    "unbacked_absence_claims",
)
_NOT_DETECTORS = frozenset({
    "health_report", "usage", "recommendation_letter_template",
    "page_limit_check", "check_kanji_ratio",
})


def _reported_something(result) -> bool:
    if not isinstance(result, dict):
        return False
    if any(isinstance(result.get(k), list) and result[k] for k in _FINDING_LISTS):
        return True
    for key, value in result.items():
        if isinstance(value, int) and value > 0 and (
            key.startswith("total_") or key.endswith(("_count", "_matches"))
        ):
            return True
    score = result.get("score")
    return isinstance(score, (int, float)) and score < 10


def audit(documents: list[dict]) -> list[dict]:
    """Report, per check, how often it applied and how often it said anything."""
    import inspect

    loaded = [(d, d["path"].read_text(encoding="utf-8")) for d in documents]
    rows = []
    for name in sorted(n for n in dir(gw) if n.startswith("grant_writing_")):
        short = name[len("grant_writing_"):]
        if short in _NOT_DETECTORS:
            continue
        fn = getattr(gw, name)
        signature = inspect.signature(fn)
        parameters = list(signature.parameters)
        if not parameters or parameters[0] not in ("text", "text_or_path"):
            continue
        applied = reported = 0
        error = ""
        for document, text in loaded:
            try:
                if "program" in signature.parameters:
                    result = fn(text, program=document["program"])
                else:
                    result = fn(text)
            except Exception as exc:                      # noqa: BLE001
                error = f"{type(exc).__name__}: {exc}"[:60]
                break
            if result.get("applicable", True):
                applied += 1
            if _reported_something(result):
                reported += 1
        rows.append({"check": short, "applied": applied,
                     "reported": reported, "error": error})
    return rows


def render_audit(rows: list[dict]) -> str:
    lines = ["%-46s %7s %8s" % ("check", "applied", "reported")]
    for row in rows:
        note = row["error"] or ("  never reported" if row["reported"] == 0 else "")
        lines.append("%-46s %7d %8d%s"
                     % (row["check"], row["applied"], row["reported"], note))
    silent = [r["check"] for r in rows if not r["error"] and r["reported"] == 0]
    lines.append("")
    lines.append("silent on every document: %d" % len(silent))
    for check in silent:
        lines.append("  " + check)
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-baseline", action="store_true",
                        help="record the current counts as the expected ones")
    parser.add_argument("--audit", action="store_true",
                        help="report which checks never say anything")
    args = parser.parse_args()

    manifest = manifest_path()
    if manifest is None:
        print(f"set {MANIFEST_ENV} to a corpus manifest; nothing to sweep")
        return 2

    documents = load_corpus(manifest)
    if args.audit:
        print(render_audit(audit(documents)))
        return 0

    results = sweep(documents)
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
