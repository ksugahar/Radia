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
     "outcome": "adopted", "program": "kaken_oss"},
    {"label": "current-draft",
     "paths": ["draft/purpose.tex", "draft/abilities.tex"],
     "pdf": "draft/proposal.pdf", "program": "kaken_oss"}
  ]
}
```

Relative paths resolve against the manifest's directory. ``outcome`` is
recorded but never scored: four measurements have found no relationship
between these checks and adoption.

Usage::

    python validation_test/grant_writing/sweep.py
    python validation_test/grant_writing/sweep.py --compare-outcomes
    python validation_test/grant_writing/sweep.py --write-baseline
"""

from __future__ import annotations

import argparse
import collections
import hashlib
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


def _resolve_path(raw: str, manifest: pathlib.Path) -> pathlib.Path:
    path = pathlib.Path(raw)
    if not path.is_absolute():
        path = manifest.parent / path
    return path


def _resolve_file_list(
    entry: dict,
    key: str,
    manifest: pathlib.Path,
) -> list[pathlib.Path]:
    raw_paths = entry.get(key, [])
    if not isinstance(raw_paths, list):
        raise TypeError(
            f"corpus document {entry.get('label', '?')!r} field {key!r} "
            "must be a list"
        )
    paths = [_resolve_path(raw, manifest) for raw in raw_paths]
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(f"corpus provenance file not found: {path}")
    return paths


def load_corpus(manifest: pathlib.Path) -> list[dict]:
    data = json.loads(manifest.read_text(encoding="utf-8"))
    documents = []
    for entry in data.get("documents", []):
        has_path = "path" in entry
        has_paths = "paths" in entry
        if has_path == has_paths:
            raise ValueError(
                f"corpus document {entry.get('label', '?')!r} must name exactly "
                "one of 'path' or 'paths'"
            )
        raw_paths = [entry["path"]] if has_path else entry["paths"]
        if not isinstance(raw_paths, list) or not raw_paths:
            raise ValueError(
                f"corpus document {entry.get('label', '?')!r} has no source paths"
            )
        paths = [_resolve_path(raw, manifest) for raw in raw_paths]
        for path in paths:
            if not path.is_file():
                raise FileNotFoundError(f"corpus document not found: {path}")
        pdf = entry.get("pdf")
        if pdf:
            pdf = _resolve_path(pdf, manifest)
            if not pdf.is_file():
                raise FileNotFoundError(f"compiled proposal not found: {pdf}")
        source_paths = _resolve_file_list(entry, "source_paths", manifest)
        outcome_evidence = _resolve_file_list(
            entry,
            "outcome_evidence",
            manifest,
        )
        outcome_basis = entry.get("outcome_basis", "")
        if not isinstance(outcome_basis, str):
            raise TypeError(
                f"corpus document {entry.get('label', '?')!r} field "
                "'outcome_basis' must be a string"
            )
        documents.append({
            "label": entry["label"],
            "paths": paths,
            "pdf": pdf,
            "outcome": entry.get("outcome", "unknown"),
            "outcome_basis": outcome_basis,
            "outcome_evidence": outcome_evidence,
            "program": entry.get("program", "generic"),
            "source_paths": source_paths,
        })
    if not documents:
        raise ValueError(f"corpus manifest lists no documents: {manifest}")
    return documents


def read_document_text(document: dict) -> str:
    """Read one text snapshot or an ordered set of live proposal sources."""
    return "\n\n".join(
        path.read_text(encoding="utf-8") for path in document["paths"]
    )


def measure(document: dict) -> dict:
    text = read_document_text(document)
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
        "source_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
    }
    if pages is not None:
        measured["pages"] = pages
    return measured


def sweep(documents: list[dict]) -> dict:
    return {d["label"]: measure(d) for d in documents}


def render(results: dict) -> str:
    lines = [f"{'document':<28} {'prose':>8} {'findings':>7}"]
    for label, row in results.items():
        lines.append(
            f"{label:<28} {row['prose_chars']:>8d} {row['finding_count']:>7d}"
        )
    totals: collections.Counter = collections.Counter()
    for row in results.values():
        totals.update(row["patterns"])
    lines.append("")
    lines.append(
        f"patterns: {len(totals)}   total findings: {sum(totals.values())}"
    )
    for pattern, count in totals.most_common():
        lines.append(f"  {pattern:<52} {count:>3d}")
    return "\n".join(lines)


def compare_outcomes(documents: list[dict], program: str | None = None) -> dict:
    """Describe detector findings by recorded outcome without predicting it.

    Counts are normalized by prose length and also reduced to document-level
    prevalence. This prevents a repeated placeholder in one budget table from
    masquerading as evidence across several applications. The comparison is a
    corpus audit aid, not an adoption model: year, scheme, reviewer panel, and
    scientific maturity remain uncontrolled.
    """
    selected = [
        document
        for document in documents
        if document["outcome"] in {"adopted", "rejected"}
        and (program is None or document["program"] == program)
    ]
    measured = sweep(selected)
    groups = {}
    for outcome in ("adopted", "rejected"):
        group = [d for d in selected if d["outcome"] == outcome]
        total_prose = sum(measured[d["label"]]["prose_chars"] for d in group)
        total_findings = sum(
            measured[d["label"]]["finding_count"] for d in group
        )
        prevalence: collections.Counter = collections.Counter()
        occurrences: collections.Counter = collections.Counter()
        for document in group:
            patterns = measured[document["label"]]["patterns"]
            occurrences.update(patterns)
            prevalence.update(patterns.keys())
        groups[outcome] = {
            "document_count": len(group),
            "total_prose_chars": total_prose,
            "total_findings": total_findings,
            "findings_per_10000_chars": (
                round(total_findings * 10000 / total_prose, 2)
                if total_prose else 0.0
            ),
            "pattern_document_prevalence": dict(prevalence.most_common()),
            "pattern_occurrences": dict(occurrences.most_common()),
        }
    return {
        "program": program or "all",
        "groups": groups,
        "interpretation": (
            "descriptive corpus audit only; do not infer that a detector "
            "finding caused adoption or rejection"
        ),
    }


def render_outcome_comparison(comparisons: list[dict]) -> str:
    lines = [
        "Descriptive only: outcome is not a score and no causal inference is valid."
    ]
    for comparison in comparisons:
        lines.extend(["", f"scope: {comparison['program']}"])
        for outcome in ("adopted", "rejected"):
            row = comparison["groups"][outcome]
            lines.append(
                f"  {outcome:<8} n={row['document_count']:<2d} "
                f"findings/10k chars={row['findings_per_10000_chars']:.2f}"
            )
            common = list(row["pattern_document_prevalence"].items())[:8]
            if common:
                lines.append(
                    "    documents with pattern: "
                    + ", ".join(f"{name}={count}" for name, count in common)
                )
    return "\n".join(lines)


# Result shapes differ per check: some return risks, some issues, some a bare
# count. A check that never reports anything on any real proposal is either
# correctly quiet or quietly broken, and only the inventory tells them apart.
_FINDING_LISTS = (
    "risks", "issues", "findings", "over_threshold_examples", "variants",
    "misuses", "undefined", "undefined_acronyms", "violations",
    "weak_expressions", "top_fixes",
    "unbacked_absence_claims",
)
_FINDING_COUNTS = frozenset({
    "issue_count", "missing_count", "over_threshold_count", "risk_count",
    "total_findings", "total_matches", "total_weak_expressions",
    "undefined_count", "violation_count",
})
_NOT_DETECTORS = frozenset({
    "argument_evidence_map", "health_report", "usage",
    "recommendation_letter_template",
    "page_limit_check", "check_kanji_ratio",
})
_PROGRAM_ONLY = {
    "kaken_oss_platform_check": frozenset({"kaken_oss", "kaken_oss_platform"}),
    "kddi_digital_check": frozenset({"kddi_digital"}),
    "kddi_power_electronics_focus_check": frozenset({"kddi_digital"}),
}


def _reported_something(result) -> bool:
    if not isinstance(result, dict):
        return False
    if any(isinstance(result.get(k), list) and result[k] for k in _FINDING_LISTS):
        return True
    if any(result.get(key, 0) > 0 for key in _FINDING_COUNTS):
        return True
    score = result.get("score")
    return isinstance(score, (int, float)) and score < 10


def _eligible_for_program(check: str, program: str) -> bool:
    allowed = _PROGRAM_ONLY.get(check)
    return allowed is None or program in allowed


def audit(documents: list[dict]) -> list[dict]:
    """Report, per check, how often it applied and how often it said anything."""
    import inspect

    loaded = [(d, read_document_text(d)) for d in documents]
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
        eligible = applied = reported = 0
        error = ""
        for document, text in loaded:
            if not _eligible_for_program(short, document["program"]):
                continue
            eligible += 1
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
        rows.append({"check": short, "eligible": eligible, "applied": applied,
                     "reported": reported, "error": error})
    return rows


def render_audit(rows: list[dict]) -> str:
    lines = [
        f"{'check':<46} {'eligible':>8} {'applied':>7} {'reported':>8}"
    ]
    for row in rows:
        if row["error"]:
            note = row["error"]
        elif row["eligible"] == 0:
            note = "  no eligible documents"
        elif row["reported"] == 0:
            note = "  never reported"
        else:
            note = ""
        lines.append(
            f"{row['check']:<46} {row['eligible']:>8d} {row['applied']:>7d} "
            f"{row['reported']:>8d}{note}"
        )
    silent = [
        r["check"] for r in rows
        if not r["error"] and r["eligible"] > 0 and r["reported"] == 0
    ]
    lines.append("")
    lines.append(f"silent on every document: {len(silent)}")
    for check in silent:
        lines.append("  " + check)
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-baseline", action="store_true",
                        help="record the current counts as the expected ones")
    parser.add_argument("--audit", action="store_true",
                        help="report which checks never say anything")
    parser.add_argument(
        "--compare-outcomes",
        action="store_true",
        help="describe findings by recorded outcome without scoring adoption",
    )
    args = parser.parse_args()

    manifest = manifest_path()
    if manifest is None:
        print(f"set {MANIFEST_ENV} to a corpus manifest; nothing to sweep")
        return 2

    documents = load_corpus(manifest)
    if args.compare_outcomes:
        comparisons = [compare_outcomes(documents)]
        if any(d["program"] == "kaken_generic" for d in documents):
            comparisons.append(compare_outcomes(documents, "kaken_generic"))
        print(render_outcome_comparison(comparisons))
        return 0
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
