"""Lock the adjudicated finding counts on a corpus of real proposals.

Each count in the baseline was read and judged: every finding behind it was
inspected against its excerpt and kept only if a reader would agree the tool
had found a real defect. A change that moves a count is either a fix worth
re-baselining or a false positive coming back, and this test forces that
question to be asked rather than discovered a session later.

Skipped unless ``GRANT_WRITING_CORPUS`` names a manifest. The corpus lives
outside the repository: real proposals belong to their authors, several are
colleagues' work, and this repository is public.
"""

from __future__ import annotations

import json

import pytest

from validation_test.grant_writing import sweep


@pytest.fixture(scope="module")
def corpus():
    manifest = sweep.manifest_path()
    if manifest is None:
        pytest.skip(f"{sweep.MANIFEST_ENV} not set; no proposal corpus available")
    baseline = manifest.parent / sweep.BASELINE_NAME
    if not baseline.is_file():
        pytest.skip(
            f"no baseline beside {manifest.name}; run sweep.py --write-baseline"
        )
    return {
        "documents": sweep.load_corpus(manifest),
        "baseline": json.loads(baseline.read_text(encoding="utf-8")),
    }


def test_every_baselined_document_is_still_measured(corpus):
    measured = {d["label"] for d in corpus["documents"]}
    missing = sorted(set(corpus["baseline"]) - measured)

    assert not missing, f"baselined documents absent from the manifest: {missing}"


def test_finding_counts_match_the_adjudicated_baseline(corpus):
    drift = []
    for document in corpus["documents"]:
        expected = corpus["baseline"].get(document["label"])
        if expected is None:
            continue
        actual = sweep.measure(document)
        if actual["finding_count"] != expected["finding_count"]:
            drift.append(
                "%s: %d findings, baseline %d"
                % (document["label"], actual["finding_count"],
                   expected["finding_count"])
            )

    assert not drift, (
        "finding counts moved; adjudicate each change, then re-baseline with "
        "sweep.py --write-baseline:\n  " + "\n  ".join(drift)
    )


def test_no_new_finding_pattern_appears(corpus):
    known = {
        pattern
        for row in corpus["baseline"].values()
        for pattern in row["patterns"]
    }
    appeared = []
    for document in corpus["documents"]:
        for pattern, count in sweep.measure(document)["patterns"].items():
            if pattern not in known:
                appeared.append(f"{document['label']}: {pattern} ({count})")

    assert not appeared, (
        "a finding pattern not present when the corpus was adjudicated:\n  "
        + "\n  ".join(appeared)
    )


def test_no_field_outgrows_its_page_allowance(corpus):
    """The one defect that gets a proposal returned before it is read."""
    over = []
    for document in corpus["documents"]:
        if not document.get("pdf"):
            continue
        pages = sweep.measure(document).get("pages") or {}
        for field, (used, allowed) in pages.items():
            if used > allowed:
                over.append(f"{document['label']}: {field} {used}/{allowed}")

    assert not over, "a field runs past its allowance:\n  " + "\n  ".join(over)


def test_page_usage_matches_the_baseline(corpus):
    drift = []
    for document in corpus["documents"]:
        expected = (corpus["baseline"].get(document["label"]) or {}).get("pages")
        if not expected:
            continue
        actual = sweep.measure(document).get("pages") or {}
        for field, allowance in expected.items():
            if actual.get(field) != allowance:
                drift.append(
                    "%s: %s now %s, baseline %s"
                    % (document["label"], field, actual.get(field), allowance)
                )

    assert not drift, (
        "page usage moved; rebuild the PDF and re-baseline if intended:\n  "
        + "\n  ".join(drift)
    )
