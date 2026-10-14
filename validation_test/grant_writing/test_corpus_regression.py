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
import pathlib

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


def test_sources_match_the_adjudicated_baseline(corpus):
    drift = []
    for document in corpus["documents"]:
        expected = corpus["baseline"].get(document["label"])
        if expected is None:
            continue
        actual_hash = sweep.measure(document)["source_sha256"]
        expected_hash = expected.get("source_sha256")
        if actual_hash != expected_hash:
            drift.append(document["label"])

    assert not drift, (
        "proposal source changed since its findings were adjudicated; inspect "
        "the sweep, then re-baseline with sweep.py --write-baseline:\n  "
        + "\n  ".join(drift)
    )


def test_finding_counts_match_the_adjudicated_baseline(corpus):
    drift = []
    for document in corpus["documents"]:
        expected = corpus["baseline"].get(document["label"])
        if expected is None:
            continue
        actual = sweep.measure(document)
        if actual["finding_count"] != expected["finding_count"]:
            drift.append(
                f"{document['label']}: {actual['finding_count']} findings, "
                f"baseline {expected['finding_count']}"
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
                    f"{document['label']}: {field} now {actual.get(field)}, "
                    f"baseline {allowance}"
                )

    assert not drift, (
        "page usage moved; rebuild the PDF and re-baseline if intended:\n  "
        + "\n  ".join(drift)
    )


def test_manifest_can_join_multiple_live_sources(tmp_path: pathlib.Path):
    (tmp_path / "purpose.tex").write_text("研究目的", encoding="utf-8")
    (tmp_path / "plan.tex").write_text("研究計画", encoding="utf-8")
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps({
            "documents": [{
                "label": "live-draft",
                "paths": ["purpose.tex", "plan.tex"],
            }],
        }),
        encoding="utf-8",
    )

    document = sweep.load_corpus(manifest)[0]

    assert sweep.read_document_text(document) == "研究目的\n\n研究計画"


def test_manifest_validates_private_provenance_files(tmp_path: pathlib.Path):
    (tmp_path / "proposal.txt").write_text("研究目的", encoding="utf-8")
    (tmp_path / "submitted.pdf").write_bytes(b"submission")
    (tmp_path / "result.pdf").write_bytes(b"result")
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps({
            "documents": [{
                "label": "adopted",
                "path": "proposal.txt",
                "outcome": "adopted",
                "source_paths": ["submitted.pdf"],
                "outcome_basis": "award notice",
                "outcome_evidence": ["result.pdf"],
            }],
        }),
        encoding="utf-8",
    )

    document = sweep.load_corpus(manifest)[0]

    assert document["source_paths"] == [tmp_path / "submitted.pdf"]
    assert document["outcome_evidence"] == [tmp_path / "result.pdf"]
    assert document["outcome_basis"] == "award notice"


@pytest.mark.parametrize(
    "entry",
    [
        {"label": "neither"},
        {"label": "both", "path": "one.tex", "paths": ["two.tex"]},
        {"label": "empty", "paths": []},
    ],
)
def test_manifest_requires_exactly_one_nonempty_source_form(
    tmp_path: pathlib.Path,
    entry: dict,
):
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps({"documents": [entry]}),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="exactly one|no source paths"):
        sweep.load_corpus(manifest)


def test_audit_does_not_read_measurement_metadata_as_a_finding():
    result = {
        "applicable": False,
        "score": None,
        "statement_count": 2,
        "statements": [{"text": "中心の問い"}],
        "risks": [],
    }

    assert not sweep._reported_something(result)
    assert sweep._reported_something({"violation_count": 1})


def test_audit_respects_program_specific_checks():
    assert sweep._eligible_for_program("kaken_oss_platform_check", "kaken_oss")
    assert not sweep._eligible_for_program(
        "kaken_oss_platform_check", "kaken_generic"
    )
    assert not sweep._eligible_for_program("kaken_oss_platform_check", "generic")
    assert sweep._eligible_for_program("kddi_digital_check", "kddi_digital")
    assert not sweep._eligible_for_program("kddi_digital_check", "kaken_oss")
    assert sweep._eligible_for_program("analyze_sentences", "generic")


def test_outcome_comparison_normalizes_and_counts_document_prevalence(monkeypatch):
    documents = [
        {"label": "a", "outcome": "adopted", "program": "kaken_generic"},
        {"label": "r", "outcome": "rejected", "program": "kaken_generic"},
        {"label": "u", "outcome": "unsubmitted", "program": "kaken_generic"},
    ]
    measured = {
        "a": {
            "prose_chars": 1000,
            "finding_count": 1,
            "patterns": {"check/rule": 1},
        },
        "r": {
            "prose_chars": 2000,
            "finding_count": 4,
            "patterns": {"check/rule": 4},
        },
    }
    monkeypatch.setattr(sweep, "sweep", lambda selected: measured)

    comparison = sweep.compare_outcomes(documents, "kaken_generic")

    assert comparison["groups"]["adopted"]["findings_per_10000_chars"] == 10.0
    assert comparison["groups"]["rejected"]["findings_per_10000_chars"] == 20.0
    assert comparison["groups"]["rejected"]["pattern_document_prevalence"] == {
        "check/rule": 1
    }
