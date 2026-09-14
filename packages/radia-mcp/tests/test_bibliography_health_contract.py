"""Read-only bibliography diagnostics must not turn failures into success."""
import pytest

from radia_mcp.bibliography.plans import T11_health_report as health
from radia_mcp.bibliography.plans.T5_dedupe import bibliography_dedupe
from radia_mcp.bibliography.plans.T7_lint import bibliography_lint
from radia_mcp.bibliography.plans.T9_self_citation_ratio import bibliography_self_citation_ratio
from radia_mcp.bibliography.plans.T12_year_distribution import bibliography_year_distribution


@pytest.mark.parametrize("mode", ["missing", "directory", "invalid_utf8", "nonbib", "malformed"])
def test_health_does_not_score_unreadable_or_invalid_input(tmp_path, mode):
    path = tmp_path / "references.bib"
    if mode == "directory":
        path.mkdir()
    elif mode != "missing":
        path.write_bytes({"invalid_utf8": b"\xff", "nonbib": b"not bibliography",
                          "malformed": b"@article{old,title={unfinished"}[mode])
    result = health.bibliography_health_report(str(path))
    assert "VERDICT: UNAVAILABLE" in result
    assert "/100" not in result and "GOOD" not in result


@pytest.mark.parametrize("axis", ["lint", "dedupe"])
@pytest.mark.parametrize("mode", ["exception", "empty", "none", "error"])
def test_health_propagates_detector_failure(monkeypatch, axis, mode):
    monkeypatch.setattr(health, "bibliography_lint", lambda _: "bibliography_lint: test\nPASS")
    monkeypatch.setattr(health, "bibliography_dedupe", lambda _: "bibliography_dedupe: test\nPASS")
    def failure(_):
        if mode == "exception":
            raise RuntimeError("detector failed")
        return {"empty": "", "none": None, "error": "Error: failed"}[mode]
    monkeypatch.setattr(health, "bibliography_" + axis, failure)
    result = health.bibliography_health_report("unused")
    assert "UNAVAILABLE" in result and "/100" not in result and axis in result


@pytest.mark.parametrize("kind,author,editor,missing", [
    ("article", "", "Editor, A", True), ("book", "", "Editor, A", False),
    ("book", "", "", True), ("article", "Author, A", "", False)])
def test_only_nonempty_book_editor_substitutes_for_author(tmp_path, kind, author, editor, missing):
    path = tmp_path / "references.bib"
    path.write_text(f"@{kind}{{author2024magnetic,author={{{author}}},editor={{{editor}}},"
                    "title={Magnetic fields},year=2024,journal={Journal},publisher={Publisher}}", encoding="utf-8")
    result = bibliography_lint(str(path))
    assert ("missing required field 'author'" in result) is missing


def test_duplicate_keys_are_high_severity_even_with_distinct_papers(tmp_path):
    path = tmp_path / "references.bib"
    path.write_text("@misc{same,title={Magnetic theory}}\n@misc{same,title={Acoustic experiments}}", encoding="utf-8")
    result = bibliography_dedupe(str(path))
    assert "[HIGH] Citation-key collision" in result and "PASS" not in result


def test_doi_duplicates_share_url_normalization(tmp_path):
    path = tmp_path / "references.bib"
    path.write_text("@misc{one,doi={10.1234/a#b?c%d}}\n"
                    "@misc{two,doi={HTTPS://DX.DOI.ORG/10.1234/a%23b%3Fc%25d}}", encoding="utf-8")
    result = bibliography_dedupe(str(path))
    assert "[HIGH] DOI duplicates" in result and "PASS" not in result


def test_japanese_titles_are_not_erased_before_comparison(tmp_path):
    path = tmp_path / "references.bib"
    path.write_text("@misc{one,title={電磁界解析と最適化}}\n@misc{two,title={電磁界解析と最適化}}", encoding="utf-8")
    assert "[MEDIUM]" in bibliography_dedupe(str(path))


def test_healthy_file_still_scores_success(tmp_path):
    path = tmp_path / "references.bib"
    path.write_text("@misc{doe2024magnetic,author={Doe, Jane},title={Magnetic fields},year=2024}", encoding="utf-8")
    result = health.bibliography_health_report(str(path))
    assert "100/100" in result and "VERDICT: GOOD" in result


@pytest.mark.parametrize("matched,total,severity", [(0, 4, "PASS"), (1, 4, "PASS"),
    (3, 10, "[MEDIUM]"), (4, 10, "[MEDIUM]"), (5, 10, "[HIGH]")])
def test_self_citation_threshold_order_and_boundaries(tmp_path, matched, total, severity):
    path = tmp_path / "references.bib"
    path.write_text("\n".join(f"@misc{{key{i},author={{{'Doe' if i < matched else 'Roe'}, Jane}}}}"
                               for i in range(total)), encoding="utf-8")
    result = bibliography_self_citation_ratio(str(path), "Doe")
    assert f"matches: {matched}/{total}" in result and severity in result
    assert "rejection" not in result and "Surname-match screening only" in result


@pytest.mark.parametrize("author,needle", [
    ("{Research and Development Group} AND Doe, Jane", "Research and Development Group"),
    ("Müller, Hans", "Muller"), ("Ludwig van Beethoven", "van Beethoven"),
    ("山田, 太郎", "山田")])
def test_self_citation_preserves_name_boundaries(tmp_path, author, needle):
    path = tmp_path / "references.bib"
    path.write_text(f"@misc{{one,author={{{author}}}}}", encoding="utf-8")
    assert "matches: 1/1" in bibliography_self_citation_ratio(str(path), needle)


@pytest.mark.parametrize("needle", ["", " ", "{}"])
def test_empty_self_citation_target_is_error(tmp_path, needle):
    path = tmp_path / "references.bib"
    path.write_text("@misc{one,author={Doe, Jane}}", encoding="utf-8")
    assert bibliography_self_citation_ratio(str(path), needle).startswith("Error:")


def test_year_distribution_excludes_future_and_malformed_years(tmp_path):
    import datetime
    now = datetime.date.today().year
    path = tmp_path / "references.bib"
    path.write_text("\n".join(f"@misc{{key{i},year={{{year}}}}}" for i, year in enumerate(
        [str(now), str(now - 20), str(now + 1), "2024abc", "0000", ""])), encoding="utf-8")
    result = bibliography_year_distribution(str(path))
    assert "total dated entries: 2" in result
    assert "in last 5 years:  1  (50%)" in result
    assert "unknown/invalid years: 3; future years: 1" in result
    assert "INCOMPLETE" in result and "[MEDIUM]" not in result


def test_future_only_years_have_no_recency_score(tmp_path):
    path = tmp_path / "references.bib"
    path.write_text("@misc{one,year={9999}}", encoding="utf-8")
    result = bibliography_year_distribution(str(path))
    assert "no valid non-future publication years" in result
    assert "future years: 1" in result and "in last 5" not in result
