from radia_mcp.presentation.cross_lint import presentation_check_notation_variants


def _findings_of_type(result: dict, finding_type: str) -> list[dict]:
    return [
        finding
        for finding in result["findings"]
        if finding["type"] == finding_type
    ]


def test_case_check_ignores_lowercase_host_name_in_url() -> None:
    result = presentation_check_notation_variants(
        "The source is on GitHub: https://github.com/example/project"
    )

    assert _findings_of_type(result, "case_variants") == []


def test_case_check_still_reports_prose_variants() -> None:
    result = presentation_check_notation_variants(
        "NGSolve provides the solver, but ngsolve is misspelled here."
    )

    findings = _findings_of_type(result, "case_variants")
    assert len(findings) == 1
    assert findings[0]["examples"][0]["variants"] == ["NGSolve", "ngsolve"]


def test_hyphen_variants_require_context_review() -> None:
    result = presentation_check_notation_variants(
        "The full-band model covers the full band."
    )

    findings = _findings_of_type(result, "hyphen_variants")
    assert len(findings) == 1
    assert findings[0]["requires_context_review"] is True
