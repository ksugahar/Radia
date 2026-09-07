"""Explicit language and genre selection for the unified writing surface."""

from typing import Literal


def paper_writing_review_route(
    document_type: Literal["paper", "grant_proposal"],
    language: Literal["ja", "en"],
    venue_or_program: str = "",
) -> dict:
    """Choose genre/language-specific review criteria before reviewing a draft.

    This is routing guidance, not a scientific-quality score or submission gate.
    Venue/program requirements must be supplied and checked independently.
    """
    if document_type not in {"paper", "grant_proposal"} or language not in {"ja", "en"}:
        raise ValueError("Declare document_type=paper|grant_proposal and language=ja|en")
    result = {
        "document_type": document_type, "language": language,
        "venue_or_program": venue_or_program,
        "venue_requirements_verified": False,
        "shared_checks": ["notation consistency", "source traceability", "figure readability"],
        "aggregate_score": None,
        "policy": "Never average scores across languages or genres. Diagnostics are not acceptance decisions.",
    }
    if document_type == "paper":
        result.update(
            objective="Make the scientific claim, method, evidence and limitations reproducible.",
            criteria=["definitions and assumptions", "claim-evidence correspondence",
                      "method reproducibility", "prior-work comparison", "limitations"],
            language_focus=(
                ["subject-predicate clarity", "modifier scope", "term introduction", "paragraph logic"]
                if language == "ja" else
                ["explicit subjects and claims", "clause density", "nominalisation", "tense and paragraph logic"]
            ),
            tools=["paper_writing_bilingual_readability_check", "paper_writing_check_imrad_balance",
                   "paper_writing_check_abstract_no_math_no_citation"],
            excluded_scoring=["grant_writing_japanese_readability_score", "grant funding persuasiveness"],
        )
    else:
        result.update(
            objective="Explain why the proposed question matters and how the team can answer it.",
            criteria=["academic question and significance", "reviewer-visible motivation",
                      "preliminary evidence and feasibility", "deliverables and risks",
                      "schedule and budget", "program-specific review requirements"],
            tools=(["grant_writing_japanese_genre_contract", "grant_writing_japanese_readability_score"]
                   if language == "ja" else []),
            tool_arguments=({"document_type": "grant_proposal"} if language == "ja" else {}),
            language_scoring=("Japanese grant diagnostic" if language == "ja"
                              else "No validated English grant score; perform criterion-based review."),
            excluded_scoring=["paper IMRaD score", "completed-paper submission verdict"],
        )
    return result
