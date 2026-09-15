"""Text-only regressions for the public sentence-ending diagnostic (P19)."""

from importlib import import_module

import pytest

# Avoid static traversal into unrelated optional PDF dependencies in tools.py.
check_endings = import_module(
    "radia_mcp.paper_writing.tools"
).paper_writing_check_sentence_ending_variety


@pytest.mark.parametrize("ending", [
    "である", "だ", "ある", "ない", "された", "示した", "考えられる",
    "した", "する", "される", "られる", "思われる", "示す", "得た", "得る",
])
def test_specific_sentence_endings_are_not_shadowed(ending):
    result = check_endings(f"対象は{ending}。")
    assert result["total_sentences"] == 1
    assert result["ending_histogram"] == {ending: 1}
    assert result["critical_runs"] == []


def test_decimal_values_do_not_split_sentences_or_shift_run_indices():
    result = check_endings(
        "値は1.25である。値は2.50である。値は3.75である。結果を示した。"
    )
    assert result["total_sentences"] == 4
    assert result["ending_histogram"] == {"である": 3, "示した": 1}
    assert result["critical_runs"] == [{
        "ending": "である", "start_sentence_index": 0, "length": 3,
        "sample": "値は1.25である。",
    }]
    assert result["shannon_entropy_bits"] == pytest.approx(0.811, abs=0.001)


def test_figure_abbreviation_does_not_add_a_sentence():
    result = check_endings("Fig. 3の結果を示した。Eq. 2の結果を示した。結果を示した。")
    assert result["total_sentences"] == 3
    assert result["ending_histogram"] == {"示した": 3}
    assert result["critical_runs"][0]["length"] == 3


@pytest.mark.parametrize("middle", ["独立文です。", "Independent text. "])
def test_unclassified_sentence_breaks_a_run_without_becoming_critical(middle):
    result = check_endings(f"例である。例である。{middle}例である。例である。")
    assert result["total_sentences"] == 5
    assert result["critical_runs"] == []


def test_configured_threshold_controls_run_detection():
    text = "例である。例である。"
    assert check_endings(text)["critical_runs"] == []
    assert check_endings(text, consecutive_threshold=2)["critical_runs"][0]["length"] == 2
