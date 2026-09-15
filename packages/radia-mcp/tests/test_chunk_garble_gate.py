"""Unit test for the ingestion quality gate: chunk_garble_fraction.

Guards the Zienkiewicz-class regression -- a born-digital PDF with a corrupt
CID/ToUnicode text layer extracts non-empty but meaningless text. The ingestion
path (literature_build_vector_index) uses chunk_garble_fraction to FLAG such
books for re-OCR instead of silently storing garbage. Portable: no NAS, no OCR,
no embedding model -- pure string heuristic.
"""
from radia_mcp.common.chroma_retriever import (
    chunk_garble_fraction, _text_readability)


# Clean, real prose -> high readability (alnum + whitespace dominated).
CLEAN = [
    {"text": "The finite element method discretizes the domain into elements."},
    {"text": "Shape functions interpolate the field over each element."},
    {"text": "Assembly yields a sparse global stiffness matrix to solve."},
]
# Corrupt text layer: symbol/control-char soup -> readability ~0 (the visible
# glyphs render fine in a viewer, but get_text() returns this).
GARBLE = [
    {"text": "§¶†‡•◊▮╳№—…‹›«»¿¡×÷±∓§¶†‡•◊▮╳№—…‹›«»¿¡×÷±∓"},
    {"text": "◊▮╳№—…‹›«»¿¡×÷±∓§¶†‡•◊▮╳№—…‹›«»¿¡×÷±∓§¶†‡"},
    {"text": "»¿¡×÷±∓§¶†‡•◊▮╳№—…‹›«»¿¡×÷±∓§¶†‡•◊▮╳№—…‹›"},
]


def test_readability_separates_clean_from_garble():
    assert _text_readability(CLEAN[0]["text"]) > 0.80
    assert _text_readability(GARBLE[0]["text"]) < 0.10


def test_clean_corpus_is_not_flagged():
    assert chunk_garble_fraction(CLEAN) == 0.0


def test_corrupt_text_layer_is_fully_flagged():
    assert chunk_garble_fraction(GARBLE) == 1.0


def test_mixed_fraction_is_proportional():
    # 3 garble out of 6 -> 0.5; above the 0.40 ingestion-gate threshold.
    frac = chunk_garble_fraction(CLEAN + GARBLE)
    assert abs(frac - 0.5) < 1e-9
    assert frac >= 0.40


def test_empty_input_is_zero():
    assert chunk_garble_fraction([]) == 0.0


def test_missing_text_key_treated_as_garble():
    # Defensive: a chunk dict without "text" reads as empty -> garble.
    assert chunk_garble_fraction([{}]) == 1.0
