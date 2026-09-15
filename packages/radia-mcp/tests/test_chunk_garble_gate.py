"""Unit test for the ingestion quality gate: chunk_garble_fraction.

Guards the Zienkiewicz-class regression -- a born-digital PDF with a corrupt
CID/ToUnicode text layer extracts non-empty but meaningless text. The ingestion
path (literature_build_vector_index) uses chunk_garble_fraction to FLAG such
books for re-OCR instead of silently storing garbage. Portable: no NAS, no OCR,
no embedding model -- pure string heuristic.
"""
from radia_mcp.common.chroma_retriever import (
    chunk_garble_fraction, _text_readability)
import pytest


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


@pytest.mark.parametrize("damage", ["\x00", "\x01", "\x1f", "\x7f", "\x85", "\ufffd", "(cid:123)"])
def test_cjk_cannot_mask_explicit_text_layer_damage(damage):
    text = "有限要素法による電磁界解析と渦電流損失の評価。" * 20 + damage
    assert _text_readability(text) == 0.0
    assert chunk_garble_fraction([{"text": text}]) == 1.0


def test_clean_japanese_and_layout_whitespace_remain_readable():
    text = "有限要素法による電磁界解析。\t磁束密度の評価。\r\n渦電流損失。\f"
    assert _text_readability(text) > 0.8
    assert chunk_garble_fraction([{"text": text}]) == 0.0


def test_peec_catalog_withholds_damaged_extracted_metadata(monkeypatch):
    from radia_mcp.peec import bibliography_index_knowledge as catalog

    entries = [dict(filename="source.pdf", title="正常な題目", authors="著者",
                    abstract_excerpt="正常に見える本文" * 20 + "\x01", page_count=3)]
    monkeypatch.setattr(catalog, "CATALOG_ENTRIES", entries)
    result = catalog.get_bibliography_index()
    assert "source.pdf" in result and "Re-OCR" in result
    assert "\x01" not in result and "**Abstract**" not in result


def test_shipped_peec_catalog_never_returns_control_character_damage():
    from radia_mcp.peec.bibliography_index_knowledge import get_bibliography_index

    result = get_bibliography_index()
    assert not any((ord(c) < 32 and c not in "\t\n\r\f") or 127 <= ord(c) <= 159 for c in result)
