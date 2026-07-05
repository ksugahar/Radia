"""Unit tests for the multilingual additions to chroma_retriever.

Tests the pure-Python helpers (filename language heuristic + chapter
regex) that do NOT require chromadb / sentence-transformers / pymupdf
to be installed. The CI runner does not have those optional deps.

Ported behavior origin: wjc9011/COMSOL_Multiphysics_MCP fork
  src/knowledge/pdf_processor.py::detect_filename_language
  src/knowledge/pdf_processor.py::PDFProcessor.CHAPTER_PATTERNS
"""

import pytest

from radia_mcp.common.chroma_retriever import (
    detect_filename_language,
    find_chapters,
    CHAPTER_PATTERNS,
)


# ----- detect_filename_language -----------------------------------------

@pytest.mark.parametrize("filename, expected", [
    # English -- only ASCII letters
    ("MaxwellBook.pdf", "en"),
    ("ACDC_Module_Users_Guide.pdf", "en"),
    ("hollaus_2016_eddy_current.pdf", "en"),

    # Japanese -- hiragana / katakana present
    ("ヒステリシス測定.pdf", "ja"),
    ("第1章 電磁界の基礎.pdf", "ja"),  # hiragana の missing but kanji+CJK punct
    ("有限要素法.pdf", "zh"),          # kanji only -> classified zh (no kana)
    ("カタカナ表記.pdf", "ja"),
    ("ひらがな文.pdf", "ja"),

    # Chinese-ish: only CJK ideographs, no kana
    ("电磁场理论.pdf", "zh"),
    ("有限元方法.pdf", "zh"),

    # Edge cases
    ("", None),
    ("   .pdf", None),
    ("123.pdf", None),                  # digits only
])
def test_detect_filename_language(filename, expected):
    assert detect_filename_language(filename) == expected


def test_detect_filename_language_mixed_kanji_kana_is_ja():
    """If both hiragana/katakana AND kanji appear, the kana wins
    (kanji alone is ambiguous between ja and zh)."""
    assert detect_filename_language("論文の概要.pdf") == "ja"
    assert detect_filename_language("テストデータ.pdf") == "ja"


def test_detect_filename_language_strips_path():
    """Heuristic should look at the basename stem, not the directory."""
    assert detect_filename_language(
        "C:/temp/english_notes.pdf") == "en"
    # NOTE: full-path with Chinese filename
    assert detect_filename_language(
        "C:/temp/测试案例.pdf") == "zh"


# ----- find_chapters ----------------------------------------------------

def test_find_chapters_english():
    text = "\n".join([
        "Chapter 1: Introduction",
        "Some intro text here.",
        "More intro text.",
        "Chapter 2: Maxwell Equations",
        "Faraday's law states...",
    ])
    chapters = find_chapters(text)
    titles = [c[0] for c in chapters]
    assert "Introduction" in titles
    assert "Maxwell Equations" in titles


def test_find_chapters_japanese_dai_n_sho():
    """Japanese 第N章 pattern detection (ported from COMSOL fork)."""
    text = "\n".join([
        "第1章 電磁界の基礎",
        "本書では電磁界を扱う。",
        "第2章 有限要素法",
        "差分方程式は次のように離散化される。",
    ])
    chapters = find_chapters(text)
    titles = [c[0] for c in chapters]
    # Japanese chapter titles should be extracted as the trailing text
    assert "電磁界の基礎" in titles
    assert "有限要素法" in titles


def test_find_chapters_japanese_n_sho_omitted_prefix():
    """'第' prefix-omitted form: '1章', '2章'."""
    text = "\n".join([
        "1章 概要",
        "概要を述べる。",
        "2章 詳細",
        "詳細を述べる。",
    ])
    chapters = find_chapters(text)
    titles = [c[0] for c in chapters]
    assert "概要" in titles
    assert "詳細" in titles


def test_find_chapters_chinese_kanji_numeral():
    """Chinese kanji-numeral chapter form: 第一章, 第二章."""
    text = "\n".join([
        "第一章 引言",
        "本书介绍电磁场。",
        "第二章 麦克斯韦方程",
        "麦克斯韦方程组如下。",
    ])
    chapters = find_chapters(text)
    titles = [c[0] for c in chapters]
    # Kanji-numeral chapter titles -> trailing text captured
    assert "引言" in titles
    assert "麦克斯韦方程" in titles


def test_find_chapters_no_headings_yields_introduction():
    """Text with no recognized chapter headings collapses into one
    'Introduction' chapter spanning the whole text."""
    text = "Plain paragraph with no chapter heading at all."
    chapters = find_chapters(text)
    assert len(chapters) == 1
    assert chapters[0][0] == "Introduction"


def test_chapter_patterns_constant_is_list():
    """CHAPTER_PATTERNS is the public regex list; callers may iterate
    it directly if they need custom chunking logic."""
    assert isinstance(CHAPTER_PATTERNS, list)
    assert len(CHAPTER_PATTERNS) >= 7   # 3 EN + 3 JA + 1 ZH


# ----- common/__init__.py re-exports ------------------------------------

def test_common_init_reexports():
    """The new helpers must be importable from radia_mcp.common
    directly (not just from .chroma_retriever)."""
    from radia_mcp.common import (
        ChromaRetriever,
        extract_pdf_chunks,
        detect_filename_language,
        find_chapters,
        CHAPTER_PATTERNS,
        DEFAULT_EMBEDDING_MODEL,
    )
    assert detect_filename_language("test.pdf") == "en"
    assert callable(extract_pdf_chunks)
    assert callable(find_chapters)
    assert isinstance(DEFAULT_EMBEDDING_MODEL, str)
