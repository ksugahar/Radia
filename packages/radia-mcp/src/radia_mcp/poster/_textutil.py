"""Shared LaTeX text extraction helpers for the poster sub-skill.

Stripping LaTeX down to readable text is needed by many lints
(word_count, line_length, caption_self_contained, …). Keep this in one
place so all the plans share the same definition of "body text".
"""
from __future__ import annotations

import pathlib
import re


_COMMENT_RE = re.compile(r"(?<!\\)%.*$", re.MULTILINE)
_INPUT_INCLUDE_RE = re.compile(r"\\(input|include|InputIfFileExists)\{([^}]+)\}")
_INCLUDEGRAPHICS_ARG = re.compile(r"\\includegraphics\b[^{]*\{[^}]*\}")
_BEGIN_END_BLOCK = re.compile(
    r"\\begin\{(tikzpicture|tabular|equation|equation\*|align|align\*|verbatim|comment|"
    r"tcolorbox\*?|pgfplots|axis|figure)\b[^\n]*\}.*?\\end\{\1\}",
    re.DOTALL,
)
_COMMAND_WITH_ARG = re.compile(r"\\[A-Za-z@]+\*?\s*(\[[^\]]*\])*\s*\{([^{}]*)\}")
_BARE_COMMAND = re.compile(r"\\[A-Za-z@]+\*?")
_BRACES = re.compile(r"[{}]")
_MULTI_SPACE = re.compile(r"[ \t]+")
_BLANK_LINES = re.compile(r"\n{3,}")


def read_tex(tex_path: str) -> tuple[pathlib.Path, str]:
    """Read a .tex file as UTF-8, returning the resolved path and source.

    Falls back to cp932 (Shift-JIS) for legacy Japanese files so the
    poster sub-skill works with files saved by Windows PowerPoint or
    historical platex tooling.
    """
    p = pathlib.Path(tex_path)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        raise FileNotFoundError(p)
    try:
        return p, p.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return p, p.read_text(encoding="cp932")


def strip_latex(tex: str) -> str:
    """Return a best-effort plain-text rendering of LaTeX source.

    The output is *good enough* for word-count and line-length lints —
    not a real LaTeX parser. It strips:

    - line comments
    - heavy environments (tikzpicture, tabular, equation, …)
    - ``\\command{arg}`` keeping ``arg``
    - bare ``\\command`` tokens
    - braces, multi-spaces, runs of blank lines
    """
    s = _COMMENT_RE.sub("", tex)
    # Drop heavy environments wholesale.
    s = _BEGIN_END_BLOCK.sub("", s)
    # Drop the includegraphics call entirely (caption is separate).
    s = _INCLUDEGRAPHICS_ARG.sub("", s)
    # Unwrap \cmd{arg} -> arg, iteratively for nested args.
    for _ in range(4):
        new = _COMMAND_WITH_ARG.sub(lambda m: m.group(2), s)
        if new == s:
            break
        s = new
    s = _BARE_COMMAND.sub("", s)
    s = _BRACES.sub("", s)
    s = _MULTI_SPACE.sub(" ", s)
    s = _BLANK_LINES.sub("\n\n", s)
    return s.strip()


def count_japanese_chars(text: str) -> int:
    """Count CJK characters (Hiragana / Katakana / Kanji) in ``text``."""
    return sum(
        1 for ch in text
        if "぀" <= ch <= "ヿ" or "一" <= ch <= "鿿"
    )


def count_words_mixed(text: str) -> tuple[int, int]:
    """Return ``(english_word_count, japanese_char_count)``.

    Japanese is character-counted (no word boundaries); English is
    whitespace-tokenized after CJK characters are removed.
    """
    jp = count_japanese_chars(text)
    en_text = re.sub(r"[぀-ヿ一-鿿、。・「」『』（）]+", " ", text)
    en_words = len([w for w in en_text.split() if any(c.isalpha() for c in w)])
    return en_words, jp


def split_lines(text: str, min_len: int = 4) -> list[str]:
    """Split into non-trivial lines (≥ ``min_len`` characters)."""
    return [ln.strip() for ln in text.splitlines() if len(ln.strip()) >= min_len]
