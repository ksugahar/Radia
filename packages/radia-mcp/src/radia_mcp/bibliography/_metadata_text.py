"""Render plain external metadata, not arbitrary TeX, as BibTeX text fields."""

import html
import re


def metadata_text(value: str) -> str:
    """Decode entities before removing simple emphasis; preserve no hidden TeX."""
    for _ in range(3):
        decoded = html.unescape(value)
        if decoded == value:
            break
        value = decoded
    value = re.sub(r"</?(?:jats:)?(?:i|b|em|strong|italic|bold)\s*>", "", value,
                   flags=re.IGNORECASE)
    if re.search(r"</?[A-Za-z][^>]*>|\\|\$", value):
        raise ValueError("TeX/math or unsupported metadata markup requires manual verification")
    return value.strip()


def bibtex_text(value: str) -> str:
    """Escape literal TeX special characters once at the external-data boundary."""
    replacements = {"&": r"\&", "%": r"\%", "_": r"\_", "#": r"\#",
                    "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}",
                    "^": r"\textasciicircum{}"}
    return "".join(replacements.get(char, char) for char in metadata_text(value))
