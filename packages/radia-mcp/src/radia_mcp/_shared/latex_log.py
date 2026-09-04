"""LaTeX log scanning shared by the paper-writing and presentation lints.

TeX reports an overfull box from four different sites, and only one of them
says "in paragraph at lines".  A scanner keyed on that phrase silently counts
zero for a table row that does not fit, for a float TeX had to break out
("detected at line"), and for anything raised while the output routine runs
("has occurred while \\output is active", which carries no line number at
all).  For a check whose stated target is zero overfull boxes, undercounting
reads exactly like a pass, so the scanner has to cover every site.

Both `paper_writing_check_overfull_hbox` and
`presentation_check_overfull_hbox` carried the same one-site regex.  They
share this module instead so the two cannot drift apart again.
"""
from __future__ import annotations

import re

# The severity group keeps whatever TeX printed ("12.34pt too wide",
# "5.0pt too high") verbatim; callers display it rather than parse it.
_OVERFULL_RE = re.compile(
    r"Overfull \\(?P<box>[hv])box \((?P<severity>[^)]+)\)"
    r"(?:"
    r"\s+in (?P<site>paragraph|alignment) at lines (?P<start>\d+)(?:--(?P<end>\d+))?"
    r"|\s+detected at line (?P<at>\d+)"
    r"|\s+has occurred while \\output is active"
    r")"
)


def scan_overfull(text: str) -> list[dict]:
    """Return one record per overfull box found in a LaTeX log.

    Each record carries ``box`` ("hbox"/"vbox"), ``severity`` as TeX printed
    it, ``site`` ("paragraph"/"alignment"/"detected"/"output") and ``lines``
    ("10-12", "30", or "" when TeX reported no line).
    """
    records: list[dict] = []
    for m in _OVERFULL_RE.finditer(text):
        if m.group("site"):
            site = m.group("site")
            start, end = m.group("start"), m.group("end")
            lines = f"{start}-{end}" if end else start
        elif m.group("at"):
            site = "detected"
            lines = m.group("at")
        else:
            # \output is active: TeX names the page, never a source line.
            site = "output"
            lines = ""
        records.append({
            "box": m.group("box") + "box",
            "severity": m.group("severity"),
            "site": site,
            "lines": lines,
        })
    return records


def summarize_overfull(text: str, max_details: int = 20) -> dict:
    """Count overfull boxes by type and return the first ``max_details``.

    ``overfull_count`` is every overfull box, not only the horizontal ones:
    a beamer frame that overflows its slide and a journal float that
    overflows the text height both report ``\\vbox``, and both are defects.
    ``hbox_count`` / ``vbox_count`` keep the breakdown visible.
    """
    records = scan_overfull(text)
    return {
        "overfull_count": len(records),
        "hbox_count": sum(1 for r in records if r["box"] == "hbox"),
        "vbox_count": sum(1 for r in records if r["box"] == "vbox"),
        "overfull_details": records[:max_details],
    }
