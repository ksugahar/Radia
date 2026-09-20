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

Two properties are load bearing and each has a regression test.  The site
suffix is OPTIONAL, because TeX also prints a bare header whose context
follows on later lines; requiring a known site would drop those.  And the
separators are ``\\s+`` rather than single spaces, because TeX wraps the
source line number onto the next log line, which would otherwise match no
alternative at all.
"""
from __future__ import annotations

import re

# The severity group keeps whatever TeX printed ("12.34pt too wide",
# "5.0pt too high") verbatim; callers display it rather than parse it.
_OVERFULL_RE = re.compile(
    r"Overfull \\(?P<box>[hv])box \((?P<severity>[^)]+)\)"
    r"(?:"
    r"\s+in (?P<site>paragraph|alignment) at\s+lines\s+(?P<start>\d+)(?:--(?P<end>\d+))?"
    r"|\s+detected at\s+line\s+(?P<at>\d+)"
    r"|\s+has occurred while \\output is active"
    r")?"
)


def scan_overfull(text: str) -> list[dict]:
    """Return one record per overfull box found in a LaTeX log.

    Each record carries ``box`` ("hbox"/"vbox"), ``severity`` as TeX printed
    it, ``site`` ("paragraph"/"alignment"/"detected"/"output"/"") and ``lines``
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
        elif m.group(0).rstrip().endswith("active"):
            # \output is active: TeX names the page, never a source line.
            site = "output"
            lines = ""
        else:
            # Bare header; TeX puts the context on following lines.
            site = ""
            lines = ""
        records.append({
            "box": m.group("box") + "box",
            "severity": m.group("severity"),
            "site": site,
            "lines": lines,
        })
    return records


def summarize_overfull(text: str, max_details: int = 20,
                       boxes: tuple[str, ...] | None = None) -> dict:
    """Count overfull boxes by type and return the first ``max_details``.

    ``boxes`` restricts which box kinds ``overfull_count`` and
    ``overfull_details`` describe.  ``paper_writing_check_overfull_hbox``
    passes ``("hbox",)`` to keep its published count meaning the horizontal
    boxes it is named for; a caller that wants both passes nothing.  A beamer
    frame overflowing its slide reports ``\\vbox`` and is a real defect, so
    ``hbox_count`` / ``vbox_count`` always report the full breakdown.
    """
    records = scan_overfull(text)
    selected = records if boxes is None else [r for r in records if r["box"] in boxes]
    return {
        "overfull_count": len(selected),
        "hbox_count": sum(1 for r in records if r["box"] == "hbox"),
        "vbox_count": sum(1 for r in records if r["box"] == "vbox"),
        "overfull_details": selected[:max_details],
    }
