"""Tier 4 — poster_qr_audit.

Audits a poster's QR-code usage against the BetterPosters 2025-05 best
practices: prominent size, concrete label, no end-of-poster ghetto.
"""
from __future__ import annotations

import re
import urllib.error
import urllib.request

from .._textutil import read_tex


_QRCODE_RE = re.compile(
    r"\\qrcode\b(?:\[(?P<opts>[^\]]*)\])?\{(?P<url>[^}]+)\}"
)
_BAD_LABEL_PHRASES = (
    "scan for more", "scan to discover more", "more info",
    "詳細はこちら", "詳しくは",
)
_GOOD_LABEL_HINTS = (
    "paper", "preprint", "data", "supplementary", "code",
    "論文", "プレプリント", "データ", "資料", "コード",
)


def _check_url_reachable(url: str, timeout: float = 5.0) -> tuple[bool, str]:
    """Try HEAD request; return (ok, hint)."""
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            code = resp.status
            return (200 <= code < 400, f"HTTP {code}")
    except urllib.error.HTTPError as e:
        return (False, f"HTTP {e.code}")
    except (urllib.error.URLError, ValueError, TimeoutError, OSError) as e:
        return (False, f"unreachable: {type(e).__name__}")


def poster_qr_audit(tex_path: str, check_url: bool = False) -> str:
    """Audit QR code(s) in a poster for prominence + labeling + URL reachability.

    Checks (per BetterPosters 2025-05):

    - QR size ≥ 2.5cm (1 inch) — anything smaller goes ignored
    - Label is concrete (mentions paper/data/code/preprint), not "Scan
      for more"
    - QR is not jammed into the very bottom edge without surrounding
      whitespace (heuristic: label exists and isn't trivial)
    - Optional: HTTP HEAD reaches the URL (off by default to keep the
      tool offline-friendly)

    Args:
        tex_path: Path to the poster ``.tex`` source.
        check_url: If True, issue a HEAD request to each QR URL.

    Returns:
        Per-QR audit report.
    """
    try:
        path, tex = read_tex(tex_path)
    except FileNotFoundError as e:
        return f"Error: file not found: {e}"

    qrs = list(_QRCODE_RE.finditer(tex))
    if not qrs:
        return f"poster_qr_audit: {path}\n  no \\qrcode found — consider poster_qr_inject"

    lines = [f"poster_qr_audit: {path}", f"  {len(qrs)} QR code(s) found"]
    fail = 0
    for i, m in enumerate(qrs, 1):
        url = m.group("url").strip()
        opts = (m.group("opts") or "").lower()
        size_cm = None
        size_m = re.search(r"height\s*=\s*(\d+(?:\.\d+)?)\s*cm", opts)
        if size_m:
            size_cm = float(size_m.group(1))
        # Find label: look forward up to 200 chars for a fontsize or text node.
        tail = tex[m.end():m.end() + 300]
        label_m = re.search(r"\\fontsize\{\d+\}\{\d+\}\\selectfont\s+([^\\\n]+)", tail)
        label = label_m.group(1).strip() if label_m else ""
        problems = []
        if size_cm is not None and size_cm < 2.5:
            problems.append(f"size {size_cm}cm < 2.5cm minimum")
        if size_cm is None:
            problems.append("size not specified (defaults are usually too small)")
        if not label:
            problems.append("no label adjacent to QR — viewer cannot tell what it does")
        else:
            low = label.lower()
            if any(p in low for p in _BAD_LABEL_PHRASES):
                problems.append(f"vague label: {label!r} — be concrete")
            elif not any(p in low for p in _GOOD_LABEL_HINTS):
                problems.append(f"label {label!r} lacks concrete deliverable hint")
        if check_url:
            ok, hint = _check_url_reachable(url)
            if not ok:
                problems.append(f"URL probe failed: {hint}")
        verdict = "OK " if not problems else "FAIL"
        if problems:
            fail += 1
        lines.append(f"  [{verdict}] QR {i}: url={url[:60]}{'…' if len(url) > 60 else ''}")
        for p in problems:
            lines.append(f"      - {p}")
    if fail:
        lines.append(f"FAIL — {fail}/{len(qrs)} QR(s) need fixes")
    else:
        lines.append("PASS — QR(s) meet BetterPosters prominence + labeling")
    return "\n".join(lines)
