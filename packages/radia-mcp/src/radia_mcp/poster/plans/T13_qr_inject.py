"""Tier 4 — poster_qr_inject.

Injects a labeled QR code into a poster .tex. Uses the LaTeX
``qrcode`` package (no external image needed).

Source: BetterPosters blog (2025-05) — QR codes get ~50% scan rate when
prominent and labeled with concrete value ("Scan for full paper"), not
"Scan for more".
"""
from __future__ import annotations

import pathlib
import re

from .._textutil import read_tex


_QRCODE_PKG_LINE = "\\usepackage[hyperlinks=false]{qrcode}\n"


def poster_qr_inject(tex_path: str,
                     url: str,
                     label: str = "Scan for full paper",
                     size_cm: float = 4.0,
                     position: str = "bottom-right",
                     in_place: bool = False) -> str:
    """Inject a labeled QR code into a poster .tex.

    Adds ``\\usepackage{qrcode}`` if missing, then inserts a
    ``\\qrcode[height=Xcm]{URL}`` block with a label text directly above
    the closing ``\\end{document}``.

    Args:
        tex_path: Path to the poster ``.tex`` source.
        url: Target URL the QR encodes. HTTPS recommended.
        label: One-line label printed below the QR. Should be concrete
               ("Scan for full paper") not vague ("Scan for more").
        size_cm: QR size in cm (Animate Your Science recommends ≥2.5cm
                 for conference posters; 4cm is comfortable).
        position: Reserved for future use ("bottom-right" only for now).
        in_place: If True, modify ``tex_path`` in place. Otherwise return
                  the modified source as a string.

    Returns:
        The modified ``.tex`` source (or a status line if ``in_place``).
    """
    try:
        path, tex = read_tex(tex_path)
    except FileNotFoundError as e:
        return f"Error: file not found: {e}"

    if size_cm < 2.5:
        size_cm = 2.5  # enforce Animate Your Science minimum

    new_tex = tex
    if "qrcode" not in new_tex:
        # Insert package right after the documentclass line.
        new_tex = re.sub(
            r"(\\documentclass[^\n]*\n)",
            lambda m: m.group(1) + _QRCODE_PKG_LINE,
            new_tex,
            count=1,
        )

    qr_block = (
        "\n% --- QR code (auto-injected by poster_qr_inject) --------------\n"
        "\\begin{tikzpicture}[remember picture,overlay]\n"
        f"  \\node[anchor=south east, inner sep=20pt] at (current page.south east) {{\n"
        f"    \\begin{{tabular}}{{c}}\n"
        f"      \\qrcode[height={size_cm}cm]{{{url}}} \\\\[4pt]\n"
        f"      \\fontsize{{22}}{{26}}\\selectfont {label}\n"
        f"    \\end{{tabular}}\n"
        f"  }};\n"
        "\\end{tikzpicture}\n"
    )

    if "\\begin{tikzpicture}[remember picture,overlay]" not in new_tex:
        # Ensure tikz is available.
        if "\\usepackage{tikz}" not in new_tex:
            new_tex = re.sub(
                r"(\\documentclass[^\n]*\n)",
                lambda m: m.group(1) + "\\usepackage{tikz}\n",
                new_tex,
                count=1,
            )

    # Insert before \end{document}
    new_tex = new_tex.replace("\\end{document}", qr_block + "\n\\end{document}", 1)

    if in_place:
        path.write_text(new_tex, encoding="utf-8")
        return f"Injected QR into: {path}\n  url={url}\n  label={label}\n  size={size_cm}cm"
    return new_tex
