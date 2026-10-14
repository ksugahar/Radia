"""Color parsing + WCAG / CVD utilities for poster accessibility lints.

Implements the WebAIM-formulaic relative-luminance + contrast-ratio
calculation (WCAG 2.1) and a deuteranopia-style channel mix so we can
flag hazardous color pairs without external dependencies.
"""
from __future__ import annotations

import re


_DEFINECOLOR_RE = re.compile(
    r"\\definecolor\{(?P<name>[^}]+)\}\{(?P<model>RGB|rgb|gray|HTML|cmyk)\}\{(?P<spec>[^}]+)\}"
)


def parse_definecolor(tex: str) -> dict[str, tuple[float, float, float]]:
    """Extract ``\\definecolor`` declarations as ``{name: (r,g,b)}`` in 0-1.

    Supported models: RGB (0-255), rgb (0-1), HTML (hex), gray (0-1), cmyk.
    Unknown models are silently ignored.
    """
    out: dict[str, tuple[float, float, float]] = {}
    for m in _DEFINECOLOR_RE.finditer(tex):
        name = m.group("name")
        model = m.group("model")
        spec = m.group("spec").strip()
        try:
            if model == "RGB":
                r, g, b = (float(x) for x in spec.split(","))
                out[name] = (r / 255, g / 255, b / 255)
            elif model == "rgb":
                r, g, b = (float(x) for x in spec.split(","))
                out[name] = (r, g, b)
            elif model == "HTML":
                v = int(spec, 16)
                out[name] = ((v >> 16) / 255, ((v >> 8) & 0xff) / 255, (v & 0xff) / 255)
            elif model == "gray":
                v = float(spec)
                out[name] = (v, v, v)
            elif model == "cmyk":
                c, mg, y, k = (float(x) for x in spec.split(","))
                r = (1 - c) * (1 - k)
                g = (1 - mg) * (1 - k)
                b = (1 - y) * (1 - k)
                out[name] = (r, g, b)
        except ValueError:
            continue
    return out


def _channel_lin(c: float) -> float:
    return c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4


def relative_luminance(rgb: tuple[float, float, float]) -> float:
    """WCAG 2.1 relative luminance (0-1)."""
    r, g, b = rgb
    return 0.2126 * _channel_lin(r) + 0.7152 * _channel_lin(g) + 0.0722 * _channel_lin(b)


def contrast_ratio(rgb1, rgb2) -> float:
    """WCAG 2.1 contrast ratio, range [1, 21]."""
    l1 = relative_luminance(rgb1)
    l2 = relative_luminance(rgb2)
    hi, lo = max(l1, l2), min(l1, l2)
    return (hi + 0.05) / (lo + 0.05)


def simulate_deuteranopia(rgb: tuple[float, float, float]) -> tuple[float, float, float]:
    """Approximate deuteranopia (green-blind) channel mix.

    Uses the Brettel-Viénot matrix for the 0,0,1 confusion axis as
    typically published in colorblind-simulator code. Adequate for
    flagging hazardous pairs; not a substitute for clinical simulation.
    """
    r, g, b = rgb
    rr = 0.625 * r + 0.375 * g
    gg = 0.700 * r + 0.300 * g
    bb = 0.000 * r + 0.300 * g + 0.700 * b
    return (max(0, min(1, rr)), max(0, min(1, gg)), max(0, min(1, bb)))


def color_distance(rgb1, rgb2) -> float:
    """Simple Euclidean distance in sRGB; good enough for flagging."""
    r1, g1, b1 = rgb1
    r2, g2, b2 = rgb2
    return ((r1 - r2) ** 2 + (g1 - g2) ** 2 + (b1 - b2) ** 2) ** 0.5
