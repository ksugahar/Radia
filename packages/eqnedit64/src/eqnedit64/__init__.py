"""Python and native-command surfaces for Eqnedit64."""
from __future__ import annotations

from .eqnedit_core import (
    Equation,
    MAX_NESTING_DEPTH,
    SvgStyle,
    compose_document,
    math_font_loaded,
    normalize_paste,
    palette_categories,
    palettes,
    symbol_commands,
    symbol_palette_count,
    tex_normalize,
    tex_to_mathml,
    tex_to_svg,
)
from .api import backend_path, copy_equation, export_powerpoint_math, render_equation, web_asset

__all__ = [
    "Equation",
    "MAX_NESTING_DEPTH",
    "SvgStyle",
    "backend_path",
    "compose_document",
    "copy_equation",
    "export_powerpoint_math",
    "math_font_loaded",
    "normalize_paste",
    "palette_categories",
    "palettes",
    "render_equation",
    "symbol_commands",
    "symbol_palette_count",
    "tex_normalize",
    "tex_to_mathml",
    "tex_to_svg",
    "web_asset",
]

__version__ = "3.1.4"
