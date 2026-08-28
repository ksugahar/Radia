"""The supported equation API is TeX-only; legacy binary formats stay gone."""
from __future__ import annotations

from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]


def test_retired_codec_sources_are_not_present() -> None:
    retired = [
        "mtef_parser.cpp",
        "mtef2tex.cpp",
        "tex2mtef.cpp",
        "eqnedt64_main.cpp",
        "eq_window.cpp",
    ]
    source = ROOT / "src" / "ext" / "equation"
    assert not [name for name in retired if (source / name).exists()]


def test_public_module_has_no_mtef_or_eqn_entry_points() -> None:
    equation = pytest.importorskip("radia.equation")
    retired = {
        "tex_to_mtef",
        "mtef_to_tex",
        "mtef_to_latex",
        "mtef_to_mathml",
        "mtef_to_omml",
        "mtef_to_rtf",
        "mtef_to_svg",
        "read_eqn",
        "write_eqn",
    }
    assert retired.isdisjoint(dir(equation))
    assert retired.isdisjoint(dir(__import__("radia._equation", fromlist=["*"])))


def test_tex_to_native_office_remains_supported() -> None:
    equation = pytest.importorskip("radia.equation")
    assert "<m:f>" in equation.tex_to_omml(r"\frac{a}{b}")
    assert "<mfrac>" in equation.tex_to_mathml(r"\frac{a}{b}")
