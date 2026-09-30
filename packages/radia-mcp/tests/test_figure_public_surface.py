"""The documented TikZ guide is routed through the existing paper server."""
import tomllib
from pathlib import Path
from radia_mcp.figure import register
from radia_mcp.figure.tools import _PROFILES, figure_tikz_recipe, figure_matlab2tikz_recipe


def test_tikz_is_registered_with_figure_tools():
    names = []
    class Registry:
        def tool(self):
            return lambda fn: names.append(fn.__name__) or fn
    register(Registry())
    assert names.count("figure_tikz_recipe") == 1


def test_figure_uses_paper_server_not_a_retired_entry_point():
    data = tomllib.loads((Path(__file__).parents[1] / "pyproject.toml").read_text(encoding="utf-8"))
    scripts = data["project"]["scripts"]
    assert "mcp-server-paper-writing" in scripts
    assert "mcp-server-figure" not in scripts


def test_tikz_recipe_has_primary_sources_and_drawing_guidance():
    recipe = figure_tikz_recipe()
    for text in ("PGFPlots", "No title inside the picture", "https://tikz.dev/", "https://ctan.org/pkg/pgfplots"):
        assert text in recipe


def test_tikz_recipe_uses_current_target_width():
    recipe = figure_tikz_recipe("pgfplots", "paper_double_column")
    assert f"width={_PROFILES['paper_double_column']['embed_width_cm']:.2f}cm" in recipe


def test_matlab_section_reuses_the_existing_recipe():
    assert figure_matlab2tikz_recipe(target="paper_single_column").strip() == figure_tikz_recipe("matlab2tikz")

