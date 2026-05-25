"""MCP server: radia_mcp.chart2d.

Exposes 22 chart-type MCP tools (`chart2d_line`, `chart2d_loglog`, ...)
plus a `chart2d_catalog` introspection tool.  Every tool accepts a
`return_mode` parameter:

    return_mode='recipe'  -> returns Python code string (default).
                              The user pastes the snippet into a
                              script and runs it locally.  Fast, no
                              matplotlib dependency on the MCP server.
    return_mode='image'   -> server-side rendering with matplotlib +
                              radia_mcp.graph profile, save to a
                              temp PDF/PNG, return the PNG bytes via
                              MCP Image content type.  Inline preview
                              in Claude Desktop / claude.ai.
    return_mode='both'    -> returns a dict {'recipe': str, 'image':
                              Image, 'path': str} -- preview + script.

Usage:
    mcp-server-chart2d              # start (stdio transport)
    mcp-server-chart2d --selftest   # smoke check all 22 types
"""

from __future__ import annotations

import sys
import tempfile
import os
from pathlib import Path
from typing import Any, Literal, Optional

from mcp.server.fastmcp import FastMCP, Image

from ..common import register_status_tool
from . import _charts as _ch
from . import _render as _rd
from . import _recipe as _rc

mcp = FastMCP("mcp-server-chart2d")


# ============================================================
# Catalog tool
# ============================================================


@mcp.tool()
def chart2d_catalog(group: str = "all") -> dict:
    """List the 22 2D chart types and their groupings.

    Args:
        group: 'all' (default) | 'line' | 'dist' | 'field' | 'point'

    Returns dict with:
        n_types: total
        groups: {group: [chart_name, ...]}
        types: [{name, group, description, aliases, ...}]
    """
    g = _ch.chart_groups()
    if group != "all":
        if group not in g:
            return {"error": f"Unknown group {group!r}. "
                              f"Available: all, {', '.join(g)}"}
        filtered = {group: g[group]}
    else:
        filtered = g
    types = []
    for name, spec in _ch.CHARTS.items():
        if group != "all" and spec.group != group:
            continue
        types.append({
            "name": spec.name,
            "group": spec.group,
            "description": spec.description,
            "aliases": list(spec.aliases),
            "mpl_call": spec.mpl_call,
            "x_units_hint": spec.x_units_hint,
            "y_units_hint": spec.y_units_hint,
            "notes": spec.notes,
        })
    return {
        "n_types": len(types),
        "groups": filtered,
        "types": types,
        "available_profiles": [
            "ieee_single_column", "ieee_double_column",
            "ieej_single_column", "ieej_double_column",
            "igte_digest_double", "igte_digest_single",
        ],
        "return_modes": ["recipe", "image", "both"],
        "hint":
            "Call chart2d_<name>(x=..., y=..., "
            "return_mode='recipe'|'image'|'both').  See radia_mcp.graph"
            " quality_rules for style policy.",
    }


# ============================================================
# Dispatch helper used by every chart tool
# ============================================================


def _dispatch(chart_type: str, return_mode: str, render_kwargs: dict,
              xlabel: str, ylabel: str,
              profile: str, rel_width: float,
              legend_labels=None) -> Any:
    """Common entry point shared by all 22 chart tools."""
    if return_mode not in ("recipe", "image", "both"):
        return {"error": f"return_mode must be 'recipe', 'image', or "
                          f"'both'; got {return_mode!r}"}

    out: dict[str, Any] = {}

    if return_mode in ("recipe", "both"):
        out["recipe"] = _rc.chart_recipe(
            chart_type, profile=profile, rel_width=rel_width,
            xlabel=xlabel or r"$f$ (Hz)",
            ylabel=ylabel or r"$|Z|$ ($\Omega$)",
        )

    if return_mode in ("image", "both"):
        # Render to a temp PDF/PNG via radia_mcp.graph profile.
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        try:
            fig, axes = _rd.render_chart(
                chart_type,
                profile=profile, rel_width=rel_width,
                xlabel=xlabel, ylabel=ylabel,
                legend_labels=legend_labels,
                **render_kwargs,
            )
            # Save PNG to temp file, return bytes via MCP Image type.
            with tempfile.NamedTemporaryFile(
                suffix=".png", delete=False) as tf:
                png_path = tf.name
            fig.savefig(png_path, dpi=300, bbox_inches="tight")
            plt.close(fig)
            with open(png_path, "rb") as f:
                png_bytes = f.read()
            try:
                os.unlink(png_path)
            except OSError:
                pass
            out["image"] = Image(data=png_bytes, format="png")
        except Exception as e:
            out["error"] = (
                f"chart2d_{chart_type} render failed: "
                f"{type(e).__name__}: {e}"
            )

    if return_mode == "recipe":
        return {"recipe": out["recipe"]}
    if return_mode == "image":
        return out.get("image") or {"error": out.get("error",
                                                       "no image")}
    return out


# ============================================================
# 22 chart tools — minimal docstrings; full guidance via chart2d_catalog
# ============================================================
#
# Implementation pattern: each tool accepts the type-specific kwargs +
# the universal ones (xlabel/ylabel/profile/rel_width/return_mode),
# packages render_kwargs, and forwards to _dispatch.
# Default profile is ieee_single_column (88.9 mm @ 10 pt) -- the most
# common lab use case.


# --- LINE / TIME-SERIES (8) ---

@mcp.tool()
def chart2d_line(x: list = None, y: list = None,
                  labels: list = None,
                  xlabel: str = "", ylabel: str = "",
                  profile: str = "ieee_single_column",
                  rel_width: float = 1.0,
                  return_mode: str = "recipe") -> Any:
    """Standard line plot. y can be 1D or 2D (multiple traces)."""
    return _dispatch("line", return_mode,
                     {"x": x, "y": y}, xlabel, ylabel, profile,
                     rel_width, legend_labels=labels)


@mcp.tool()
def chart2d_loglog(x: list = None, y: list = None,
                    label: str = "",
                    xlabel: str = "", ylabel: str = "",
                    profile: str = "ieee_single_column",
                    rel_width: float = 1.0,
                    return_mode: str = "recipe") -> Any:
    """Log-log plot.  Straight line of slope alpha = power-law y=k*x^alpha."""
    return _dispatch("loglog", return_mode,
                     {"x": x, "y": y}, xlabel, ylabel, profile,
                     rel_width, legend_labels=[label] if label else None)


@mcp.tool()
def chart2d_semilogx(x: list = None, y: list = None,
                      label: str = "",
                      xlabel: str = "", ylabel: str = "",
                      profile: str = "ieee_single_column",
                      rel_width: float = 1.0,
                      return_mode: str = "recipe") -> Any:
    """Linear y, log10 x.  Bode-magnitude / frequency-domain default."""
    return _dispatch("semilogx", return_mode,
                     {"x": x, "y": y}, xlabel, ylabel, profile,
                     rel_width, legend_labels=[label] if label else None)


@mcp.tool()
def chart2d_semilogy(x: list = None, y: list = None,
                      label: str = "",
                      xlabel: str = "", ylabel: str = "",
                      profile: str = "ieee_single_column",
                      rel_width: float = 1.0,
                      return_mode: str = "recipe") -> Any:
    """Log10 y.  Decaying / growing quantities (relaxation, noise tails)."""
    return _dispatch("semilogy", return_mode,
                     {"x": x, "y": y}, xlabel, ylabel, profile,
                     rel_width, legend_labels=[label] if label else None)


@mcp.tool()
def chart2d_step(x: list = None, y: list = None,
                  where: str = "post",
                  label: str = "",
                  xlabel: str = "", ylabel: str = "",
                  profile: str = "ieee_single_column",
                  rel_width: float = 1.0,
                  return_mode: str = "recipe") -> Any:
    """Stair-step plot.  PWM / discrete-time / step histogram."""
    return _dispatch("step", return_mode,
                     {"x": x, "y": y, "where": where}, xlabel, ylabel,
                     profile, rel_width,
                     legend_labels=[label] if label else None)


@mcp.tool()
def chart2d_errorbar(x: list = None, y: list = None,
                      yerr: list = None, xerr: list = None,
                      fmt: str = "o", capsize: float = 3.0,
                      label: str = "",
                      xlabel: str = "", ylabel: str = "",
                      profile: str = "ieee_single_column",
                      rel_width: float = 1.0,
                      return_mode: str = "recipe") -> Any:
    """Line + error bars for measured data with uncertainty."""
    return _dispatch("errorbar", return_mode,
                     {"x": x, "y": y, "yerr": yerr, "xerr": xerr,
                      "fmt": fmt, "capsize": capsize},
                     xlabel, ylabel, profile, rel_width,
                     legend_labels=[label] if label else None)


@mcp.tool()
def chart2d_fill_between(x: list = None, y: list = None,
                          y2: list = None, alpha: float = 0.3,
                          label: str = "",
                          xlabel: str = "", ylabel: str = "",
                          profile: str = "ieee_single_column",
                          rel_width: float = 1.0,
                          return_mode: str = "recipe") -> Any:
    """Filled band between y and y2 (defaults to 0).  Confidence / hysteresis."""
    return _dispatch("fill_between", return_mode,
                     {"x": x, "y": y, "y2": y2, "alpha": alpha},
                     xlabel, ylabel, profile, rel_width,
                     legend_labels=[label] if label else None)


@mcp.tool()
def chart2d_bode(freq: list = None, magnitude: list = None,
                  phase: list = None,
                  mag_db: bool = True, phase_unit: str = "deg",
                  label: str = "",
                  xlabel: str = "", ylabel: str = "",
                  profile: str = "ieee_single_column",
                  rel_width: float = 1.0,
                  return_mode: str = "recipe") -> Any:
    """Bode pair: magnitude (top) + phase (bottom), shared frequency axis."""
    return _dispatch("bode", return_mode,
                     {"freq": freq, "magnitude": magnitude,
                      "phase": phase, "mag_db": mag_db,
                      "phase_unit": phase_unit},
                     xlabel, ylabel, profile, rel_width,
                     legend_labels=[label] if label else None)


# --- DISTRIBUTION / STATISTICS (5) ---

@mcp.tool()
def chart2d_histogram(data: list = None, bins: int = 30,
                       density: bool = False,
                       xlabel: str = "", ylabel: str = "",
                       profile: str = "ieee_single_column",
                       rel_width: float = 1.0,
                       return_mode: str = "recipe") -> Any:
    """1D histogram."""
    return _dispatch("histogram", return_mode,
                     {"data": data, "bins": bins, "density": density},
                     xlabel, ylabel, profile, rel_width)


@mcp.tool()
def chart2d_bar(categories: list = None, heights: list = None,
                 width: float = 0.8,
                 xlabel: str = "", ylabel: str = "",
                 profile: str = "ieee_single_column",
                 rel_width: float = 1.0,
                 return_mode: str = "recipe") -> Any:
    """Vertical bar chart."""
    return _dispatch("bar", return_mode,
                     {"categories": categories, "heights": heights,
                      "width": width},
                     xlabel, ylabel, profile, rel_width)


@mcp.tool()
def chart2d_box(data: list = None, labels: list = None,
                 showfliers: bool = True,
                 xlabel: str = "", ylabel: str = "",
                 profile: str = "ieee_single_column",
                 rel_width: float = 1.0,
                 return_mode: str = "recipe") -> Any:
    """Box-and-whisker plot."""
    return _dispatch("box", return_mode,
                     {"data": data, "labels": labels,
                      "showfliers": showfliers},
                     xlabel, ylabel, profile, rel_width)


@mcp.tool()
def chart2d_violin(data: list = None, labels: list = None,
                    showmedians: bool = True,
                    xlabel: str = "", ylabel: str = "",
                    profile: str = "ieee_single_column",
                    rel_width: float = 1.0,
                    return_mode: str = "recipe") -> Any:
    """Violin plot (KDE-as-fill)."""
    return _dispatch("violin", return_mode,
                     {"data": data, "labels": labels,
                      "showmedians": showmedians},
                     xlabel, ylabel, profile, rel_width)


@mcp.tool()
def chart2d_ecdf(data: list = None,
                  xlabel: str = "", ylabel: str = "",
                  profile: str = "ieee_single_column",
                  rel_width: float = 1.0,
                  return_mode: str = "recipe") -> Any:
    """Empirical cumulative distribution."""
    return _dispatch("ecdf", return_mode,
                     {"data": data},
                     xlabel, ylabel, profile, rel_width)


# --- SCIENTIFIC 2D FIELDS (7) ---

@mcp.tool()
def chart2d_contour(X=None, Y=None, Z=None, levels=None,
                     xlabel: str = "", ylabel: str = "",
                     profile: str = "ieee_single_column",
                     rel_width: float = 1.0,
                     return_mode: str = "recipe") -> Any:
    """Contour lines of Z(X, Y)."""
    return _dispatch("contour", return_mode,
                     {"X": X, "Y": Y, "Z": Z, "levels": levels},
                     xlabel, ylabel, profile, rel_width)


@mcp.tool()
def chart2d_contourf(X=None, Y=None, Z=None, levels=None,
                      cmap: str = "viridis",
                      xlabel: str = "", ylabel: str = "",
                      profile: str = "ieee_single_column",
                      rel_width: float = 1.0,
                      return_mode: str = "recipe") -> Any:
    """Filled contours.  Default viridis = perceptually uniform + CVD safe."""
    return _dispatch("contourf", return_mode,
                     {"X": X, "Y": Y, "Z": Z, "levels": levels,
                      "cmap": cmap},
                     xlabel, ylabel, profile, rel_width)


@mcp.tool()
def chart2d_pcolormesh(X=None, Y=None, Z=None,
                        cmap: str = "viridis",
                        shading: str = "auto",
                        xlabel: str = "", ylabel: str = "",
                        profile: str = "ieee_single_column",
                        rel_width: float = 1.0,
                        return_mode: str = "recipe") -> Any:
    """Irregular-grid heatmap."""
    return _dispatch("pcolormesh", return_mode,
                     {"X": X, "Y": Y, "Z": Z, "cmap": cmap,
                      "shading": shading},
                     xlabel, ylabel, profile, rel_width)


@mcp.tool()
def chart2d_quiver(X=None, Y=None, U=None, V=None,
                    scale: float = None,
                    scale_units: str = "xy",
                    xlabel: str = "", ylabel: str = "",
                    profile: str = "ieee_single_column",
                    rel_width: float = 1.0,
                    return_mode: str = "recipe") -> Any:
    """2D vector field as arrows."""
    return _dispatch("quiver", return_mode,
                     {"X": X, "Y": Y, "U": U, "V": V,
                      "scale": scale, "scale_units": scale_units},
                     xlabel, ylabel, profile, rel_width)


@mcp.tool()
def chart2d_streamplot(X=None, Y=None, U=None, V=None,
                        density: float = 1.0, color=None,
                        xlabel: str = "", ylabel: str = "",
                        profile: str = "ieee_single_column",
                        rel_width: float = 1.0,
                        return_mode: str = "recipe") -> Any:
    """Field lines of a 2D vector field on a REGULAR grid."""
    return _dispatch("streamplot", return_mode,
                     {"X": X, "Y": Y, "U": U, "V": V,
                      "density": density, "color": color},
                     xlabel, ylabel, profile, rel_width)


@mcp.tool()
def chart2d_imshow(data=None, cmap: str = "viridis",
                    origin: str = "lower", extent=None,
                    aspect: str = "auto",
                    xlabel: str = "", ylabel: str = "",
                    profile: str = "ieee_single_column",
                    rel_width: float = 1.0,
                    return_mode: str = "recipe") -> Any:
    """2D array as image / heatmap."""
    return _dispatch("imshow", return_mode,
                     {"data": data, "cmap": cmap, "origin": origin,
                      "extent": extent, "aspect": aspect},
                     xlabel, ylabel, profile, rel_width)


@mcp.tool()
def chart2d_polar(theta: list = None, r: list = None,
                   profile: str = "ieee_single_column",
                   rel_width: float = 1.0,
                   return_mode: str = "recipe") -> Any:
    """Polar / radial plot."""
    return _dispatch("polar", return_mode,
                     {"theta": theta, "r": r},
                     "", "", profile, rel_width)


# --- POINT (2) ---

@mcp.tool()
def chart2d_scatter(x: list = None, y: list = None,
                     s=None, c=None,
                     alpha: float = 0.7, cmap: str = "viridis",
                     label: str = "",
                     xlabel: str = "", ylabel: str = "",
                     profile: str = "ieee_single_column",
                     rel_width: float = 1.0,
                     return_mode: str = "recipe") -> Any:
    """Scatter plot of (x, y) points."""
    return _dispatch("scatter", return_mode,
                     {"x": x, "y": y, "s": s, "c": c,
                      "alpha": alpha, "cmap": cmap},
                     xlabel, ylabel, profile, rel_width,
                     legend_labels=[label] if label else None)


@mcp.tool()
def chart2d_phase(z=None, label: str = "",
                   xlabel: str = "Re", ylabel: str = "Im",
                   profile: str = "ieee_single_column",
                   rel_width: float = 1.0,
                   return_mode: str = "recipe") -> Any:
    """Complex-plane (Re vs Im) -- Nyquist / impedance locus / root locus.
    Aliases: nyquist, complex_plane, argand."""
    return _dispatch("phase", return_mode,
                     {"z": z},
                     xlabel, ylabel, profile, rel_width,
                     legend_labels=[label] if label else None)


# ============================================================
# Self-introspection
# ============================================================

register_status_tool(
    mcp,
    server_name="mcp-server-chart2d",
    description=(
        "22 paper-quality 2D charts as MCP tools.  Recipe + MCP-Image "
        "dual return.  Inherits radia_mcp.graph profile + gate stack "
        "(10 pt @ 8 cm, Okabe-Ito CVD palette, Type-42 PDF embed)."
    ),
    subpackage="radia_mcp.chart2d",
    related_servers=["graph", "mathematica"],
    optional_deps=["matplotlib"],
)


# ============================================================
# Entry point
# ============================================================


def main():
    if "--selftest" in sys.argv:
        print("chart2d MCP server self-test:")
        cat = chart2d_catalog("all")
        print(f"  Catalog: {cat['n_types']} chart types in "
              f"{len(cat['groups'])} groups")
        for grp_name, grp_charts in cat["groups"].items():
            print(f"    [{grp_name}] {', '.join(grp_charts)}")

        # Smoke test: every chart type's recipe is non-empty Python
        # (server-side render check is skipped here -- it would import
        # matplotlib AND require numpy arrays.  The render path is
        # exercised by tests/test_chart2d.py).
        print(f"  Recipe smoke test ({cat['n_types']} types):")
        for spec_name in _ch.CHARTS:
            try:
                recipe = _rc.chart_recipe(spec_name)
                assert "paper_figure(" in recipe
                assert "emit_paper_figure(" in recipe
                print(f"    {spec_name:18s} -> {len(recipe):5d} chars")
            except Exception as e:
                print(f"    {spec_name:18s} -> FAIL: "
                      f"{type(e).__name__}: {e}")
                sys.exit(1)

        # Alias resolution
        print(f"  Aliases:")
        for alias, target in _ch.ALIASES.items():
            print(f"    {alias:18s} -> {target}")
        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
