"""radia_mcp.chart2d — 22 paper-quality 2D charts as MCP tools.

Companion to radia_mcp.graph (which handles styling / profiles /
gates).  chart2d is the "data -> rendered 2D chart" layer: the AI
calls `chart2d_line(x=..., y=..., return_mode='image'|'recipe')` and
gets either a rendered PNG (MCP Image content type, viewable inline
in Claude Desktop) OR a copy-paste Python recipe.

22 chart types organised in 4 groups:

    LINE/TIME-SERIES (8):
        chart2d_line, chart2d_loglog, chart2d_semilogx,
        chart2d_semilogy, chart2d_step, chart2d_errorbar,
        chart2d_fill_between, chart2d_bode

    DISTRIBUTION/STATS (5):
        chart2d_histogram, chart2d_bar, chart2d_box,
        chart2d_violin, chart2d_ecdf

    SCIENTIFIC 2D FIELDS (7):
        chart2d_contour, chart2d_contourf, chart2d_pcolormesh,
        chart2d_quiver, chart2d_streamplot, chart2d_imshow,
        chart2d_polar

    POINT (2):
        chart2d_scatter, chart2d_phase

Every chart inherits radia_mcp.graph's profile system (10 pt @ 8 cm
absolute font, Okabe-Ito CVD-safe palette, Type-42 PDF embed) AND
goes through emit_paper_figure()'s 5-gate stack on save (no in-figure
title, no legend overlap, axes-area fraction floor, font embedding
verification).

Survey-driven design (2026-05-26): patterns absorbed from
StacklokLabs/plotting-mcp, antvis/mcp-server-chart (4.1k★),
LindseyyyLi/MCP-Server, isaacwasserman/mcp-vegalite-server,
arshlibruh/plotly-mcp-cursor.  The "recipe + image" dual return
generalises the lessons that no single MCP plotting server combines.
"""

from ._charts import CHARTS, ALIASES, get_spec, chart_groups, ChartSpec  # noqa: F401
from ._render import render_chart  # noqa: F401
from ._recipe import chart_recipe  # noqa: F401
