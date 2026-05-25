"""Render a chart spec to (fig, axes) using radia_mcp.graph profiles.

This is the matplotlib-backed implementation of the 22 chart types in
_charts.py.  Every renderer:
  1. calls radia_mcp.graph.paper_figure(profile, ...) so the figure
     inherits the journal profile (10 pt @ 8 cm, Okabe-Ito CVD-safe
     palette, Type-42 font embed, etc.)
  2. delegates to matplotlib's per-chart-type method (ax.plot,
     ax.contour, ...)
  3. sets axis labels with units-in-parentheses (lab convention)
  4. does NOT call set_title() (lab rule: titles go in LaTeX caption)
  5. does NOT call savefig() -- the caller decides whether to save,
     measure, gate, etc.

The renderer is intentionally minimal: it produces the figure but
applies NO gate.  Callers (server.py) thread the result through
emit_paper_figure() if they want save + Type-42 verify + legend
overlap check etc.
"""
from __future__ import annotations

from typing import Any, Optional


def render_chart(
    chart_type: str,
    *,
    profile: str = "ieee_single_column",
    rel_width: float = 1.0,
    xlabel: str = "",
    ylabel: str = "",
    legend_labels=None,
    **chart_kwargs,
):
    """Build a paper-quality (fig, axes) for the requested chart type.

    Args:
        chart_type: name or alias from radia_mcp.chart2d._charts.CHARTS.
        profile: radia_mcp.graph PaperProfile name (default IEEE 1col).
        rel_width: fraction of profile width.
        xlabel: x-axis label.  Lab convention: units in parentheses,
            e.g. 'f (Hz)' not 'f [Hz]'.
        ylabel: y-axis label, same convention.
        legend_labels: optional list of per-trace labels (line plots
            with multiple traces).  When set, ax.legend(loc='best',
            frameon=False) is called.  emit_paper_figure() will then
            reject if the legend overlaps any data line.
        **chart_kwargs: type-specific kwargs (see _charts.py per-spec
            `special_kwargs`).

    Returns:
        (fig, axes) -- axes is an ndarray (1x1 for single panel,
        2x1 for Bode pair) per radia_mcp.graph.paper_figure
        convention.
    """
    from radia_mcp.graph import paper_figure
    from ._charts import get_spec

    spec = get_spec(chart_type)

    # The Bode pair is a 2x1 panel; everything else is 1x1.
    if spec.name == "bode":
        return _render_bode(profile=profile, rel_width=rel_width,
                            xlabel=xlabel, ylabel=ylabel,
                            legend_labels=legend_labels, **chart_kwargs)

    # Polar needs projection='polar' subplot kwarg, which paper_figure
    # does not currently expose.  Build manually for that case so the
    # rest of the line / scatter / field renderers stay uniform.
    if spec.name == "polar":
        return _render_polar(profile=profile, rel_width=rel_width,
                             **chart_kwargs)

    fig, axes = paper_figure(profile, nrows=1, ncols=1,
                              rel_width=rel_width)
    ax = axes[0, 0]

    # Dispatch to per-type matplotlib call
    _dispatch_chart(ax, spec, chart_kwargs, legend_labels)

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    if legend_labels:
        # Lab rule: frameless legend; emit_paper_figure() detects
        # overlap with data.  AI should call find_best_legend_loc
        # before emit if it cares.
        ax.legend(loc="best", frameon=False)

    return fig, axes


def _dispatch_chart(ax, spec, kw, legend_labels):
    """Call the appropriate matplotlib method on ax."""
    name = spec.name
    if name == "line":
        x = kw.get("x")
        y = kw["y"]
        # Multi-trace support: y can be 2D (n_trace, n_point)
        import numpy as np
        y_arr = np.asarray(y)
        if y_arr.ndim == 2:
            for i, y_i in enumerate(y_arr):
                lbl = (legend_labels[i] if legend_labels and
                       i < len(legend_labels) else None)
                ax.plot(x if x is not None else range(len(y_i)),
                        y_i, label=lbl)
        else:
            ax.plot(x if x is not None else range(len(y_arr)),
                    y_arr,
                    label=legend_labels[0] if legend_labels else None)
        return

    if name in ("loglog", "semilogx", "semilogy"):
        x = kw.get("x")
        y = kw["y"]
        getattr(ax, name)(x, y,
                          label=legend_labels[0] if legend_labels else None)
        return

    if name == "step":
        x = kw.get("x")
        y = kw["y"]
        where = kw.get("where", "post")
        ax.step(x, y, where=where,
                label=legend_labels[0] if legend_labels else None)
        return

    if name == "errorbar":
        ax.errorbar(
            kw.get("x"), kw["y"],
            yerr=kw.get("yerr"), xerr=kw.get("xerr"),
            capsize=kw.get("capsize", 3.0),
            fmt=kw.get("fmt", "o"),
            label=legend_labels[0] if legend_labels else None,
        )
        return

    if name == "fill_between":
        x = kw["x"]
        y1 = kw["y"]
        y2 = kw.get("y2", 0)
        ax.fill_between(x, y1, y2, alpha=kw.get("alpha", 0.3),
                        label=legend_labels[0] if legend_labels else None)
        return

    if name == "histogram":
        data = kw["data"]
        ax.hist(data, bins=kw.get("bins", 30),
                density=kw.get("density", False))
        return

    if name == "bar":
        cats = kw["categories"]
        h = kw["heights"]
        w = kw.get("width", 0.8)
        ax.bar(cats, h, width=w)
        return

    if name == "box":
        data = kw["data"]
        ax.boxplot(data, tick_labels=kw.get("labels"),
                   showfliers=kw.get("showfliers", True))
        return

    if name == "violin":
        data = kw["data"]
        parts = ax.violinplot(data,
                               showmedians=kw.get("showmedians", True))
        labels = kw.get("labels")
        if labels:
            ax.set_xticks(range(1, len(labels) + 1))
            ax.set_xticklabels(labels)
        return

    if name == "ecdf":
        data = kw["data"]
        if hasattr(ax, "ecdf"):
            ax.ecdf(data)
        else:
            # matplotlib < 3.8 fallback
            import numpy as np
            sorted_d = np.sort(np.asarray(data))
            n = len(sorted_d)
            ax.step(sorted_d, (np.arange(n) + 1) / n, where="post")
        return

    if name == "contour":
        X, Y, Z = kw["X"], kw["Y"], kw["Z"]
        cs = ax.contour(X, Y, Z, levels=kw.get("levels"))
        ax.clabel(cs, inline=True, fontsize=8)
        return

    if name == "contourf":
        X, Y, Z = kw["X"], kw["Y"], kw["Z"]
        cs = ax.contourf(X, Y, Z, levels=kw.get("levels"),
                          cmap=kw.get("cmap", "viridis"))
        ax.figure.colorbar(cs, ax=ax)
        return

    if name == "pcolormesh":
        X, Y, Z = kw["X"], kw["Y"], kw["Z"]
        m = ax.pcolormesh(X, Y, Z, cmap=kw.get("cmap", "viridis"),
                           shading=kw.get("shading", "auto"))
        ax.figure.colorbar(m, ax=ax)
        return

    if name == "quiver":
        ax.quiver(kw["X"], kw["Y"], kw["U"], kw["V"],
                  scale=kw.get("scale"),
                  scale_units=kw.get("scale_units", "xy"))
        return

    if name == "streamplot":
        sp = ax.streamplot(kw["X"], kw["Y"], kw["U"], kw["V"],
                            density=kw.get("density", 1.0),
                            color=kw.get("color"))
        return

    if name == "imshow":
        im = ax.imshow(kw["data"],
                       cmap=kw.get("cmap", "viridis"),
                       origin=kw.get("origin", "lower"),
                       extent=kw.get("extent"),
                       aspect=kw.get("aspect", "auto"))
        ax.figure.colorbar(im, ax=ax)
        return

    if name == "scatter":
        sc = ax.scatter(
            kw["x"], kw["y"],
            s=kw.get("s"), c=kw.get("c"),
            alpha=kw.get("alpha", 0.7),
            cmap=kw.get("cmap", "viridis"),
            label=legend_labels[0] if legend_labels else None,
        )
        return

    if name == "phase":
        # Complex-plane plot: Re vs Im
        import numpy as np
        z = kw["z"]
        if isinstance(z, tuple) and len(z) == 2:
            re, im = np.asarray(z[0]), np.asarray(z[1])
        else:
            z_arr = np.asarray(z, dtype=complex)
            re, im = z_arr.real, z_arr.imag
        ax.plot(re, im,
                label=legend_labels[0] if legend_labels else None)
        ax.set_aspect("equal")
        ax.axhline(0, color="0.7", lw=0.5, zorder=0)
        ax.axvline(0, color="0.7", lw=0.5, zorder=0)
        return

    raise NotImplementedError(
        f"chart type {name!r} dispatch not implemented"
    )


# ============================================================
# Special-case renderers (multi-panel / non-rectangular axes)
# ============================================================


def _render_bode(*, profile, rel_width, xlabel, ylabel,
                  legend_labels, freq, magnitude, phase,
                  mag_db=True, phase_unit="deg"):
    """Bode pair: magnitude (top) + phase (bottom), shared frequency axis.

    Top:    semilogx of |H| (in dB if mag_db else linear).
    Bottom: semilogx of arg(H) (deg or rad).
    """
    from radia_mcp.graph import paper_figure
    import numpy as np

    fig, axes = paper_figure(profile, nrows=2, ncols=1,
                              rel_width=rel_width, sharex=True,
                              panel_labels=False)
    ax_mag = axes[0, 0]
    ax_phs = axes[1, 0]

    f = np.asarray(freq)
    M = np.asarray(magnitude)
    P = np.asarray(phase)

    ax_mag.semilogx(f, M,
                    label=legend_labels[0] if legend_labels else None)
    ax_phs.semilogx(f, P,
                    label=legend_labels[0] if legend_labels else None)

    # Mag y label: |H| (dB) or |H|
    if not ylabel:
        ax_mag.set_ylabel(r"$|H|$ (dB)" if mag_db else r"$|H|$")
    else:
        ax_mag.set_ylabel(ylabel)

    # Phase y label: arg(H) (deg) or arg(H) (rad)
    pu = "deg" if phase_unit.lower().startswith("deg") else "rad"
    ax_phs.set_ylabel(
        r"$\arg(H)$ (" + ("$\\degree$" if pu == "deg" else "rad") + ")"
    )

    # X label only on the bottom; lab convention to suppress top.
    ax_phs.set_xlabel(xlabel if xlabel else r"$f$ (Hz)")
    # Hide top panel's x tick labels (sharex keeps the tick marks).
    for tl in ax_mag.get_xticklabels():
        tl.set_visible(False)

    if legend_labels:
        ax_mag.legend(loc="best", frameon=False)

    return fig, axes


def _render_polar(*, profile, rel_width, theta, r):
    """Polar plot: not a paper_figure() one-shot because polar axes need
    projection='polar' subplot kwarg.  Build manually."""
    import matplotlib.pyplot as plt
    from radia_mcp.graph import paper_figure
    import numpy as np

    # Set up rcParams via a throwaway paper_figure call, then close it.
    fig0, _ = paper_figure(profile, rel_width=rel_width)
    figsize = fig0.get_size_inches()
    plt.close(fig0)

    fig = plt.figure(figsize=figsize, dpi=600)
    ax = fig.add_subplot(111, projection="polar")
    ax.plot(theta, r)
    # Return axes as 2D ndarray for caller consistency
    import numpy as _np
    axes = _np.array([[ax]])
    return fig, axes
