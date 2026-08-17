"""radia_mcp.figure._builders -- one-call publication-figure builders.

Each builder returns ``(fig, axes)`` authored via :func:`lab_figure`
(Times New Roman, 10 pt at the embed width, Okabe-Ito CVD-safe colors).
Pass the result to :func:`save_lab_figure`, which runs the fail-loud gates
and emits the ``\\includegraphics[width=<cm>]`` snippet.

These distil the figures hand-rolled for the CEFC 2026 oral deck so the
next deck does not re-derive them:

    scaling_loglog  -- O(N log N) memory/time vs DOF, dense curves end at OOM
    grouped_bars    -- solver-comparison bars (None = a failed case omitted)
    convergence     -- iteration curves (loop content, residual, ...)
    quiver_pair     -- two vector-field panels (aligned vs loop magnetization)
    bh_curve        -- monotonic B-H magnetization curve
"""
from __future__ import annotations

# numpy is an optional runtime dep (not installed in CI minimal env).
# Import lazily inside each function so the module-level import succeeds
# without numpy (radia_mcp servers must be importable without heavy deps).
try:
    import numpy as np
    _HAS_NUMPY = True
except ImportError:  # pragma: no cover
    np = None  # type: ignore[assignment]
    _HAS_NUMPY = False

from ._paper_figure import _get_mpl, OKABE_ITO
from ._lab_api import lab_figure, legend_no_overlap


def _cycle(keys, colors):
    return colors or {k: OKABE_ITO[i % len(OKABE_ITO)] for i, k in enumerate(keys)}


def scaling_loglog(x, series, ylabel, *, xlabel="DOF", embed_width_cm=7.6,
                   aspect=0.74, ram_line_gb=None, colors=None, markers=None):
    """Log-log scaling: ``x`` vs each ``series[label]``.  ``None``/NaN in a
    series = the run stopped (e.g. out of memory) -> the point is dropped so
    the curve visibly ENDS.  Optional horizontal RAM-ceiling line."""
    _get_mpl()
    fig, ax = lab_figure(embed_width_cm, aspect=aspect)
    x = np.asarray(x, float)
    keys = list(series.keys())
    cyc = _cycle(keys, colors)
    mk = markers or {k: "os^vD"[i % 5] for i, k in enumerate(keys)}
    for k in keys:
        y = np.array([np.nan if v is None else v for v in series[k]], float)
        g = ~np.isnan(y)
        ax.loglog(x[g], y[g], mk[k] + "-", color=cyc[k], label=k, ms=4.5, lw=1.3)
    if ram_line_gb:
        ax.axhline(ram_line_gb, color="0.45", ls="--", lw=0.9)
        ax.text(x[0] * 1.05, ram_line_gb * 1.15, f"{ram_line_gb:g} GB RAM",
                fontsize=10, color="0.35")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.55)
    legend_no_overlap(ax, fontsize=10)
    return fig, ax


def grouped_bars(group_labels, series, ylabel, *, embed_width_cm=8.4,
                 aspect=0.62, log=True, colors=None, width=0.26):
    """Grouped bars.  ``None`` in a series omits that bar (a failed case).
    Returns ``(fig, ax)`` so the caller can add FAIL text / speedup arrows."""
    _get_mpl()
    fig, ax = lab_figure(embed_width_cm, aspect=aspect)
    keys = list(series.keys())
    n = len(keys)
    x = np.arange(len(group_labels))
    cyc = _cycle(keys, colors)
    for i, k in enumerate(keys):
        vals = series[k]
        xs = [x[j] + (i - (n - 1) / 2.0) * width for j, v in enumerate(vals) if v is not None]
        ys = [v for v in vals if v is not None]
        ax.bar(xs, ys, width, label=k, color=cyc[k], edgecolor="black", linewidth=0.5, zorder=3)
    if log:
        ax.set_yscale("log")
    ax.set_ylabel(ylabel)
    ax.set_xticks(x)
    ax.set_xticklabels(group_labels)
    ax.grid(axis="y", which="both", ls=":", lw=0.5, alpha=0.6, zorder=0)
    legend_no_overlap(ax, fontsize=10)
    return fig, ax


def convergence(curves, ylabel, *, xlabel="iteration", embed_width_cm=7.6,
                aspect=0.66, colors=None, markers=None, ylim=None, logy=False):
    """One axes; each ``curves[label]`` plotted vs its 1-based index."""
    _get_mpl()
    fig, ax = lab_figure(embed_width_cm, aspect=aspect)
    keys = list(curves.keys())
    cyc = _cycle(keys, colors)
    mk = markers or {k: "os^vD"[i % 5] for i, k in enumerate(keys)}
    for k in keys:
        y = np.asarray(curves[k], float)
        kx = np.arange(1, len(y) + 1)
        (ax.semilogy if logy else ax.plot)(kx, y, mk[k] + "-", ms=3, color=cyc[k], label=k)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if ylim:
        ax.set_ylim(*ylim)
    ax.grid(alpha=0.3)
    legend_no_overlap(ax)
    return fig, ax


def quiver_pair(X, Y, vec_list, *, tags=("(a)", "(b)"), embed_width_cm=7.0,
                aspect=0.56, xlabel="x (mm)", ylabel="y (mm)", scale=11):
    """Two side-by-side quiver panels (e.g. aligned M vs circulating loop M).
    ``vec_list`` is a list of ``(U, V)`` arrays sampled at ``(X, Y)``."""
    _get_mpl()
    fig, ax = lab_figure(embed_width_cm, aspect=aspect, ncols=2)
    axes = np.atleast_1d(ax)
    for a, (U, V), tag in zip(axes, vec_list, tags):
        a.quiver(X, Y, U, V, np.hypot(U, V), cmap="viridis", scale=scale, width=0.010)
        a.set_aspect("equal")
        a.set_xlabel(xlabel)
        a.text(0.03, 0.97, tag, transform=a.transAxes, fontweight="bold", va="top")
    axes[0].set_ylabel(ylabel)
    return fig, ax


def bh_curve(H, B, *, embed_width_cm=6.0, aspect=0.62, xlabel="H (kA/m)", ylabel="B (T)"):
    """A monotonic B-H magnetization curve."""
    _get_mpl()
    fig, ax = lab_figure(embed_width_cm, aspect=aspect)
    ax.plot(H, B, "-", color=OKABE_ITO[2], lw=1.6)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, ls=":", lw=0.5, alpha=0.55)
    ax.set_xlim(0, None)
    ax.set_ylim(0, None)
    return fig, ax
