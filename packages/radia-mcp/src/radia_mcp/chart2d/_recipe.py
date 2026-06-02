"""Generate copy-pasteable Python recipes for each chart type.

Symmetric with _render.py: every renderer has a recipe-generator that
emits the equivalent Python code as a string.  The recipe always
ends in emit_paper_figure() so the user gets the same gate-stack
guarantees they would have gotten from server-side rendering.
"""
from __future__ import annotations

from textwrap import dedent

from ._charts import get_spec


def chart_recipe(
    chart_type: str,
    *,
    profile: str = "ieee_single_column",
    rel_width: float = 1.0,
    xlabel: str = r"$f$ (Hz)",
    ylabel: str = r"$|Z|$ ($\Omega$)",
) -> str:
    """Return a self-contained Python recipe for the chart type."""
    spec = get_spec(chart_type)
    head = (
        f"# {spec.description}\n"
        f"# group = {spec.group}\n"
        f"# profile = {profile}  (rel_width = {rel_width})\n"
        f"# Lab style: 10 pt @ 8 cm, Okabe-Ito CVD palette, no title in"
        f" figure, units in PARENTHESES (lab convention).\n\n"
        "from radia_mcp.figure import paper_figure, emit_paper_figure\n"
    )
    if spec.notes:
        head += f"\n# {spec.notes}\n"

    if spec.name == "bode":
        return head + _recipe_bode(profile, rel_width)
    if spec.name == "polar":
        return head + _recipe_polar(profile, rel_width)
    body = _recipe_generic(spec, profile, rel_width, xlabel, ylabel)
    return head + body


def _recipe_generic(spec, profile, rel_width, xlabel, ylabel) -> str:
    """One-axes recipe for the 20 non-Bode-non-polar types."""
    plot_call = _emit_plot_call(spec)
    return dedent(f"""
        import numpy as np

        # --- your data ---
        x = np.linspace(0, 10, 100)
        y = x ** 1.5   # replace with real data

        # --- build paper-quality figure ---
        fig, axes = paper_figure(
            profile={profile!r},
            nrows=1, ncols=1,
            rel_width={rel_width},
        )
        ax = axes[0, 0]

        # --- chart call ---
        {plot_call}

        ax.set_xlabel(r{xlabel!r})
        ax.set_ylabel(r{ylabel!r})
        # ax.legend(loc='best', frameon=False)   # if multiple traces

        # --- emit with all 5 gates ---
        emit_paper_figure(
            fig, 'out', {profile!r},
            min_axes_fraction=0.72,
            on_fail='raise',     # 'raise' | 'warn' | 'auto_tighten'
        )
    """).lstrip()


def _emit_plot_call(spec) -> str:
    """One matplotlib line for the dispatch call."""
    n = spec.name
    if n == "line":
        return "ax.plot(x, y, label='data')"
    if n == "loglog":
        return "ax.loglog(x, y, label='data')"
    if n == "semilogx":
        return "ax.semilogx(x, y, label='data')"
    if n == "semilogy":
        return "ax.semilogy(x, y, label='data')"
    if n == "step":
        return "ax.step(x, y, where='post')"
    if n == "errorbar":
        return ("yerr = 0.1 * y\n        "
                "ax.errorbar(x, y, yerr=yerr, fmt='o', capsize=3)")
    if n == "fill_between":
        return ("y_low, y_high = y - 0.1*y, y + 0.1*y\n        "
                "ax.fill_between(x, y_low, y_high, alpha=0.3, "
                "label='+/- 10%')")
    if n == "histogram":
        return ("data = np.random.normal(size=500)\n        "
                "ax.hist(data, bins=30, density=True)")
    if n == "bar":
        return ("categories = ['A','B','C','D']\n        "
                "heights = [3.5, 5.1, 2.7, 4.4]\n        "
                "ax.bar(categories, heights)")
    if n == "box":
        return ("data = [np.random.normal(size=80) for _ in range(4)]\n"
                "        labels = ['group1','group2','group3','group4']\n"
                "        ax.boxplot(data, tick_labels=labels)")
    if n == "violin":
        return ("data = [np.random.normal(size=120) for _ in range(4)]\n"
                "        ax.violinplot(data, showmedians=True)\n"
                "        ax.set_xticks([1,2,3,4])\n"
                "        ax.set_xticklabels(['g1','g2','g3','g4'])")
    if n == "ecdf":
        return ("data = np.random.normal(size=500)\n        "
                "ax.ecdf(data)")
    if n in ("contour", "contourf"):
        suffix = ", cmap='viridis'" if n == "contourf" else ""
        return ("X, Y = np.meshgrid(np.linspace(-2, 2, 50), "
                "np.linspace(-2, 2, 50))\n        "
                "Z = np.exp(-(X**2 + Y**2))\n        "
                f"cs = ax.{n}(X, Y, Z, levels=10{suffix})\n        "
                + ("ax.clabel(cs, inline=True, fontsize=8)"
                   if n == "contour"
                   else "fig.colorbar(cs, ax=ax)"))
    if n == "pcolormesh":
        return ("X, Y = np.meshgrid(np.linspace(-2, 2, 50), "
                "np.linspace(-2, 2, 50))\n        "
                "Z = np.exp(-(X**2 + Y**2))\n        "
                "m = ax.pcolormesh(X, Y, Z, cmap='viridis', "
                "shading='auto')\n        "
                "fig.colorbar(m, ax=ax)")
    if n == "quiver":
        return ("X, Y = np.meshgrid(np.linspace(-2, 2, 20), "
                "np.linspace(-2, 2, 20))\n        "
                "U, V = -Y, X    # rotation field\n        "
                "ax.quiver(X, Y, U, V, scale=30, scale_units='xy')")
    if n == "streamplot":
        return ("X = np.linspace(-2, 2, 30)\n        "
                "Y = np.linspace(-2, 2, 30)\n        "
                "XX, YY = np.meshgrid(X, Y)\n        "
                "U, V = -YY, XX\n        "
                "ax.streamplot(X, Y, U, V, density=1.2)")
    if n == "imshow":
        return ("Z = np.random.rand(50, 50)\n        "
                "im = ax.imshow(Z, cmap='viridis', origin='lower', "
                "extent=(-2, 2, -2, 2), aspect='auto')\n        "
                "fig.colorbar(im, ax=ax)")
    if n == "scatter":
        return ("ax.scatter(x, y, s=20, alpha=0.7)")
    if n == "phase":
        return ("z = (1 + 1j*x) / (1 - 0.3j*x)   # impedance locus\n"
                "        ax.plot(z.real, z.imag, marker='o', "
                "markersize=3)\n        "
                "ax.set_aspect('equal')\n        "
                "ax.axhline(0, color='0.7', lw=0.5, zorder=0)\n        "
                "ax.axvline(0, color='0.7', lw=0.5, zorder=0)")
    return f"# TODO: dispatch for {n!r}"


def _recipe_bode(profile, rel_width) -> str:
    return dedent(f"""
        import numpy as np

        # --- your H(jw) ---
        freq = np.logspace(0, 6, 200)                        # 1 Hz - 1 MHz
        omega = 2 * np.pi * freq
        H = 1.0 / (1.0 + 1j * omega / (2*np.pi*1e3))         # 1st-order LP
        mag_dB = 20 * np.log10(np.abs(H))
        phase_deg = np.angle(H, deg=True)

        # --- build 2x1 Bode pair ---
        fig, axes = paper_figure(
            profile={profile!r}, nrows=2, ncols=1, sharex=True,
            rel_width={rel_width},
            panel_labels=False,   # Bode is 'one figure', not (a)(b)
        )
        ax_mag, ax_phs = axes[0, 0], axes[1, 0]
        ax_mag.semilogx(freq, mag_dB)
        ax_phs.semilogx(freq, phase_deg)
        ax_mag.set_ylabel(r'$|H|$ (dB)')
        ax_phs.set_ylabel(r'$\\arg(H)$ ($\\degree$)')
        ax_phs.set_xlabel(r'$f$ (Hz)')
        # hide top panel's x tick text
        for tl in ax_mag.get_xticklabels():
            tl.set_visible(False)

        emit_paper_figure(
            fig, 'bode_out', {profile!r},
            min_axes_fraction=0.72, on_fail='raise',
        )
    """).lstrip()


def _recipe_polar(profile, rel_width) -> str:
    return dedent(f"""
        import numpy as np
        import matplotlib.pyplot as plt
        from radia_mcp.figure import paper_figure, emit_paper_figure

        # paper_figure() doesn't expose projection='polar' yet, so we
        # snapshot the profile rcParams + figsize via a throwaway call,
        # then build a polar axes manually:
        fig0, _ = paper_figure({profile!r}, rel_width={rel_width})
        figsize = fig0.get_size_inches()
        plt.close(fig0)

        fig = plt.figure(figsize=figsize, dpi=600)
        ax = fig.add_subplot(111, projection='polar')

        theta = np.linspace(0, 2*np.pi, 360)
        r = 1 + 0.5 * np.cos(3 * theta)        # 3-lobed rose
        ax.plot(theta, r)
        ax.set_theta_zero_location('N')         # 0 deg at top (compass)
        ax.set_theta_direction(-1)              # clockwise

        emit_paper_figure(
            fig, 'polar_out', {profile!r},
            min_axes_fraction=0.50,             # polar wastes corners
            on_fail='warn',                      # polar is naturally loose
        )
    """).lstrip()
