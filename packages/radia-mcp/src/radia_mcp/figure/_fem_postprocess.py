"""FEM post-processing plot recipes (Sugahara Lab).

Promoted on 2026-05-26 from a lab flux-line MATLAB workflow (2020-01).
The originals were a small MATLAB pipeline that drove an axisymmetric FE
field-probe (``getb`` / ``geta``-style) via ``ode45`` to trace a
charged-particle / flux-line trajectory through an axisymmetric magnet,
overlay the trajectory on the FEM mesh, and inspect the field along a
sample line with element-boundary tick marks.

Three plot recipes are exposed here in matplotlib form, plus one
trajectory-integration helper, all lab-style-conforming (Times New
Roman, 10 pt body text, units in parentheses).  They are designed to
be composed with ``paper_figure`` / ``paper_figure_8cm``:

  - ``plot_flux_line_trajectory(ax, traj_x, traj_y, mesh_*, source_xy)``
      Overlays the FEM mesh wireframe (blue thin lines) with the black
      trajectory and a red filled circle at the source.  Square aspect
      with equal data scales.  Used in the original ``flux_line.png``.

  - ``plot_field_probe_line(ax, y, Bx, By, element_change_indices)``
      Plots ``Bx`` and ``By`` versus ``y`` along a sample line, with
      vertical dotted lines at element-boundary crossings.  This
      visualises the PIECEWISE-LINEAR nature of a first-order FEM
      solution (the field is constant per element when computed from
      ``grad A_z``).  Used in the original ``B_in_elements.png``.

  - ``plot_field_vs_arclength(ax, s, Bx, By, Az=None)``
      Plots ``Bx(s)``, ``By(s)``, and optionally ``Az(s)`` versus the
      flux-line arc-length parameter ``s``.  Used in the original
      ``s-B.png`` and ``s-Az.png``.  The arc-length has units of m/T
      because the flux-line ODE is ``dr/ds = B``, not ``dr/dt = v``.

  - ``trace_flux_line(B_func, x0, s_span, ...)``
      Modern scipy.integrate.solve_ivp replacement for MATLAB's
      ``ode45`` driving ``dx/ds = Bx, dy/ds = By``.  Returns
      (s_array, xy_array) where xy_array.shape == (len(s), 2).

REFERENCE EQUATIONS (FROM 磁束線の方程式.docx + main.m):

  A magnetic flux line is parameterised by arc-length ``s`` (m/T) as

      dx/ds = B_x(x, y)
      dy/ds = B_y(x, y)

  Integrating from a seed point (x_0, y_0) traces the magnetic-flux
  contour through that point.  Along the trajectory, ``A_z`` is
  CONSTANT (in the ideal continuum limit) because ``B = grad x (A_z
  z_hat)`` is by construction perpendicular to grad A_z.  The
  ``A_z(s)`` plot in the original deck is therefore a diagnostic for
  the FEM piecewise-linear error: the steps in ``A_z(s)`` are jumps
  across element boundaries.

  For Hamiltonian systems (charged-particle tracking with H = (p^2 -
  qA)^2 / 2m), the lab also documents SYMPLECTIC INTEGRATORS:

      explicit (first order):   x_{n+1} = x_n + dt (dH/dp)_{x_n, p_n}
                                p_{n+1} = p_n - dt (dH/dx)_{x_{n+1}, p_n}

      implicit (first order):   x_{n+1} = x_n + dt (dH/dp)_{x_{n+1}, p_n}
                                p_{n+1} = p_n - dt (dH/dx)_{x_{n+1}, p_n}

  Both preserve the symplectic 2-form dx ∧ dp = dq ∧ dp exactly,
  which gives long-term energy conservation in cyclic orbits where
  generic RK4/ode45 drifts.  The 2020-01-06 lab note compared the
  two for circular-orbit tracking in a uniform field.

These recipes use **matplotlib only** — no scipy.integrate dependency
at import time (``trace_flux_line`` imports lazily so the rest of the
module loads on a stripped-down scientific Python stack).
"""

from __future__ import annotations

from typing import Optional, Sequence, Callable, Tuple

import numpy as np


# Lab default colour pair for (Bx, By) components.
# Red = x-component (horizontal), Blue = y-component (vertical) —
# matches the colour convention in the lab MATLAB scripts and is
# distinguishable in greyscale print (red is darker).
_BX_COLOR = "#D55E00"   # Okabe-Ito vermillion (CVD-safe red)
_BY_COLOR = "#0072B2"   # Okabe-Ito blue (CVD-safe blue)


# ============================================================
# trace_flux_line — scipy.integrate.solve_ivp wrapper
# ============================================================


def trace_flux_line(B_func: Callable,
                    x0: Sequence[float],
                    s_span: Tuple[float, float] = (0.0, 10000.0),
                    *,
                    method: str = "RK45",
                    max_step: Optional[float] = None,
                    rtol: float = 1e-3,
                    atol: float = 1e-6,
                    dense_output: bool = True,
                    n_eval: Optional[int] = None,
                    ):
    """Integrate the magnetic flux line dx/ds = B(x,y) from x0.

    Modern scipy replacement for the lab's MATLAB ``ode45(@velocity,
    [0,10000], [5,0])`` flux-line tracing pattern.

    Args:
        B_func: callable ``(x, y) -> (Bx, By)``.  Both inputs may be
            scalars or 1-D arrays; the integrator passes scalars.
            In the lab MATLAB, this was ``mo_getb(x, y)`` driven by
            an open FEMM solution.  Use any equivalent in Radia /
            NGSolve land (``rad.Fld`` on a Radia object, NGSolve
            ``CoefficientFunction`` evaluation, etc.).
        x0: ``(x0, y0)`` seed point.  Units must match what B_func
            returns: if x in mm and B in T then s has units of mm/T.
        s_span: ``(s_start, s_end)`` arc-length interval (m/T or mm/T
            depending on caller convention).
        method: solve_ivp method name (RK45 / RK23 / Radau / BDF /
            LSODA).  RK45 is the closest analogue to MATLAB ``ode45``.
        max_step: max step in arc-length units.  Use a small value
            (e.g. 1.0) to avoid the integrator skipping past small
            mesh elements.  Default None lets solve_ivp choose.
        rtol, atol: relative / absolute tolerances.
        dense_output: produce a dense interpolant.
        n_eval: if not None, evaluate the dense solution at this many
            equally-spaced points in s_span (returned as the 2nd
            element).  If None, returns the solver's chosen sample
            points.

    Returns:
        (s, xy) where s has shape (N,) and xy has shape (N, 2)
        (columns: x, y).

    Notes:
        * The lab MATLAB ``velocity.m`` scaled by 1e6 so that ode45
          would take reasonably-sized steps; this Python wrapper
          assumes the caller has already chosen units that give an
          arc-length range of O(1) to O(1e4).  If you see a stuck
          integrator, try scaling ``s_span`` and ``max_step``.
        * For closed flux lines (lab's circular case), one period of
          the orbit returns to (x0, y0); after N periods you can see
          ``A_z(s)`` drift due to per-element interpolation error
          — this is the diagnostic plot from ``s-Az.png``.

    Example:
        >>> import radia as rad
        >>> obj = rad.ObjBckg(lambda p: [0.0, 0.1, 0.0])
        >>> def B(x, y):
        ...     bx, by, _ = rad.Fld(obj, 'b', [x, y, 0.0])
        ...     return (bx, by)
        >>> s, xy = trace_flux_line(B, x0=(0.005, 0.0),
        ...                          s_span=(0.0, 1.0), max_step=0.01)
    """
    try:
        from scipy.integrate import solve_ivp
    except ImportError as e:
        raise ImportError(
            "trace_flux_line requires scipy: pip install scipy"
        ) from e

    def rhs(s, xy):
        bx, by = B_func(xy[0], xy[1])
        return [bx, by]

    kw = dict(method=method, rtol=rtol, atol=atol,
              dense_output=dense_output)
    if max_step is not None:
        kw["max_step"] = float(max_step)

    sol = solve_ivp(rhs, list(s_span), list(x0), **kw)

    if not sol.success:
        raise RuntimeError(
            f"trace_flux_line: solve_ivp failed: {sol.message}"
        )

    if n_eval is not None and dense_output:
        s_eval = np.linspace(s_span[0], s_span[1], int(n_eval))
        xy_eval = sol.sol(s_eval).T
        return s_eval, xy_eval

    return sol.t, sol.y.T


# ============================================================
# plot_flux_line_trajectory — trajectory + mesh + source marker
# ============================================================


def plot_flux_line_trajectory(ax,
                              traj_x: Sequence[float],
                              traj_y: Sequence[float],
                              *,
                              mesh_nodes: Optional[np.ndarray] = None,
                              mesh_elements: Optional[np.ndarray] = None,
                              source_xy: Optional[Sequence[float]] = None,
                              traj_color: str = "k",
                              traj_lw: float = 1.0,
                              traj_label: Optional[str] = "trajectory",
                              mesh_color: str = "#1f77b4",
                              mesh_lw: float = 0.3,
                              mesh_alpha: float = 0.6,
                              source_color: str = "#D55E00",
                              source_size: float = 30.0,
                              source_label: Optional[str] = None,
                              equal_aspect: bool = True,
                              ) -> dict:
    """Plot a flux-line trajectory over the FEM mesh + source marker.

    Reproduces the layout of the lab's ``flux_line.png``: thin blue mesh
    wireframe
    underneath, black trajectory line on top, red filled dot at the
    source.

    Args:
        ax: matplotlib Axes.
        traj_x, traj_y: trajectory point arrays (1-D, same length).
        mesh_nodes: shape ``(n_nodes, 2)`` or ``(2, n_nodes)`` array
            of (x, y) node coordinates.  Optional — pass None to skip
            the mesh wireframe.
        mesh_elements: shape ``(n_elem, 3)`` array of 0-indexed triangle
            connectivity referencing ``mesh_nodes``.  4-node quads also
            accepted (shape ``(n_elem, 4)``).
        source_xy: ``(x, y)`` source location for the red marker.
            None to skip.
        traj_color, traj_lw, traj_label: trajectory styling.
        mesh_color, mesh_lw, mesh_alpha: mesh wireframe styling.
        source_color, source_size: source-marker styling.
        source_label: legend label for the source (None = no entry).
        equal_aspect: True (default) calls ``ax.set_aspect('equal')``
            so the trajectory's circular orbit looks circular.

    Returns:
        dict with the matplotlib artists actually drawn (lines for
        mesh, line for trajectory, scatter for source).  Useful for
        post-tweaks (e.g. ``out['trajectory'].set_zorder(...)``).

    Lab convention (from the original .m + .png pair):
        - mesh thin blue lines (mesh_color='#1f77b4', lw=0.3)
        - trajectory black (traj_color='k', lw=1.0)
        - source red filled circle (source_color='#D55E00', size=30)
        - axes equal, no in-figure title.  Caller writes a caption.
    """
    artists = {"mesh_lines": [], "trajectory": None, "source": None}

    # --- Mesh wireframe (drawn FIRST so trajectory overlays it) ---
    if mesh_nodes is not None and mesh_elements is not None:
        nodes = np.asarray(mesh_nodes, dtype=float)
        elems = np.asarray(mesh_elements, dtype=int)
        # Accept (2, N) or (N, 2)
        if nodes.shape[0] == 2 and nodes.shape[1] != 2:
            nodes = nodes.T
        # Accept (3, M) or (M, 3) for triangles, similarly for quads
        if elems.ndim == 2 and elems.shape[0] in (3, 4) and elems.shape[1] not in (3, 4):
            elems = elems.T
        for tri in elems:
            # Close the polygon by repeating first vertex
            idx = list(tri) + [tri[0]]
            ln = ax.plot(nodes[idx, 0], nodes[idx, 1],
                          color=mesh_color, lw=mesh_lw, alpha=mesh_alpha,
                          zorder=1)
            artists["mesh_lines"].extend(ln)

    # --- Trajectory (black on top of mesh) ---
    tx = np.asarray(traj_x, dtype=float)
    ty = np.asarray(traj_y, dtype=float)
    line, = ax.plot(tx, ty,
                    color=traj_color, lw=traj_lw,
                    label=traj_label, zorder=3)
    artists["trajectory"] = line

    # --- Source marker (top of stack) ---
    if source_xy is not None:
        sx, sy = float(source_xy[0]), float(source_xy[1])
        sc = ax.scatter([sx], [sy], s=source_size,
                        c=source_color, zorder=5,
                        label=source_label)
        artists["source"] = sc

    if equal_aspect:
        ax.set_aspect("equal", adjustable="box")

    return artists


# ============================================================
# plot_field_probe_line — Bx, By vs y with element-boundary ticks
# ============================================================


def plot_field_probe_line(ax,
                          y: Sequence[float],
                          Bx: Sequence[float],
                          By: Sequence[float],
                          *,
                          element_change_indices: Optional[Sequence[int]] = None,
                          xlabel: str = r"$y$ (mm)",
                          ylabel: str = r"$B$ (T)",
                          bx_label: str = r"$B_x$",
                          by_label: str = r"$B_y$",
                          marker_size: float = 3.0,
                          line_width: float = 0.5,
                          boundary_color: str = "k",
                          boundary_ls: str = ":",
                          boundary_lw: float = 0.4,
                          legend: bool = True,
                          legend_loc: str = "best",
                          ) -> dict:
    """Plot Bx(y), By(y) along a sample line with element-boundary ticks.

    Reproduces the lab's ``B_in_elements.png``: red Bx, blue By, both with
    small markers, and vertical dotted lines wherever the probe
    crossed an element boundary.

    The vertical lines make visible the PIECEWISE-LINEAR / element-wise
    nature of a first-order FEM B-field (which is the curl of a
    piecewise-linear A_z, hence piecewise-constant per element).
    Reviewers find this much more informative than a smooth
    interpolated line that hides the discretisation.

    Args:
        ax: matplotlib Axes.
        y: 1-D array of probe positions (length N).
        Bx, By: 1-D field-component arrays at each probe point.
        element_change_indices: optional 1-D int array of indices into
            ``y`` where the probe enters a new FEM element.  If None,
            no vertical boundary lines are drawn.  Compute via
            ``np.where(np.abs(np.diff(elem_id_array)) > 0)[0]``.
        xlabel, ylabel: axis labels.
        bx_label, by_label: legend labels.
        marker_size: matplotlib marker size in pt (default 3.0).
        line_width: line width in pt (default 0.5; thin so the markers
            stand out).
        boundary_color, boundary_ls, boundary_lw: styling for the
            element-boundary vertical lines.
        legend: True (default) to draw a legend.
        legend_loc: passed to ``ax.legend(loc=...)``.

    Returns:
        dict with the drawn artists: ``{"bx_line": ..., "by_line": ...,
        "boundaries": [list of Line2D]}``.

    Lab default colours:
        Bx = vermillion (#D55E00), By = blue (#0072B2) — Okabe-Ito
        CVD-safe; distinguishable in greyscale.
    """
    y = np.asarray(y, dtype=float)
    bx_arr = np.asarray(Bx, dtype=float)
    by_arr = np.asarray(By, dtype=float)
    if y.shape != bx_arr.shape or y.shape != by_arr.shape:
        raise ValueError(
            f"plot_field_probe_line: y/Bx/By shapes must match, got "
            f"{y.shape}, {bx_arr.shape}, {by_arr.shape}"
        )

    # Element-boundary verticals first so they sit underneath the data
    boundary_lines = []
    if element_change_indices is not None:
        ymin, ymax = ax.get_ylim()  # may be (0, 1) until data plotted
        for idx in element_change_indices:
            xv = float(y[int(idx)])
            ln = ax.axvline(xv,
                            color=boundary_color, ls=boundary_ls,
                            lw=boundary_lw, alpha=0.6, zorder=1)
            boundary_lines.append(ln)

    # Field curves
    bx_line, = ax.plot(y, bx_arr,
                        color=_BX_COLOR, marker=".",
                        markersize=marker_size,
                        lw=line_width,
                        label=bx_label, zorder=3)
    by_line, = ax.plot(y, by_arr,
                        color=_BY_COLOR, marker=".",
                        markersize=marker_size,
                        lw=line_width,
                        label=by_label, zorder=3)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if legend:
        ax.legend(loc=legend_loc, frameon=False)

    return {"bx_line": bx_line, "by_line": by_line,
            "boundaries": boundary_lines}


# ============================================================
# plot_field_vs_arclength — Bx(s), By(s), Az(s)
# ============================================================


def plot_field_vs_arclength(ax,
                            s: Sequence[float],
                            Bx: Optional[Sequence[float]] = None,
                            By: Optional[Sequence[float]] = None,
                            Az: Optional[Sequence[float]] = None,
                            *,
                            s_unit: str = "m/T",
                            field_unit: str = "T",
                            az_unit: str = r"Wb/m$^2$",
                            line_width: float = 1.0,
                            bx_label: str = r"$B_x$",
                            by_label: str = r"$B_y$",
                            az_label: str = r"$A_z$",
                            ax_az: Optional["object"] = None,
                            legend: bool = True,
                            legend_loc: str = "best",
                            ) -> dict:
    """Plot B-field components (and optionally A_z) versus arc-length s.

    Reproduces the lab's ``s-B.png`` and ``s-Az.png``.  Either pass just
    (Bx, By) to get the ``s-B.png`` style or just Az to get the ``s-Az.png``
    style; pass all three to get a twin-axis plot with B on the
    primary y-axis and A_z on a twin y-axis.

    Args:
        ax: primary matplotlib Axes.
        s: arc-length values (1-D, length N).
        Bx, By, Az: component arrays at each s.  Any may be None.
        s_unit: shown in xlabel as ``$s$ ({s_unit})``.  Default ``m/T``
            because the flux-line ODE is ``dr/ds = B`` (arc-length is
            length per Tesla).
        field_unit: shown as ``$B$ ({field_unit})``.
        az_unit: shown as ``$A_z$ ({az_unit})`` — default Wb/m^2 matches
            the lab MATLAB plot.
        line_width: in pt.
        bx_label, by_label, az_label: legend labels.
        ax_az: pre-existing twin-y Axes to use for A_z.  If None and
            both B* and Az are supplied, ``ax.twinx()`` creates one
            automatically.
        legend: True (default) draws a legend.

    Returns:
        dict ``{"bx_line": ..., "by_line": ..., "az_line": ...,
        "ax_az": ...}`` with whatever artists were drawn (others None).

    A_z constancy diagnostic:
        For an IDEAL flux-line, A_z is CONSTANT along the curve.  The
        FEM piecewise-linear A_z + piecewise-constant B = grad A_z
        cause A_z(s) to drift in small steps; the size of the drift
        is the per-element interpolation error.  This is the original
        lab diagnostic from ``s-Az.png``.
    """
    s = np.asarray(s, dtype=float)
    out = {"bx_line": None, "by_line": None, "az_line": None, "ax_az": None}

    has_b = (Bx is not None) or (By is not None)
    has_az = Az is not None

    # B-components on primary axis
    if Bx is not None:
        bx_arr = np.asarray(Bx, dtype=float)
        if bx_arr.shape != s.shape:
            raise ValueError(
                f"plot_field_vs_arclength: Bx shape {bx_arr.shape} != "
                f"s shape {s.shape}"
            )
        out["bx_line"], = ax.plot(s, bx_arr,
                                   color=_BX_COLOR, lw=line_width,
                                   label=bx_label)
    if By is not None:
        by_arr = np.asarray(By, dtype=float)
        if by_arr.shape != s.shape:
            raise ValueError(
                f"plot_field_vs_arclength: By shape {by_arr.shape} != "
                f"s shape {s.shape}"
            )
        out["by_line"], = ax.plot(s, by_arr,
                                   color=_BY_COLOR, lw=line_width,
                                   label=by_label)

    if has_b:
        ax.set_xlabel(rf"$s$ ({s_unit})")
        ax.set_ylabel(rf"$B$ ({field_unit})")

    # A_z on twin axis (or solo on primary if no B at all)
    if has_az:
        az_arr = np.asarray(Az, dtype=float)
        if az_arr.shape != s.shape:
            raise ValueError(
                f"plot_field_vs_arclength: Az shape {az_arr.shape} != "
                f"s shape {s.shape}"
            )
        if has_b:
            ax_az = ax_az if ax_az is not None else ax.twinx()
            out["ax_az"] = ax_az
            out["az_line"], = ax_az.plot(s, az_arr,
                                          color="k", lw=line_width,
                                          label=az_label)
            ax_az.set_ylabel(rf"$A_z$ ({az_unit})")
        else:
            out["az_line"], = ax.plot(s, az_arr,
                                       color="k", lw=line_width,
                                       label=az_label)
            ax.set_xlabel(rf"$s$ ({s_unit})")
            ax.set_ylabel(rf"$A_z$ ({az_unit})")

    if legend:
        # Combine handles across primary + twin if both exist
        handles, labels = ax.get_legend_handles_labels()
        if out["ax_az"] is not None:
            h2, l2 = out["ax_az"].get_legend_handles_labels()
            handles += h2
            labels += l2
        if handles:
            ax.legend(handles, labels, loc=legend_loc, frameon=False)

    return out


# ============================================================
# Knowledge tool string — for MCP server consumption
# ============================================================


_FLUX_LINE_KNOWLEDGE = """\
[flux_line_recipe]

LAB FLUX-LINE TRACING + VISUALIZATION RECIPE (2020-01-06, Sugahara Lab)

Promoted from a lab flux-line MATLAB workflow — MATLAB driving an
axisymmetric FE field-probe to trace flux lines through a 2D magnet and
inspect the field along sample lines / along the trajectory.  All three
plot recipes are now reproducible in matplotlib via radia_mcp.figure:

  trace_flux_line(B_func, x0, s_span, ...)
  plot_flux_line_trajectory(ax, traj_x, traj_y, mesh_*, source_xy)
  plot_field_probe_line(ax, y, Bx, By, element_change_indices)
  plot_field_vs_arclength(ax, s, Bx, By, Az)

EQUATIONS (磁束線の方程式.docx, equation form):

  dx/ds = B_x(x, y)
  dy/ds = B_y(x, y)

Arc-length s has units of [length / Tesla].  Integrating this ODE
from a seed point (x0, y0) traces the magnetic-flux contour through
that point.  In the IDEAL continuum limit, A_z is CONSTANT along the
trajectory (because B = curl(A_z * z_hat) is perpendicular to grad
A_z by construction).  The piecewise-linear A_z used by first-order
FEM causes A_z(s) to drift in small steps -- those steps are the
per-element interpolation error.  Plotting A_z(s) is the lab's
canonical diagnostic for first-order-FEM accuracy along a flux line.

WHY THIS MATTERS:

  1. Validating a FEM solution: a real flux line should close on
     itself.  Trace one through your magnet, plot the trajectory; if
     it spirals outward or fails to close, the FEM mesh is too coarse
     OR you have a sign error in B_y.

  2. Visualising piecewise-linear FEM: the original B_in_elements.png
     plot shows Bx and By with vertical dotted lines at element
     boundaries.  This makes the discretisation HONEST in the figure
     -- reviewers can see that the field is piecewise-constant per
     element rather than a fictitious smooth curve.

  3. Charged-particle tracking: the same dx/ds = B integration scheme
     traces electron / proton trajectories through static magnets,
     useful for accelerator beam dynamics.  For non-zero rest mass
     and time-dependent fields, switch to a SYMPLECTIC INTEGRATOR
     (see below) to preserve energy over many orbits.

SYMPLECTIC INTEGRATORS (陽的_陰的シンプレクティック.jpg):

  For Hamiltonian H(x, p), 1st-order symplectic Euler:

      explicit:   x_{n+1} = x_n + dt (dH/dp)_{x_n, p_n}
                  p_{n+1} = p_n - dt (dH/dx)_{x_{n+1}, p_n}

      implicit:   x_{n+1} = x_n + dt (dH/dp)_{x_{n+1}, p_n}
                  p_{n+1} = p_n - dt (dH/dx)_{x_{n+1}, p_n}

  Both preserve the symplectic 2-form dx ^ dp exactly.  For circular
  orbits in a uniform B field, the symplectic scheme keeps the energy
  CONSTANT to the precision of the dt step, while a generic RK4 /
  ode45 drifts ~O(dt^4 * N_steps).  Use a symplectic integrator
  whenever you need long-term phase-space accuracy (many-orbit
  particle tracking, magnetic-mirror trapping, etc.).

  scipy.integrate.solve_ivp does NOT provide a symplectic option;
  for that use one of: SciML / DifferentialEquations.jl, scipy-symp,
  or a hand-rolled stepper.  trace_flux_line() here uses RK45 because
  the dx/ds = B problem is NOT Hamiltonian -- it's a 2D vector-field
  integration and RK45 is appropriate.

LAB COLOR CONVENTION (Bx vs By):

  Bx (horizontal) = vermillion #D55E00 (Okabe-Ito CVD-safe red)
  By (vertical)   = blue       #0072B2 (Okabe-Ito CVD-safe blue)

  These match the lab MATLAB scripts (which used pure red / pure blue)
  but are upgraded to the CVD-safe Okabe-Ito palette so they pass
  the lab's colorblind-safe gate.  Greyscale-distinguishable too.

UNITS GOTCHAS:

  * The lab MATLAB scaled the rhs by 1e6 in velocity.m so ode45
    would take sensibly-sized steps.  Choose your units (m vs mm)
    so that s_span * |B| ~ length range of interest.  If your seed
    is at radius 5 mm and B ~ 1e-8 T (small toroidal field) then
    s_span = (0, 1e7) m/T to trace one orbit, with max_step ~ 1e3.

  * The s-axis unit is m/T (or mm/T), NOT m -- always label as such.

  * A_z has units Wb/m on the s-axis plot — same as flux per
    unit length along the symmetry axis.

EXAMPLE: trace a flux line through a Radia/NGSolve magnet, plot all 3
diagnostic figures using paper_figure_8cm(panels=2):

    import radia as rad
    import numpy as np
    from radia_mcp.figure import (
        paper_figure_8cm, emit_paper_figure,
        trace_flux_line,
        plot_flux_line_trajectory,
        plot_field_vs_arclength,
    )

    # ... build your Radia model 'magnet' here ...
    def B_xy(x, y):
        bx, by, _ = rad.Fld(magnet, 'b', [x, y, 0.0])
        return (bx, by)

    # Trace one flux line
    s, xy = trace_flux_line(B_xy, x0=(0.005, 0.0),
                             s_span=(0.0, 1.0e6),
                             max_step=1e3, n_eval=2000)
    Az = [rad.Fld(magnet, 'a', [x, y, 0.0])[2] for x, y in xy]
    Bx = [rad.Fld(magnet, 'b', [x, y, 0.0])[0] for x, y in xy]
    By = [rad.Fld(magnet, 'b', [x, y, 0.0])[1] for x, y in xy]

    # Two panels at the 8 cm column: left = trajectory, right = s vs B
    fig, axes = paper_figure_8cm(panels=2, sharey=False)

    plot_flux_line_trajectory(axes[0, 0], xy[:, 0], xy[:, 1],
                              source_xy=(0.0, 0.0))
    axes[0, 0].set_xlabel(r'$x$ (m)')
    axes[0, 0].set_ylabel(r'$y$ (m)')

    plot_field_vs_arclength(axes[0, 1], s, Bx, By, Az=Az)

    emit_paper_figure(fig, 'flux_line_diagnostic',
                       'ieee_single_column')
"""


def get_flux_line_knowledge() -> str:
    """Return the lab flux-line tracing + visualization recipe text."""
    return _FLUX_LINE_KNOWLEDGE
