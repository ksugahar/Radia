# -*- coding: utf-8 -*-
r"""Parametric **EM-device archetypes** for build123d -- ready-made magnet/coil/yoke generators.

Composes the generic operations in :mod:`radia_mcp.build123d.modeling` (annular wedge, tube, arrays)
into the recurring building blocks of magnetostatic devices: permanent-magnet primitives, a segmented
**Halbach ring**, a **C-core** electromagnet yoke, and a **solenoid** winding bundle.  Each returns a
clean, labelled, Netgen-meshable build123d solid / Compound wired for the Radia / NGSolve pipeline.

MAGNETIZATION LABEL CONVENTION.  A permanent-magnet region carries its easy-axis (magnetization)
direction in its ``.label`` as an in-plane angle in degrees: ``"{name}_{index:02d}_M{angle:.1f}"`` (the
``_M<deg>`` suffix).  The solver side recovers it with :func:`parse_magnetization` and sets
``M = Br * (cos angle, sin angle, 0)``.  This is the labels-drive-magnetization idiom the lab's
``run_pipeline_multi`` carries through to Gmsh physical groups -- so a Halbach array's per-segment easy
axes survive the geometry -> mesh -> solver hand-off with no side-channel.

These archetypes replace the abridged snippets in the ``examples_lab_patterns`` documentation with
tested generators (the Halbach segment that was ``pass`` in the docs is now a real wedge with the
Mallinson easy-axis angle baked into each segment label).
"""
from __future__ import annotations

import math
import re

from build123d import (BuildLine, BuildSketch, Box, Compound, Cylinder, Plane, Polyline, Pos, Rot,
                       extrude, make_face)

from .modeling import annular_segment, assembly, polar_array, tube

__all__ = ["magnetization_tag", "parse_magnetization", "magnetization_map", "cylindrical_magnet",
           "block_magnet", "halbach_ring", "c_core", "solenoid", "pole_tip", "multipole_yoke",
           "h_dipole", "helmholtz_pair", "cos_theta_dipole"]


def magnetization_tag(name, index, angle_deg):
    r"""Build a PM region label encoding the in-plane magnetization angle (to 1e-3 deg):
    ``{name}_{index:02d}_M{deg:.3f}``."""
    return f"{name}_{index:02d}_M{angle_deg % 360.0:.3f}"


def parse_magnetization(label):
    r"""Recover the in-plane magnetization angle [deg] from a label's ``_M<deg>`` suffix (or ``None``)."""
    m = re.search(r"_M(-?\d+\.?\d*)", label or "")
    return float(m.group(1)) if m else None


def cylindrical_magnet(radius, h, m_angle_deg=0.0, name="magnet", index=0):
    r"""A solid cylindrical permanent magnet (z-centred), easy-axis angle ``m_angle_deg`` in the xy
    plane encoded in the label (see the magnetization-label convention)."""
    mag = Cylinder(radius=radius, height=h).solid()
    mag.label = magnetization_tag(name, index, m_angle_deg)
    return mag


def block_magnet(length, width, h, m_angle_deg=0.0, name="magnet", index=0):
    r"""A rectangular block permanent magnet (z-centred), easy-axis angle ``m_angle_deg`` encoded in
    the label."""
    mag = Box(length, width, h).solid()
    mag.label = magnetization_tag(name, index, m_angle_deg)
    return mag


def halbach_ring(r_in, r_out, h, n_segments, pole_pairs=1, name="halbach", gap_deg=0.0):
    r"""A **segmented Halbach permanent-magnet ring** (``r_in < r < r_out``, height ``h``,
    ``n_segments`` wedges).  Each segment is an :func:`~radia_mcp.build123d.modeling.annular_segment`
    whose label encodes the Halbach easy-axis angle

        ``alpha_k = (pole_pairs + 1) * theta_k``          (Mallinson / Halbach cylinder)

    at the segment centre ``theta_k`` -- so ``pole_pairs=1`` is the classic dipole array (uniform
    transverse field inside the bore, ``alpha`` advancing at twice the mechanical angle),
    ``pole_pairs=2`` a quadrupole, etc.  ``gap_deg`` leaves an inter-segment air gap.  Returns a
    labelled :class:`~build123d.Compound`; feed the per-segment ``_M<deg>`` labels to the solver via
    :func:`parse_magnetization`.
    """
    if n_segments < 3:
        raise ValueError("a Halbach ring needs >= 3 segments")
    seg_span = 360.0 / n_segments - gap_deg
    if seg_span <= 0.0:
        raise ValueError("gap_deg too large for n_segments")
    children = []
    for k in range(n_segments):
        a0 = k * 360.0 / n_segments
        seg = annular_segment(r_in, r_out, h, a0, a0 + seg_span)
        theta_c = a0 + seg_span / 2.0
        seg.label = magnetization_tag(name, k, (pole_pairs + 1) * theta_c)
        children.append(seg)
    return Compound(children=children, label=name)


def c_core(width, height, depth, leg, gap, name="core"):
    r"""A **C-core electromagnet yoke**: a rectangular window frame (outer ``width`` x ``height`` x
    ``depth``, wall thickness ``leg``) opened on the top centre by an air ``gap`` between two pole
    faces -- the canonical attractive-force / dipole yoke.  Returns a single labelled iron Solid;
    place a coil (see :func:`solenoid`) around a leg and a target across the gap.
    """
    if 2 * leg >= min(width, height):
        raise ValueError("leg too thick for the frame")
    if gap >= width - 2 * leg:
        raise ValueError("gap wider than the window")
    frame = Box(width, height, depth) - Box(width - 2 * leg, height - 2 * leg, depth + 2)
    core = (frame - Pos(0, height / 2, 0) * Box(gap, leg * 2.2, depth + 2)).solid()
    core.label = name
    return core


def solenoid(r_in, r_out, h, name="coil"):
    r"""A **solenoid winding bundle** as a hollow-cylinder (tube) current region, z-centred, labelled
    as a coil.  The azimuthal current density is the solver's to set; this is the meshable conductor
    volume (mean-turn bundle), the simplest realistic inductor / IH-coil region."""
    coil = tube(r_in, r_out, h, label=name)
    return coil


def magnetization_map(parts, Br=1.0):
    r"""Close the geometry -> field loop: for every PM region (a label carrying an ``_M<deg>`` suffix)
    return ``{label: (Mx, My)}`` with ``M = Br (cos angle, sin angle)`` -- the magnetization vector the
    solver assigns to that region.  ``parts`` may be a :class:`~build123d.Compound` (its ``.children``
    are scanned) or any iterable of labelled solids.  Non-magnet regions (no ``_M`` tag) are skipped.

    Example (Halbach -> per-segment easy axis -> solver):
        hb = halbach_ring(40, 55, 20, 12)
        Mmap = magnetization_map(hb, Br=1.2)        # {"hb_00_M30.000": (Mx, My), ...}
    """
    children = parts.children if isinstance(parts, Compound) else list(parts)
    out = {}
    for c in children:
        ang = parse_magnetization(getattr(c, "label", "") or "")
        if ang is not None:
            r = math.radians(ang)
            out[c.label] = (Br * math.cos(r), Br * math.sin(r))
    return out


def pole_tip(base_width, tip_width, height, depth, name="pole"):
    r"""A **trapezoidal pole piece** (base ``base_width`` at z=0 tapering to ``tip_width`` at z=height,
    ``depth`` along y), z-base at 0 pointing +z.  The building block of a shaped magnet pole; place it
    on a yoke (see :func:`multipole_yoke`, :func:`h_dipole`) with the tip toward the bore."""
    with BuildSketch(Plane.XZ) as sk:
        with BuildLine():
            Polyline((-base_width / 2, 0), (base_width / 2, 0),
                     (tip_width / 2, height), (-tip_width / 2, height), close=True)
        make_face()
    p = (Pos(0, -depth / 2, 0) * extrude(sk.sketch, amount=depth)).solid()
    p.label = name
    return p


def multipole_yoke(n_poles, r_bore, pole_len, pole_width, yoke_thickness, depth, name="magnet"):
    r"""An **n-pole iron yoke**: ``n_poles`` radial pole bars pointing inward to a bore of radius
    ``r_bore`` (length ``pole_len``, tangential width ``pole_width``), closed by an outer return ring
    (thickness ``yoke_thickness``).  ``n_poles=2`` is a dipole, ``4`` a quadrupole, ``6`` a sextupole.
    Returns a labelled multi-region :class:`~build123d.Compound` (``pole_kk`` + ``yoke``)."""
    if n_poles < 2:
        raise ValueError("need >= 2 poles")
    pole = (Pos(r_bore + pole_len / 2, 0, 0) * Box(pole_len, pole_width, depth)).solid()
    pole.label = "pole"
    poles = polar_array(pole, n_poles, 360.0, label="pole")
    ring = tube(r_bore + pole_len, r_bore + pole_len + yoke_thickness, depth, label="yoke")
    return assembly(*poles.children, ring, label=name)


def h_dipole(width, height, depth, leg, pole_width, gap, name="yoke"):
    r"""An **H-frame dipole yoke**: a rectangular window frame (outer ``width`` x ``height`` x
    ``depth``, wall ``leg``) with two poles protruding from the top and bottom of the window toward a
    central ``gap`` (the field region), pole tangential width ``pole_width``.  The classic
    accelerator/laboratory dipole.  Returns one fused iron Solid."""
    win_h = height - 2 * leg
    pole_h = (win_h - gap) / 2.0
    if pole_h <= 0:
        raise ValueError("gap too large for the window height")
    frame = Box(width, height, depth) - Box(width - 2 * leg, win_h, depth + 2)
    top = Pos(0, gap / 2 + pole_h / 2, 0) * Box(pole_width, pole_h, depth)
    bot = Pos(0, -(gap / 2 + pole_h / 2), 0) * Box(pole_width, pole_h, depth)
    yoke = (frame + top + bot).solid()
    yoke.label = name
    return yoke


def helmholtz_pair(r_in, r_out, h, separation, name="coil"):
    r"""A **Helmholtz coil pair**: two coaxial solenoid (tube) windings on the z-axis, separated by
    ``separation`` (for the classic uniform-centre field, set ``separation`` = mean coil radius).
    Returns a labelled :class:`~build123d.Compound` (``coil_0`` at -sep/2, ``coil_1`` at +sep/2)."""
    c1 = (Pos(0, 0, -separation / 2) * tube(r_in, r_out, h)).solid(); c1.label = "coil_0"
    c2 = (Pos(0, 0, separation / 2) * tube(r_in, r_out, h)).solid(); c2.label = "coil_1"
    return assembly(c1, c2, label=name)


def cos_theta_dipole(radius, conductor_w, conductor_h, length, n_per_half, name="coil"):
    r"""A **cos-theta dipole winding** approximated by ``2*n_per_half`` axial conductor bars on a shell
    of radius ``radius`` (cross-section ``conductor_w`` x ``conductor_h``, axial ``length``).  The bars
    are placed at **arcsin-spaced** angles so that ``sin(theta)`` is uniform -- i.e. the bar density
    is ``~ cos(theta)``, the current layout that produces a pure dipole field.  Two polarity groups
    (``{name}_go_kk`` for theta in (-90,90), ``{name}_ret_kk`` on the opposite side); the solver drives
    +I / -I.  Returns a labelled :class:`~build123d.Compound`."""
    if n_per_half < 2:
        raise ValueError("need >= 2 bars per half")
    bars = []
    for i in range(n_per_half):
        s = 2.0 * (i + 0.5) / n_per_half - 1.0                       # uniform in (-1, 1)
        th = math.degrees(math.asin(s))                             # arcsin spacing -> density ~ cos(theta)
        for off, grp in ((0.0, "go"), (180.0, "ret")):
            bar = (Rot(0, 0, th + off) * Pos(radius, 0, 0) * Box(conductor_w, conductor_h, length)).solid()
            bar.label = f"{name}_{grp}_{i:02d}"
            bars.append(bar)
    return Compound(children=bars, label=name)
