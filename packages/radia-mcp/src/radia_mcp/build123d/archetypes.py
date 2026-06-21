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

from build123d import Box, Compound, Cylinder, Pos

from .modeling import annular_segment, tube

__all__ = ["magnetization_tag", "parse_magnetization", "cylindrical_magnet", "block_magnet",
           "halbach_ring", "c_core", "solenoid"]


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
