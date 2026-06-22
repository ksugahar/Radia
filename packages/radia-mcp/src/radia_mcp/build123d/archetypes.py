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

import itertools
import math
import re

from build123d import (BuildLine, BuildSketch, Box, Circle, Compound, Cylinder, Helix, Plane, Polyline,
                       Pos, Rectangle, RegularPolygon, Rot, Spline, extrude, make_face, sweep)

from .modeling import annular_segment, assembly, polar_array, tube

__all__ = ["magnetization_tag", "parse_magnetization", "magnetization_map", "cylindrical_magnet",
           "block_magnet", "halbach_ring", "c_core", "solenoid", "pole_tip", "multipole_yoke",
           "h_dipole", "helmholtz_pair", "cos_theta_dipole", "e_core", "slotted_stator", "spm_rotor",
           "litz_packing_radius", "litz_fill_factor", "litz_wire", "hierarchical_litz",
           "rectangular_litz"]


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


def e_core(width, height, depth, leg_width, back_thickness, name="iron_core"):
    r"""An **E-core** (transformer / inductor): a rectangular block (``width`` x ``height`` x ``depth``)
    with two windows cut out, leaving a back spine of ``back_thickness`` and three legs of ``leg_width``
    (centre + two outer).  Pair with two windows for windings around the centre leg.  Returns one iron
    Solid.  (For an EI-core, add a separate bar across the open face.)"""
    ww = (width - 3 * leg_width) / 2.0
    wh = height - back_thickness
    if ww <= 0 or wh <= 0:
        raise ValueError("legs/back too thick for the given width/height")
    win = Box(ww, wh, depth + 2)
    y_win = -height / 2 + back_thickness + wh / 2
    x1 = -width / 2 + leg_width + ww / 2
    x2 = width / 2 - leg_width - ww / 2
    core = (Box(width, height, depth) - Pos(x1, y_win, 0) * win - Pos(x2, y_win, 0) * win).solid()
    core.label = name
    return core


def slotted_stator(r_bore, r_yoke, n_slots, slot_depth, slot_span_deg, h, name="stator"):
    r"""A **slotted stator lamination**: an annular ring (bore ``r_bore``, yoke ``r_yoke``, height ``h``)
    with ``n_slots`` radial slots opening to the bore (depth ``slot_depth``, angular width
    ``slot_span_deg``), leaving the teeth between them.  The single-region iron core a motor / actuator
    winding sits in; feed it to the AGE rotating-machine solver (see
    ``radia_ngsolve.airgap_motor_workflow``).  Returns one iron Solid."""
    if not (0 < r_bore < r_yoke):
        raise ValueError("require 0 < r_bore < r_yoke")
    if slot_depth >= r_yoke - r_bore:
        raise ValueError("slot deeper than the stator radial thickness")
    stator = tube(r_bore, r_yoke, h, label=name)
    slot = annular_segment(r_bore - 0.5, r_bore + slot_depth, h, -slot_span_deg / 2, slot_span_deg / 2)
    slots = polar_array(slot, n_slots, 360.0)
    core = (stator - slots).solid()
    core.label = name
    return core


def spm_rotor(r_shaft, r_rotor, n_poles, magnet_thickness, magnet_span_deg, h, name="magnet"):
    r"""A **surface-PM (SPM) rotor**: an iron hub (``r_shaft < r < r_rotor``) carrying ``n_poles``
    surface magnets (annular segments, radial thickness ``magnet_thickness``, angular width
    ``magnet_span_deg``) whose easy axes are RADIAL and ALTERNATE N/S pole-to-pole (encoded in each
    magnet label via the magnetization convention).  Returns a labelled multi-region Compound
    (``rotor_iron`` + ``{name}_kk_M<deg>``); ``magnetization_map`` turns the labels into M vectors.
    The rotor half of a PMSM you can drop into the AGE solver."""
    if not (0 <= r_shaft < r_rotor):
        raise ValueError("require 0 <= r_shaft < r_rotor")
    hub = tube(r_shaft, r_rotor, h, label="rotor_iron")
    children = [hub]
    for k in range(n_poles):
        c = k * 360.0 / n_poles
        seg = annular_segment(r_rotor, r_rotor + magnet_thickness, h,
                              c - magnet_span_deg / 2, c + magnet_span_deg / 2)
        seg.label = magnetization_tag(name, k, c + (180.0 if k % 2 else 0.0))   # radial, alternating
        children.append(seg)
    return assembly(*children, label="spm_rotor")


def litz_packing_radius(n_strands, strand_radius):
    r"""Bundle radius for ``n_strands`` strands packed touching on a SINGLE LAYER (ring):
    ``R = strand_radius / sin(pi / n_strands)`` (neighbour centres ``2 R sin(pi/n)`` apart =
    ``2 strand_radius``).  Use for the single-served-layer Litz idealisation (``n_strands >= 3``)."""
    if n_strands < 3:
        raise ValueError("single-layer packing needs >= 3 strands")
    return strand_radius / math.sin(math.pi / n_strands)


def litz_fill_factor(n_strands, strand_radius, bundle_radius):
    r"""Copper **fill factor** of a round Litz bundle: ``n_strands * (strand_radius / bundle_radius)**2``
    (total strand copper area / bundle envelope area), where ``bundle_radius`` is the bundle's OUTER
    (envelope) radius.  For the single-layer touching ring the envelope is the strand-centre ring radius
    plus one strand radius, ``litz_packing_radius(n, rs) + rs``, which gives the closed form
    ``n sin(pi/n)**2 / (1 + sin(pi/n))**2``; real served / multi-layer Litz sits around 0.4-0.5 once
    insulation, twist take-up and voids are counted.  A pure number in ``(0, 1)`` -- the headline metric
    the rectangular / profiled constructions exist to raise."""
    if bundle_radius <= 0:
        raise ValueError("bundle_radius must be > 0")
    return n_strands * (strand_radius / bundle_radius) ** 2


def _sweep_section(path, section_2d):
    r"""Sweep a 2D ``section_2d`` (a build123d face/sketch such as ``Circle(r)``, ``Rectangle(w, h)`` or
    ``RegularPolygon(r, n)``) along ``path`` (a wire / edge), orienting the section normal to the path
    tangent at the start, and return the swept :class:`~build123d.Solid`."""
    sec = Plane(origin=path @ 0.0, z_dir=path % 0.0) * section_2d
    return sweep(sec, path=path).solid()


def _superposed_centerline(levels, indices, length, n):
    r"""Sampled centreline (``n+1`` points over axial ``length``) of one hierarchical-Litz strand: the
    **additive superposition** of one circular orbit per level.  ``levels`` is ``[(count, radius, pitch),
    ...]`` and ``indices`` the strand's per-level index; level ``l`` contributes radius ``radius_l`` at
    angular rate ``2*pi/pitch_l`` with phase ``2*pi*index_l/count_l``."""
    pts = []
    for s in range(n + 1):
        z = length * s / n
        x = y = 0.0
        for (count, radius, pitch), idx in zip(levels, indices):
            ph = 2.0 * math.pi * z / pitch + 2.0 * math.pi * idx / count
            x += radius * math.cos(ph)
            y += radius * math.sin(ph)
        pts.append((x, y, z))
    return pts


def litz_wire(n_strands, strand_radius, bundle_radius, length, pitch, name="litz", strand_section=None):
    r"""A **Litz wire**: ``n_strands`` insulated strands twisted together on a bundle of radius
    ``bundle_radius`` with a twist ``pitch`` (axial length per turn), over an axial ``length``.  Each
    strand is a conductor swept along its own helix; the strands are the same helix rotated by
    ``2*pi*i/n_strands`` (a phase-shifted twist), so the model is ONE swept helical strand replicated by
    :func:`~radia_mcp.build123d.modeling.polar_array`.  By default the strand section is a circle of
    radius ``strand_radius``; pass ``strand_section`` (any build123d 2D face, e.g. ``Rectangle(w, h)`` for
    a flat / edge wire or ``RegularPolygon(r, 6)`` for a compacted hexagonal strand) to override -- real
    Litz is not always round.  Returns a labelled multi-region :class:`~build123d.Compound`
    (``{name}_kk`` per strand) -- each strand a separate conductor region for AC-loss (skin / proximity)
    analysis (PEEC / FE), the whole point of Litz wire.  Use :func:`litz_packing_radius` for the
    single-layer touching ``bundle_radius`` and :func:`litz_fill_factor` for the copper fill.
    """
    if n_strands < 1:
        raise ValueError("n_strands must be >= 1")
    helix = Helix(pitch=pitch, height=length, radius=bundle_radius)
    sec2d = Circle(strand_radius) if strand_section is None else strand_section
    strand = _sweep_section(helix, sec2d)
    strand.label = "strand"
    return polar_array(strand, n_strands, 360.0, label=name)


def hierarchical_litz(levels, strand_radius, length, name="litz", n_axial=160, strand_section=None):
    r"""A **hierarchical (multi-level) Litz** cable -- a *coiled coil* / bundle-of-bundles, the real
    construction of high-strand-count Litz (e.g. ``5x5x5``).  ``levels`` is ``[(count, radius, pitch),
    ...]`` ordered OUTER cabling -> INNER strand orbit; the cable carries ``prod(count)`` strands, one for
    every combination of per-level indices.  Each strand's centreline is the **additive superposition**
    (:func:`_superposed_centerline`) of one circular orbit per level -- the analytic idealisation in which
    the orbits add in the lab frame (not carried / Frenet), matching the closed-form ``sum r_l cos(.)``
    centreline.  The section defaults to ``Circle(strand_radius)``; override with ``strand_section``.
    Returns a labelled :class:`~build123d.Compound` (``{name}_ii_jj_..`` per strand), each a separate
    conductor region.  ``n_axial`` is the centreline sampling (raise for many turns of the finest pitch).
    """
    if not levels:
        raise ValueError("levels must be non-empty")
    sec2d = Circle(strand_radius) if strand_section is None else strand_section
    strands = []
    for combo in itertools.product(*[range(int(c)) for (c, _, _) in levels]):
        pts = _superposed_centerline(levels, combo, length, n_axial)
        strand = _sweep_section(Spline(*pts), sec2d)
        strand.label = name + "_" + "_".join(f"{i:02d}" for i in combo)
        strands.append(strand)
    return assembly(*strands, label=name)


def rectangular_litz(nx, ny, strand_radius, pitch, length, twist_pitch=None, name="litz",
                     strand_section=None, n_axial=160):
    r"""A **rectangular / profiled Litz bundle** -- ``nx x ny`` strands on a rectangular grid (centre
    spacing ``pitch``), the slot-fill ("compacted" / rectangular) Litz idealisation that fills a
    rectangular winding window far better than a round bundle.  With ``twist_pitch`` the whole rectangular
    bundle is twisted about its axis (each strand follows a helix at its grid radius
    ``hypot(gx, gy)``); ``twist_pitch=None`` leaves the strands straight and parallel.  The section
    defaults to ``Circle(strand_radius)``; override with ``strand_section`` (e.g. ``Rectangle`` for a
    foil-like compacted strand).  Returns a labelled :class:`~build123d.Compound` (``{name}_ii_jj`` per
    strand), each a separate conductor region.  Pair with :func:`litz_fill_factor` (use the rectangle
    ``nx*pitch x ny*pitch`` as the envelope) to report the slot-fill gain over a round bundle.
    """
    if nx < 1 or ny < 1:
        raise ValueError("nx, ny must be >= 1")
    sec2d = Circle(strand_radius) if strand_section is None else strand_section
    gx = [(i - (nx - 1) / 2.0) * pitch for i in range(nx)]
    gy = [(j - (ny - 1) / 2.0) * pitch for j in range(ny)]
    strands = []
    for i in range(nx):
        for j in range(ny):
            rad = math.hypot(gx[i], gy[j])
            if twist_pitch is None or rad < 1e-9:
                sec = Plane(origin=(gx[i], gy[j], 0.0), z_dir=(0.0, 0.0, 1.0)) * sec2d
                strand = extrude(sec, amount=length).solid()
            else:
                ph0 = math.atan2(gy[j], gx[i])
                pts = [(rad * math.cos(2.0 * math.pi * (length * s / n_axial) / twist_pitch + ph0),
                        rad * math.sin(2.0 * math.pi * (length * s / n_axial) / twist_pitch + ph0),
                        length * s / n_axial) for s in range(n_axial + 1)]
                strand = _sweep_section(Spline(*pts), sec2d)
            strand.label = f"{name}_{i:02d}_{j:02d}"
            strands.append(strand)
    return assembly(*strands, label=name)
