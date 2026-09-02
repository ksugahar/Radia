"""Post-hoc Kelvin Periodic identification on an existing NGSolve mesh.

This module provides a **public API** for adding the Kelvin point-pair
identification on a NGSolve mesh that already has `kelvin_int` and
`kelvin_ext` boundary labels but no `Identifications` section yet.

Why this exists
---------------
There are three upstream paths that produce a Kelvin-ready `.vol`:

1. **Cubit C++ exporter** (`ExportNetgenCommand.cpp`): writes point-pair
   identifications during `export netgen` if every kelvin_int
   vertex matches a kelvin_ext vertex within tolerance after the rigid
   Kelvin translation.  All-or-nothing: a single unmatched vertex
   skips the C++ identification entirely.

2. **NGSolve OCC** + `Identify(IdentificationType.PERIODIC)`: must be
   called AFTER `Glue` and BEFORE `Mesh()`.  Bakes both point pairs
   and segment / surface-element pairs into the mesh.

3. **This module** (`add_kelvin_identification`): works on an
   ALREADY-LOADED `ngsolve.Mesh` -- no Cubit, no OCC, no re-export.
   Useful when:
     - the Cubit exporter skipped due to a single tolerance miss
       (typical 1/8-octant geometry, copy-mesh sub-pixel imprecision)
     - you received a `.vol` from elsewhere without identifications
     - you want to compare a "with-Kelvin" and a "without-Kelvin" run
       on the same mesh without re-running Cubit
     - you're scripting from a Jupyter notebook with NGSolve only
       (no Cubit license available)

Method
------
Walks `mesh.Elements(BND)`, collects vertex indices on `kelvin_int`
(boundary name) and `kelvin_ext`, builds 1-based netgen `PointIndex`
arrays, translates the interior set by the Kelvin offset (auto-detected
as `mean(kelvin_ext) - mean(kelvin_int)` if not provided), then solves
a Hungarian-algorithm assignment (`scipy.optimize.linear_sum_assignment`)
to find the optimal vertex pairing.  Falls back to a greedy nearest-
neighbour if scipy is unavailable.

Each surviving pair is registered via
``ngmesh.AddPointIdentification(p1, p2, identnr=1,
type=IdentificationType.PERIODIC)``.

NGSolve's `Periodic(H1(mesh, ...))` FES then deduces high-order edge /
face DOF coupling from the point pairs + mesh topology automatically;
verified 2026-04-25 on em_elf_quarter at p=2 (slaved=6085 DOFs,
boundary `Set(1)`-integral ratio kelvin_int/kelvin_ext = 1.0).

Relation to `calc_common.py`
----------------------------
This is the public, importable version of the helpers in
`radia.panels.calc_common`:
- `add_periodic_kelvin` (panel-internal, calls
  `_add_kelvin_point_idents_manual` after skip-if-existing check)
- `_add_kelvin_point_idents_manual` (the matching algorithm)
- `detect_kelvin_offset` (offset auto-detection)

This module re-exports them as `add_kelvin_identification(...)` and
`detect_kelvin_offset(...)` with cleaner names + a richer return dict.
The panel internals continue to use the `calc_common` versions for
backward compatibility (those are now thin shims calling here).

Public API
----------
- ``add_kelvin_identification(mesh, kelvin_offset=None, ...) -> dict``
- ``detect_kelvin_offset(mesh) -> tuple[float, float, float]``
- ``has_kelvin_identification(mesh) -> bool``

Example
-------
>>> from ngsolve import Mesh
>>> from radia.kelvin_identify_ngsolve import add_kelvin_identification
>>>
>>> mesh = Mesh("ih_fem_kelvin_demo.vol")  # has kelvin_int/kelvin_ext but
>>>                                         # NO Identifications section
>>> info = add_kelvin_identification(mesh)
>>> print(f"Added {info['n_pairs']} Kelvin pairs, max_dist={info['max_dist']:.2e}")
>>> # Now mesh.GetIdentifications() returns the new pairs and you can
>>> # build a Periodic FES on top.
"""
from __future__ import annotations

import math
from typing import Optional, Tuple

import numpy as np


__all__ = [
    "add_kelvin_identification",
    "detect_kelvin_offset",
    "has_kelvin_identification",
]


def detect_kelvin_offset(mesh) -> Tuple[float, float, float]:
    """Detect Kelvin translation offset from kelvin_int / kelvin_ext BND.

    The offset is the rigid translation mapping interior Kelvin sphere
    centre onto exterior Kelvin sphere centre.  Cubit's `copy mesh
    surface ... source vertex ... target vertex` and OCC's
    `Identify(IdentificationType.PERIODIC)` both produce 1:1 vertex
    correspondence; therefore
    ``mean(kelvin_ext) - mean(kelvin_int)`` is exact regardless of
    partial-octant geometry.

    Args:
        mesh: ngsolve.Mesh

    Returns:
        (x, y, z) offset.  Returns (0, 0, 0) if either boundary is
        missing from the mesh (i.e. nothing to detect).
    """
    from ngsolve import BND

    boundaries = mesh.GetBoundaries()
    if "kelvin_int" not in boundaries or "kelvin_ext" not in boundaries:
        return (0.0, 0.0, 0.0)

    int_pts: list = []
    ext_pts: list = []
    for el in mesh.Elements(BND):
        bc = boundaries[el.index]
        if bc == "kelvin_int":
            int_pts.extend(mesh.vertices[v.nr].point for v in el.vertices)
        elif bc == "kelvin_ext":
            ext_pts.extend(mesh.vertices[v.nr].point for v in el.vertices)

    if not int_pts or not ext_pts:
        return (0.0, 0.0, 0.0)

    int_mean = np.mean(int_pts, axis=0)
    ext_mean = np.mean(ext_pts, axis=0)
    return (float(ext_mean[0] - int_mean[0]),
            float(ext_mean[1] - int_mean[1]),
            float(ext_mean[2] - int_mean[2]))


def _boundary_vertex_ids(mesh, boundary: str) -> set[int]:
    """Return 1-based Netgen point IDs carried by one named boundary."""
    boundaries = mesh.GetBoundaries()
    vertices: set[int] = set()
    for element in mesh.Elements(__import__("ngsolve").BND):
        if boundaries[element.index] == boundary:
            vertices.update(int(vertex.nr) + 1 for vertex in element.vertices)
    return vertices


def has_kelvin_identification(mesh) -> bool:
    """Return whether the Kelvin boundaries, rather than any boundary, match.

    A rotational FFAG sector and a Kelvin exterior may both use NGSolve
    ``PERIODIC`` point pairs.  Treating a non-empty global identification list
    as Kelvin evidence makes a sector-only mesh look Kelvin-ready, then causes
    the reduced-potential solvers to apply an invalid exterior constraint.  A
    valid Kelvin identification must therefore pair every ``kelvin_int``
    boundary vertex with a ``kelvin_ext`` boundary vertex.
    """
    try:
        inner = _boundary_vertex_ids(mesh, "kelvin_int")
        outer = _boundary_vertex_ids(mesh, "kelvin_ext")
        if not inner or not outer:
            return False
        pairs = mesh.ngmesh.GetIdentifications()
        matched_inner: set[int] = set()
        matched_outer: set[int] = set()
        for first, second in pairs:
            # ``GetIdentifications`` returns PointId objects on current
            # Netgen, whose public numeric member is ``nr`` rather than
            # Python's ``__int__`` protocol.
            first = int(first.nr)
            second = int(second.nr)
            if first in inner and second in outer:
                matched_inner.add(first)
                matched_outer.add(second)
            elif second in inner and first in outer:
                matched_inner.add(second)
                matched_outer.add(first)
        return matched_inner == inner and matched_outer == outer
    except Exception:
        return False


def add_kelvin_identification(
    mesh,
    kelvin_offset: Optional[Tuple[float, float, float]] = None,
    *,
    inner_bnd: str = "kelvin_int",
    outer_bnd: str = "kelvin_ext",
    idnr: int | None = None,
    point_tolerance: Optional[float] = None,
    skip_if_existing: bool = True,
) -> dict:
    """Add Kelvin Periodic point-pair identification to an existing mesh.

    Walks the BND surface elements, collects vertices on `inner_bnd`
    and `outer_bnd`, translates the inner set by `kelvin_offset` (auto-
    detected if None), pairs them via Hungarian assignment
    (`scipy.optimize.linear_sum_assignment`), and registers each pair
    with ``ngmesh.AddPointIdentification(p1, p2, identnr,
    type=IdentificationType.PERIODIC)``.

    NGSolve's `Periodic(H1(mesh, ...))` then derives high-order DOF
    coupling from the point pairs + mesh topology automatically; no
    explicit segment / surface-element identification is added by this
    helper (verified 2026-04-25 on em_elf_quarter at p=2, slaved=6085).

    Args:
        mesh: ngsolve.Mesh with `inner_bnd` and `outer_bnd` boundary
            labels (and matching 1:1 vertex count on the two sides).
        kelvin_offset: (x, y, z) rigid translation from interior to
            exterior Kelvin sphere.  Pass None to auto-detect via
            ``detect_kelvin_offset(mesh)`` (recommended).
        inner_bnd: boundary name for the interior Kelvin sphere
            (typically the air sphere outer surface).
        outer_bnd: boundary name for the exterior Kelvin shell.
        idnr: identification number for ``AddPointIdentification``.  It
            defaults to 1 for a mesh without other periodic constraints.  A
            mesh that already has non-Kelvin identifications (for example an
            FFAG azimuthal sector) must pass a distinct, explicit number.
        point_tolerance: max distance for vertex-pair acceptance [m].
            Default = max(5% of |kelvin_offset|, 1 mm).  Pairs farther
            than this are dropped (counted as `n_unmatched`).
        skip_if_existing: if True, return early only when the Kelvin boundary
            pairs already exist.  Unrelated rotational or mirror pairs do not
            suppress the Kelvin operation.

    Returns:
        dict with keys:
            ``n_pairs``: number of identification pairs registered
            ``n_unmatched``: vertices that exceeded `point_tolerance`
            ``kelvin_offset``: offset used (auto-detected or provided)
            ``max_dist``: largest residual distance after translation [m]
            ``point_tolerance``: tolerance used [m]
            ``was_pre_existing``: True if `skip_if_existing` returned early

    Raises:
        ValueError: if neither inner_bnd nor outer_bnd is in
            mesh.GetBoundaries() (i.e. the mesh has no Kelvin labels).
            This is a configuration error -- you cannot add Kelvin
            identification to a mesh without Kelvin labels.

    Notes:
        - This function does NOT verify that the kelvin shell mesh is
          a 1:1 image of the air outer surface.  If the two surfaces
          have different vertex counts the Hungarian pairing will leave
          unmatched vertices (counted in `n_unmatched`); the FES will
          still be correct over the matched subset but the open-
          boundary truncation will be partial.  Check
          ``info["n_unmatched"] == 0`` to confirm full coverage.
        - The `idnr` argument matters only if you plan to add multiple
          independent identification sets to the same mesh (rare).
          Default `1` is conventional.
    """
    from netgen.meshing import IdentificationType

    # --- Validate boundaries -------------------------------------------------
    boundaries = mesh.GetBoundaries()
    have_inner = inner_bnd in boundaries
    have_outer = outer_bnd in boundaries
    if not (have_inner and have_outer):
        missing = [n for n, ok in [(inner_bnd, have_inner),
                                   (outer_bnd, have_outer)] if not ok]
        raise ValueError(
            f"Mesh does not have the required Kelvin boundary labels "
            f"(missing: {missing}; available: {sorted(set(boundaries))}). "
            f"Cannot add Kelvin identification to a mesh without "
            f"kelvin_int / kelvin_ext labels.  Either rebuild the mesh "
            f"with add_kelvin_cubit / OCC Identify, or pass the correct "
            f"inner_bnd / outer_bnd names.")

    # --- Skip-if-existing pre-flight ----------------------------------------
    kelvin_already_identified = has_kelvin_identification(mesh)
    if skip_if_existing and kelvin_already_identified:
        return {
            "n_pairs": 0,
            "n_unmatched": 0,
            "kelvin_offset": kelvin_offset or (0.0, 0.0, 0.0),
            "max_dist": 0.0,
            "point_tolerance": 0.0,
            "was_pre_existing": True,
        }
    existing_pairs = tuple(mesh.ngmesh.GetIdentifications())
    if idnr is None:
        if existing_pairs and not kelvin_already_identified:
            raise ValueError(
                "mesh already has non-Kelvin point identifications; pass an "
                "explicit distinct idnr when adding the Kelvin constraint")
        idnr = 1
    idnr = int(idnr)
    if idnr < 1:
        raise ValueError("Kelvin idnr must be a positive integer")

    # --- Auto-detect kelvin_offset if not provided --------------------------
    if kelvin_offset is None:
        kelvin_offset = detect_kelvin_offset(mesh)
        if kelvin_offset == (0.0, 0.0, 0.0):
            raise ValueError(
                "detect_kelvin_offset(mesh) returned (0, 0, 0).  Either "
                "the kelvin_int / kelvin_ext boundaries have no surface "
                "elements, or the two spheres are at the same centre "
                "(degenerate Kelvin geometry).  Pass kelvin_offset "
                "explicitly to bypass auto-detection.")

    # --- Tolerance default ---------------------------------------------------
    if point_tolerance is None:
        off_len = math.sqrt(sum(c * c for c in kelvin_offset))
        point_tolerance = max(0.05 * off_len, 1e-3)

    # --- Collect vertex PointIndex sets on each side ------------------------
    from ngsolve import BND
    inner_pids = _boundary_vertex_ids(mesh, inner_bnd)
    outer_pids = _boundary_vertex_ids(mesh, outer_bnd)

    if not inner_pids or not outer_pids:
        return {
            "n_pairs": 0,
            "n_unmatched": 0,
            "kelvin_offset": tuple(kelvin_offset),
            "max_dist": 0.0,
            "point_tolerance": point_tolerance,
            "was_pre_existing": False,
        }

    inner_pids_l = sorted(inner_pids)
    outer_pids_l = sorted(outer_pids)
    inner_pos = np.array([mesh.vertices[p - 1].point for p in inner_pids_l])
    outer_pos = np.array([mesh.vertices[p - 1].point for p in outer_pids_l])
    offset_v = np.array(kelvin_offset, dtype=float)

    # --- Translate + pair (Hungarian, greedy fallback) ----------------------
    inner_pos_t = inner_pos + offset_v
    try:
        from scipy.optimize import linear_sum_assignment
        cost = np.sqrt(np.sum(
            (inner_pos_t[:, None, :] - outer_pos[None, :, :]) ** 2, axis=2))
        row_ind, col_ind = linear_sum_assignment(cost)
        pairs: list = []
        max_dist = 0.0
        n_unmatched = 0
        for i, j in zip(row_ind, col_ind):
            d = float(cost[i, j])
            if d <= point_tolerance:
                pairs.append((inner_pids_l[i], outer_pids_l[j]))
                max_dist = max(max_dist, d)
            else:
                n_unmatched += 1
    except ImportError:
        used: set = set()
        pairs = []
        max_dist = 0.0
        n_unmatched = 0
        for i, ipt in enumerate(inner_pos_t):
            d2_all = np.sum((outer_pos - ipt) ** 2, axis=1)
            if used:
                d2_all = d2_all.copy()
                for j in used:
                    d2_all[j] = float("inf")
            j_best = int(np.argmin(d2_all))
            d_best = float(np.sqrt(d2_all[j_best]))
            if d_best <= point_tolerance:
                pairs.append((inner_pids_l[i], outer_pids_l[j_best]))
                used.add(j_best)
                max_dist = max(max_dist, d_best)
            else:
                n_unmatched += 1

    # --- Register identifications -------------------------------------------
    ngmesh = mesh.ngmesh
    for p1, p2 in pairs:
        ngmesh.AddPointIdentification(p1, p2, identnr=idnr,
                                      type=IdentificationType.PERIODIC)

    return {
        "n_pairs": len(pairs),
        "n_unmatched": n_unmatched,
        "kelvin_offset": tuple(float(c) for c in kelvin_offset),
        "max_dist": float(max_dist),
        "point_tolerance": float(point_tolerance),
        "was_pre_existing": False,
    }
