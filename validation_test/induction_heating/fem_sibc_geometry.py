"""Meshes and coils for the 3D FEM-SIBC induction-heating validations.

The workpiece is a hole in the air mesh (its surface is the ``sibc``
boundary); the air is a sphere of radius ``r_air`` with ``outer`` A = 0.
A ring region around the part, invariant under rotation about z, is meshed
finer, so turning the part (world frame) or the coil (body frame) sees the
same mesh quality.  Meshes are curved to order 2 before saving; the .vol
keeps the curved elements.

``air_mesh`` takes a function that builds the part, and builds it afresh:
meshing a part shape changes its state, so re-using one shape object gives
a different mesh the second time (measured: 19507 against 19442 elements
for the same inputs).
"""
from __future__ import annotations

import math

import numpy as np


def cylinder_part(R, H, *, bore=0.0, cross_hole=0.0, hole_maxh=None):
    """Cylinder of radius R, height H centred on the origin, axis z; an
    optional axial through bore and an optional through cross hole along
    x (both radii).  ``hole_maxh`` refines the hole walls and so the edges
    where they meet the outer surface, at which the field is singular."""
    from netgen.occ import Axes, Cylinder, Pnt, X, Z
    part = Cylinder(Axes(Pnt(0, 0, -H / 2), Z), r=R, h=H)
    if bore:
        part = part - Cylinder(Axes(Pnt(0, 0, -H), Z), r=bore, h=2 * H)
    if cross_hole:
        hole = Cylinder(Pnt(-2 * R, 0, 0), X, r=cross_hole, h=4 * R)
        if hole_maxh:
            hole.faces.maxh = hole_maxh
        part = part - hole
    return part


def air_mesh(make_part, *, angle=0.0, r_air, ring_r, ring_h, maxh_part,
             maxh_ring, maxh_far, order=2):
    """Curved air mesh with the part ``make_part()`` (turned by ``angle``
    about z) as a hole: boundaries ``sibc`` (the part) and ``outer``."""
    from netgen.occ import Axes, Axis, Cylinder, Glue, OCCGeometry, Pnt, \
        Sphere, Z
    from ngsolve import Mesh
    if not callable(make_part):
        raise TypeError("air_mesh needs a function that builds the part; a "
                        "shape that was meshed before gives another mesh")
    part = make_part()
    if angle:
        part = part.Rotate(Axis((0, 0, 0), Z), math.degrees(angle))
    part.faces.name = "sibc"
    for face in part.faces:
        # keep a finer size set on the hole walls
        face.maxh = min(face.maxh, maxh_part)
    ring = Cylinder(Axes(Pnt(0, 0, -ring_h / 2), Z), r=ring_r, h=ring_h)
    near = ring - part
    near.mat("air")
    near.maxh = maxh_ring
    sphere = Sphere(Pnt(0, 0, 0), r_air)
    sphere.faces.name = "outer"
    far = sphere - ring
    far.mat("air")
    mesh = Mesh(OCCGeometry(Glue([far, near])).GenerateMesh(maxh=maxh_far))
    if order > 1:
        mesh.Curve(order)
    names = set(mesh.GetBoundaries())
    if not {"sibc", "outer"} <= names:
        raise RuntimeError(f"air mesh boundaries {sorted(names)} lack "
                           "'sibc' or 'outer'")
    return mesh


def loop_paths(center, normal, radius, *, n_seg=720, angle=0.0):
    """One circular filament loop, turned by ``angle`` about z; returns the
    ``(p1, p2)`` segment list."""
    c = np.asarray(center, float)
    n = np.asarray(normal, float) / np.linalg.norm(normal)
    a = np.cross(n, [0, 0, 1.0]) if abs(n[2]) < 0.9 else np.cross(n, [1.0, 0, 0])
    a /= np.linalg.norm(a)
    b = np.cross(n, a)
    s = np.linspace(0, 2 * np.pi, n_seg + 1)
    pts = c + radius * (np.cos(s)[:, None] * a + np.sin(s)[:, None] * b)
    cz, sz = math.cos(angle), math.sin(angle)
    pts = pts @ np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]]).T
    return [(pts[i], pts[i + 1]) for i in range(n_seg)]
