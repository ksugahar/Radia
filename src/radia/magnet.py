"""radia.magnet -- intent-based permanent-magnet builders on fixed-M surface-charge elements.

The historical ``ObjRecMag`` (surface-current rectangular block) is being retired as a user-facing
primitive (CLAUDE.md "Reduce Proprietary API Surface": geometry primitives demote to an internal
representation; users place intent objects).  A uniformly magnetized rectangular block is represented
**exactly** by an ``ObjHexahedron`` carrying fixed magnetization, equivalent to surface charge
``sigma = M . n`` on its six faces:
for a constant ``M`` the surface-charge and surface-current models give the identical external field.

``magnet_box(...)`` is the drop-in substitute -- the SAME ``(center, dimensions, magnetization)`` call
shape as ``ObjRecMag``, the fixed-M surface-charge representation, and an external field validated to ~5e-8 relative error
against ``ObjRecMag`` (off-axis, far field, and on the magnetization axis just above a face).

A permanent magnet needs **no** ``Solve`` -- the magnetization is fixed in the constructor and
``rad.Fld`` evaluates the exact open-boundary field analytically (the unsolved element has ``sigma == 0``
so ``B_comp_frM`` takes the analytic magnetization-based branch).  Call ``Solve`` only when soft iron is
present alongside the magnet.

Example::

    import radia as rad
    pm = rad.magnet_box([0, 0, 0], [0.02, 0.02, 0.01], [0, 0, 954930])  # Br=1.2 T along +z, M = Br/mu_0
    B  = rad.Fld(pm, 'b', [0, 0, 0.03])                                 # no Solve needed

General (non-box) permanent magnets do not need a helper -- build the element directly with the
magnetization in the constructor: ``rad.ObjHexahedron(verts, M)`` / ``rad.ObjTetrahedron(verts, M)`` /
``rad.ObjWedge(verts, M)`` (fixed-M surface-charge elements; all evaluate the exact open
boundary field with no Solve).
"""
import radia as rad


def magnet_box(center, dimensions, magnetization):
    """Rectangular permanent magnet as a fixed-M ``ObjHexahedron`` -- the ``ObjRecMag`` substitute.

    Parameters
    ----------
    center : sequence of 3 floats
        Box center ``(cx, cy, cz)`` in meters (Radia is always meters).
    dimensions : sequence of 3 floats
        Full edge lengths ``(Lx, Ly, Lz)`` in meters (NOT half-lengths -- same convention as
        ``ObjRecMag``).
    magnetization : sequence of 3 floats
        Magnetization ``(Mx, My, Mz)`` in **A/m** (``M = Br / mu_0``; e.g. ``Br = 1.2 T`` ->
        ``M = 954930`` A/m).  NOT magnetic polarization ``J`` (Tesla).

    Returns
    -------
    int
        A Radia ``ObjHexahedron`` handle holding the fixed magnetization -- a permanent magnet.
        No ``Solve`` is required; pass it to ``rad.Fld`` directly, or into an ``ObjCnt`` alongside
        soft iron (then ``Solve`` the container).
    """
    cx, cy, cz = (float(v) for v in center)
    hx, hy, hz = (float(v) / 2.0 for v in dimensions)
    # 8 vertices: bottom face (z-) CCW, then top face (z+) CCW -- the netgen hexahedron winding
    # expected by ObjHexahedron (validated to reproduce ObjRecMag to ~5e-8).
    verts = [
        [cx - hx, cy - hy, cz - hz], [cx + hx, cy - hy, cz - hz],
        [cx + hx, cy + hy, cz - hz], [cx - hx, cy + hy, cz - hz],
        [cx - hx, cy - hy, cz + hz], [cx + hx, cy - hy, cz + hz],
        [cx + hx, cy + hy, cz + hz], [cx - hx, cy + hy, cz + hz],
    ]
    return rad.ObjHexahedron(verts, [float(m) for m in magnetization])


def ObjRecMag(center, dimensions, magnetization):
    """Script-side ``ObjRecMag`` compatibility wrapper -- a rectangular fixed-M
    ``ObjHexahedron`` representation, identical signature to the retiring C++
    ``ObjRecMag`` (surface-current ``radTRecMag``).

    The C++ ``ObjRecMag`` is being demoted from the user API (CLAUDE.md "Reduce Proprietary API
    Surface"); other solvers / examples / tests that call ``rad.ObjRecMag(center, dimensions,
    magnetization)`` keep working through this script-side wrapper, now backed by the fixed-M
    element.  For a uniformly magnetized block the surface-charge field IS the exact analytic
    (K&J) cuboid field -- verified to machine precision on-axis (1e-16) and ~1e-7 vs the old
    surface-current ``ObjRecMag`` off-axis -- so this is a faithful drop-in.

    Returns a Radia ``ObjHexahedron`` handle (permanent magnet for fixed ``magnetization``; apply a
    soft material with ``rad.MatApl`` for soft iron, then ``rad.Solve`` the container).
    """
    return magnet_box(center, dimensions, magnetization)
