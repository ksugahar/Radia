"""radia.vim._shapes -- mesh-less-SHAPE intent constructors for soft iron -> HDiv-VIM.

The mesh-less soft-iron capability (build the iron from a simple SHAPE, no external NGSolve mesh) is KEPT
by these intent constructors.  Each internally STRUCTURE-MESHES the shape into a subdivided structured hex
mesh and registers it via :func:`radia.vim.MeshSoftIron`, so ``rad.Solve`` auto-routes it to the FEEC HDiv-VIM
(BDM1).  This is the API-compatible replacement for the retired ``ObjHexahedron + MatApl(MatLin) +
rad.Solve`` mesh-less surface-charge route.

Why a constructor and not transparent ``ObjHexahedron`` interception: Radia exposes no per-element vertex
getter, so a mesh-less handle's geometry is not recoverable at ``rad.Solve`` time -- the geometry must be
given up front.  These constructors take the shape explicitly (the intent-based user layer of the
"Reduce Proprietary API Surface" direction), build the mesh, and hand back a normal ``vim.MeshSoftIron``
container.

The returned container IS the subdivided mesh's Radia hex elements (built by ``vim.MeshSoftIron``), so
``rad.Solve(cont)`` -> HDiv-VIM BDM1 and ``rad.Fld(cont, ...)`` reflects the resolved per-sub-element M --
NO write-back plumbing is needed (the container already is the sub-mesh).  Exactly one of ``mu_r`` (linear)
or ``bh_table`` (nonlinear ``[[H,B],...]``) per iron.  The caller opens ``with ng.TaskManager():``.
"""
import numpy as np


def _structured_hex_mesh(mapping, nsub):
    """A structured all-hex NGSolve mesh of the unit cube pushed through ``mapping(X,Y,Z)`` (X,Y,Z in
    [0,1]).  ``nsub`` = subdivisions per dimension (int -> nxnxn, or a 3-tuple (nx,ny,nz))."""
    from ngsolve.meshes import MakeStructured3DMesh
    n = (int(nsub), int(nsub), int(nsub)) if np.isscalar(nsub) else tuple(int(v) for v in nsub)
    if min(n) < 1:
        raise ValueError("soft-iron shape: nsub must be >= 1 per dimension (got %r)" % (nsub,))
    return MakeStructured3DMesh(hexes=True, nx=n[0], ny=n[1], nz=n[2], mapping=mapping)


def soft_iron_box(center, size, mu_r=None, bh_table=None, nsub=4, material_filter=None, verbose=False):
    """Axis-aligned soft-iron BOX -> subdivided structured hex mesh -> HDiv-VIM.

    center = (cx,cy,cz), size = (sx,sy,sz) in METERS; the box spans ``center +- size/2``.  nsub = the
    subdivision count per dimension (int, or (nx,ny,nz)) -- the BDM1 mesh resolution; nsub=4 (=64 hexes)
    gives the converged cube demag ~1/3.  Exactly one of ``mu_r`` (linear, > 1) or ``bh_table`` (nonlinear
    [[H,B],...]).  Returns a :func:`radia.vim.MeshSoftIron` container (a mesh-backed hex iron): ``rad.Solve(cont)``
    -> HDiv-VIM BDM1, ``rad.Fld(cont, ...)`` reflects the resolved M.  (Caller opens ``with ng.TaskManager():``.)"""
    from ._radsolve import soft_iron_from_mesh as _mesh_soft_iron_impl
    cx, cy, cz = (float(v) for v in center)
    sx, sy, sz = (float(v) for v in size)
    if min(sx, sy, sz) <= 0.0:
        raise ValueError("soft_iron_box: size must be positive in every dimension (got %r)" % (size,))
    mesh = _structured_hex_mesh(
        lambda X, Y, Z: (cx + sx * (X - 0.5), cy + sy * (Y - 0.5), cz + sz * (Z - 0.5)), nsub)
    return _mesh_soft_iron_impl(mesh, mu_r=mu_r, bh_table=bh_table,
                                material_filter=material_filter, verbose=verbose)


# Trilinear (8-node) hex shape functions in the reference-cube corner order that MATCHES the CHEXA vertex
# order used by rad.ObjHexahedron (bottom z- face 0-3 CCW, top z+ face 4-7 CCW): corner k sits at the
# reference point (rx,ry,rz) below.
_HEX_REF = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                     [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)


def soft_iron_hex(vertices, mu_r=None, bh_table=None, nsub=4, material_filter=None, verbose=False):
    """General (possibly non-axis-aligned) soft-iron HEXAHEDRON from 8 vertices -> subdivided structured
    hex mesh (trilinear map of the reference cube onto the 8 corners) -> HDiv-VIM.

    vertices = 8 corners in the ``rad.ObjHexahedron`` CHEXA order (bottom face 0-3 CCW, top face 4-7 CCW),
    in METERS.  nsub = subdivisions per dimension.  Exactly one of ``mu_r`` / ``bh_table``.  Returns a
    :func:`radia.vim.MeshSoftIron` container (rad.Solve -> HDiv-VIM BDM1)."""
    from ._radsolve import soft_iron_from_mesh as _mesh_soft_iron_impl
    V = np.asarray(vertices, float)
    if V.shape != (8, 3):
        raise ValueError("soft_iron_hex: vertices must be 8x3 (CHEXA order); got shape %r" % (V.shape,))

    def tri(X, Y, Z):
        x = np.asarray(X, float); y = np.asarray(Y, float); z = np.asarray(Z, float)
        P = np.zeros(x.shape + (3,))
        for k in range(8):
            rx, ry, rz = _HEX_REF[k]
            w = (x if rx else 1 - x) * (y if ry else 1 - y) * (z if rz else 1 - z)
            P += w[..., None] * V[k]
        return (P[..., 0], P[..., 1], P[..., 2])

    mesh = _structured_hex_mesh(tri, nsub)
    return _mesh_soft_iron_impl(mesh, mu_r=mu_r, bh_table=bh_table,
                                material_filter=material_filter, verbose=verbose)


# --------------------------------------------------------------------------------------------------
# PERMANENT-MAGNET shape constructors (the "Magnet(...)" half of the intent-based user layer).
#
# A uniform permanent magnet is a FIXED-M body whose field is EXACT analytically (surface-current /
# surface-charge closed form) -- it needs NO mesh and NO demag solve, unlike soft iron.  So these
# constructors return a SINGLE Radia element (not a subdivided mesh): rad.ObjRecMag for an axis-aligned
# box, rad.ObjHexahedron for a general hexahedron.  They exist for API SYMMETRY with soft_iron_box/hex
# and for the composition below.
#
# PM + SOFT IRON: a permanent magnet is a SOURCE.  Place an analytic magnet element alongside a soft_iron_box in
# one container and rad.Solve auto-routes the (registered) iron to the HDiv-VIM with the magnet's field as
# the applied H_ext --
#     iron = radia.vim.soft_iron_box(center=..., size=..., mu_r=1000)
#     mag  = radia.vim.magnet_box(center=..., size=..., M=(0, 0, Br/MU0))
#     res  = rad.Solve(rad.ObjCnt([iron, mag]))     # HDiv-VIM solves the iron in the magnet's field
# (Br/MU0 converts a remanence Br [T] to the magnetization M [A/m]; MU0 = 4e-7*pi.)
# A spatially distributed given M uses vim.MagnetizationSource(pm_mesh, M_given),
# whose independent HDiv space supplies a native C++ field CF to the iron solve.


def magnet_box(center, size, M):
    """Axis-aligned uniform PERMANENT-MAGNET box (fixed magnetization ``M`` in A/m) -- analytic, no
    mesh/solve.  center = (cx,cy,cz), size = (sx,sy,sz) in METERS (the box spans ``center +- size/2``).
    Returns a single ``rad.ObjRecMag`` (exact surface-current field).  Compose with soft_iron_box in a
    container so the magnet drives the iron as a fixed-M source (see the module note)."""
    import radia as rad
    c = [float(v) for v in center]
    s = [float(v) for v in size]
    m = [float(v) for v in M]
    if len(c) != 3 or len(s) != 3 or len(m) != 3:
        raise ValueError("magnet_box: center, size, M must each be length-3 (got %r, %r, %r)" % (center, size, M))
    if min(s) <= 0.0:
        raise ValueError("magnet_box: size must be positive in every dimension (got %r)" % (size,))
    return rad.ObjRecMag(c, s, m)


def magnet_hex(vertices, M):
    """General (8-vertex) uniform PERMANENT-MAGNET hexahedron (fixed ``M`` in A/m) -- analytic, no
    mesh/solve.  vertices = 8 corners in the ``rad.ObjHexahedron`` CHEXA order, in METERS.  Returns a
    single ``rad.ObjHexahedron`` (exact surface-charge field)."""
    import radia as rad
    V = np.asarray(vertices, float)
    if V.shape != (8, 3):
        raise ValueError("magnet_hex: vertices must be 8x3 (CHEXA order); got shape %r" % (V.shape,))
    m = [float(v) for v in M]
    if len(m) != 3:
        raise ValueError("magnet_hex: M must be length-3 (A/m); got %r" % (M,))
    return rad.ObjHexahedron([[float(c) for c in v] for v in V], m)
