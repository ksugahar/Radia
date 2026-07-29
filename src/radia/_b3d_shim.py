"""build123d-compatible API shim built directly on pythonocc (OCP).

Why this module exists
======================
build123d's ``__init__.py`` unconditionally imports ``objects_curve``
which does ``import sympy`` for the ``TangentArc`` class.  That single
import pulls in 419 sympy submodules and costs 5-8 s on a cold Python
process.  build123d also wraps every OCP ``gp_Pnt`` / ``TopoDS_*`` in
its own ``Vector`` / ``Face`` / ``Edge`` Python class, which adds ~10 s
of Vector iteration overhead on the PEEC pipeline.

The PEEC pipeline only uses ~14 names from build123d (Vector, Plane,
import_step, section, GeomType, Compound, Solid/Face/Edge/Wire
accessors, BBox, PositionMode, export_brep).  Re-implementing those
14 names directly against OCP (which build123d itself wraps) skips
both the sympy import and the Vector wrapper layer, dropping the
``filaments_from_step`` cold-run cost from ~22 s to ~6 s on the
3-turn pancake sample.

Drop-in replacement
===================
The shim exposes the same names with the same signatures and the
same access patterns as build123d:

    from radia._b3d_shim import (Vector, Plane, import_step, section,
                                   GeomType, Compound, PositionMode,
                                   Face, Edge, export_brep)

so caller code does not need any other edits beyond switching the
import line.  See ``src/radia/coil_from_cad.py`` for the canonical
use site.

Scope
=====
This is a READ-ONLY CAD inspection shim.  build123d's authoring
features (BuildPart, Sketch, Wire constructors, etc.) are NOT
implemented — for STEP authoring continue to use build123d directly.
``coil_profile_occ.py`` is the only other ``from build123d import ...``
call site in radia, and it does authoring (Circle / Rectangle / loft);
that module stays on build123d.

Coordinates throughout are CAD units (typically mm); callers scale to
meters via ``cad_units_per_meter``.

Reference
---------
* OCP API: ``C:/Program Files/Python312/Lib/site-packages/OCP/`` (pythonocc-core 7.8.1.1)
* build123d source (for API parity): ``C:/Program Files/Python312/Lib/site-packages/build123d/``
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, IntEnum
from typing import List, Optional, Tuple

# OCP imports are deferred to function/method bodies so that
# ``import radia._b3d_shim`` itself stays under 0.1 s.  Importing all
# of OCP at module level adds ~0.6 s.  See coil_from_cad.py top-of-
# function ``from radia._b3d_shim import ...`` pattern.


# ============================================================
# GeomType enum (parity with build123d.build_enums.GeomType)
# ============================================================

class GeomType(Enum):
    """Surface / curve geometry type, parity with build123d.

    build123d's ``GeomType`` is a string-valued ``Enum``.  We mirror
    that and additionally override ``__eq__`` / ``__hash__`` so a
    shim ``GeomType.X`` compares equal to build123d's ``GeomType.X``
    by ``.name``.  This is essential because ``coil_topology.py``
    and ``coil_from_cad.py`` may receive Faces from EITHER the shim
    or build123d directly (tests typically use ``build123d.import_step``
    which returns build123d Face objects whose ``.geom_type`` is
    build123d's GeomType class).
    """
    PLANE = "PLANE"
    CYLINDER = "CYLINDER"
    CONE = "CONE"
    SPHERE = "SPHERE"
    TORUS = "TORUS"
    BEZIER = "BEZIER"
    BSPLINE = "BSPLINE"
    REVOLUTION = "REVOLUTION"
    EXTRUSION = "EXTRUSION"
    OFFSET = "OFFSET"
    OTHER = "OTHER"
    LINE = "LINE"
    CIRCLE = "CIRCLE"
    ELLIPSE = "ELLIPSE"
    HYPERBOLA = "HYPERBOLA"
    PARABOLA = "PARABOLA"

    def __eq__(self, other):
        # Cross-Enum-class equality by name.  Allows tests that load
        # via ``build123d.import_step`` (returning build123d Face
        # whose ``.geom_type`` is build123d's GeomType enum) to
        # compare against the shim's enum: build123d_GeomType.PLANE
        # == shim_GeomType.PLANE iff names match.
        if hasattr(other, "name") and isinstance(other.name, str):
            return self.name == other.name
        return NotImplemented

    def __hash__(self):
        return hash(self.name)


def _surface_geom_type(face_or_topods) -> GeomType:
    """Map an OCP ``TopoDS_Face`` to ``GeomType``."""
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.GeomAbs import (
        GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_Cone,
        GeomAbs_Sphere, GeomAbs_Torus, GeomAbs_BezierSurface,
        GeomAbs_BSplineSurface, GeomAbs_SurfaceOfRevolution,
        GeomAbs_SurfaceOfExtrusion, GeomAbs_OffsetSurface,
        GeomAbs_OtherSurface,
    )
    _MAP = {
        GeomAbs_Plane: GeomType.PLANE,
        GeomAbs_Cylinder: GeomType.CYLINDER,
        GeomAbs_Cone: GeomType.CONE,
        GeomAbs_Sphere: GeomType.SPHERE,
        GeomAbs_Torus: GeomType.TORUS,
        GeomAbs_BezierSurface: GeomType.BEZIER,
        GeomAbs_BSplineSurface: GeomType.BSPLINE,
        GeomAbs_SurfaceOfRevolution: GeomType.REVOLUTION,
        GeomAbs_SurfaceOfExtrusion: GeomType.EXTRUSION,
        GeomAbs_OffsetSurface: GeomType.OFFSET,
        GeomAbs_OtherSurface: GeomType.OTHER,
    }
    return _MAP.get(BRepAdaptor_Surface(face_or_topods).GetType(),
                     GeomType.OTHER)


def _curve_geom_type(edge_topods) -> GeomType:
    """Map an OCP ``TopoDS_Edge`` to ``GeomType``."""
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    from OCP.GeomAbs import (
        GeomAbs_Line, GeomAbs_Circle, GeomAbs_Ellipse,
        GeomAbs_Hyperbola, GeomAbs_Parabola,
        GeomAbs_BezierCurve, GeomAbs_BSplineCurve,
        GeomAbs_OtherCurve,
    )
    _MAP = {
        GeomAbs_Line: GeomType.LINE,
        GeomAbs_Circle: GeomType.CIRCLE,
        GeomAbs_Ellipse: GeomType.ELLIPSE,
        GeomAbs_Hyperbola: GeomType.HYPERBOLA,
        GeomAbs_Parabola: GeomType.PARABOLA,
        GeomAbs_BezierCurve: GeomType.BEZIER,
        GeomAbs_BSplineCurve: GeomType.BSPLINE,
        GeomAbs_OtherCurve: GeomType.OTHER,
    }
    return _MAP.get(BRepAdaptor_Curve(edge_topods).GetType(),
                     GeomType.OTHER)


# ============================================================
# PositionMode enum (parity with build123d.build_enums.PositionMode)
# ============================================================

class PositionMode(IntEnum):
    LENGTH = 0
    PARAMETER = 1


# ============================================================
# Vector — wrapping gp_Pnt / gp_Vec / gp_Dir uniformly
# ============================================================

class Vector:
    """3D vector / point, parity with build123d.geometry.Vector.

    build123d's ``Vector`` is used for both points and direction
    vectors interchangeably.  We follow that convention.  Internally
    we store three floats; OCP ``gp_*`` objects are constructed on
    demand via :meth:`to_pnt` / :meth:`to_dir` / :meth:`to_vec`.
    """
    __slots__ = ("X", "Y", "Z")

    def __init__(self, *args, **kwargs):
        # build123d's Vector(x, y, z) is the primary form used in
        # coil_from_cad.py.  Also support Vector(pnt) and
        # Vector(other_vector) and the (X=, Y=, Z=) keyword form.
        if len(args) == 3:
            self.X = float(args[0])
            self.Y = float(args[1])
            self.Z = float(args[2])
        elif len(args) == 1:
            other = args[0]
            if hasattr(other, "X") and not callable(other.X):
                self.X = float(other.X)
                self.Y = float(other.Y)
                self.Z = float(other.Z)
            elif hasattr(other, "X") and callable(other.X):  # gp_Pnt.X()
                self.X = float(other.X())
                self.Y = float(other.Y())
                self.Z = float(other.Z())
            else:
                t = tuple(other)
                self.X = float(t[0])
                self.Y = float(t[1])
                self.Z = float(t[2])
        elif len(args) == 0 and kwargs:
            self.X = float(kwargs.get("X", 0.0))
            self.Y = float(kwargs.get("Y", 0.0))
            self.Z = float(kwargs.get("Z", 0.0))
        else:
            self.X = 0.0
            self.Y = 0.0
            self.Z = 0.0

    # ---- arithmetic (build123d uses .length on differences) ----

    def __sub__(self, other: "Vector") -> "Vector":
        return Vector(self.X - other.X, self.Y - other.Y, self.Z - other.Z)

    def __add__(self, other: "Vector") -> "Vector":
        return Vector(self.X + other.X, self.Y + other.Y, self.Z + other.Z)

    def __mul__(self, k: float) -> "Vector":
        return Vector(self.X * k, self.Y * k, self.Z * k)

    __rmul__ = __mul__

    def __neg__(self) -> "Vector":
        return Vector(-self.X, -self.Y, -self.Z)

    def __iter__(self):
        return iter((self.X, self.Y, self.Z))

    def __repr__(self):
        return f"Vector({self.X:.6g}, {self.Y:.6g}, {self.Z:.6g})"

    @property
    def length(self) -> float:
        return (self.X ** 2 + self.Y ** 2 + self.Z ** 2) ** 0.5

    def dot(self, other: "Vector") -> float:
        return self.X * other.X + self.Y * other.Y + self.Z * other.Z

    def cross(self, other: "Vector") -> "Vector":
        return Vector(
            self.Y * other.Z - self.Z * other.Y,
            self.Z * other.X - self.X * other.Z,
            self.X * other.Y - self.Y * other.X,
        )

    # ---- OCP bridges ----

    def to_pnt(self):
        from OCP.gp import gp_Pnt
        return gp_Pnt(self.X, self.Y, self.Z)

    def to_dir(self):
        from OCP.gp import gp_Dir
        return gp_Dir(self.X, self.Y, self.Z)

    def to_vec(self):
        from OCP.gp import gp_Vec
        return gp_Vec(self.X, self.Y, self.Z)


def _pnt_to_vector(pnt) -> Vector:
    """``gp_Pnt`` (or ``gp_Vec``/``gp_Dir``) → :class:`Vector`."""
    return Vector(float(pnt.X()), float(pnt.Y()), float(pnt.Z()))


# ============================================================
# BBox
# ============================================================

@dataclass
class BBox:
    min: Vector
    max: Vector


# ============================================================
# Plane — keyword constructor (origin=, z_dir=)
# ============================================================

class Plane:
    """Plane defined by origin + z-axis direction.

    Parity with build123d.geometry.Plane(origin=Vector, z_dir=Vector).
    Internally builds an OCP ``gp_Pln`` via ``gp_Ax3`` lazily.
    """
    __slots__ = ("origin", "z_dir")

    def __init__(self, origin: Vector, z_dir: Vector):
        self.origin = origin
        self.z_dir = z_dir

    def to_gp_pln(self):
        from OCP.gp import gp_Pln, gp_Ax3
        ax3 = gp_Ax3(self.origin.to_pnt(), self.z_dir.to_dir())
        return gp_Pln(ax3)


# ============================================================
# Shape wrappers (Face / Edge / Wire / Solid / Compound)
# ============================================================

class _ShapeBase:
    """Base class holding the underlying OCP ``TopoDS_Shape``.

    ``self.wrapped`` is the OCP object — exposed to mimic
    build123d's escape-hatch property of the same name (which the
    caller code uses to bridge into raw OCP operations).

    No __slots__: subclasses (Solid) attach lazy caches
    (_bb_cache, _faces_cache, ...) dynamically.  __dict__ overhead
    is negligible compared to the OCP TopoDS_Shape held internally.
    """

    def __init__(self, topods):
        self.wrapped = topods
        self.label = ""  # build123d exposes .label even when empty


class Edge(_ShapeBase):
    """``TopoDS_Edge`` wrapper, parity with build123d.topology.Edge."""

    @property
    def geom_type(self) -> GeomType:
        return _curve_geom_type(self.wrapped)

    @property
    def length(self) -> float:
        from OCP.BRepAdaptor import BRepAdaptor_Curve
        from OCP.GCPnts import GCPnts_AbscissaPoint
        adaptor = BRepAdaptor_Curve(self.wrapped)
        return float(GCPnts_AbscissaPoint.Length_s(adaptor))

    def center(self) -> Vector:
        # build123d uses mass-properties centroid.  For a curve
        # the BRepGProp_LinearProperties gives mass = length and
        # centre of mass = arc-length midpoint analogue.
        from OCP.BRepGProp import BRepGProp
        from OCP.GProp import GProp_GProps
        gp = GProp_GProps()
        BRepGProp.LinearProperties_s(self.wrapped, gp)
        return _pnt_to_vector(gp.CentreOfMass())

    def start_point(self) -> Vector:
        from OCP.BRepAdaptor import BRepAdaptor_Curve
        ad = BRepAdaptor_Curve(self.wrapped)
        return _pnt_to_vector(ad.Value(ad.FirstParameter()))

    def end_point(self) -> Vector:
        from OCP.BRepAdaptor import BRepAdaptor_Curve
        ad = BRepAdaptor_Curve(self.wrapped)
        return _pnt_to_vector(ad.Value(ad.LastParameter()))

    def position_at(self, distance: float,
                    position_mode: PositionMode = PositionMode.PARAMETER) -> Vector:
        """Position along the curve.

        Parity with build123d.topology.Mixin1D.position_at:
        * PositionMode.PARAMETER (default) -- distance in (0, 1]
          is treated as a NORMALISED parameter linearly mapped to
          [FirstParameter, LastParameter]; other values are passed
          through as the raw curve parameter.
        * PositionMode.LENGTH -- distance is treated as arc length
          from the curve start.
        """
        from OCP.BRepAdaptor import BRepAdaptor_Curve
        ad = BRepAdaptor_Curve(self.wrapped)
        if position_mode == PositionMode.LENGTH:
            from OCP.GCPnts import GCPnts_AbscissaPoint
            total = GCPnts_AbscissaPoint.Length_s(ad)
            target = (float(distance) * total
                      if 0.0 <= float(distance) <= 1.0
                      else float(distance))
            ap = GCPnts_AbscissaPoint(ad, target, ad.FirstParameter())
            return _pnt_to_vector(ad.Value(ap.Parameter()))
        # PARAMETER mode (default, matches build123d)
        u_first = ad.FirstParameter()
        u_last = ad.LastParameter()
        d = float(distance)
        if 0.0 < d <= 1.0:
            u = u_first + (u_last - u_first) * d
        else:
            u = u_first + d
        return _pnt_to_vector(ad.Value(u))

    def __matmul__(self, t: float) -> Vector:
        """``edge @ t`` -> ``edge.position_at(t)`` (build123d shorthand).

        Default is PositionMode.PARAMETER, normalised t in (0, 1] maps
        linearly to [FirstParameter, LastParameter].  This is the form
        used in ``coil_from_cad.py`` for spine sampling:
        ``spine @ 0.5`` -> midpoint, ``spine @ t`` for t in [0, 1] ->
        uniform parametric stations along the spine.
        """
        return self.position_at(t, PositionMode.PARAMETER)

    @property
    def is_closed(self) -> bool:
        from OCP.BRep import BRep_Tool
        # An edge is closed iff its first and last vertices coincide
        try:
            from OCP.TopoDS import TopoDS_Vertex, TopoDS
            from OCP.TopAbs import TopAbs_VERTEX, TopAbs_FORWARD, TopAbs_REVERSED
            from OCP.TopExp import TopExp
            v1 = TopExp.FirstVertex_s(self.wrapped, True)
            v2 = TopExp.LastVertex_s(self.wrapped, True)
            return v1.IsSame(v2)
        except Exception:
            return False

    @property
    def radius(self) -> float:
        """Circle radius (only meaningful when ``geom_type == CIRCLE``)."""
        from OCP.BRepAdaptor import BRepAdaptor_Curve
        ad = BRepAdaptor_Curve(self.wrapped)
        return float(ad.Circle().Radius())

    def Location(self) -> Vector:
        """For a circular edge, the centre of the circle."""
        from OCP.BRepAdaptor import BRepAdaptor_Curve
        ad = BRepAdaptor_Curve(self.wrapped)
        return _pnt_to_vector(ad.Circle().Location())

    def Axis(self):
        """For a circular edge, the axis of the circle.

        Returns an :class:`_AxisProxy` with ``.Direction()``
        returning a :class:`Vector` (matches the build123d/OCP
        accessor pattern used in coil_from_cad.py).
        """
        from OCP.BRepAdaptor import BRepAdaptor_Curve
        ad = BRepAdaptor_Curve(self.wrapped)
        return _AxisProxy(ad.Circle().Axis())


class _AxisProxy:
    """Wraps gp_Ax1 / gp_Ax2 with a .Direction() → Vector and
    .Location() → Vector accessor, mimicking what coil_from_cad.py
    expects from build123d's circle/cylinder Axis() return."""

    __slots__ = ("_ax",)

    def __init__(self, gp_axis):
        self._ax = gp_axis

    def Direction(self) -> Vector:
        d = self._ax.Direction()
        return Vector(d.X(), d.Y(), d.Z())

    def Location(self) -> Vector:
        p = self._ax.Location()
        return Vector(p.X(), p.Y(), p.Z())


def _unique_subshapes(root_topods, abs_type, wrap_cls, downcast):
    """Iterate over unique sub-shapes of the given TopAbs_* type.

    OCP's ``TopExp_Explorer`` visits each EDGE once per containing
    face -- shared edges between adjacent faces appear multiple times.
    build123d's ``.edges()`` accessor dedupes; we mirror that via
    ``TopExp.MapShapes_s`` which builds an indexed map keyed by
    geometric identity (each TopoDS_Shape stored once)."""
    from OCP.TopExp import TopExp
    from OCP.TopTools import TopTools_IndexedMapOfShape
    m = TopTools_IndexedMapOfShape()
    TopExp.MapShapes_s(root_topods, abs_type, m)
    out = []
    for i in range(1, m.Extent() + 1):
        out.append(wrap_cls(downcast(m.FindKey(i))))
    return out


class Wire(_ShapeBase):
    """``TopoDS_Wire`` wrapper."""

    def edges(self) -> List[Edge]:
        from OCP.TopAbs import TopAbs_EDGE
        from OCP.TopoDS import TopoDS
        return _unique_subshapes(self.wrapped, TopAbs_EDGE, Edge, TopoDS.Edge_s)


class Face(_ShapeBase):
    """``TopoDS_Face`` wrapper, parity with build123d.topology.Face."""

    @property
    def geom_type(self) -> GeomType:
        return _surface_geom_type(self.wrapped)

    @property
    def area(self) -> float:
        from OCP.BRepGProp import BRepGProp
        from OCP.GProp import GProp_GProps
        gp = GProp_GProps()
        BRepGProp.SurfaceProperties_s(self.wrapped, gp)
        return float(gp.Mass())

    def center(self) -> Vector:
        from OCP.BRepGProp import BRepGProp
        from OCP.GProp import GProp_GProps
        gp = GProp_GProps()
        BRepGProp.SurfaceProperties_s(self.wrapped, gp)
        return _pnt_to_vector(gp.CentreOfMass())

    def normal_at(self, point: Optional[Vector] = None) -> Vector:
        """Outward normal at a point on the face.

        If ``point`` is None, evaluate at the centre of the (u, v)
        parameter range — matches build123d's default behaviour.
        """
        from OCP.BRepAdaptor import BRepAdaptor_Surface
        from OCP.gp import gp_Pnt, gp_Vec
        from OCP.GeomLProp import GeomLProp_SLProps
        from OCP.BRep import BRep_Tool
        from OCP.BRepTools import BRepTools
        from OCP.TopAbs import TopAbs_REVERSED

        surf = BRep_Tool.Surface_s(self.wrapped)
        umin, umax, vmin, vmax = BRepTools.UVBounds_s(self.wrapped)
        if point is None:
            u = 0.5 * (umin + umax)
            v = 0.5 * (vmin + vmax)
        else:
            # Project the requested point to (u, v).
            from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
            proj = GeomAPI_ProjectPointOnSurf(point.to_pnt(), surf)
            u, v = proj.LowerDistanceParameters()

        props = GeomLProp_SLProps(surf, u, v, 1, 1e-9)
        if not props.IsNormalDefined():
            return Vector(0.0, 0.0, 0.0)
        n = props.Normal()
        nv = Vector(n.X(), n.Y(), n.Z())
        # Respect the face's orientation: a REVERSED face flips the normal
        if self.wrapped.Orientation() == TopAbs_REVERSED:
            nv = -nv
        return nv

    def outer_wire(self) -> Wire:
        from OCP.BRepTools import BRepTools
        return Wire(BRepTools.OuterWire_s(self.wrapped))

    def bounding_box(self) -> BBox:
        return _bounding_box(self.wrapped)

    # build123d exposes Axis() on TORUS / CYLINDER / REVOLUTION faces
    # to access the underlying axis.  coil_from_cad.py uses this for
    # cylinder lead detection (line 484-489) and revolution faces.
    def Axis(self):
        from OCP.BRepAdaptor import BRepAdaptor_Surface
        from OCP.GeomAbs import (
            GeomAbs_Cylinder, GeomAbs_Cone, GeomAbs_Sphere,
            GeomAbs_Torus, GeomAbs_SurfaceOfRevolution,
        )
        ad = BRepAdaptor_Surface(self.wrapped)
        t = ad.GetType()
        if t == GeomAbs_Cylinder:
            return _AxisProxy(ad.Cylinder().Axis())
        if t == GeomAbs_Cone:
            return _AxisProxy(ad.Cone().Axis())
        if t == GeomAbs_Sphere:
            return _AxisProxy(ad.Sphere().Position().Axis())
        if t == GeomAbs_Torus:
            return _AxisProxy(ad.Torus().Axis())
        if t == GeomAbs_SurfaceOfRevolution:
            return _AxisProxy(ad.AxeOfRevolution())
        raise TypeError(f"Face.Axis() not defined for geom_type={t}")


def _bounding_box(topods) -> BBox:
    """Optimal (tight) bbox via Newton-iteration on the geometry.

    Matches build123d's default ``optimal=True``.  The looser
    triangulation-based ``BRepBndLib.Add_s`` returns coordinates
    inflated by 1.5-3 mm for curved geometries (the torus radius
    on the IH samples), which breaks the
    ``_check_filaments_cover_solid_bbox`` test downstream.
    """
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    bb = Bnd_Box()
    BRepBndLib.AddOptimal_s(topods, bb, True, True)
    xmin, ymin, zmin, xmax, ymax, zmax = bb.Get()
    return BBox(min=Vector(xmin, ymin, zmin), max=Vector(xmax, ymax, zmax))


class Solid(_ShapeBase):
    """``TopoDS_Solid`` (or generic shape) wrapper.

    ``import_step`` returns ``Solid`` for a single-solid STEP and
    ``Compound`` for a multi-solid one — caller code in
    coil_from_cad.py branches via ``isinstance(solid, Compound)``
    to extract member solids.

    Lazy caches on first access (the underlying TopoDS_Shape never
    changes after import_step / section, so these are safe to cache
    for the wrapper's lifetime):

    * ``bounding_box()``: optimal-bbox is the dominant ~0.5 s call on
      large coil solids and was being re-computed 6 × in the original
      coil_from_cad.py path because multiple validation predicates
      ask for it independently.  Cached on first call.
    * ``faces()`` / ``edges()`` / ``solids()``: TopExp deduplication
      is fast (~0.01 s) but is also called from several predicates.
    """

    def faces(self) -> List[Face]:
        if getattr(self, "_faces_cache", None) is not None:
            return self._faces_cache
        from OCP.TopAbs import TopAbs_FACE
        from OCP.TopoDS import TopoDS
        self._faces_cache = _unique_subshapes(
            self.wrapped, TopAbs_FACE, Face, TopoDS.Face_s)
        return self._faces_cache

    def edges(self) -> List[Edge]:
        if getattr(self, "_edges_cache", None) is not None:
            return self._edges_cache
        from OCP.TopAbs import TopAbs_EDGE
        from OCP.TopoDS import TopoDS
        self._edges_cache = _unique_subshapes(
            self.wrapped, TopAbs_EDGE, Edge, TopoDS.Edge_s)
        return self._edges_cache

    def solids(self) -> List["Solid"]:
        if getattr(self, "_solids_cache", None) is not None:
            return self._solids_cache
        from OCP.TopAbs import TopAbs_SOLID
        from OCP.TopoDS import TopoDS
        self._solids_cache = _unique_subshapes(
            self.wrapped, TopAbs_SOLID, Solid, TopoDS.Solid_s)
        return self._solids_cache

    def bounding_box(self) -> BBox:
        if getattr(self, "_bb_cache", None) is not None:
            return self._bb_cache
        self._bb_cache = _bounding_box(self.wrapped)
        return self._bb_cache


class Compound(Solid):
    """``TopoDS_Compound`` wrapper.

    Inherits .faces() / .edges() / .solids() from Solid.  Callers
    use ``isinstance(shape, Compound)`` to decide whether to call
    .solids() and pick one — see coil_from_cad.py:1746.
    """

    @property
    def children(self):
        """Return the immediate child shapes (mimicking build123d's
        Part/Assembly tree).  For an auto-created Compound this is
        typically an empty tuple — coil_from_cad.py:2833 has the
        comment noting this case."""
        # Walk the first-level TopExp iterator (Solid + Shell + Face)
        from OCP.TopoDS import TopoDS_Iterator
        out = []
        it = TopoDS_Iterator(self.wrapped)
        while it.More():
            out.append(_wrap_topods(it.Value()))
            it.Next()
        return out


def _wrap_topods(topods):
    """Wrap an OCP ``TopoDS_Shape`` in the most specific shim class."""
    from OCP.TopAbs import (TopAbs_COMPOUND, TopAbs_SOLID,
                              TopAbs_FACE, TopAbs_WIRE, TopAbs_EDGE)
    from OCP.TopoDS import TopoDS
    st = topods.ShapeType()
    if st == TopAbs_COMPOUND:
        return Compound(TopoDS.Compound_s(topods))
    if st == TopAbs_SOLID:
        return Solid(TopoDS.Solid_s(topods))
    if st == TopAbs_FACE:
        return Face(TopoDS.Face_s(topods))
    if st == TopAbs_WIRE:
        return Wire(TopoDS.Wire_s(topods))
    if st == TopAbs_EDGE:
        return Edge(TopoDS.Edge_s(topods))
    # Fall back to generic Solid wrapper — gives .faces() etc.
    return Solid(topods)


# ============================================================
# import_step / section / export_brep
# ============================================================

def import_step(file_path: str):
    """Load a STEP file → :class:`Solid` (or :class:`Compound`).

    Parity with build123d.importers.import_step.  Uses OCP
    ``STEPControl_Reader`` directly — no XCAF label support yet
    (coil_from_cad.py does not depend on labels; the user-facing
    coil import in coil_profile_occ.py still uses build123d's
    labelled importer if labels are needed).
    """
    from OCP.STEPControl import STEPControl_Reader
    from OCP.IFSelect import IFSelect_RetDone
    reader = STEPControl_Reader()
    status = reader.ReadFile(str(file_path))
    if status != IFSelect_RetDone:
        raise IOError(f"STEP read failed (status={status}): {file_path}")
    reader.TransferRoots()
    shape = reader.OneShape()
    return _wrap_topods(shape)


def section(solid, section_by: Plane):
    """Section a solid with a plane → :class:`Compound` of the cut FACES.

    Parity with build123d.operations_part.section(obj=solid,
    section_by=Plane(...)).  build123d builds an oversized planar
    ``Face`` covering the bounding box and computes the BOOLEAN
    INTERSECTION of the solid with that face -- the result is the
    set of cross-section FACES on the cutting plane, NOT just the
    cut edges.  OCP's lower-level ``BRepAlgoAPI_Section`` only
    returns the curves; we use ``BRepAlgoAPI_Common`` (boolean
    intersect) against a bounded planar face built large enough to
    cover the solid's bounding box.

    ``solid`` may be a shim :class:`Solid`/:class:`Compound` (read
    .wrapped) or a raw OCP TopoDS_*.
    """
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace

    src = solid.wrapped if hasattr(solid, "wrapped") else solid

    # Size the planar face from the solid's bounding box (with margin).
    bb = _bounding_box(src)
    diag = (
        (bb.max.X - bb.min.X) ** 2
        + (bb.max.Y - bb.min.Y) ** 2
        + (bb.max.Z - bb.min.Z) ** 2
    ) ** 0.5
    # 4 x diagonal is generous; build123d uses 2 * max(|min|, |max|) + diag.
    half = max(diag * 2.0, 1.0)

    pln = section_by.to_gp_pln()
    # BRepBuilderAPI_MakeFace(gp_Pln, umin, umax, vmin, vmax) makes a
    # bounded planar face in the plane's UV (X,Y of the plane).
    face_builder = BRepBuilderAPI_MakeFace(pln, -half, half, -half, half)
    if not face_builder.IsDone():
        return None
    plane_face = face_builder.Face()

    common = BRepAlgoAPI_Common(src, plane_face)
    common.Build()
    if not common.IsDone():
        return None
    return _wrap_topods(common.Shape())


def export_brep(shape, path: str) -> None:
    """Write a shape to a BRep file."""
    from OCP.BRepTools import BRepTools
    src = shape.wrapped if hasattr(shape, "wrapped") else shape
    BRepTools.Write_s(src, str(path))


# ============================================================
# section_by= alias for callers that pass it positionally
# ============================================================
# build123d's section() takes ``section_by`` as the second kwarg.
# coil_from_cad.py mostly uses ``section(solid, plane)`` positionally
# — our :func:`section` already supports that pattern.

__all__ = [
    "GeomType",
    "PositionMode",
    "Vector",
    "Plane",
    "BBox",
    "Edge",
    "Wire",
    "Face",
    "Solid",
    "Compound",
    "import_step",
    "section",
    "export_brep",
]
