"""Curved Netgen ``.vol`` export from a Cubit model in an external process.

This module replaces the native Cubit plugin exporter (``export netgen``).
One CPython process imports Cubit's Python API and Netgen, rebuilds the
linear Cubit mesh as a ``netgen.meshing.Mesh`` with face and edge
descriptors and surface UV parameters, attaches a ``NetgenGeometry``
subclass that projects onto the Cubit (ACIS) geometry, curves the mesh with
``Mesh.Curve(order)`` (plus ``Mesh.GeometricRefit()`` for order >= 3) and
writes the ``.vol`` and its companion ``.vol.json``.

Run it as a subprocess::

    python -I -B -m radia_mcp.cubit.vol_export --cubit-bin <Cubit bin dir> \
        --journal model.jou --out model.vol --order 2 --result status.json

Requirements: a Netgen runtime whose ``NetgenGeometry`` trampoline passes
points and geometry info by reference and forwards ``edgenr`` to the edge
callbacks, and that provides ``Mesh.GeometricRefit``.  The worker checks for
``GeometricRefit`` and fails loudly when it is missing.

Projection route.  Cubit's Python ``Surface`` has no UV-seeded closest-point
query (the native exporter used ACIS ``closest_point_uv_guess``), so every
surface projection uses ``Surface.closest_point_trimmed`` followed by
``Surface.u_v_from_position``.  This is the native exporter's own fallback
route; it is recorded as ``projection_route`` in the diagnostics.

All lengths are in the Cubit model's units.  No unit conversion is applied.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from radia_mcp.cubit import periodic as _periodic

# Thresholds in model length units, identical to the native exporter.
PROJECT_REJECT_THRESHOLD = 1.5
EDGE_REJECT_THRESHOLD = 1.5
PROJECTION_ROUTE = "closest_point_trimmed"

_VOLUME_KINDS = (("hex", 8), ("tet", 4), ("wedge", 6), ("pyramid", 5))
_VOLUME_FACES = {
    "tet": ((0, 1, 2), (0, 1, 3), (1, 2, 3), (2, 0, 3)),
    "hex": ((0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5),
            (2, 3, 7, 6), (3, 0, 4, 7)),
    "wedge": ((0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)),
    "pyramid": ((0, 1, 2, 3), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)),
}


class ExportError(RuntimeError):
    """The model cannot be exported faithfully; nothing is substituted."""


# ---------------------------------------------------------------------------
# Cubit runtime
# ---------------------------------------------------------------------------

def import_cubit(cubit_bin: str | os.PathLike[str]):
    """Import and initialise Cubit's Python API in batch/nographics mode.

    ``-noinitfile`` keeps the user's ``$HOME/.cubit`` and ``$CWD/.cubit`` (and
    anything they play) out of the export.
    """
    bin_dir = Path(cubit_bin)
    if not (bin_dir / "cubit").is_dir():
        raise ExportError(f"Cubit Python package not found under {bin_dir}")
    if str(bin_dir) not in sys.path:
        sys.path.append(str(bin_dir))
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(str(bin_dir))
    import cubit  # noqa: PLC0415  (Cubit is an optional, process-local runtime)

    cubit.init(["cubit", "-nographics", "-nojournal", "-batch", "-noinitfile"])
    return cubit


# ---------------------------------------------------------------------------
# Mesh extraction (port of MeshData.cpp)
# ---------------------------------------------------------------------------

@dataclass
class VolumeElement:
    block: int
    kind: str
    conn: tuple[int, ...]


@dataclass
class Sideset:
    id: int
    name: str
    faces: list[tuple[int, ...]] = field(default_factory=list)


@dataclass
class Nodeset:
    id: int
    name: str
    node_ids: list[int]


@dataclass
class CubitMesh:
    node_ids: list[int]
    coords: dict[int, tuple[float, float, float]]
    block_ids: list[int]
    elements: list[VolumeElement]
    sidesets: list[Sideset]
    nodesets: list[Nodeset]
    warnings: list[str]


def _connectivity(cubit, kind: str, eid: int, nv: int) -> tuple[int, ...]:
    conn = tuple(cubit.get_connectivity(kind, eid))
    if len(conn) < nv:
        raise ExportError(f"{kind} {eid} has {len(conn)} nodes, expected {nv}")
    return conn[:nv]


def extract_mesh(cubit) -> CubitMesh:
    """Collect nodes, block elements, sidesets and nodesets from Cubit."""
    node_ids = sorted(cubit.parse_cubit_list("node", "all"))
    if not node_ids:
        raise ExportError("No mesh nodes found. Nothing to export.")
    coords = {nid: tuple(cubit.get_nodal_coordinates(nid)) for nid in node_ids}

    warnings: list[str] = []
    block_ids: list[int] = []
    elements: list[VolumeElement] = []
    seen: set[tuple[int, ...]] = set()
    duplicates = 0
    for bid in sorted(cubit.parse_cubit_list("block", "all")):
        if bid not in block_ids:
            block_ids.append(bid)
        for kind, nv in _VOLUME_KINDS:
            # get_block_<kind>s lists only elements added directly; the block
            # query also resolves elements owned through block volumes.
            for eid in cubit.parse_cubit_list(kind, f"in block {bid}"):
                conn = _connectivity(cubit, kind, eid, nv)
                key = tuple(sorted(conn))
                if key in seen:
                    duplicates += 1
                    continue
                seen.add(key)
                elements.append(VolumeElement(bid, kind, conn))
    if duplicates:
        warnings.append(f"Skipped {duplicates} duplicate elements across blocks.")

    exported = {kind: 0 for kind, _ in _VOLUME_KINDS}
    for element in elements:
        exported[element.kind] += 1
    totals = {"tet": cubit.get_tet_count(), "hex": cubit.get_hex_count(),
              "wedge": cubit.get_wedge_count(), "pyramid": cubit.get_pyramid_count()}
    if not elements:
        raise ExportError(
            f"No block contains volume elements; the model has {sum(totals.values())} "
            "meshed volume elements outside every block. Add them to a block "
            "before exporting, e.g. 'block 1 add volume all'.")
    missing = [f"{totals[k] - exported[k]} of {totals[k]} {k}"
               for k in totals if totals[k] > exported[k]]
    if missing:
        warnings.append("Volume elements outside every block are not exported: "
                        + ", ".join(missing))

    sidesets: list[Sideset] = []
    for sid_id in sorted(cubit.parse_cubit_list("sideset", "all")):
        sideset = Sideset(sid_id, cubit.get_exodus_entity_name("sideset", sid_id))
        seen_faces: set[tuple[int, ...]] = set()

        def append(kind: str, fid: int, nv: int) -> None:
            conn = _connectivity(cubit, kind, fid, nv)
            key = tuple(sorted(conn))
            if key not in seen_faces:
                seen_faces.add(key)
                sideset.faces.append(conn)

        for surface in cubit.get_sideset_surfaces(sid_id):
            for tid in cubit.parse_cubit_list("tri", f"in surface {surface}"):
                append("tri", tid, 3)
            for qid in cubit.parse_cubit_list("face", f"in surface {surface}"):
                append("face", qid, 4)
        for tid in cubit.get_sideset_tris(sid_id):
            append("tri", tid, 3)
        for qid in cubit.get_sideset_quads(sid_id):
            append("face", qid, 4)
        sidesets.append(sideset)

    nodesets = [Nodeset(ns, cubit.get_exodus_entity_name("nodeset", ns),
                        list(cubit.parse_cubit_list("node", f"in nodeset {ns}")))
                for ns in sorted(cubit.parse_cubit_list("nodeset", "all"))]
    return CubitMesh(node_ids, coords, block_ids, elements, sidesets, nodesets, warnings)


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def _dist(a, b) -> float:
    dx, dy, dz = a[0] - b[0], a[1] - b[1], a[2] - b[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _periods(cubit, sid: int) -> tuple[float, float]:
    """Parameter periods (0 = not periodic) of a Cubit surface.

    Cubit's Python API reports only whether a surface is periodic, so each
    parameter direction is tested geometrically: it is periodic when the two
    ends of its parameter range map to the same points.
    """
    if not cubit.is_periodic("surface", sid):
        return 0.0, 0.0
    surface = cubit.surface(sid)
    (u0, u1), (v0, v1) = surface.get_param_range_U(), surface.get_param_range_V()
    lo, hi = surface.bounding_box()[0:3], surface.bounding_box()[3:6]
    tol = 1e-9 * max(1.0, _dist(lo, hi))
    samples = [i / 6.0 for i in range(7)]

    def closes(along_u: bool) -> bool:
        for t in samples:
            if along_u:
                a = surface.position_from_u_v(u0, v0 + t * (v1 - v0))
                b = surface.position_from_u_v(u1, v0 + t * (v1 - v0))
            else:
                a = surface.position_from_u_v(u0 + t * (u1 - u0), v0)
                b = surface.position_from_u_v(u0 + t * (u1 - u0), v1)
            if _dist(a, b) > tol:
                return False
        return True

    up = (u1 - u0) if u1 > u0 and closes(True) else 0.0
    vp = (v1 - v0) if v1 > v0 and closes(False) else 0.0
    if up == 0.0 and vp == 0.0:
        raise ExportError(f"Cubit reports surface {sid} periodic, but no parameter "
                          "direction closes on itself")
    return up, vp


def _has_pole(cubit, sid: int) -> bool:
    """True when one side of the parameter rectangle collapses to a point."""
    surface = cubit.surface(sid)
    (u0, u1), (v0, v1) = surface.get_param_range_U(), surface.get_param_range_V()
    lo, hi = surface.bounding_box()[0:3], surface.bounding_box()[3:6]
    tol = 1e-9 * max(1.0, _dist(lo, hi))
    samples = [i / 6.0 for i in range(7)]
    sides = [[surface.position_from_u_v(u, v0 + t * (v1 - v0)) for t in samples] for u in (u0, u1)]
    sides += [[surface.position_from_u_v(u0 + t * (u1 - u0), v) for t in samples] for v in (v0, v1)]
    return any(max(_dist(points[0], q) for q in points) <= tol for points in sides)


def _seam(u_ref: float, u: float, period: float) -> float:
    if period <= 0.0:
        return u
    while u - u_ref > 0.5 * period:
        u -= period
    while u - u_ref < -0.5 * period:
        u += period
    return u


# ---------------------------------------------------------------------------
# Netgen geometry backed by Cubit (port of callbackgeom.cpp + NetgenCurver.cpp)
# ---------------------------------------------------------------------------

def _make_geometry_class(ngm):
    class CubitGeometry(ngm.NetgenGeometry):
        """Netgen geometry callbacks that project onto Cubit surfaces and curves."""

        def __init__(self, cubit, fd_to_surface, edge_bindings, polar, bboxes, pad):
            super().__init__()
            self.cubit = cubit
            self.fd_to_surface = fd_to_surface
            self.edge_bindings = edge_bindings
            self.polar = polar
            self.bboxes = bboxes
            self.pad = pad
            # same bounds as lo[i] - pad and hi[i] + pad, computed once
            self.padded_bboxes = {sid: ([v - pad for v in lo], [v + pad for v in hi])
                                  for sid, (lo, hi) in bboxes.items()}
            self._surfaces = {}
            self._curves = {}
            self.diag = {"projection_route": PROJECTION_ROUTE,
                         "project_calls": 0, "project_rejects": 0,
                         "project_bbox_identity": 0, "polar_skips": 0,
                         "project_max_disp": 0.0, "project_max_disp_surface": -1,
                         "edge_calls": 0, "edge_rejects": 0, "edge_failures": 0,
                         "edge_max_disp": 0.0, "between_edge_calls": 0,
                         "between_edge_projection": 0, "between_edge_arc": 0,
                         "normal_calls": 0, "tangent_calls": 0}

        # cubit.surface()/cubit.curve() look the entity up on every call
        # (about 1 ms each on a meshed model), so keep one handle per entity.
        def _surface(self, sid):
            handle = self._surfaces.get(sid)
            if handle is None:
                handle = self._surfaces[sid] = self.cubit.surface(sid)
            return handle

        def _curve_object(self, cid):
            handle = self._curves.get(cid)
            if handle is None:
                handle = self._curves[cid] = self.cubit.curve(cid)
            return handle

        # -- surface projection -------------------------------------------
        def _project(self, surfnr, x):
            self.diag["project_calls"] += 1
            sid = self.fd_to_surface.get(surfnr)
            if sid is None:
                return x, None
            if sid in self.polar:
                self.diag["polar_skips"] += 1
                return x, None
            lo, hi = self.padded_bboxes[sid]
            x0, x1, x2 = x
            if (x0 < lo[0] or x0 > hi[0] or x1 < lo[1] or x1 > hi[1]
                    or x2 < lo[2] or x2 > hi[2]):
                self.diag["project_bbox_identity"] += 1
                return x, None
            surface = self._surfaces.get(sid) or self._surface(sid)
            closest = tuple(surface.closest_point_trimmed([x0, x1, x2]))
            disp = _dist(closest, x)
            if not math.isfinite(disp):
                raise ExportError(f"Non-finite projection on Cubit surface {sid}")
            if disp > self.diag["project_max_disp"]:
                self.diag["project_max_disp"] = disp
                self.diag["project_max_disp_surface"] = sid
            if disp > PROJECT_REJECT_THRESHOLD:
                self.diag["project_rejects"] += 1
                return x, tuple(surface.u_v_from_position(list(x)))
            return closest, tuple(surface.u_v_from_position(list(closest)))

        def _apply(self, surfnr, p, gi):
            x = (p[0], p[1], p[2])
            q, uv = self._project(surfnr, x)
            p[0], p[1], p[2] = q
            if gi is not None:
                gi.trignum = surfnr
                if uv is not None:
                    gi.u, gi.v = uv
            return True

        def ProjectPoint(self, surfnr, p):
            gi = ngm.PointGeomInfo()
            self._apply(surfnr, p, gi)
            return gi

        def ProjectPointGI(self, surfnr, p, gi):
            return self._apply(surfnr, p, gi)

        def CalcPointGeomInfo(self, surfnr, gi, p):
            _, uv = self._project(surfnr, (p[0], p[1], p[2]))
            gi.trignum = surfnr
            if uv is not None:
                gi.u, gi.v = uv
            return True

        def PointBetween(self, p1, p2, s, surfnr, gi1, gi2, newp, newgi):
            for i in range(3):
                newp[i] = p1[i] + s * (p2[i] - p1[i])
            if surfnr >= 0:
                newgi.u = gi1.u + s * (gi2.u - gi1.u)
                newgi.v = gi1.v + s * (gi2.v - gi1.v)
                self._apply(surfnr, newp, newgi)

        def GetNormal(self, surfnr, p, gi):
            self.diag["normal_calls"] += 1
            sid = self.fd_to_surface.get(surfnr)
            if sid is None:
                return ngm.Vec(0.0, 0.0, 1.0)
            n = self._surface(sid).normal_at([p[0], p[1], p[2]])
            return ngm.Vec(n[0], n[1], n[2])

        # -- curve projection ---------------------------------------------
        def _curve(self, edgenr, s1, s2):
            binding = self.edge_bindings.get(edgenr)
            if binding is None:
                raise ExportError(f"Unknown Netgen geometric edge {edgenr}")
            curve_id, f1, f2 = binding
            if {f1, f2} != {s1, s2}:
                raise ExportError(f"Surface mismatch for Netgen geometric edge {edgenr}: "
                                  f"({s1},{s2}) vs ({f1},{f2})")
            return curve_id, self._curve_object(curve_id)

        def ProjectPointEdge(self, s1, s2, p, gi, edgenr):
            self.diag["edge_calls"] += 1
            curve_id, curve = self._curve(edgenr, s1, s2)
            x = (p[0], p[1], p[2])
            closest = tuple(curve.closest_point_trimmed(list(x)))
            disp = _dist(closest, x)
            if not math.isfinite(disp):
                self.diag["edge_failures"] += 1
                raise ExportError(f"Non-finite projection for Netgen geometric edge {edgenr}")
            self.diag["edge_max_disp"] = max(self.diag["edge_max_disp"], disp)
            if disp > EDGE_REJECT_THRESHOLD:
                self.diag["edge_rejects"] += 1
                self.diag["edge_failures"] += 1
                raise ExportError("Curve projection displacement exceeds unchanged limit "
                                  f"for Cubit curve {curve_id} edge {edgenr}")
            p[0], p[1], p[2] = closest
            if gi is not None:
                _, uv = self._project(s1, closest)
                gi.dist = 0.0
                if uv is not None:
                    gi.u, gi.v = uv

        def PointBetweenEdge(self, p1, p2, s, s1, s2, a1, a2, newp, newgi, edgenr):
            self.diag["between_edge_calls"] += 1
            q = self._between_edge(edgenr, s1, s2, (p1[0], p1[1], p1[2]),
                                   (p2[0], p2[1], p2[2]), s)
            newp[0], newp[1], newp[2] = q
            if s1 >= 0:
                _, uv = self._project(s1, q)
                if uv is not None:
                    newgi.u, newgi.v = uv
            newgi.dist = a1.dist + s * (a2.dist - a1.dist)

        def _between_edge(self, edgenr, s1, s2, x1, x2, s):
            curve_id, curve = self._curve(edgenr, s1, s2)

            def fail(reason):
                self.diag["edge_failures"] += 1
                raise ExportError(f"{reason} for Cubit curve {curve_id} edge {edgenr}")

            length = curve.length()
            if not math.isfinite(length) or length <= 0.0:
                fail("Invalid curve length")
            u1 = curve.u_from_position(list(x1))
            u2 = curve.u_from_position(list(x2))
            periodic = curve.is_periodic()
            period = curve.end_param() - curve.start_param() if periodic else 0.0
            if periodic and period <= 0.0:
                fail("Invalid curve period")

            def unwrap(reference, u):
                return u + period * round((reference - u) / period)

            u2_local = unwrap(u1, u2) if periodic else u2
            lm = tuple(x1[i] + s * (x2[i] - x1[i]) for i in range(3))
            cp = tuple(curve.closest_point_trimmed(list(lm)))
            u_cp = curve.u_from_position(list(cp))
            u_expected = u1 + s * (u2_local - u1)
            u_cp_local = unwrap(u_expected, u_cp) if periodic else u_cp
            umin, umax = min(u1, u2_local), max(u1, u2_local)
            tol = 0.01 * abs(umax - umin)
            chord = _dist(x1, x2)
            endpoint_tol = 1.0e-8 * chord + 1.0e-15
            stuck = 1.0e-8 < s < 1.0 - 1.0e-8 and (
                _dist(cp, x1) <= endpoint_tol or _dist(cp, x2) <= endpoint_tol)
            if not stuck and umin - tol <= u_cp_local <= umax + tol:
                self.diag["between_edge_projection"] += 1
                return cp
            if periodic:
                out = tuple(curve.position_from_u(u_expected))
            else:
                segment = curve.length_from_u(u1, u2_local)
                out = tuple(curve.point_from_arc_length(u1, s * segment))
            dd = _dist(out, lm)
            if not math.isfinite(dd) or dd > chord:
                fail("Arc interpolation exceeds unchanged chord limit")
            self.diag["between_edge_arc"] += 1
            return out

        def GetTangent(self, p, s1, s2, egi, edgenr):
            self.diag["tangent_calls"] += 1
            _, curve = self._curve(edgenr, s1, s2)
            t = curve.tangent([p[0], p[1], p[2]])
            norm = math.sqrt(t[0] ** 2 + t[1] ** 2 + t[2] ** 2)
            if not math.isfinite(norm) or norm <= 0.0:
                self.diag["edge_failures"] += 1
                raise ExportError(f"Invalid tangent for Netgen geometric edge {edgenr}")
            return ngm.Vec(t[0] / norm, t[1] / norm, t[2] / norm)

    return CubitGeometry


# ---------------------------------------------------------------------------
# Netgen mesh construction and export (port of NetgenCurver + ExportNetgenCommand)
# ---------------------------------------------------------------------------

_GEOM_ACIS = 13  # Netgen Mesh::GEOM_ACIS, the geometry type of the native exporter


def _explicit_surface_label(cubit, sid: int) -> str:
    label = cubit.get_entity_name("surface", sid)
    if label in (f"Surface {sid}", f"Surface_{sid}", f"surface_{sid}"):
        return ""
    return label


def _boundary_name(cubit, bc: int, sideset_label: dict[int, str]) -> str:
    """Sideset name, else the Cubit surface name, else ``surface_<id>``; no spaces.

    ``bc`` is a Cubit surface id, or ``-sideset_id`` for a free sideset face set.
    """
    label = sideset_label.get(bc)
    if label is None:
        label = cubit.get_entity_name("surface", bc) or f"surface_{bc}"
    return label.replace(" ", "_")


def _mesh_area(points: list[tuple[float, float, float]]) -> float:
    import numpy as np  # noqa: PLC0415

    p = np.asarray(points, dtype=float)
    if len(p) == 3:
        return 0.5 * float(np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0])))
    x, w = np.polynomial.legendre.leggauss(8)
    x, w = 0.5 * (x + 1.0), 0.5 * w
    area = 0.0
    for u, wu in zip(x, w):
        for v, wv in zip(x, w):
            du = -(1 - v) * p[0] + (1 - v) * p[1] + v * p[2] - v * p[3]
            dv = -(1 - u) * p[0] - u * p[1] + u * p[2] + (1 - u) * p[3]
            area += wu * wv * float(np.linalg.norm(np.cross(du, dv)))
    return area


def _serialize(mesh, near: Path) -> bytes:
    """Save through Netgen and return the bytes with LF line ends only."""
    temporary = near.with_name(near.name + f".{os.getpid()}.tmp.vol")
    mesh.Save(str(temporary))
    try:
        data = temporary.read_bytes().replace(b"\r\n", b"\n")
    finally:
        temporary.unlink()
    if b"\r" in data:
        raise ExportError("Netgen .vol serialization contains a carriage return")
    return data


def export_vol(out_path: str | os.PathLike[str], order: int, *, cubit,
               journal: str | os.PathLike[str] | None = None,
               periodic: dict[str, Any] | None = None) -> dict[str, Any]:
    """Export the current Cubit model (optionally after playing a journal).

    ``periodic`` holds the nine explicit ``periodic_*`` contract fields (see
    :mod:`radia_mcp.cubit.periodic`); the output must not exist yet.

    Returns the companion JSON record, which is also written next to the
    ``.vol`` as ``<out_path>.json``.  Raises :class:`ExportError` when the
    model cannot be exported faithfully.
    """
    import netgen.meshing as ngm  # noqa: PLC0415

    if not 1 <= int(order) <= 5:
        raise ExportError("Order must be 1-5.")
    order = int(order)
    if not hasattr(ngm.Mesh, "GeometricRefit"):
        raise ExportError("This Netgen runtime lacks Mesh.GeometricRefit")
    t_start = time.perf_counter()
    try:
        contract = _periodic.parse_contract(periodic or {})
        if contract is not None:
            _periodic.require_fresh(Path(out_path).resolve())
    except _periodic.PeriodicError as error:
        raise ExportError(f"Named periodic contract rejected: {error}") from error
    if journal is not None:
        journal = Path(journal).resolve()
        if not journal.is_file():
            raise ExportError(f"Journal not found: {journal}")
        cubit.cmd(f'playback "{journal.as_posix()}"')
        if cubit.get_error_count():
            raise ExportError(f"Cubit reported {cubit.get_error_count()} errors playing {journal}")

    md = extract_mesh(cubit)
    out_path = Path(out_path).resolve()

    # -- domains and face ownership ---------------------------------------
    block_to_domain = {bid: i + 1 for i, bid in enumerate(md.block_ids)}
    ndomains = len(md.block_ids)
    vol_to_domain: dict[int, int] = {}
    for bid in md.block_ids:
        for vid in cubit.get_block_volumes(bid):
            vol_to_domain[vid] = block_to_domain[bid]
    face_domains: dict[tuple[int, ...], set[int]] = {}
    for element in md.elements:
        for local in _VOLUME_FACES[element.kind]:
            key = tuple(sorted(element.conn[i] for i in local))
            face_domains.setdefault(key, set()).add(block_to_domain[element.block])

    sideset_label: dict[int, str] = {}
    for ss in md.sidesets:
        label = ss.name or f"sideset_{ss.id}"
        sideset_label[-ss.id] = label
        for sid in cubit.get_sideset_surfaces(ss.id):
            sideset_label[sid] = label

    mesh = ngm.Mesh(dim=3)
    mesh.geomtype = _GEOM_ACIS  # writes the surface element UV columns
    pid = {nid: mesh.Add(ngm.MeshPoint(ngm.Pnt(*md.coords[nid]))) for nid in md.node_ids}

    # -- face descriptors, one per meshed Cubit surface -----------------------
    surfaces = list(cubit.get_entities("surface"))
    fd_of_surface: dict[int, int] = {}
    fd_bc: list[int] = []  # original BC property per descriptor (surface id or -sideset)
    fd_domains: list[tuple[int, int]] = []
    surface_faces: dict[int, tuple[list[tuple[int, ...]], list[tuple[int, ...]]]] = {}
    seam_descriptors: list[int] = []
    seam_surfaces: list[int] = []
    for sid in surfaces:
        tris = [_connectivity(cubit, "tri", t, 3)
                for t in cubit.parse_cubit_list("tri", f"in surface {sid}")]
        quads = [_connectivity(cubit, "face", q, 4)
                 for q in cubit.parse_cubit_list("face", f"in surface {sid}")]
        if not tris and not quads:
            continue
        parents = list(cubit.get_relatives("surface", sid, "volume"))
        din = vol_to_domain.get(parents[0], 0) if parents else 0
        dout = vol_to_domain.get(parents[1], 0) if len(parents) >= 2 else 0
        fd = len(fd_bc) + 1
        if din > 0 and din == dout:
            # A webcut seam inside one material: it supports curving and the
            # conforming mesh, but is not an FE boundary.  It is removed after
            # curving, unless it carries a user label (then refuse).
            label = sideset_label.get(sid) or _explicit_surface_label(cubit, sid)
            if label:
                raise ExportError(
                    f"Same-material internal surface {sid} is explicitly labelled '{label}'. "
                    "Netgen cannot encode it as a boundary with identical DomainIn/DomainOut. "
                    "Remove that surface from the sideset, or model the intended interface "
                    "with distinct domains.")
            seam_descriptors.append(fd)
            seam_surfaces.append(sid)
        fd_of_surface[sid] = fd
        fd_bc.append(sid)
        fd_domains.append((din, dout))
        surface_faces[sid] = (tris, quads)

    # -- canonical per-node UV per surface and bounding boxes -----------------
    polar = {sid for sid in fd_of_surface
             if cubit.surface(sid).is_planar() and _has_pole(cubit, sid)}
    canonical_uv: dict[tuple[int, int], tuple[float, float]] = {}
    bboxes: dict[int, tuple[list[float], list[float]]] = {}
    for sid, (tris, quads) in surface_faces.items():
        surface = cubit.surface(sid)
        lo, hi = [math.inf] * 3, [-math.inf] * 3
        for conn in tris + quads:
            for nid in conn:
                if (sid, nid) not in canonical_uv:
                    x = md.coords[nid]
                    canonical_uv[(sid, nid)] = tuple(surface.u_v_from_position(list(x)))
                    for i in range(3):
                        lo[i], hi[i] = min(lo[i], x[i]), max(hi[i], x[i])
        bboxes[sid] = (lo, hi)
    global_diam = max(_dist(lo, hi) for lo, hi in bboxes.values()) if bboxes else 0.0
    pad = 0.05 * global_diam

    added_faces: set[tuple[int, ...]] = set()
    seam_shifts = 0
    for i, (din, dout) in enumerate(fd_domains):
        descriptor = ngm.FaceDescriptor(surfnr=i + 1, domin=din, domout=dout, bc=i + 1)
        descriptor.tlosurf = i + 1
        mesh.Add(descriptor)
    for sid, (tris, quads) in surface_faces.items():
        fd = fd_of_surface[sid]
        up, vp = _periods(cubit, sid)
        for conn in tris + quads:
            uv = [canonical_uv[(sid, nid)] for nid in conn]
            us = [uv[0][0]] + [_seam(uv[0][0], u, up) for u, _ in uv[1:]]
            vs = [uv[0][1]] + [_seam(uv[0][1], v, vp) for _, v in uv[1:]]
            seam_shifts += sum(1 for a, b in zip(us, (u for u, _ in uv)) if a != b)
            seam_shifts += sum(1 for a, b in zip(vs, (v for _, v in uv)) if a != b)
            mesh.Add(ngm.Element2D(fd, [pid[n] for n in conn], uv=list(zip(us, vs))))
            added_faces.add(tuple(sorted(conn)))

    # -- free sideset faces (no CAD surface) ----------------------------------
    free_faces = 0
    for ss in md.sidesets:
        pending: dict[tuple[int, int], list[tuple[int, ...]]] = {}
        for conn in ss.faces:
            key = tuple(sorted(conn))
            if key in added_faces:
                continue
            adjacent = sorted(face_domains.get(key, ()))
            if not adjacent:
                raise ExportError(f"Free sideset {ss.id} contains a face with no adjacent volume element.")
            if len(adjacent) > 2:
                raise ExportError(f"Free sideset {ss.id} contains a non-manifold face adjacent to "
                                  f"{len(adjacent)} domains.")
            pending.setdefault((adjacent[0], adjacent[1] if len(adjacent) == 2 else 0), []).append(conn)
        for (din, dout), faces in sorted(pending.items()):
            fd = len(fd_bc) + 1
            fd_bc.append(-ss.id)
            fd_domains.append((din, dout))
            descriptor = ngm.FaceDescriptor(surfnr=fd, domin=din, domout=dout, bc=fd)
            descriptor.tlosurf = fd
            mesh.Add(descriptor)
            for conn in faces:
                mesh.Add(ngm.Element2D(fd, [pid[n] for n in conn], uv=[(0.0, 0.0)] * len(conn)))
                added_faces.add(tuple(sorted(conn)))
                free_faces += 1
    for fd, (din, dout) in enumerate(fd_domains, start=1):
        if fd_bc[fd - 1] < 0 and (not 1 <= din <= ndomains or not 0 <= dout <= ndomains or din == dout):
            raise ExportError(f"Free sideset descriptor {-fd_bc[fd - 1]} has invalid domain ownership "
                              f"({din} -> {dout}) for a {ndomains}-domain mesh.")

    # -- segments on curves shared by two meshed surfaces ---------------------
    curve_to_edge: dict[int, int] = {}
    edge_bindings: dict[int, tuple[int, int, int]] = {}
    edgenr = 1
    for cid in cubit.parse_cubit_list("curve", "all"):
        parents = list(cubit.parse_cubit_list("surface", f"in curve {cid}"))
        if len(parents) < 2:
            continue
        fd1, fd2 = fd_of_surface.get(parents[0]), fd_of_surface.get(parents[1])
        if fd1 is None or fd2 is None:
            continue
        ed = ngm.EdgeDescriptor()
        ed.edgenr = edgenr
        ed.surfnr = (fd1, fd2)
        ed.index = fd2
        descriptor = mesh.Add(ed)
        segments = 0
        edge_ids = list(cubit.parse_cubit_list("edge", f"in curve {cid}"))
        curve = cubit.curve(cid)
        length, u_start = curve.length(), curve.start_param()
        for eid in edge_ids:
            conn = _connectivity(cubit, "edge", eid, 2)
            dist = [0.0, 0.0]
            if length > 0:
                dist = [curve.length_from_u(u_start, curve.u_from_position(list(md.coords[n])))
                        / length for n in conn[:2]]
            mesh.Add(ngm.Element1D([pid[conn[0]], pid[conn[1]]], index=descriptor, dist=dist))
            segments += 1
        if segments != len(edge_ids):
            raise ExportError(f"Cubit curve {cid} exported only {segments}/{len(edge_ids)} segments")
        if segments:
            curve_to_edge[cid] = descriptor
        if edgenr in edge_bindings:
            raise ExportError(f"Ambiguous Netgen geometric edge {edgenr}")
        edge_bindings[edgenr] = (cid, fd1, fd2)
        edgenr += 1

    # -- volume elements (Cubit vertex order, as the native exporter) ---------
    for element in md.elements:
        mesh.Add(ngm.Element3D(block_to_domain[element.block], [pid[n] for n in element.conn]))
    for bid in md.block_ids:
        mesh.SetMaterial(block_to_domain[bid], cubit.get_exodus_entity_name("block", bid)
                         or f"volume_{bid}")

    # -- curving ----------------------------------------------------------------
    mesh.UpdateTopology()
    geometry_class = _make_geometry_class(ngm)
    fd_to_surface = {fd: sid for sid, fd in fd_of_surface.items()}
    geometry = geometry_class(cubit, fd_to_surface, edge_bindings, polar, bboxes, pad)
    refit: dict[str, Any] | None = None
    if order >= 2:
        mesh.SetGeometry(geometry)
        try:
            mesh.Curve(order)
            if order >= 3:
                refit = dict(mesh.GeometricRefit())
        finally:
            mesh.SetGeometry(ngm.NetgenGeometry())
    removed_seam_elements = 0
    if seam_descriptors:
        removed_seam_elements = mesh.DeleteFaceDescriptors(seam_descriptors)
        dropped = set(seam_descriptors)
        fd_bc = [bc for fd, bc in enumerate(fd_bc, start=1) if fd not in dropped]
        fd_domains = [d for fd, d in enumerate(fd_domains, start=1) if fd not in dropped]

    # -- boundary names ---------------------------------------------------------
    for fd, bc in enumerate(fd_bc, start=1):
        mesh.SetBCName(fd - 1, _boundary_name(cubit, bc, sideset_label))

    names = {bid: cubit.get_exodus_entity_name("block", bid).lower() for bid in md.block_ids}
    air = next((block_to_domain[b] for b, n in names.items() if n == "air"), 0)
    kelvin = next((block_to_domain[b] for b, n in names.items() if n == "kelvin"), 0)
    kelvin_counts = {"int": 0, "ext": 0}
    if air and kelvin:
        for fd, (din, dout) in enumerate(fd_domains, start=1):
            if fd_bc[fd - 1] in sideset_label:
                continue
            if {din, dout} == {air, kelvin}:
                mesh.SetBCName(fd - 1, "kelvin_int")
                kelvin_counts["int"] += 1
            elif {din, dout} == {kelvin, 0}:
                mesh.SetBCName(fd - 1, "kelvin_ext")
                kelvin_counts["ext"] += 1

    # -- curve sidesets (BBND names) and nodesets (CD3 names) -----------------
    edge_labels: dict[int, str] = {}
    for ss in md.sidesets:
        label = (ss.name or f"sideset_{ss.id}").replace(" ", "_")
        for cid in cubit.get_sideset_curves(ss.id):
            if label == "default":
                raise ExportError(f"Curve sideset {ss.id} uses reserved BBND name 'default'.")
            if cid not in curve_to_edge:
                raise ExportError(f"Curve sideset {ss.id} '{label}': curve {cid} cannot be "
                                  "represented as a complete BBND edge in the exported volume mesh.")
            previous = edge_labels.setdefault(curve_to_edge[cid], label)
            if previous != label:
                raise ExportError(f"Curve {cid} has conflicting BBND sideset names "
                                  f"'{previous}' and '{label}'.")
    for descriptor, label in edge_labels.items():
        mesh.SetCD2Name(descriptor, label)
    # CD3 numbers are 1-based in Mesh::SetCD3Name and Element0D.index; a
    # repeated nodeset name reuses its number, as Mesh::AddCD3Name does.
    cd3_index: dict[str, int] = {}
    for ns in md.nodesets:
        if not ns.name or not ns.node_ids:
            continue
        cd3 = cd3_index.get(ns.name)
        if cd3 is None:
            cd3 = cd3_index[ns.name] = len(cd3_index) + 1
            mesh.SetCD3Name(cd3, ns.name)
        for nid in ns.node_ids:
            if nid in pid:
                mesh.Add(ngm.Element0D(pid[nid], index=cd3))

    curved_pyramids = sum(1 for el in mesh.Elements3D()
                          if len(el.vertices) == 5 and el.curved)

    periodic_result = None
    if contract is not None:
        try:
            periodic_result = _periodic.apply_identification(
                mesh, contract, lambda m: _serialize(m, out_path).decode("utf-8"))
        except _periodic.PeriodicError as error:
            raise ExportError(f"Named periodic export rejected: {error}") from error
    vol_bytes = _serialize(mesh, out_path)

    # -- companion JSON -----------------------------------------------------------
    materials: dict[str, float] = {}
    mesh_only: list[str] = []
    for bid in md.block_ids:
        name = cubit.get_exodus_entity_name("block", bid) or f"volume_{bid}"
        volumes = list(cubit.parse_cubit_list("volume", f"in block {bid}"))
        if not volumes:
            mesh_only.append(name)
            continue
        materials[name] = sum(cubit.get_volume_volume(v) for v in volumes)
    boundaries: dict[str, float] = {}
    for fd, bc in enumerate(fd_bc, start=1):
        name = mesh.GetBCName(fd - 1)
        if bc < 0:
            area = sum(_mesh_area([md.coords[n] for n in conn])
                       for ss in md.sidesets if ss.id == -bc for conn in ss.faces
                       if tuple(sorted(conn)) in added_faces and conn and
                       face_domains.get(tuple(sorted(conn))))
        else:
            area = cubit.get_surface_area(bc)
        boundaries[name] = boundaries.get(name, 0.0) + area
    edges: dict[str, float] = {}
    support = {bc for bc in fd_bc if bc > 0} | set(seam_surfaces)
    for cid in cubit.parse_cubit_list("curve", "all"):
        parents = list(cubit.parse_cubit_list("surface", f"in curve {cid}"))
        if len(parents) < 2 or parents[0] not in support or parents[1] not in support:
            continue
        if cubit.parse_cubit_list("edge", f"in curve {cid}"):
            edges[f"curve_{cid}"] = cubit.get_curve_length(cid)

    diag = dict(geometry.diag)
    diag.update({"error": "", "polar_surfaces": sorted(polar), "seam_uv_shifts": seam_shifts,
                 "free_sideset_faces": free_faces,
                 "same_material_seams": {"descriptors": len(seam_descriptors),
                                         "surface_elements": removed_seam_elements}, "curved_pyramids": curved_pyramids,
                 "geometric_refit": refit, "kelvin_autodetect": kelvin_counts,
                 "warnings": md.warnings})
    record: dict[str, Any] = {"export_passed": True}
    if periodic_result is not None:
        record["periodic_identification"] = _periodic.companion_record(
            periodic_result, {fd: bc for fd, bc in enumerate(fd_bc, start=1)})
    record.update({
        "curving": diag,
        "materials": materials,
        "mesh_only_materials": mesh_only,
        "boundaries": dict(sorted(boundaries.items())),
        "edges": edges,
        "n_elements": len(md.elements),
        "n_points": len(md.node_ids),
        "order": order,
        "export_time_s": round(time.perf_counter() - t_start, 3),
    })
    record_bytes = (json.dumps(record, indent=2) + "\n").encode("utf-8")
    if contract is not None:
        try:
            _periodic.publish_bundle(out_path, vol_bytes, record_bytes)
        except (_periodic.PeriodicError, OSError) as error:
            raise ExportError(f"Named periodic bundle publication failed: {error}. Failed outputs "
                              "are unqualified; remove only your own leftovers before retrying.") from error
    else:
        out_path.write_bytes(vol_bytes)
        Path(str(out_path) + ".json").write_bytes(record_bytes)
    return record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cubit-bin", default=os.environ.get("RADIA_CUBIT_BIN"),
                        help="Cubit bin directory (default: $RADIA_CUBIT_BIN)")
    parser.add_argument("--journal", required=True, help="Cubit journal that builds and meshes the model")
    parser.add_argument("--out", required=True, help="output .vol path")
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--result", help="status JSON path (written on success and failure)")
    for name in _periodic.FIELDS:
        parser.add_argument("--" + name.replace("_", "-"), dest=name,
                            help="named periodic contract field (all nine or none)")
    args = parser.parse_args(argv)
    status: dict[str, Any]
    try:
        if not args.cubit_bin:
            raise ExportError("Pass --cubit-bin or set RADIA_CUBIT_BIN")
        cubit = import_cubit(args.cubit_bin)
        record = export_vol(args.out, args.order, cubit=cubit, journal=args.journal,
                            periodic={k: getattr(args, k) for k in _periodic.FIELDS})
        status = {"status": "ok", "vol": str(Path(args.out).resolve()), "record": record}
        code = 0
    except Exception as error:  # noqa: BLE001  (reported to the caller as JSON)
        status = {"status": "error", "error_type": type(error).__name__, "error": str(error)}
        code = 1
    if args.result:
        Path(args.result).write_text(json.dumps(status, indent=2) + "\n",
                                     encoding="utf-8", newline="\n")
    print(json.dumps({k: status[k] for k in status if k != "record"}))
    return code


if __name__ == "__main__":
    raise SystemExit(main())
