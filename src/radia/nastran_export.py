"""Write a Netgen/NGSolve mesh as a first- or second-order Nastran BDF.

The input is a checked ``.vol`` file or an in-memory mesh (``ngsolve.Mesh``
or ``netgen.meshing.Mesh``).  Second-order nodes are the curved geometry map
of each element evaluated at its reference edge midpoints, so a curved
``.vol`` (Cubit plugin export or ``mesh.Curve(p)``) keeps its geometry.

The BDF is a derived mesh for other tools; the solver boundary stays the
checked ``.vol``.

Cards written:
    GRID*            every vertex, then every mid-edge node (large field)
    CTETRA, CHEXA, CPENTA, CPYRAM   volume elements (4/10, 8/20, 6/15, 5/13)
    CTRIA3/6, CQUAD4/8              boundary faces, or the elements of a 2-D mesh
    PSOLID / PSHELL                 one property per material / boundary name

Material and boundary names go into comments beside their property IDs.

Nastran has no face or centre nodes: CHEXA20, CPENTA15, CPYRAM13 and CQUAD8
carry the vertices and mid-edge nodes only.  A curved order-2 Netgen hex also
has face and cell terms, so its CHEXA20 volume differs slightly from the .vol
(about 1e-4 relative on a coarse swept cylinder); tets and triangles are exact.

Usage:
    from radia.nastran_export import write_nastran_bdf
    write_nastran_bdf("model.vol", "model.bdf", order=2)

    python -m radia.nastran_export model.vol model.bdf --order 2
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np

__all__ = ["write_nastran_bdf"]

# Mid-edge node order per Nastran card, as vertex pairs in card-local
# numbering (G1 = 0).
_BDF_EDGES = {
    "TET": ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)),
    "HEX": ((0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 5), (2, 6), (3, 7),
            (4, 5), (5, 6), (6, 7), (7, 4)),
    "PRISM": ((0, 1), (1, 2), (2, 0), (0, 3), (1, 4), (2, 5),
              (3, 4), (4, 5), (5, 3)),
    "PYRAMID": ((0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4)),
    "TRIG": ((0, 1), (1, 2), (2, 0)),
    "QUAD": ((0, 1), (1, 2), (2, 3), (3, 0)),
}

_CARD = {"TET": "CTETRA", "HEX": "CHEXA", "PRISM": "CPENTA",
         "PYRAMID": "CPYRAM", "TRIG": ("CTRIA3", "CTRIA6"),
         "QUAD": ("CQUAD4", "CQUAD8")}

# Corner frame (origin, a, b, c) whose determinant is positive for a
# correctly oriented solid, and the reflection that flips a negative one.
# Each reflection is a symmetry of the reference cell, so the set of edges
# (and so the mid-edge nodes) is unchanged by it.
_FRAME = {"TET": (0, 1, 2, 3), "HEX": (0, 1, 3, 4),
          "PRISM": (0, 1, 2, 3), "PYRAMID": (0, 1, 3, 4)}
_MIRROR = {"TET": (0, 2, 1, 3), "HEX": (4, 5, 6, 7, 0, 1, 2, 3),
           "PRISM": (3, 4, 5, 0, 1, 2), "PYRAMID": (0, 3, 2, 1, 4)}

# A shared mid-edge node is evaluated from every element that holds the
# edge; the copies must agree to this fraction of the mesh size.
_SHARED_NODE_TOLERANCE = 1e-10


def _et_name(el) -> str:
    return str(el.type).split(".")[-1]


def _as_ngsolve_mesh(mesh):
    import ngsolve
    if isinstance(mesh, (str, os.PathLike)):
        path = Path(mesh)
        if not path.is_file():
            raise FileNotFoundError(f"mesh file not found: {path}")
        return ngsolve.Mesh(str(path))
    if isinstance(mesh, ngsolve.Mesh):
        return mesh
    from netgen.meshing import Mesh as NetgenMesh
    if isinstance(mesh, NetgenMesh):
        return ngsolve.Mesh(mesh)
    raise TypeError(f"expected a .vol path, ngsolve.Mesh or netgen Mesh, not {type(mesh)!r}")


def _reference_vertices(et_name: str):
    from ngsolve import ET
    from ngsolve.fem import ElementTopology
    verts = [tuple(float(c) for c in v)
             for v in ElementTopology(getattr(ET, et_name)).vertices]
    if et_name == "PYRAMID":
        verts[4] = (0.0, 0.0, 1.0)  # the topology lists the apex at 1 - 1e-12
    return verts


def _mid_edge_coordinates(mesh, vb, elements, edges_of):
    """Evaluate the element map at reference edge midpoints.

    Returns {(gv_a, gv_b): xyz} with gv_a < gv_b, and the largest
    disagreement between elements that share an edge.
    """
    from ngsolve import ET, IntegrationRule, CF, x, y, z

    by_type = {}
    for el in elements:
        by_type.setdefault(_et_name(el), []).append(el)
    if not by_type:
        return {}, 0.0

    # Point 0 of each rule is reference vertex 0: it checks that the mapped
    # points are returned in element order.
    rules, local_edges = {}, {}
    for name in by_type:
        ref = _reference_vertices(name)
        dim = len(ref[0])
        edges = edges_of(name)
        pts = [ref[0]] + [tuple(0.5 * (ref[a][k] + ref[b][k]) for k in range(dim))
                          for a, b in edges]
        rules[getattr(ET, name)] = IntegrationRule(pts, [0.0] * len(pts))
        local_edges[name] = edges
    mapped = mesh.MapToAllElements(rules, vb)
    coords = np.asarray(CF((x, y, z))(mapped), dtype=float).reshape(-1, 3)

    entries, row = [], 0
    for el in mesh.Elements(vb):
        name = _et_name(el)
        n_pts = 1 + len(local_edges[name])
        block = coords[row:row + n_pts]
        row += n_pts
        gv = [v.nr for v in el.vertices]
        v0 = np.asarray(_point3(mesh.vertices[gv[0]].point))
        if np.linalg.norm(block[0] - v0) > 1e-9 * (1.0 + np.linalg.norm(v0)):
            raise ValueError(
                f"the map of element {el.nr} does not pass through its vertex "
                f"{gv[0]} (off by {np.linalg.norm(block[0] - v0):.3e}); a deformed "
                "mesh cannot be written, its GRID points would not match")
        entries.append((gv, local_edges[name], block[1:]))
    if row != coords.shape[0]:
        raise RuntimeError(f"mapped {coords.shape[0]} points, consumed {row}")
    return _merge_mid_edges(entries)


def _merge_mid_edges(entries):
    """One node per global edge from (vertices, local edges, xyz rows) entries.

    Returns {(gv_a, gv_b): xyz} with gv_a < gv_b, and the largest distance
    between two elements' copies of the same edge node.
    """
    out, mismatch = {}, 0.0
    for gv, edges, rows in entries:
        for (a, b), xyz in zip(edges, rows):
            key = (min(gv[a], gv[b]), max(gv[a], gv[b]))
            prev = out.get(key)
            if prev is None:
                out[key] = np.asarray(xyz, dtype=float)
            else:
                mismatch = max(mismatch, float(np.linalg.norm(prev - np.asarray(xyz))))
    return out, mismatch


def _point3(p):
    return (float(p[0]), float(p[1]) if len(p) > 1 else 0.0,
            float(p[2]) if len(p) > 2 else 0.0)


def _frame_det(pts, frame):
    o, a, b, c = (np.asarray(pts[i]) for i in frame)
    return float(np.linalg.det(np.stack([a - o, b - o, c - o])))


def _fields(name: str, values) -> list[str]:
    """Small-field card lines: 8 data fields per line, '+' continuations."""
    values = list(values)
    lines, first = [], True
    while values or first:
        chunk, values = values[:8], values[8:]
        head = f"{name:<8}" if first else "+" + " " * 7
        line = head + "".join(f"{int(v):8d}" for v in chunk)
        if values:
            line = f"{line:<72}+"
        lines.append(line)
        first = False
    return lines


def write_nastran_bdf(mesh, path, *, order: int | None = None,
                      boundaries: bool = True, overwrite: bool = False) -> dict:
    """Write ``mesh`` to ``path`` as a Nastran BDF and return a summary.

    Args:
        mesh: ``.vol`` path, ``ngsolve.Mesh`` or ``netgen.meshing.Mesh``.
        path: output ``.bdf``; an existing file is refused unless ``overwrite``.
        order: 1 or 2.  Default: 2 for a curved mesh, else 1.  Order 2 needs a
            curved mesh (curve order >= 2); a straight mesh is refused rather
            than written with straight mid-edge nodes.
        boundaries: also write the boundary faces of a 3-D mesh as shell
            elements, one PSHELL per boundary name.  Interface faces are
            boundary faces too in a ``.vol``.
        overwrite: replace an existing output file.

    Raises:
        ValueError: unsupported order or element type, an order-2 request on a
            straight mesh, an inverted element, or shared mid-edge nodes that
            disagree.
        FileExistsError: ``path`` exists and ``overwrite`` is False.
    """
    from ngsolve import VOL, BND

    m = _as_ngsolve_mesh(mesh)
    curve_order = int(m.GetCurveOrder())
    if order is None:
        order = 2 if curve_order >= 2 else 1
    if order not in (1, 2):
        raise ValueError(f"Nastran BDF order must be 1 or 2, not {order}")
    if order == 2 and curve_order < 2:
        raise ValueError(
            "order 2 needs a curved mesh: this mesh has curve order "
            f"{curve_order}. Export a curved .vol or call mesh.Curve(2) first")
    out_path = Path(path)
    if out_path.exists() and not overwrite:
        raise FileExistsError(f"{out_path} exists; pass overwrite=True to replace it")

    vol_elems = list(m.Elements(VOL))
    shell_elems = list(m.Elements(BND)) if (boundaries and m.dim == 3) else []
    if m.dim == 2:
        vol_elems, shell_elems = [], list(m.Elements(VOL))
    supported = set(_BDF_EDGES)
    for el in vol_elems + shell_elems:
        if _et_name(el) not in supported:
            raise ValueError(f"element type {_et_name(el)} has no Nastran card here")

    vertices = [_point3(v.point) for v in m.vertices]
    lo, hi = np.min(vertices, axis=0), np.max(vertices, axis=0)
    size = float(np.linalg.norm(hi - lo)) or 1.0

    mid = {}
    mismatch = 0.0
    if order == 2:
        def edges_of(name):
            return _BDF_EDGES[name]
        if vol_elems:
            got, mm = _mid_edge_coordinates(m, VOL, vol_elems, edges_of)
            mid.update(got)
            mismatch = max(mismatch, mm)
        if shell_elems:
            got, mm = _mid_edge_coordinates(
                m, BND if m.dim == 3 else VOL, shell_elems, edges_of)
            for key, xyz in got.items():
                prev = mid.get(key)
                if prev is None:
                    mid[key] = xyz
                else:
                    mismatch = max(mismatch, float(np.linalg.norm(prev - xyz)))
            mismatch = max(mismatch, mm)
        if mismatch > _SHARED_NODE_TOLERANCE * size:
            raise ValueError(
                f"shared mid-edge nodes disagree by {mismatch:.3e} "
                f"(mesh size {size:.3e}); the element maps are not continuous")

    nv = len(vertices)
    mid_id = {key: nv + 1 + i for i, key in enumerate(sorted(mid))}

    solid_names = sorted({el.mat for el in vol_elems})
    if m.dim == 2:
        shell_names = sorted({el.mat for el in shell_elems})
    else:
        shell_names = sorted({m.GetBoundaries()[el.index] for el in shell_elems})
    pid_solid = {n: i + 1 for i, n in enumerate(solid_names)}
    pid_shell = {n: len(solid_names) + i + 1 for i, n in enumerate(shell_names)}

    lines = ["$ Nastran BDF written by radia.nastran_export",
             f"$ element order {order}, mesh curve order {curve_order}",
             "BEGIN BULK"]
    for i, (px, py, pz) in enumerate(vertices, start=1):
        lines += [f"GRID*   {i:16d}{0:16d}{px:16.9e}{py:16.9e}+",
                  f"*       {pz:16.9e}"]
    for key in sorted(mid):
        px, py, pz = mid[key]
        lines += [f"GRID*   {mid_id[key]:16d}{0:16d}{px:16.9e}{py:16.9e}+",
                  f"*       {pz:16.9e}"]

    counts, flipped, eid = {}, 0, 0

    def card_nodes(name, gv):
        nodes = [g + 1 for g in gv]
        if order == 2:
            nodes += [mid_id[(min(gv[a], gv[b]), max(gv[a], gv[b]))]
                      for a, b in _BDF_EDGES[name]]
        return nodes

    for el in vol_elems:
        name = _et_name(el)
        gv = [v.nr for v in el.vertices]
        det = _frame_det([vertices[g] for g in gv], _FRAME[name])
        if det < 0:
            gv = [gv[i] for i in _MIRROR[name]]
            flipped += 1
            det = _frame_det([vertices[g] for g in gv], _FRAME[name])
        if not det > 0:
            raise ValueError(f"volume element {el.nr} ({name}) is degenerate "
                             f"(corner frame determinant {det:.3e})")
        eid += 1
        card = _CARD[name]
        lines += _fields(card, [eid, pid_solid[el.mat]] + card_nodes(name, gv))
        counts[card] = counts.get(card, 0) + 1

    for el in shell_elems:
        name = _et_name(el)
        gv = [v.nr for v in el.vertices]
        label = el.mat if m.dim == 2 else m.GetBoundaries()[el.index]
        eid += 1
        card = _CARD[name][order - 1]
        lines += _fields(card, [eid, pid_shell[label]] + card_nodes(name, gv))
        counts[card] = counts.get(card, 0) + 1

    for n in solid_names:
        lines.append(f"$ material {n}")
        lines += _fields("PSOLID", [pid_solid[n], pid_solid[n]])
    for n in shell_names:
        lines.append(f"$ {'material' if m.dim == 2 else 'boundary'} {n}")
        lines += _fields("PSHELL", [pid_shell[n], pid_shell[n]])
    lines.append("ENDDATA")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_name(out_path.name + ".tmp")
    tmp.write_text("\n".join(lines) + "\n", encoding="ascii", newline="\n")
    os.replace(tmp, out_path)

    return {"path": str(out_path), "order": order, "curve_order": curve_order,
            "n_grid": nv + len(mid), "n_vertices": nv, "n_mid_edge": len(mid),
            "elements": counts, "flipped_to_positive": flipped,
            "shared_node_mismatch": mismatch,
            "properties": {"solid": pid_solid, "shell": pid_shell}}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m radia.nastran_export",
                                 description=__doc__.split("\n")[0])
    ap.add_argument("mesh", help="input .vol")
    ap.add_argument("bdf", help="output .bdf")
    ap.add_argument("--order", type=int, choices=(1, 2))
    ap.add_argument("--no-boundaries", action="store_true",
                    help="write only the volume elements of a 3-D mesh")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args(argv)
    try:
        info = write_nastran_bdf(args.mesh, args.bdf, order=args.order,
                                 boundaries=not args.no_boundaries,
                                 overwrite=args.overwrite)
    except (ValueError, FileExistsError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(info, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
