"""Standalone writer of a static magnetic solver input in the .jcf container format.

Writes from mesh data in memory with the standard library and NumPy only; no solver,
pre-processor, license or other Python package is needed. Only the parameters known to be
read by the solver are written (static analysis: linear materials, external field, linear
permanent magnets, nonlinear B-H steel with lamination factor, coils with a flow-in face or a
direction vector, symmetry boundaries; hexahedra, tetrahedra, prisms and pyramids).

Provenance: derived from the template-free builder of the author's private solver-interface
package (MIT), reduced to writing. Additions here: a magnet whose direction is not +z gets its
own coordinate system referenced by ``local_coordinate_id`` (the solver reads the
magnetisation along that system's z axis and ignores ``axis2``), and a coil's flow-in face may
be given as face node sets (a Cubit sideset).
"""
from __future__ import annotations

import collections
import os
import struct
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path
from xml.sax.saxutils import escape, quoteattr

import numpy as np

TET, HEX, PRISM, PYRAMID = 3, 1, 5, 13
UNITS = {"m": 1.0, "mm": 1e-3, "cm": 1e-2, "um": 1e-6}
_NODES_PER = {1: 8, 3: 4, 5: 6, 13: 5}
_TEXT_TO_BIN = {7: 3, 3: 5, 1: 6, 5: 8, 13: 9, 11: 2}

# Local face k (1-based) of a face group entry (element, k), as local node indices.
LOCAL_FACES = {
    HEX: [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)],
    PRISM: [(0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
    TET: [(1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)],
    PYRAMID: [(2, 3, 4), (3, 0, 4), (0, 1, 4), (1, 2, 4), (0, 1, 2, 3)],
}


class JcfError(ValueError):
    pass


# ----------------------------------------------------------------------------- mesh
@dataclass
class Mesh:
    node_ids: np.ndarray    # (N,) int, sorted
    xyz: np.ndarray         # (N, 3) float, metres
    elem_ids: np.ndarray    # (E,) int
    code: np.ndarray        # (E,) element code
    part: np.ndarray        # (E,) part id (== material_id)
    conn: np.ndarray        # (E, 8) node ids, 0-padded

    @property
    def nn(self) -> np.ndarray:
        lut = np.zeros(16, np.int64)
        for k, v in _NODES_PER.items():
            lut[k] = v
        return lut[self.code]


def mesh_from_data(nodes: dict, elements, unit: str = "mm") -> Mesh:
    """nodes {id: (x, y, z)} in ``unit``; elements [(id, element code, part, node ids)]."""
    if unit not in UNITS:
        raise JcfError(f"unit must be one of {sorted(UNITS)}")
    if not elements:
        raise JcfError("no volume elements")
    used = sorted({n for e in elements for n in e[3]})
    missing = [n for n in used if n not in nodes]
    if missing:
        raise JcfError(f"elements refer to missing nodes, e.g. {missing[:5]}")
    nid = np.array(sorted(nodes), np.int64)
    xyz = np.array([nodes[i] for i in nid], float) * UNITS[unit]
    ne = len(elements)
    eid = np.empty(ne, np.int64); code = np.empty(ne, np.int64); part = np.empty(ne, np.int64)
    conn = np.zeros((ne, 8), np.int64)
    for k, (i, c, pid, g) in enumerate(elements):
        if int(c) not in _NODES_PER or len(g) != _NODES_PER[int(c)]:
            raise JcfError(f"element {i}: code {c} with {len(g)} nodes is not supported")
        eid[k], code[k], part[k] = i, c, pid
        conn[k, :len(g)] = g
    if len(set(eid.tolist())) != ne:
        raise JcfError("duplicate element ids")
    return Mesh(nid, xyz, eid, code, part, conn)


def _node_index(m: Mesh, ids) -> np.ndarray:
    order = np.argsort(m.node_ids, kind="stable")
    pos = np.searchsorted(m.node_ids, ids, sorter=order)
    return order[np.clip(pos, 0, len(order) - 1)]


def element_volumes(m: Mesh) -> np.ndarray:
    """Exact volume of each element (outward faces fanned from their centroids)."""
    faces = {TET: [[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]],
             PRISM: [[0, 2, 1], [3, 4, 5], [0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5]],
             HEX: [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]],
             PYRAMID: [[0, 3, 2, 1], [0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]]}
    tet = lambda a, b, c, d: np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0
    out = np.zeros(len(m.elem_ids))
    for code in np.unique(m.code):
        sel = np.nonzero(m.code == code)[0]
        p = m.xyz[_node_index(m, m.conn[sel, :_NODES_PER[int(code)]])]
        c = p.mean(axis=1)
        v = np.zeros(len(sel))
        for f in faces[int(code)]:
            fc = p[:, f].mean(axis=1)
            for i in range(len(f)):
                v += tet(c, fc, p[:, f[i]], p[:, f[(i + 1) % len(f)]])
        out[sel] = v
    return out


def _write_mesh_bin(m: Mesh) -> bytes:
    hdr = b"m".ljust(4, b"\x00") + b"\x00" * 40
    out = [b"MJ\x00\x00", struct.pack("<4i", 1, len(m.node_ids), len(m.elem_ids), (len(hdr) - 40) // 4), hdr]
    nd = np.zeros(len(m.node_ids), np.dtype([("id", "<i4"), ("xyz", "<f8", 3)]))
    nd["id"], nd["xyz"] = m.node_ids, m.xyz
    out.append(nd.tobytes())
    nn = m.nn
    rows = [np.array([m.elem_ids[i], _TEXT_TO_BIN[int(m.code[i])], 1, m.part[i], *m.conn[i, :nn[i]]], "<i4")
            for i in range(len(m.elem_ids))]
    out.append(np.concatenate(rows).tobytes() if rows else b"")
    return b"".join(out)


def _write_entity_groups(groups) -> bytes:
    lines = ["%10d%10d" % (len(groups), 1)]
    for key, typ, gid, data in groups:
        data = np.asarray(data, np.int64).reshape(len(data), -1)
        cnt, width = data.shape
        lines.append("#")
        lines.append("%10d%10d%10d%10d%10d" % (key, typ, cnt, width, gid))
        if width == 1:
            flat = data.ravel()
            for i in range(0, len(flat), 8):
                lines.append("".join("%10d" % x for x in flat[i:i + 8]))
        else:
            lines.extend("".join("%10d" % x for x in row) for row in data)
    return ("\n".join(lines) + "\n").encode("ascii")


# ----------------------------------------------------------------------------- faces
def _faces(m: Mesh):
    """Every local face: (element index, local face 1-based, sorted node ids padded with -1)."""
    E, K, KEY = [], [], []
    for code, fl in LOCAL_FACES.items():
        sel = np.nonzero(m.code == code)[0]
        if not len(sel):
            continue
        conn = m.conn[sel]
        for k, f in enumerate(fl, 1):
            key = np.sort(conn[:, list(f)], axis=1)
            if key.shape[1] == 3:
                key = np.hstack([np.full((len(sel), 1), -1, key.dtype), key])
            E.append(sel); K.append(np.full(len(sel), k)); KEY.append(key)
    return np.concatenate(E), np.concatenate(K), np.vstack(KEY)


def _outward_normal(m: Mesh, ei: int, kk: int) -> np.ndarray:
    f = LOCAL_FACES[int(m.code[ei])][kk - 1]
    P = m.xyz[_node_index(m, m.conn[ei, :int(m.nn[ei])])]
    Q = P[list(f)]
    n = np.cross(Q[1] - Q[0], Q[2] - Q[0]); n /= np.linalg.norm(n)
    return -n if n @ (Q.mean(0) - P.mean(0)) < 0 else n


def _unit(v, what) -> np.ndarray:
    """``v`` normalized; scaled by its largest component first so that no finite input
    overflows or underflows the norm."""
    v = np.asarray(v, float)
    if v.shape != (3,) or not np.all(np.isfinite(v)) or not np.any(v):
        raise JcfError(f"{what} must be a finite nonzero 3-vector")
    v = v / np.abs(v).max()
    return v / np.linalg.norm(v)


def _flowin_faces_by_part(m: Mesh, coil_parts: set, fp: int, normal) -> np.ndarray:
    """Faces between part ``fp`` and the coil's other parts, owned by the other part's element,
    whose outward normal points along ``normal``."""
    v = _unit(normal, "flowin normal")
    e, k, key = _faces(m)
    _u, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    inv = inv.ravel()
    shared = cnt[inv] == 2
    part = m.part[e]
    from_keys = np.unique(inv[shared & (part == fp)])
    own = shared & np.isin(part, list(coil_parts - {fp})) & np.isin(inv, from_keys)
    rows = [(int(m.elem_ids[e[i]]), int(k[i])) for i in np.nonzero(own)[0]
            if _outward_normal(m, int(e[i]), int(k[i])) @ v > 0.5]
    if not rows:
        raise JcfError(f"no flow-in faces between part {fp} and the coil along {list(normal)}")
    return np.array(rows, np.int64)


def _flowin_faces_by_nodes(m: Mesh, coil_parts: set, face_nodes, normal, label) -> np.ndarray:
    """Coil element faces made of the given node sets whose outward normal points along
    ``normal``. On a face inside the coil this is the element the current flows out of."""
    v = _unit(normal, f"flowin normal of {label!r}")
    e, k, keys = _faces(m)
    inside = np.isin(m.part[e], list(coil_parts))
    by_key = collections.defaultdict(list)
    for i in np.flatnonzero(inside):
        by_key[tuple(int(n) for n in keys[i] if n != -1)].append(i)
    rows = []
    for nodes in face_nodes:
        cand = [i for i in by_key.get(tuple(sorted(int(n) for n in nodes)), ())
                if _outward_normal(m, int(e[i]), int(k[i])) @ v > 0.5]
        if len(cand) != 1:
            raise JcfError(f"flow-in face {sorted(nodes)} of {label!r}: {len(cand)} coil element faces "
                           f"with outward normal along {list(normal)}")
        rows.append((int(m.elem_ids[e[cand[0]]]), int(k[cand[0]])))
    if not rows or len(set(rows)) != len(rows):
        raise JcfError(f"flow-in faces of {label!r} must be nonempty and unique")
    return np.array(rows, np.int64)


def _exterior_faces(m: Mesh, face_nodes=None, label="") -> np.ndarray:
    e, k, keys = _faces(m)
    _u, inv, cnt = np.unique(keys, axis=0, return_inverse=True, return_counts=True)
    ext = np.flatnonzero(cnt[inv.ravel()] == 1)
    if face_nodes is None:
        return np.stack([m.elem_ids[e[ext]], k[ext]], axis=1).astype(np.int64)
    faces = {tuple(int(n) for n in keys[i] if n != -1): (int(m.elem_ids[e[i]]), int(k[i])) for i in ext}
    rows = []
    for nodes in face_nodes:
        key = tuple(sorted(int(n) for n in nodes))
        if key not in faces:
            raise JcfError(f"face {list(key)} of {label!r} is not an exterior mesh face")
        rows.append(faces[key])
    if not rows or len(set(rows)) != len(rows):
        raise JcfError(f"faces of {label!r} must be nonempty and unique")
    return np.array(rows, np.int64)


# ----------------------------------------------------------------------------- root.xml
def _fmt(v) -> str:
    if isinstance(v, np.generic):
        v = v.item()
    if isinstance(v, float):
        r = repr(v)
        return r[:-2] if r.endswith(".0") else r
    return str(v)


def _hex_table(a) -> str:
    return "".join(struct.pack(">d", float(v)).hex().upper()[::-1] for v in np.ravel(a))


def _p(name, value, unit=None, typ=None):
    u = f' unit="{unit}"' if unit is not None else ""
    t = f' type="{typ}"' if typ else ""
    return f'\t\t\t<parameter name="{name}"{u}{t}>{escape(_fmt(value))}</parameter>\n'


def _table(name, rows, column_units=""):
    rows = np.atleast_2d(np.asarray(rows, float))
    body = "".join("\t\t\t\t" + " ".join(_fmt(float(v)) for v in r) + " \n" for r in rows)
    return (f'\t\t\t<parameter rows="{rows.shape[0]}" cols="{rows.shape[1]}" encoding="0" '
            f'column_units="{column_units}" name="{name}" type="table">\n{body}\t\t\t</parameter>\n')


def _point(name, v):
    x, y, z = (_fmt(float(a)) for a in v)
    return (f'\t\t\t<parameter name="{name}" type="point" unit="">\n\t\t\t\t<x>{x}</x>\n'
            f'\t\t\t\t<y>{y}</y>\n\t\t\t\t<z>{z}</z>\n\t\t\t</parameter>\n')


CONTROL_DEFAULTS = {"circuit_divide_multiplier": 1, "flg_auto_accel": 1, "flg_minimize_distribution": 1,
                    "iccg_tolerance": 1e-08, "ralaxation_factor_type": 5, "nonlinear_max_iteration": 15,
                    "flg_io_permeance_factor": 1, "flg_io_vector_potential": 1}
_INTENSITY = {"circuit_divide_multiplier", "iccg_tolerance"}


def magnet_frame(direction):
    """(x axis, xy-plane axis) of a right-handed system whose z axis is ``direction``: the
    rotation that takes +z to it along the shortest arc. For +x: x = (0, 0, -1), y = (0, 1, 0)."""
    d = _unit(direction, "magnet direction")
    return _shortest_arc_frame(d)


def _shortest_arc_frame(d):
    """Rodrigues rotation R = I + K + K^2 / (1 + c) with K = [z x d]_x, c = d_z. Near -z,
    1 + c cancels, so it is taken as s^2 / (1 - c) with s = |z x d| from the transverse
    components. Only an antiparallel d (s^2 below the smallest normal float) uses the fixed
    half turn about x."""
    c = float(d[2])
    a = np.array([-d[1], d[0], 0.0])                     # z x d
    s2 = float(d[0] * d[0] + d[1] * d[1])
    if c < 0 and s2 < np.finfo(float).tiny:              # -z: half turn about x
        return np.array([1.0, 0, 0]), np.array([0, -1.0, 0])
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    R = np.eye(3) + K + (K @ K / (1 + c) if c >= 0 else K @ K * ((1 - c) / s2))
    clean = lambda v: np.where(np.abs(v) < 1e-15, 0.0, v)
    return clean(R[:, 0]), clean(R[:, 1])


def _cs_xml(index, x, xy, unit) -> str:
    body = "".join(f'\t\t\t<point name="{n}">\n\t\t\t\t<x>{_fmt(float(a))}</x>\n\t\t\t\t<y>{_fmt(float(b))}</y>\n'
                   f'\t\t\t\t<z>{_fmt(float(c))}</z>\n\t\t\t</point>\n'
                   for n, (a, b, c) in (("position_center", (0, 0, 0)), ("position_on_x", x), ("position_on_xy", xy)))
    if unit:
        body += f'\t\t\t<parameter name="unit">{unit}</parameter>\n'
    return f'\t\t<coordinate_system index="{index}" type="1">\n{body}\t\t</coordinate_system>\n'


def _is_plus_z(d) -> bool:
    return bool(np.allclose(_unit(d, "magnet direction"), (0, 0, 1), atol=1e-12))


def _material_category(m: dict) -> int:
    if m.get("magnet"):
        return 5
    if m.get("bh") is not None:
        return 6
    if m.get("resistivity") is not None:
        return 3
    return 4 if float(m.get("mu_r", 1)) != 1 else 1


def _elem_group_xml(key, gid, n, ftype="elements", caption=""):
    return (f'\t\t<group type="{ftype}" key="{key}" caption={quoteattr(caption)} groupId="{gid}" '
            f'size="{n}" external="1"')


def build_root_xml(spec: dict, parts: set, groups: list) -> bytes:
    mats = spec["materials"]
    ids = [int(m["part"]) for m in mats]
    if len(set(ids)) != len(ids) or set(ids) != parts:
        raise JcfError(f"materials must cover each mesh part {sorted(parts)} once; got {sorted(ids)}")
    if spec.get("study", "static") != "static":
        raise JcfError("only a static study is supported")
    u = spec.get("unit", "mm")
    lu = u if u in ("m", "mm", "cm") else None
    out = ['<?xml version="1.0" encoding="UTF-8"?>\n<jmag>\n']
    if groups:
        out.append("\t<group_data>\n")
        for key, gid, data, kind, x in groups:
            if kind in ("face", "boundary"):
                out.append(_elem_group_xml(key, gid, len(data), "faces", x[1] if kind == "boundary" else "") + " />\n")
            elif kind == "magnet":
                rows = "".join("\t\t\t\t0 0 0 \n" for _ in data)
                out.append(_elem_group_xml(key, gid, len(data)) + f'>\n\t\t\t<group_table rows="{len(data)}" '
                           f'cols="3" encoding="0" column_units="" name="elemental_magnetization_table">\n'
                           f'{rows}\t\t\t</group_table>\n\t\t</group>\n')
            else:
                out.append(_elem_group_xml(key, gid, len(data)) + " />\n")
        out.append("\t</group_data>\n")
    out.append('\t<material_data module="jiba">\n')
    for m in mats:
        name = escape(str(m.get("name") or f"part{m['part']}"), {'"': "&quot;"})
        out.append(f'\t\t<material name="{name}">\n' + _p("material_id", int(m["part"])))
        mag = m.get("magnet")
        if mag:                          # linear magnet from (H=-Hcb, B=0) to (H=0, B=Br)
            br, hcb = float(mag["br"]), float(mag["hcb"])
            out.append(_p("flg_exist_magnet", 1))
            out.append(f'\t\t\t<parameter rows="2" cols="2" encoding="1" column_units="A/m;tesla" '
                       f'name="mag_table_1" type="table">{_hex_table([[-hcb, 0.0], [0.0, br]])}</parameter>\n')
            out.append(_p("flg_magnetic_property1", 7) + _p("flg_io_material", 1) + _p("flg_use_anisotropic_magnet", 1))
        elif m.get("bh") is not None:    # nonlinear B-H rows (H [A/m], B [T]) from (0, 0)
            bh = np.asarray(m["bh"], float)
            lf = m.get("lamination_factor")
            out.append(_p("flg_anisotropy", 1) + _p("flg_magnetic_property1", 7))
            if lf is not None:
                out.append(_p("flg_use_lamination_factor", 1))
            out.append(_p("flg_io_material", 1))
            out.append(f'\t\t\t<parameter rows="{bh.shape[0]}" cols="2" encoding="1" column_units="A/m;tesla" '
                       f'name="mag_table_1" type="table">{_hex_table(bh)}</parameter>\n')
            if lf is not None:
                out.append(_p("lamination_factor", lf, unit="percent"))
        else:
            out.append(_p("permeability_1", f"({_fmt(float(m.get('mu_r', 1)))}, 0)", unit="Intensity"))
            out.append(_p("flg_magnetic_property1", 1))
            if m.get("output"):
                out.append(_p("flg_io_material", 1))
        if m.get("resistivity") is not None:
            out.append(_p("flg_conductivity_unit", 1) + _p("flg_use_conductivity", 0)
                       + _p("conductivity_1", float(m["resistivity"]), unit="ohm_m") + _p("flg_electric_property1", 0))
        out.append(_p("material_category", _material_category(m)))   # importers need it
        out.append("\t\t</material>\n")
    out.append('\t</material_data>\n\t<control_parameters module="mesh">\n\t</control_parameters>\n')
    out.append("\t<condition_data>\n")
    coils = [g for g in groups if g[3] == "coil"]
    ctl = dict(CONTROL_DEFAULTS)
    if coils:
        ctl["flg_io_current_density"] = 1
    over = spec.get("control") or {}
    unknown = set(over) - set(CONTROL_DEFAULTS)
    if unknown:
        raise JcfError(f"control keys {sorted(unknown)} are not supported")
    ctl.update(over)
    out.append('\t\t<condition name="Control">\n')
    for k, v in ctl.items():
        out.append(_p(k, v, unit="Intensity" if k in _INTENSITY else None))
    out.append("\t\t</condition>\n")
    out.append('\t\t<condition name="Step">\n' + _p("step_type", 1) + _p("num_analysis_step", 1)
               + _table("step_table", [[0, 0, 1], [1, 0, 1]]) + "\t\t</condition>\n")
    ef = spec.get("external_field")
    if ef:
        table = ef.get("table") or [[0, ef["value"], 0], [1000, ef["value"], 0]]
        out.append('\t\t<condition name="ExternalField">\n'
                   + _p("ring_radius", ef.get("ring_radius", 100000), unit="m")
                   + _p("num_division", ef.get("num_division", 120))
                   + _point("axis_vector", ef.get("axis", [0, 0, 1]))
                   + ("" if ef.get("table") else _p("function_type", 0)
                      + _p("designer_const_amplitude", ef["value"], unit="tesla"))
                   + _table("external_field_table", table, ";tesla") + "\t\t</condition>\n")
    extra_cs = []
    for i, (key, _gid, _e, _k, d) in enumerate([g for g in groups if g[3] == "magnet"], 1):
        lcs = ""
        if not _is_plus_z(d):
            idx = 2 + len(extra_cs)
            extra_cs.append(_cs_xml(idx, *magnet_frame(d), lu))
            lcs = _p("local_coordinate_id", idx)
        out.append('\t\t<condition name="PermanentMagnet">\n' + _p("pattern_type", 1)
                   + _point("axis2", d) + _p("magnetization", 1) + _p("id", i)
                   + _p("angle_from_z_axis", 0) + lcs
                   + f'\t\t\t<group type="elements" key="{key}" />\n\t\t</condition>\n')
    for i, (key, _gid, _e, _k, d) in enumerate([g for g in groups if g[3] == "steel"], 1):
        out.append('\t\t<condition name="Anisotropy">\n' + _point("axis2", d) + _p("coordinate_type", 1)
                   + _p("local_coordinate_id", 1) + _p("id", i)
                   + f'\t\t\t<group type="elements" key="{key}" />\n\t\t</condition>\n')
    for i, (key, _gid, _e, _k, (c, fkey, scale)) in enumerate(coils, 1):
        amp = float(c["current"])
        body = (_p("turn", c.get("turns", 1), unit="turn")
                + _table("current_table", [[0, amp, 0], [1000, amp, 0]], ";A") + _p("id", i))
        if fkey is None:
            o = [float(v) * scale for v in c["origin"]]
            body += (_p("flowin_type", 1) + _point("flowin_vector", c["direction"])
                     + _point("flowin_origin", o).replace('unit=""', 'unit="m"'))
        else:                            # explicit: readers differ in the default
            body += _p("flowin_type", 0)
        body += f'\t\t\t<group type="elements" key="{key}" />\n'
        if fkey is not None:
            body += f'\t\t\t<group type="faces" key="{fkey}" />\n'
        out.append('\t\t<condition name="CurrentCtrl">\n' + body + "\t\t</condition>\n")
    for key, _gid, _e, kind, extra in groups:
        if kind == "boundary":
            out.append(f'\t\t<condition name="{extra[0]}">\n\t\t\t<group type="faces" key="{key}" />\n'
                       '\t\t</condition>\n')
    out.append("\t</condition_data>\n")
    out.append("\t<coordinate_system_data>\n" + _cs_xml(1, (1, 0, 0), (0, 1, 0), lu)
               + "".join(extra_cs) + "\t</coordinate_system_data>\n")
    out.append('\t<control_parameters module="common">\n')
    com = {"max_plot_size": 1000, "version": 191, "unit_input": 1, "plot_mode": 1}
    if lu:                               # length unit of the model
        com["unit_label"] = lu
        com["unit_output"] = UNITS[lu]
    for k, v in com.items():
        out.append(_p(k, v).replace("\t\t\t", "\t\t", 1))
    out.append("\t</control_parameters>\n</jmag>\n")
    return "".join(out).encode("utf-8")


# ----------------------------------------------------------------------------- groups
_BOUNDARY = {"symmetry": "SymmetryBoundary3D"}
_MATERIAL_KEYS = {"part", "name", "mu_r", "magnet", "bh", "lamination_factor", "resistivity", "output"}


def _groups(m: Mesh, spec: dict) -> list:
    groups = []
    mats = {int(x["part"]): x for x in spec["materials"]}
    for x in spec["materials"]:
        unknown = set(x) - _MATERIAL_KEYS
        if unknown:
            raise JcfError(f"material {x.get('name')!r}: unsupported keys {sorted(unknown)}")
        if x.get("magnet") and x.get("bh") is not None:
            raise JcfError(f"material {x.get('name')!r}: give either 'bh' or 'magnet', not both")
    for mg in spec.get("magnets") or []:
        part = int(mg["part"])
        if not (mats.get(part) or {}).get("magnet"):
            raise JcfError(f"magnet part {part} needs a material with 'magnet': {{br, hcb}}")
        el = m.elem_ids[m.part == part]
        if not len(el):
            raise JcfError(f"magnet part {part} has no elements")
        d = mg.get("direction", [0, 0, 1])
        _unit(d, f"magnet direction of part {part}")
        groups.append((len(groups) + 1, len(groups) + 1, el, "magnet", [float(v) for v in d]))
    for x in spec["materials"]:          # every B-H material gets an Anisotropy group
        if x.get("bh") is None:
            continue
        bh = np.asarray(x["bh"], float)
        if bh.ndim != 2 or bh.shape[1] != 2 or len(bh) < 2 or not np.all(np.diff(bh, axis=0) > 0):
            raise JcfError(f"material {x.get('name')!r}: 'bh' must be rows (H, B), both strictly increasing")
        el = m.elem_ids[m.part == int(x["part"])]
        if not len(el):
            raise JcfError(f"B-H part {x['part']} has no elements")
        groups.append((len(groups) + 1, len(groups) + 1, el, "steel", [0, 0, 1]))
    scale = UNITS[spec.get("unit", "mm")]
    for c in spec.get("coils") or []:
        cp = {int(p) for p in c.get("parts") or []}
        if not cp or not cp <= set(mats):
            raise JcfError(f"coil parts {sorted(cp)} must be mesh parts with materials")
        if any(mats[p].get("resistivity") is None for p in cp):
            raise JcfError("every coil part needs a material 'resistivity' [ohm m]")
        if "current" not in c:
            raise JcfError("a coil needs 'current' [A]")
        el = m.elem_ids[np.isin(m.part, list(cp))]
        fi = c.get("flowin") or {}
        fkey = None
        if fi.get("part") is not None or fi.get("face_nodes") is not None:
            if fi.get("face_nodes") is not None:
                faces = _flowin_faces_by_nodes(m, cp, fi["face_nodes"], fi.get("normal"), fi.get("label", ""))
            else:
                if int(fi["part"]) not in cp:
                    raise JcfError("flowin part must be one of the coil parts")
                faces = _flowin_faces_by_part(m, cp, int(fi["part"]), fi.get("normal"))
            groups.append((len(groups) + 1, len(groups) + 1, faces, "face", None))
            fkey = len(groups)
        elif c.get("direction") is None or c.get("origin") is None:
            raise JcfError("a coil needs a flow-in face or 'direction' and 'origin'")
        groups.append((len(groups) + 1, len(groups) + 1, el, "coil", (c, fkey, scale)))
    seen = set()
    for b in spec.get("boundaries") or []:
        cname = _BOUNDARY.get(b.get("type"))
        if cname is None:
            raise JcfError(f"boundary type {b.get('type')!r} is not supported (only 'symmetry')")
        faces = _exterior_faces(m, b.get("face_nodes"), b.get("label", ""))
        rows = {tuple(r) for r in faces.tolist()}
        if seen & rows:
            raise JcfError("symmetry boundary groups overlap")
        seen |= rows
        groups.append((len(groups) + 1, len(groups) + 1, faces, "boundary", (cname, b.get("label", ""))))
    return groups


def write_jcf(mesh: Mesh, out, spec: dict) -> dict:
    """Write ``out`` (.jcf, must not exist) from ``mesh`` and a numeric spec (parts by id)."""
    out = Path(out)
    if out.suffix.lower() != ".jcf":
        raise JcfError("output must end in .jcf")
    if out.exists():
        raise JcfError(f"refusing to overwrite {out}")
    vol = element_volumes(mesh)
    if (vol <= 0).any():
        bad = mesh.elem_ids[vol <= 0][:5].tolist()
        raise JcfError(f"{int((vol <= 0).sum())} elements with non-positive volume, e.g. {bad}")
    groups = _groups(mesh, spec)
    root = build_root_xml(spec, set(mesh.part.tolist()), groups)
    out.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=out.name + ".", suffix=".part", dir=out.parent)
    os.close(fd)
    tmp = Path(tmp)
    try:
        with zipfile.ZipFile(tmp, "w") as z:
            # stored, not deflated: a deflated mesh.bin makes the solver report no mesh data
            z.writestr("mesh.bin", _write_mesh_bin(mesh), compress_type=zipfile.ZIP_STORED)
            z.writestr("root.xml", root, compress_type=zipfile.ZIP_DEFLATED)
            if groups:
                eg = [(key, 2 if kind in ("face", "boundary") else 1, gid, d) for key, gid, d, kind, _x in groups]
                z.writestr("entity_groups.dat", _write_entity_groups(eg), compress_type=zipfile.ZIP_DEFLATED)
        try:
            os.link(tmp, out)                # never replaces an existing file
        except FileExistsError:
            raise JcfError(f"refusing to overwrite {out}") from None
        except OSError:
            if os.name != "nt":
                raise
            try:
                os.rename(tmp, out)
            except FileExistsError:
                raise JcfError(f"refusing to overwrite {out}") from None
    finally:
        tmp.unlink(missing_ok=True)
    return {"ok": True, "path": str(out), "nodes": len(mesh.node_ids), "elements": len(mesh.elem_ids),
            "by_code": {int(k): v for k, v in collections.Counter(mesh.code.tolist()).items()},
            "parts": sorted(set(mesh.part.tolist())), "volume_m3": float(vol.sum()),
            "groups": [(kind, len(d)) for _k, _g, d, kind, _x in groups]}
