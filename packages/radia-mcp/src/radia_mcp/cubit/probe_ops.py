"""
probe_ops.py — SHARED Cubit-side probe queries for both session runners.

The GUI file-drop runner (bootstrap.py, executed inside Cubit's Qt/Python)
and the batch stdio runner (daemon.py, executed by Cubit's bundled
Python 3.10) both dispatch the "probe" op here, so the query surface can
never drift between the two transports again (it had: `per_volume`
existed only in bootstrap, `entities`/`labels` only in daemon).

Both runners execute as plain scripts (no package context), so this
module is imported BY PATH (the runner inserts its own directory into
``sys.path``).  Keep it importable under Python 3.10 and free of any
import beyond the stdlib; the ``cubit`` module is passed in as an
argument, never imported here.
"""

import os
import sys
import traceback
import uuid

_HERE = os.path.dirname(os.path.abspath(__file__))

_ENTITY_DUMP_CAP = 300


def trimmed_traceback(exclude_basenames=("daemon.py", "bootstrap.py",
                                         "probe_ops.py")) -> str:
    """Format the ACTIVE exception with runner-harness frames removed.

    The model debugging a failed cubit command should see the frames of
    the command/script it controls, not the runner's dispatch loop
    (MathWorks CustomStackException.trimmedBefore pattern).  Falls back
    to the full traceback when trimming would leave nothing.
    """
    etype, evalue, tb = sys.exc_info()
    if etype is None:
        return ""
    frames = traceback.extract_tb(tb)
    kept = [f for f in frames
            if os.path.basename(f.filename) not in exclude_basenames]
    if not kept:
        return "".join(traceback.format_exception(etype, evalue, tb))
    out = ["Traceback (most recent call last):\n"]
    out.extend(traceback.format_list(kept))
    out.extend(traceback.format_exception_only(etype, evalue))
    return "".join(out)


def _round3(seq, ndigits=6):
    return [round(float(v), ndigits) for v in seq]


def _bbox_fields(bb):
    """Cubit's GeomEntity.bounding_box() returns a 6-tuple
    (min_x, min_y, min_z, max_x, max_y, max_z) -- VERIFIED on Coreform
    Cubit 2025.12 (brick x 1 -> [-0.5]*3 + [0.5]*3).  Extent is derived;
    a flat face has zero extent in its normal direction."""
    lo, hi = bb[0:3], bb[3:6]
    return {
        "bbox_min": _round3(lo),
        "bbox_max": _round3(hi),
        "extent": _round3([h - l for l, h in zip(lo, hi)]),
    }


def _probe_entities(cubit_mod):
    """Per-entity geometric dump implementing the lab's Probe-Don't-Guess
    policy: every volume/surface with centroid, bbox, extent, and measure,
    so entity-classification predicates are derived from printed values
    rather than a-priori reasoning about the .jou source."""
    vol_ids = [int(v) for v in cubit_mod.parse_cubit_list("volume", "all")]
    surf_ids = [int(s) for s in cubit_mod.parse_cubit_list("surface", "all")]
    result = {
        "volume_count": len(vol_ids),
        "surface_count": len(surf_ids),
        "volumes": [],
        "surfaces": [],
    }
    for vid in vol_ids[:_ENTITY_DUMP_CAP]:
        v = cubit_mod.volume(vid)
        entry = {
            "id": vid,
            "centroid": _round3(v.centroid()),
            "volume": float(v.volume()),
        }
        entry.update(_bbox_fields(v.bounding_box()))
        result["volumes"].append(entry)
    for sid in surf_ids[:_ENTITY_DUMP_CAP]:
        s = cubit_mod.surface(sid)
        cx, cy, cz = s.center_point()   # Surface API: NOT .centroid()
        entry = {
            "id": sid,
            "center": _round3((cx, cy, cz)),
            "area": float(s.area()),
        }
        entry.update(_bbox_fields(s.bounding_box()))
        result["surfaces"].append(entry)
    if len(vol_ids) > _ENTITY_DUMP_CAP or len(surf_ids) > _ENTITY_DUMP_CAP:
        result["truncated"] = (
            f"listing capped at {_ENTITY_DUMP_CAP} entities per kind; "
            "use parse_cubit_list + per-entity calls via cubit_exec for the rest"
        )
    return result


def _probe_labels(cubit_mod):
    """Blocks/sidesets with names + convention audit (phantom/mixed-block
    detection before `export netgen`).

    VERIFIED on Coreform Cubit 2025.12: `block N add <other-kind>` on an
    already-typed block returns success but adds NOTHING (silent no-op),
    and `block N add tri ...` on a hex mesh creates NO block at all
    (phantom).  The primary value of this probe is therefore verifying
    ACTUAL block membership after every block/sideset command -- an
    intended-but-absent element set shows up here as missing counts."""
    def _ids_in(kind: str, scope: str, eid: int) -> set[int]:
        try:
            return {int(item) for item in cubit_mod.parse_cubit_list(
                kind, f"in {scope} {eid}")}
        except Exception:
            return set()

    blocks = []
    for bid in cubit_mod.get_block_id_list():
        bid = int(bid)
        volumes = [int(v) for v in cubit_mod.get_block_volumes(bid)]
        surfaces = [int(s) for s in cubit_mod.get_block_surfaces(bid)]
        # EFFECTIVE element counts.  Cubit binds elements to a block two
        # ways: DIRECT element membership (get_block_hexes/... -- e.g.
        # `block N add face in surface S`) and GEOMETRY membership
        # (`block N add volume V` -- get_block_hexes returns 0 there
        # even for a fully meshed volume, VERIFIED Cubit 2025.12
        # 2026-08-05: 125-hex brick in a volume block reported 0).  Take
        # the UNION so a block that has both forms of membership does not
        # double-count an element.
        vol_ids = {
            int(item)
            for getter in (cubit_mod.get_block_hexes,
                           cubit_mod.get_block_tets,
                           cubit_mod.get_block_wedges,
                           cubit_mod.get_block_pyramids)
            for item in getter(bid)
        }
        for vid in volumes:
            for kind in ("hex", "tet", "wedge", "pyramid"):
                vol_ids.update(_ids_in(kind, "volume", vid))
        surf_ids = {
            int(item)
            for getter in (cubit_mod.get_block_tris,
                           cubit_mod.get_block_faces)
            for item in getter(bid)
        }
        for sid in surfaces:
            for kind in ("face", "tri"):
                surf_ids.update(_ids_in(kind, "surface", sid))
        blocks.append({
            "id": bid,
            "name": str(cubit_mod.get_exodus_entity_name("block", bid)),
            "volumes": volumes,
            "surfaces": surfaces,
            "volume_elems": len(vol_ids),
            "surface_elems": len(surf_ids),
        })
    sidesets = []
    for sid in cubit_mod.get_sideset_id_list():
        sid = int(sid)
        sidesets.append({
            "id": sid,
            "name": str(cubit_mod.get_exodus_entity_name("sideset", sid)),
            "surfaces": [int(s) for s in cubit_mod.get_sideset_surfaces(sid)],
        })
    # Runners execute standalone (no package context); import the shared
    # cubit-free audit module by path.
    if _HERE not in sys.path:
        sys.path.insert(0, _HERE)
    from label_audit import audit_label_records
    audit = audit_label_records(blocks, sidesets)
    return {"blocks": blocks, "sidesets": sidesets, "audit": audit}


def _probe_quality(cubit_mod):
    """Scaled-jacobian stats over all hexes+tets: min / max / mean /
    count below 0.2 -- small, stable proxy for mesh quality without a
    full export."""
    try:
        hex_ids = list(cubit_mod.parse_cubit_list("hex", "all"))
        tet_ids = list(cubit_mod.parse_cubit_list("tet", "all"))
    except Exception:
        hex_ids, tet_ids = [], []
    vals = []
    try:
        if hex_ids:
            vals.extend(float(v) for v in cubit_mod.get_quality_values(
                "hex", [int(x) for x in hex_ids], "scaled jacobian"))
    except Exception:
        pass
    try:
        if tet_ids:
            vals.extend(float(v) for v in cubit_mod.get_quality_values(
                "tet", [int(x) for x in tet_ids], "scaled jacobian"))
    except Exception:
        pass
    if not vals:
        return {"hex_count": len(hex_ids), "tet_count": len(tet_ids),
                "min": None, "max": None, "mean": None,
                "below_0.2": 0}
    below = sum(1 for v in vals if v < 0.2)
    return {
        "hex_count": int(len(hex_ids)),
        "tet_count": int(len(tet_ids)),
        "min": round(min(vals), 6),
        "max": round(max(vals), 6),
        "mean": round(sum(vals) / len(vals), 6),
        "below_0.2": int(below),
        "below_0.2_pct": round(100 * below / len(vals), 2),
    }


def _probe_per_volume(cubit_mod):
    """One row per volume: id / meshing scheme / hex / tet / meshed."""
    try:
        vol_ids = list(cubit_mod.parse_cubit_list("volume", "all"))
    except Exception:
        vol_ids = []
    rows = []
    for vid in vol_ids:
        vid = int(vid)
        try:
            scheme = cubit_mod.get_volume_meshing_scheme(vid)
        except Exception:
            scheme = None
        try:
            n_hex = len(cubit_mod.parse_cubit_list("hex", f"in volume {vid}"))
        except Exception:
            n_hex = 0
        try:
            n_tet = len(cubit_mod.parse_cubit_list("tet", f"in volume {vid}"))
        except Exception:
            n_tet = 0
        rows.append({
            "id": vid,
            "scheme": str(scheme) if scheme is not None else None,
            "hex": int(n_hex),
            "tet": int(n_tet),
            "meshed": (n_hex + n_tet) > 0,
        })
    return rows


def op_snapshot(cubit_mod, args):
    """Hardcopy the current view to PNG -- SHARED by both runners.

    Uses the plain ``hardcopy "<path>" png`` form: the ``window <w> <h>``
    variant returns rc=1 but writes a 0-byte file in GUI mode (measured
    Coreform Cubit 2025.12, 2026-08-05).  Any requested size is reported
    back as ignored rather than silently dropped; the image renders at
    the current graphics-window size.  The written file is verified
    (existence + non-zero size) because ``cubit.cmd`` returns success
    for silently-failing hardcopies.
    """
    if len(args) < 1:
        return {"error": "snapshot requires at least 1 arg (path)"}
    out_path = str(args[0])
    parent = os.path.dirname(os.path.abspath(out_path))
    os.makedirs(parent, exist_ok=True)
    tmp_path = os.path.join(
        parent,
        f".{os.path.basename(out_path)}.{uuid.uuid4().hex}.tmp.png",
    )
    try:
        try:
            rc = cubit_mod.cmd(f'hardcopy "{tmp_path}" png')
            size = (os.path.getsize(tmp_path)
                    if os.path.isfile(tmp_path) else 0)
        except OSError:
            rc, size = False, 0
        ok = bool(rc) and size > 0
        result = {"path": out_path, "ok": ok, "bytes": int(size)}
        if len(args) >= 3:
            result["requested_size_ignored"] = [int(args[1]), int(args[2])]
        if not ok:
            result["error"] = (
                "hardcopy wrote no image (0 bytes). The graphics window is "
                "not rendering -- in batch (-nographics) mode snapshots are "
                "unavailable; in GUI mode issue a draw/display command first.")
        else:
            try:
                os.replace(tmp_path, out_path)
            except OSError as exc:
                result.update(ok=False, error=f"cannot publish snapshot: {exc}")
        return result
    finally:
        try:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
        except OSError:
            pass


def op_probe(cubit_mod, args):
    """Dispatch one probe query.  The single implementation for BOTH the
    GUI file-drop runner and the batch stdio runner."""
    query = args[0] if args else ""
    q = str(query).strip().lower()
    try:
        if q in ("volume_count", "volumes"):
            return int(cubit_mod.get_volume_count())
        if q in ("surface_count", "surfaces"):
            return int(cubit_mod.get_surface_count())
        if q in ("vertex_count", "vertices"):
            return int(cubit_mod.get_vertex_count())
        if q in ("curve_count", "curves"):
            return int(cubit_mod.get_curve_count())
        if q in ("node_count", "nodes"):
            return int(cubit_mod.get_node_count())
        if q in ("hex_count", "hexes", "hex"):
            return int(cubit_mod.get_hex_count())
        if q in ("tet_count", "tets", "tet"):
            return int(cubit_mod.get_tet_count())
        if q in ("summary",):
            return {
                "volumes": int(cubit_mod.get_volume_count()),
                "surfaces": int(cubit_mod.get_surface_count()),
                "curves": int(cubit_mod.get_curve_count()),
                "vertices": int(cubit_mod.get_vertex_count()),
                "nodes": int(cubit_mod.get_node_count()),
                "hexes": int(cubit_mod.get_hex_count()),
                "tets": int(cubit_mod.get_tet_count()),
            }
        if q in ("quality_summary", "quality"):
            return _probe_quality(cubit_mod)
        if q in ("per_volume", "volumes_detail"):
            return _probe_per_volume(cubit_mod)
        if q in ("entities", "geometry"):
            return _probe_entities(cubit_mod)
        if q in ("labels", "blocks", "block_sideset"):
            return _probe_labels(cubit_mod)
    except Exception:
        return {"error": trimmed_traceback()}
    return {"error": f"Unknown probe query: {query!r}. "
                     "Try: volume_count / surface_count / vertex_count / "
                     "curve_count / node_count / hex_count / tet_count / "
                     "summary / quality / per_volume / entities / labels"}
