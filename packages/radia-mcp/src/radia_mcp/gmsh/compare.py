"""Cross-mesh field comparison: sample two .msh files at the SAME points.

``diff_msh`` (msh_inspect) is a REGRESSION check: it compares structure
and per-view statistics, so two different meshings of the same geometry
come back as "different" even when they carry the identical analytic
field.  This module answers the other question -- *do two solvers /
discretizations agree on the field?* -- by probing both files at one
shared point cloud and reporting the pointwise difference.

Measured semantics (gmsh 4.15.2, locked by tests):

- Two meshes can live in ONE gmsh session as two models
  (``gmsh.model.add`` + ``gmsh.merge``); a view keeps the model it was
  merged with, so ``gmsh.view.probe`` on a view of model A is unaffected
  by model B being current.  New view tags appear only for the file
  merged last, which is how the two files' views are told apart.
- ``gmsh.view.probe`` returns ``([], distance)`` outside the mesh, so
  points that miss either mesh are skipped and counted rather than
  silently contributing a zero.
- ``gmsh.view.probe(step=N)`` with ``N >= NbTimeStep`` is a native
  out-of-bounds read that KILLS the child process, so multi-step views
  are rejected before the first probe.
"""

from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path
from typing import Any

from ._gmsh_subprocess import run_gmsh_json_subprocess

_SAMPLE_MODES = ("random", "grid")

_COMPARE_SCRIPT = r"""
import json
import math
import random
import sys

cfg_path, out_path = sys.argv[1], sys.argv[2]
with open(cfg_path, encoding="utf-8") as f:
    cfg = json.load(f)
result = {"ok": False, "ran": False}


def _view_name(tag):
    import gmsh
    return gmsh.option.getString("View[%d].Name" % gmsh.view.getIndex(tag))


def _n_steps(tag):
    import gmsh
    return int(gmsh.option.getNumber(
        "View[%d].NbTimeStep" % gmsh.view.getIndex(tag)))


def _pick(tags, selector, label, fallback_name):
    # Explicit name/index wins; otherwise the caller-independent
    # fallback (the name shared by both files) is used.  Never a silent
    # "first view" guess when the choice is ambiguous.
    import gmsh
    names = [_view_name(t) for t in tags]
    if isinstance(selector, int) and not isinstance(selector, bool):
        if not 0 <= selector < len(tags):
            raise RuntimeError(
                "%s index %d out of range (0..%d)"
                % (label, selector, len(tags) - 1))
        return tags[selector]
    want = selector if selector is not None else fallback_name
    if want not in names:
        raise RuntimeError(
            "%s %r not found; available: %s" % (label, want, sorted(names)))
    return tags[names.index(want)]


def _points(box, n, mode, seed):
    lo, hi = box[0], box[1]
    if mode == "grid":
        k = max(2, int(round(n ** (1.0 / 3.0))))
        pts = []
        for i in range(k):
            for j in range(k):
                for m in range(k):
                    frac = ((i + 0.5) / k, (j + 0.5) / k, (m + 0.5) / k)
                    pts.append([lo[c] + (hi[c] - lo[c]) * frac[c]
                                for c in range(3)])
        return pts
    rng = random.Random(seed)
    return [[rng.uniform(lo[c], hi[c]) for c in range(3)] for _ in range(n)]


try:
    import gmsh
    gmsh.initialize(["-noconfig"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("compare_a")
        gmsh.merge(cfg["path_a"])
        tags_a = [int(t) for t in gmsh.view.getTags()]
        bbox_a = [float(v) for v in gmsh.model.getBoundingBox(-1, -1)]

        gmsh.model.add("compare_b")
        gmsh.merge(cfg["path_b"])
        seen = set(tags_a)
        tags_b = [int(t) for t in gmsh.view.getTags() if int(t) not in seen]
        bbox_b = [float(v) for v in gmsh.model.getBoundingBox(-1, -1)]

        if not tags_a:
            raise RuntimeError("no post-processing views in a: %s"
                               % cfg["path_a"])
        if not tags_b:
            raise RuntimeError("no post-processing views in b: %s"
                               % cfg["path_b"])

        names_a = [_view_name(t) for t in tags_a]
        names_b = [_view_name(t) for t in tags_b]
        common = [n for n in names_a if n in names_b]
        fallback = None
        if cfg["view_a"] is None or cfg["view_b"] is None:
            if len(common) == 1:
                fallback = common[0]
            elif not common:
                raise RuntimeError(
                    "a and b share no view name (a: %s, b: %s); pass "
                    "view_a and view_b explicitly"
                    % (sorted(set(names_a)), sorted(set(names_b))))
            else:
                raise RuntimeError(
                    "a and b share %d view names %s; pass view_a/view_b to "
                    "choose one" % (len(common), sorted(set(common))))
        tag_a = _pick(tags_a, cfg["view_a"], "view_a", fallback)
        tag_b = _pick(tags_b, cfg["view_b"], "view_b", fallback)
        name_a, name_b = _view_name(tag_a), _view_name(tag_b)

        for tag, label, name in ((tag_a, "view_a", name_a),
                                 (tag_b, "view_b", name_b)):
            steps = _n_steps(tag)
            if steps != 1:
                raise RuntimeError(
                    "%s %r has %d time steps; compare_fields is a "
                    "single-step verb (probing a step that does not exist "
                    "crashes gmsh). Extract one step first."
                    % (label, name, steps))

        box = cfg["bbox"]
        if box is None:
            lo = [max(bbox_a[c], bbox_b[c]) for c in range(3)]
            hi = [min(bbox_a[c + 3], bbox_b[c + 3]) for c in range(3)]
            if any(hi[c] < lo[c] for c in range(3)):
                raise RuntimeError(
                    "bounding boxes do not overlap (a: %s, b: %s); pass an "
                    "explicit bbox if the models are in different frames"
                    % (bbox_a, bbox_b))
            box = [lo, hi]

        pts = _points(box, int(cfg["n_points"]), cfg["sample"],
                      int(cfg["seed"]))

        ncomp = None
        diffs = []
        kept = []
        mag_a2 = 0.0
        mag_b2 = 0.0
        peak_a = 0.0
        peak_b = 0.0
        n_skipped = 0
        for pt in pts:
            va, da = gmsh.view.probe(tag_a, pt[0], pt[1], pt[2], step=0)
            vb, db = gmsh.view.probe(tag_b, pt[0], pt[1], pt[2], step=0)
            va = [float(v) for v in va]
            vb = [float(v) for v in vb]
            if not va or not vb or float(da) > 0.0 or float(db) > 0.0:
                n_skipped += 1
                continue
            if len(va) != len(vb):
                raise RuntimeError(
                    "view_a %r has %d components but view_b %r has %d"
                    % (name_a, len(va), name_b, len(vb)))
            if ncomp is None:
                ncomp = len(va)
            # VECTOR DIFFERENCE norm -- never |a| - |b| (a magnitude
            # difference hides a rotated field completely).
            d = math.sqrt(sum((x - y) ** 2 for x, y in zip(va, vb)))
            na = math.sqrt(sum(x * x for x in va))
            nb = math.sqrt(sum(y * y for y in vb))
            diffs.append(d)
            kept.append(pt)
            mag_a2 += na * na
            mag_b2 += nb * nb
            peak_a = max(peak_a, na)
            peak_b = max(peak_b, nb)

        n_valid = len(diffs)
        if n_valid == 0:
            raise RuntimeError(
                "no sample point was found in BOTH meshes (%d probed); "
                "check the bbox and that both files carry a mesh"
                % len(pts))
        l2 = math.sqrt(math.fsum(d * d for d in diffs) / n_valid)
        worst = max(range(n_valid), key=lambda i: diffs[i])
        linf = diffs[worst]
        rms_ref = max(math.sqrt(mag_a2 / n_valid), math.sqrt(mag_b2 / n_valid))
        peak_ref = max(peak_a, peak_b)
        out_file = cfg.get("out_file")
        if out_file:
            with open(out_file, "w", encoding="utf-8") as fh:
                fh.write('View "field difference magnitude" {\n')
                for pt, d in zip(kept, diffs):
                    fh.write("SP(%.17g,%.17g,%.17g){%.17g};\n"
                             % (pt[0], pt[1], pt[2], d))
                fh.write("};\n")
        result.update({
            "ok": True, "ran": True,
            "view_a": name_a, "view_b": name_b, "ncomp": ncomp,
            "l2": l2, "linf": linf,
            "l2_rel": (l2 / rms_ref) if rms_ref > 0.0 else None,
            "linf_rel": (linf / peak_ref) if peak_ref > 0.0 else None,
            "n_valid": n_valid, "n_skipped": n_skipped,
            "n_probed": len(pts),
            "worst_point": kept[worst], "worst_diff": linf,
            "rms_a": math.sqrt(mag_a2 / n_valid),
            "rms_b": math.sqrt(mag_b2 / n_valid),
            "bbox": box, "sample": cfg["sample"], "seed": int(cfg["seed"]),
            "out_file": out_file,
        })
    finally:
        gmsh.finalize()
except Exception as exc:  # noqa: BLE001 - reported to the caller as JSON
    result["error"] = "%s: %s" % (type(exc).__name__, exc)

with open(out_path, "w", encoding="utf-8") as f:
    json.dump(result, f)
"""


def _clean_bbox(bbox: Any) -> list[list[float]]:
    """[[xmin,ymin,zmin],[xmax,ymax,zmax]] or a flat 6-list -> nested."""
    flat: list[float]
    try:
        rows = list(bbox)
    except TypeError as exc:
        raise ValueError(
            "bbox must be [[xmin,ymin,zmin],[xmax,ymax,zmax]] "
            "or a flat 6-number list") from exc
    if len(rows) == 2 and all(hasattr(r, "__len__") for r in rows):
        flat = [float(v) for r in rows for v in r]
    else:
        flat = [float(v) for v in rows]
    if len(flat) != 6 or not all(math.isfinite(v) for v in flat):
        raise ValueError(
            "bbox needs 6 finite numbers (xmin,ymin,zmin,xmax,ymax,zmax), "
            f"got {bbox!r}")
    lo, hi = flat[:3], flat[3:]
    if any(hi[c] < lo[c] for c in range(3)):
        raise ValueError(f"bbox max is below bbox min on some axis: {flat}")
    return [lo, hi]


def compare_fields(path_a: str | Path, path_b: str | Path, *,
                   view_a: str | int | None = None,
                   view_b: str | int | None = None,
                   n_points: int = 1000, seed: int = 0,
                   sample: str = "random",
                   bbox: Any = None,
                   out_file: str | Path | None = None,
                   timeout_s: float = 300.0) -> dict[str, Any]:
    """Compare one view of two DIFFERENT meshes at shared sample points.

    This is the cross-validation verb (FEM vs HDiv-VIM, coarse vs fine,
    solver A vs solver B): both files are probed at the SAME points, so
    the meshes need share nothing but the region they cover.  Use
    ``diff_msh`` instead for a before/after regression on one mesh.

    Both views are probed at time step 0 and must be single-step; the
    difference of a vector view is the norm of the VECTOR DIFFERENCE
    ``||a - b||``, never ``| |a| - |b| |`` (a magnitude difference is
    blind to a rotated field).  Points that miss either mesh are skipped
    and reported in ``n_skipped`` -- never counted as agreement.

    Args:
        path_a, path_b: the two .msh files (each carrying its own mesh
            and at least one view).
        view_a, view_b: view name or index in each file.  When omitted,
            the single view NAME shared by both files is used; an
            ambiguous or empty intersection is an error, not a guess.
        n_points: sample count.  ``sample="grid"`` rounds it to the
            nearest cube ``k**3``.
        seed: RNG seed -- the point cloud is reproducible.
        sample: ``"random"`` (uniform in the bbox) or ``"grid"``
            (cell centers of a regular lattice).
        bbox: explicit sampling box; default is the intersection of the
            two model bounding boxes.
        out_file: optional .pos point cloud of the difference magnitude,
            ready for ``gmsh_render``.

    Returns:
        ``{ok, l2, linf, l2_rel, linf_rel, n_valid, n_skipped,
        worst_point, worst_diff, view_a, view_b, ncomp, ...}`` where
        ``l2`` is the RMS of the pointwise difference and ``l2_rel`` /
        ``linf_rel`` normalize by the larger of the two fields' RMS /
        peak magnitude (``None`` when that reference is zero).
    """
    for label, path in (("path_a", path_a), ("path_b", path_b)):
        if not Path(path).is_file():
            return {"ok": False, "error": f"{label} not found: {path}"}
    if sample not in _SAMPLE_MODES:
        return {"ok": False,
                "error": f"sample must be one of {list(_SAMPLE_MODES)}, "
                         f"got {sample!r}"}
    try:
        count = int(n_points)
    except (TypeError, ValueError):
        return {"ok": False, "error": f"n_points must be an integer, "
                                      f"got {n_points!r}"}
    if count < 1:
        return {"ok": False,
                "error": f"n_points must be >= 1, got {count}"}
    try:
        clean_box = None if bbox is None else _clean_bbox(bbox)
    except ValueError as exc:
        return {"ok": False, "error": str(exc)}

    out_path = None
    if out_file is not None:
        out = Path(out_file)
        out.parent.mkdir(parents=True, exist_ok=True)
        out_path = str(out)

    cfg = {"path_a": str(path_a), "path_b": str(path_b),
           "view_a": view_a, "view_b": view_b,
           "n_points": count, "seed": int(seed), "sample": sample,
           "bbox": clean_box, "out_file": out_path}
    with tempfile.TemporaryDirectory(prefix="radia_mcp_gmsh_cmp_") as work:
        cfg_path = Path(work) / "cfg.json"
        cfg_path.write_text(json.dumps(cfg), encoding="utf-8")
        result = run_gmsh_json_subprocess(
            _COMPARE_SCRIPT, [str(cfg_path)], timeout_s=timeout_s,
            prefix="radia_mcp_gmsh_cmp_")
    result.setdefault("a", str(path_a))
    result.setdefault("b", str(path_b))
    return result
