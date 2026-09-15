"""Cubit mesh quality and volume referees, isolated from Gmsh server tooling."""
from __future__ import annotations
import math
from pathlib import Path
from typing import Any
from ._gmsh_subprocess import run_gmsh_json_subprocess

_MESH_QUALITY_SCRIPT = r"""
import json
import sys

msh_path, quadrature, threshold_s, worst_n_s, stats_s, out_path = sys.argv[1:7]
threshold = float(threshold_s)
worst_n = int(worst_n_s)
want_stats = stats_s == "1"
result = {"ok": False, "ran": False}
try:
    import numpy as np
    import gmsh
    gmsh.initialize(["-noconfig"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(msh_path)
        by_type = []
        total_neg = 0
        total_below = 0
        hist_edges = [0.0, 0.1, 0.3, 0.6, 0.9, 1.0]
        hist_total = [0] * (len(hist_edges) - 1)
        etypes, etags_all, _ = gmsh.model.mesh.getElements(3)
        for et, etags in zip(etypes, etags_all):
            name, _dim, order, _nn, _coords, _ = \
                gmsh.model.mesh.getElementProperties(int(et))
            local, weights = gmsh.model.mesh.getIntegrationPoints(
                int(et), quadrature)
            _jac, det, _pts = gmsh.model.mesh.getJacobians(int(et), local)
            det = np.asarray(det, dtype=float).reshape(len(etags),
                                                       len(weights))
            e_min = det.min(axis=1)
            e_max = det.max(axis=1)
            quality = np.asarray(
                gmsh.model.mesh.getElementQualities(etags, "minSICN"),
                dtype=float)
            inverted = (e_min <= 0.0) | (quality <= 0.0)
            neg = int(inverted.sum())
            finite_quality = quality[np.isfinite(quality)]
            # Keep min(detJ)/max(detJ) as a curvature diagnostic. It is
            # exactly 1 for every affine element, including a bad sliver,
            # so it must never be used as the shape-quality gate.
            with np.errstate(divide="ignore", invalid="ignore"):
                jac_ratio = np.where(e_max > 0.0, e_min / e_max, -1.0)
            below_mask = ((quality < threshold) | ~np.isfinite(quality)) & ~inverted
            below = int(below_mask.sum())
            hist, _ = np.histogram(np.clip(finite_quality, 0.0, 1.0),
                                   bins=hist_edges)
            worst_idx = np.argsort(
                np.nan_to_num(quality, nan=-np.inf))[:worst_n]
            tags_arr = np.asarray(etags)
            worst = [{"tag": int(tags_arr[i]),
                      "quality": (float(quality[i])
                                  if np.isfinite(quality[i]) else None),
                      "metric": "minSICN",
                      "jacobian_ratio": float(jac_ratio[i]),
                      "min_det": float(e_min[i])} for i in worst_idx]
            aniso = None
            if want_stats:
                # minSICN is a SHAPE number and says little about stretching;
                # thin-gap / lamination meshes are exactly where that gap
                # bites, so report the edge-length aspect ratio explicitly.
                e_lo = np.asarray(gmsh.model.mesh.getElementQualities(
                    etags, "minEdge"), dtype=float)
                e_hi = np.asarray(gmsh.model.mesh.getElementQualities(
                    etags, "maxEdge"), dtype=float)
                iso = np.asarray(gmsh.model.mesh.getElementQualities(
                    etags, "minIsotropy"), dtype=float)
                with np.errstate(divide="ignore", invalid="ignore"):
                    ar = np.where(e_lo > 0.0, e_hi / e_lo, np.inf)
                ar_f = ar[np.isfinite(ar)]
                aniso = {
                    "aspect_ratio": {
                        "min": float(ar_f.min()) if ar_f.size else None,
                        "mean": float(ar_f.mean()) if ar_f.size else None,
                        "p95": (float(np.percentile(ar_f, 95))
                                if ar_f.size else None),
                        "max": float(ar_f.max()) if ar_f.size else None,
                        "n_above_10": int((ar_f > 10.0).sum()),
                        "n_nonfinite": int((~np.isfinite(ar)).sum()),
                    },
                    "min_isotropy": {
                        "min": float(iso.min()) if iso.size else None,
                        "mean": float(iso.mean()) if iso.size else None,
                    },
                }
            by_type.append({
                "type": int(et), "name": str(name), "order": int(order),
                "n_elements": int(len(etags)),
                "metric": "minSICN",
                "min_quality": (float(finite_quality.min())
                                if finite_quality.size else None),
                "mean_quality": (float(finite_quality.mean())
                                 if finite_quality.size else None),
                "nonfinite_quality": int((~np.isfinite(quality)).sum()),
                "min_jacobian_ratio": float(jac_ratio.min()) if jac_ratio.size else None,
                "mean_jacobian_ratio": float(jac_ratio.mean()) if jac_ratio.size else None,
                "negative": neg,
                "below_threshold": below,
                "worst": worst,
                **(aniso or {}),
            })
            total_neg += neg
            total_below += below
            hist_total = [a + int(b) for a, b in zip(hist_total, hist)]
        result.update({
            "ran": True,
            "ok": total_neg == 0 and total_below == 0,
            "applicable": bool(by_type),
            "metric": "minSICN",
            "quadrature": quadrature,
            "threshold": threshold,
            "by_type": by_type,
            "total_negative": total_neg,
            "total_below_threshold": total_below,
            "histogram": {"edges": hist_edges, "counts": hist_total},
        })
        if want_stats and by_type:
            # ---- cost axis + anisotropy -------------------------------
            # Measured (validation_test/radia_mcp/mesh_quality_study):
            # element COUNT is not the cost -- the linear system is sized
            # by dof, and the share of INTERIOR nodes is what separated
            # the meshers at matched dof. minSICN alone reports neither.
            gmsh.model.mesh.createEdges()
            gmsh.model.mesh.createFaces()
            edge_tags, _ = gmsh.model.mesh.getAllEdges()
            tri_tags, _ = gmsh.model.mesh.getAllFaces(3)
            quad_tags, _ = gmsh.model.mesh.getAllFaces(4)
            node_tags, _, _ = gmsh.model.mesh.getNodes()
            n_elem3d = sum(bt["n_elements"] for bt in by_type)

            # boundary faces are those incident to exactly ONE 3D element
            bnd_nodes = set()
            n_bnd_faces = 0
            for et, etags in zip(etypes, etags_all):
                for fnn in (3, 4):
                    try:
                        fn = gmsh.model.mesh.getElementFaceNodes(int(et), fnn)
                    except Exception:
                        continue
                    if not len(fn):
                        continue
                    fa = np.asarray(fn, dtype=np.int64).reshape(-1, fnn)
                    key = np.sort(fa, axis=1)
                    uniq, inv_idx, counts = np.unique(
                        key, axis=0, return_inverse=True, return_counts=True)
                    once = counts == 1
                    n_bnd_faces += int(once.sum())
                    bnd_nodes.update(uniq[once].ravel().tolist())
            n_nodes = int(len(node_tags))
            n_bnd_nodes = len(bnd_nodes)
            result["mesh_stats"] = {
                "n_nodes": n_nodes,
                "n_edges": int(len(edge_tags)),
                "n_faces_tri": int(len(tri_tags)),
                "n_faces_quad": int(len(quad_tags)),
                "n_elements_3d": int(n_elem3d),
                "n_boundary_faces": n_bnd_faces,
                "n_boundary_nodes": n_bnd_nodes,
                "n_interior_nodes": n_nodes - n_bnd_nodes,
                "interior_node_fraction": (
                    (n_nodes - n_bnd_nodes) / n_nodes if n_nodes else None),
                "dof_estimate": {
                    "h1_p1": n_nodes,
                    "hcurl_lowest": int(len(edge_tags)),
                    "hdiv_lowest": int(len(tri_tags)) + int(len(quad_tags)),
                    "l2_p0": int(n_elem3d),
                },
                "note": ("dof_estimate is the TOTAL (unconstrained) dof "
                         "count of the named lowest-order space on this "
                         "mesh -- the honest cost axis for ranking meshes; "
                         "element count is not."),
            }
        if not by_type:
            result["ok"] = True
            result["note"] = "no 3D elements; quality gate not applicable"
    finally:
        gmsh.finalize()
except Exception as exc:
    result["error"] = f"{type(exc).__name__}: {exc}"
with open(out_path, "w", encoding="utf-8") as f:
    json.dump(result, f)
"""


def mesh_quality(msh_path: str | Path,
                 threshold: float = 0.1,
                 quadrature: str = "Gauss4",
                 worst_n: int = 10,
                 include_mesh_stats: bool = True,
                 timeout_s: float = 600.0) -> dict[str, Any]:
    """Gmsh minSICN shape-quality distribution for all 3D elements.

    Complements the sign-only Jacobian gate: an affine element can be
    non-inverted yet nearly degenerate. Gmsh's signed inverse condition
    number detects that shape degradation. The sampled min(detJ)/max(detJ)
    ratio is retained separately as a curvature diagnostic.

    With ``include_mesh_stats`` (default on) the report also carries the
    axes that the lab mesh-quality study found actually discriminate
    meshes -- none of which minSICN expresses:

    * ``mesh_stats.dof_estimate`` -- nodes / edges / faces, i.e. the dof
      count of H1-p1, lowest-order HCurl and lowest-order HDiv on this
      mesh. Element count is NOT the cost; the linear system is sized by
      dof, and ranking meshes by element count inverts the verdict.
    * ``mesh_stats.interior_node_fraction`` -- the share of nodes that
      are not on the boundary. A mesher that spends its nodes on the
      surface delivers less accuracy per dof.
    * per-element-type ``aspect_ratio`` (maxEdge/minEdge) and
      ``min_isotropy`` -- stretching, which a single shape number misses
      and which is exactly what thin-gap and lamination meshes exhibit.

    Set ``include_mesh_stats=False`` to skip the extra gmsh edge/face
    construction on very large meshes.
    """
    path = Path(msh_path)
    if not path.is_file():
        return {"ok": False, "ran": False, "path": str(path),
                "error": f"file not found: {path}"}
    result = run_gmsh_json_subprocess(
        _MESH_QUALITY_SCRIPT,
        [str(path), quadrature, str(float(threshold)), str(int(worst_n)),
         "1" if include_mesh_stats else "0"],
        timeout_s=timeout_s, prefix="radia_mcp_gmsh_quality_")
    result["path"] = str(path)
    return result


_MESH_VOLUME_SCRIPT = r"""
import json
import sys

msh_path, quadrature, out_path = sys.argv[1:4]
result = {"ok": False, "ran": False}
try:
    import numpy as np
    import gmsh
    gmsh.initialize(["-noconfig"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(msh_path)
        total = 0.0
        n_elem = 0
        min_det = None
        etypes, etags_all, _ = gmsh.model.mesh.getElements(3)
        for et, etags in zip(etypes, etags_all):
            local, weights = gmsh.model.mesh.getIntegrationPoints(
                int(et), quadrature)
            _jac, det, _pts = gmsh.model.mesh.getJacobians(int(et), local)
            det = np.asarray(det, dtype=float).reshape(len(etags),
                                                       len(weights))
            weights = np.asarray(weights, dtype=float)
            if not np.isfinite(det).all() or not np.isfinite(weights).all():
                raise RuntimeError("non-finite Jacobian or quadrature weight")
            total += float((det * weights).sum())
            n_elem += int(len(etags))
            m = float(det.min())
            min_det = m if min_det is None else min(min_det, m)
        valid = (n_elem > 0 and np.isfinite(total)
                 and min_det is not None and np.isfinite(min_det))
        result.update({"ran": True, "ok": bool(valid),
                       "total_volume": total, "n_elements_3d": n_elem,
                       "min_jacobian_det": min_det})
        if n_elem == 0:
            result["error"] = "no 3D elements"
        elif not valid:
            result["error"] = "non-finite integrated volume or Jacobian"
    finally:
        gmsh.finalize()
except Exception as exc:
    result["error"] = f"{type(exc).__name__}: {exc}"
with open(out_path, "w", encoding="utf-8") as f:
    json.dump(result, f)
"""


def mesh_total_volume(msh_path: str | Path,
                      quadrature: str = "Gauss4",
                      timeout_s: float = 600.0) -> dict[str, Any]:
    """Jacobian-integrated total 3D volume of a ``.msh`` file.

    The closure referee for shape-regeneration pipelines: comparing this
    against the source STL's watertight volume bounds the geometric loss of
    an ``import stl`` -> mesh -> export round trip.  Also reports
    ``min_jacobian_det`` as the inversion guard.
    """
    path = Path(msh_path)
    if not path.is_file():
        return {"ok": False, "ran": False, "path": str(path),
                "error": f"file not found: {path}"}
    if not isinstance(quadrature, str) or not quadrature.strip():
        raise ValueError("quadrature must be a nonempty string")
    try:
        timeout_s = float(timeout_s)
    except (TypeError, ValueError) as exc:
        raise ValueError("timeout_s must be a finite positive number") from exc
    if not math.isfinite(timeout_s) or timeout_s <= 0.0:
        raise ValueError("timeout_s must be a finite positive number")
    result = run_gmsh_json_subprocess(
        _MESH_VOLUME_SCRIPT, [str(path), quadrature],
        timeout_s=timeout_s, prefix="radia_mcp_gmsh_volume_")
    result["path"] = str(path)
    return result
