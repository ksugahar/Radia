"""
Mesh evaluation: CAD vs NGSolve (p=1..max_order) volume and surface area.

Called as subprocess from Cubit panel:
    python calc_mesh_eval.py --cub5 model.cub5 --max-order 5

Workflow:
  1. Open .cub5, get CAD volumes and surface areas (exact from ACIS kernel)
  2. For p = 1..max_order:
     export_NGSolveCurvedMesh(cubit, order=p) -> NGSolve Mesh
     Integrate(CF(1), mesh) for volume
     Integrate(CF(1), mesh, BND) for surface area
  3. Output error table: CAD vs NGSolve for each order

IMPORTANT: NGSolve must be imported BEFORE cubit.
Outputs JSON to stdout.
"""

import argparse
import os
import sys
import time

_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import setup_paths, setup_cubit, progress, calc_main


def _log(msg):
    progress("EVAL", msg)


def evaluate_mesh(cub5_file, max_order=5):
    """Evaluate mesh: CAD vs NGSolve for p=1..max_order.

    Returns dict with per-volume CAD values and per-order NGSolve results.
    """
    from ngsolve import Mesh as NGMesh, Integrate, CF, BND, BBND, VOL

    setup_paths()
    cubit = setup_cubit(cub5_file)
    if cubit is None:
        return {"error": "Cubit not available"}

    vol_ids = list(cubit.get_entities("volume"))
    if not vol_ids:
        return {"error": "No volumes found in model"}

    # CAD values (exact from ACIS kernel)
    volumes = []
    for vid in vol_ids:
        v = cubit.volume(vid)
        surfaces = cubit.get_relatives("volume", vid, "surface")
        cad_area = sum(cubit.surface(sid).area() for sid in surfaces)
        volumes.append({
            "id": vid,
            "name": cubit.get_entity_name("volume", vid) or f"Volume {vid}",
            "cad_volume": v.volume(),
            "cad_area": cad_area,
        })
    cad_vol_total = sum(r["cad_volume"] for r in volumes)
    cad_area_total = sum(r["cad_area"] for r in volumes)

    # CAD edge lengths: boundary curves (shared by >= 2 surfaces)
    curve_ids = list(cubit.parse_cubit_list("curve", "all"))
    cad_length_total = 0.0
    for cid in curve_ids:
        parent_surfs = list(cubit.parse_cubit_list("surface", f"in curve {cid}"))
        if len(parent_surfs) >= 2:
            cad_length_total += cubit.curve(cid).length()

    _log(f"CAD: {len(vol_ids)} volumes, V={cad_vol_total:.6e}, A={cad_area_total:.6e}, L={cad_length_total:.6e}")

    # Check if meshed
    total_elems = sum(
        len(cubit.get_volume_tets(vid)) + len(cubit.get_volume_hexes(vid))
        + len(cubit.get_volume_wedges(vid))
        for vid in vol_ids)
    if total_elems == 0:
        return {
            "volumes": volumes,
            "cad_vol_total": cad_vol_total,
            "cad_area_total": cad_area_total,
            "error": "Volumes are not meshed.",
        }

    _log(f"MESH: {total_elems} elements")

    # Evaluate p=1..max_order
    from cubit_mesh_export import extract_curved_mesh

    orders = []
    for p in range(1, max_order + 1):
        _log(f"p={p}: exporting curved mesh...")
        t0 = time.perf_counter()

        try:
            if p == 1:
                # p=1: linear mesh (no curving needed)
                ng_mesh = extract_curved_mesh(cubit, order=2)
                mesh = NGMesh(ng_mesh)
                mesh.Curve(1)  # reset to linear
            else:
                ng_mesh = extract_curved_mesh(cubit, order=p)
                mesh = NGMesh(ng_mesh)
        except Exception as e:
            orders.append({
                "order": p,
                "error": str(e),
            })
            continue

        ng_vol = Integrate(CF(1), mesh)
        ng_area = Integrate(CF(1), mesh, VOL_or_BND=BND)
        try:
            ng_length = Integrate(CF(1), mesh, BBND)
        except Exception:
            ng_length = 0.0
        n_vol = sum(1 for _ in mesh.Elements(VOL))
        n_bnd = sum(1 for _ in mesh.Elements(BND))
        t_elapsed = time.perf_counter() - t0

        vol_err = (ng_vol - cad_vol_total) / cad_vol_total * 100 if cad_vol_total else 0
        area_err = (ng_area - cad_area_total) / cad_area_total * 100 if cad_area_total else 0
        len_err = (ng_length - cad_length_total) / cad_length_total * 100 if cad_length_total else 0

        _log(f"p={p}: V_err={vol_err:+.4f}%, A_err={area_err:+.4f}%, L_err={len_err:+.4f}%, t={t_elapsed:.2f}s")

        entry = {
            "order": p,
            "ng_volume": ng_vol,
            "ng_area": ng_area,
            "ng_length": ng_length,
            "vol_error_pct": vol_err,
            "area_error_pct": area_err,
            "len_error_pct": len_err,
            "n_vol_elems": n_vol,
            "n_bnd_elems": n_bnd,
            "time": t_elapsed,
        }

        # Per-material volumes
        mats = mesh.GetMaterials()
        per_mat = []
        for i, mat in enumerate(mats):
            try:
                mv = Integrate(CF(1), mesh, definedon=mesh.Materials(mat))
                per_mat.append({"name": mat, "ng_volume": mv})
            except Exception:
                pass
        entry["materials"] = per_mat

        orders.append(entry)

    return {
        "volumes": volumes,
        "cad_vol_total": cad_vol_total,
        "cad_area_total": cad_area_total,
        "cad_length_total": cad_length_total,
        "orders": orders,
        "max_order": max_order,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Mesh evaluation: CAD vs NGSolve (p=1..N)")
    parser.add_argument("--cub5", required=True, help="Cubit .cub5 model file")
    parser.add_argument("--max-order", type=int, default=5,
                        help="Maximum curve order (1-5, default 5)")

    def run(args):
        return evaluate_mesh(args.cub5, args.max_order)

    calc_main(run, parser)


if __name__ == "__main__":
    main()
