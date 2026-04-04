"""
Mesh evaluation: export .vol for p=1..N and verify against CAD.

Called as subprocess from Cubit Mesh Evaluation:
    python calc_mesh_eval.py --cub5 temp.cub5 --max-order 5

Phase 1: Open cub5, get CAD values, export .vol for each order
Phase 2: Read each .vol with NGSolve, compute Volume/Area/Length
Phase 3: Output p-convergence table as JSON

IMPORTANT: NGSolve must be imported BEFORE cubit.
Outputs JSON to stdout (last line).
"""

import argparse
import json
import os
import sys
import tempfile
import time

_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import setup_paths, setup_cubit, progress, calc_main


def _log(msg):
    progress("EVAL", msg)


def evaluate_mesh(cub5_file, max_order=5):
    from ngsolve import Mesh as NGMesh, Integrate, CF, BND, BBND, VOL

    setup_paths()
    cubit = setup_cubit(cub5_file)
    if cubit is None:
        return {"error": "Cubit not available"}

    vol_ids = list(cubit.get_entities("volume"))
    if not vol_ids:
        return {"error": "No volumes found in model"}

    # CAD values
    cad_mats = {}
    for vid in vol_ids:
        name = cubit.get_entity_name("volume", vid) or f"Volume {vid}"
        cad_mats[name] = cubit.volume(vid).volume()
    cad_bnds = {}
    for sid in cubit.parse_cubit_list("surface", "all"):
        surfaces = cubit.get_relatives("volume", vol_ids[0], "surface") if len(vol_ids) == 1 else [sid]
        cad_bnds[f"surface_{sid}"] = cubit.surface(sid).area()
    cad_edges = {}
    for cid in cubit.parse_cubit_list("curve", "all"):
        ps = list(cubit.parse_cubit_list("surface", f"in curve {cid}"))
        if len(ps) >= 2:
            cad_edges[f"curve_{cid}"] = cubit.curve(cid).length()

    cad_vol_total = sum(cad_mats.values())
    cad_area_total = sum(cad_bnds.values())
    cad_length_total = sum(cad_edges.values())

    _log(f"CAD: V={cad_vol_total:.6e}, A={cad_area_total:.6e}, L={cad_length_total:.6e}")

    # Check meshed
    total_elems = sum(
        len(cubit.get_volume_tets(vid)) + len(cubit.get_volume_hexes(vid))
        for vid in vol_ids)
    if total_elems == 0:
        return {"error": "Volumes are not meshed",
                "cad_vol_total": cad_vol_total,
                "cad_area_total": cad_area_total,
                "cad_length_total": cad_length_total}

    # Export .vol for each order + verify
    from cubit_mesh_export import extract_curved_mesh

    tmpdir = tempfile.mkdtemp(prefix="radia_eval_")
    orders = []

    for p in range(1, max_order + 1):
        _log(f"p={p}: exporting...")
        t0 = time.perf_counter()
        vol_path = os.path.join(tmpdir, f"eval_p{p}.vol")

        try:
            if p == 1:
                ng_mesh = extract_curved_mesh(cubit, order=2)
                mesh = NGMesh(ng_mesh)
                mesh.Curve(1)
            else:
                ng_mesh = extract_curved_mesh(cubit, order=p)
                mesh = NGMesh(ng_mesh)
        except Exception as e:
            orders.append({"order": p, "error": str(e)})
            continue

        ng_vol = Integrate(CF(1), mesh)
        ng_area = Integrate(CF(1), mesh, VOL_or_BND=BND)
        try:
            ng_length = Integrate(CF(1), mesh, BBND)
        except Exception:
            ng_length = 0.0

        t_elapsed = time.perf_counter() - t0

        vol_err = (ng_vol - cad_vol_total) / cad_vol_total * 100 if cad_vol_total else 0
        area_err = (ng_area - cad_area_total) / cad_area_total * 100 if cad_area_total else 0
        len_err = (ng_length - cad_length_total) / cad_length_total * 100 if cad_length_total else 0

        _log(f"p={p}: V_err={vol_err:+.4f}%, A_err={area_err:+.4f}%, L_err={len_err:+.4f}%, t={t_elapsed:.2f}s")

        orders.append({
            "order": p,
            "ng_volume": ng_vol,
            "ng_area": ng_area,
            "ng_length": ng_length,
            "vol_error_pct": vol_err,
            "area_error_pct": area_err,
            "len_error_pct": len_err,
            "time": t_elapsed,
        })

    # Cleanup temp dir
    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)

    return {
        "cad_vol_total": cad_vol_total,
        "cad_area_total": cad_area_total,
        "cad_length_total": cad_length_total,
        "orders": orders,
        "max_order": max_order,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Mesh evaluation: export + verify p=1..N")
    parser.add_argument("--cub5", required=True, help="Cubit .cub5 model file")
    parser.add_argument("--max-order", type=int, default=5)

    def run(args):
        return evaluate_mesh(args.cub5, args.max_order)

    calc_main(run, parser)


if __name__ == "__main__":
    main()
