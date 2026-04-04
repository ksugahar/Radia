"""
Mesh evaluation: verify .vol files (p=1..N) against CAD reference JSON.

Called as subprocess from Cubit Mesh Evaluation GUI:
    python calc_mesh_eval.py --vol base_p1.vol base_p2.vol ... base_p5.vol

Each .vol has a companion .json (from C++ export) with CAD reference values.
Reads each .vol with NGSolve, computes Volume/Area/Length, outputs results.

Does NOT require Cubit — only NGSolve.
Outputs JSON to stdout (last line).
"""

import argparse
import json
import sys


def evaluate_vols(vol_paths):
    """Evaluate list of .vol files against their companion JSONs.

    Returns dict with cad_* totals and per-order NGSolve results.
    """
    from ngsolve import Mesh, Integrate, CF, BND, BBND, VOL
    import os

    if not vol_paths:
        return {"error": "No .vol files provided"}

    # Read CAD reference from first companion JSON
    json_path = vol_paths[0] + ".json"
    cad_ref = {}
    if os.path.exists(json_path):
        with open(json_path, "r") as f:
            cad_ref = json.load(f)

    cad_mats = cad_ref.get("materials", {})
    cad_bnds = cad_ref.get("boundaries", {})
    cad_edges = cad_ref.get("edges", {})
    cad_vol_total = sum(cad_mats.values())
    cad_area_total = sum(cad_bnds.values())
    cad_length_total = sum(cad_edges.values())

    orders = []
    for i, vp in enumerate(vol_paths):
        p = i + 1
        if not os.path.exists(vp):
            orders.append({"order": p, "error": f"File not found: {vp}"})
            continue

        try:
            mesh = Mesh(vp)
        except Exception as e:
            orders.append({"order": p, "error": str(e)})
            continue

        ng_vol = Integrate(CF(1), mesh)
        ng_area = Integrate(CF(1), mesh, VOL_or_BND=BND)
        try:
            ng_length = Integrate(CF(1), mesh, BBND)
        except Exception:
            ng_length = 0.0

        n_vol = sum(1 for _ in mesh.Elements(VOL))
        n_bnd = sum(1 for _ in mesh.Elements(BND))

        vol_err = (ng_vol - cad_vol_total) / cad_vol_total * 100 if cad_vol_total else 0
        area_err = (ng_area - cad_area_total) / cad_area_total * 100 if cad_area_total else 0
        len_err = (ng_length - cad_length_total) / cad_length_total * 100 if cad_length_total else 0

        orders.append({
            "order": p,
            "ng_volume": ng_vol,
            "ng_area": ng_area,
            "ng_length": ng_length,
            "vol_error_pct": vol_err,
            "area_error_pct": area_err,
            "len_error_pct": len_err,
            "n_vol_elems": n_vol,
            "n_bnd_elems": n_bnd,
        })

    return {
        "cad_vol_total": cad_vol_total,
        "cad_area_total": cad_area_total,
        "cad_length_total": cad_length_total,
        "orders": orders,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Mesh evaluation: verify .vol files against CAD JSON")
    parser.add_argument("--vol", nargs="+", required=True,
                        help=".vol files (one per order, p=1..N)")

    args = parser.parse_args()
    result = evaluate_vols(args.vol)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
