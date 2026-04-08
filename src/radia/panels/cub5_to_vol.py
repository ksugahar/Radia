"""
Convert Cubit .cub5 to Netgen .vol with high-order curving (Path B).

This is the Python reference implementation of mesh curving.
Path A (C++ `radia_export netgen` APREPRO command) must produce identical
results for tet meshes. This script is maintained for:
  (a) debugging and Path A/B comparison tests
  (b) sharing with Joachim (NGSolve team) as a demo
  (c) standalone .vol generation from .cub5 without Cubit GUI

Usage:
  python cub5_to_vol.py --cub5 model.cub5 --order 3
  python cub5_to_vol.py --cub5 model.cub5 --order 3 --output model.vol

IMPORTANT: NGSolve must be imported BEFORE cubit (DLL conflict avoidance).
"""

import argparse
import json
import os
import sys
import time


def cub5_to_vol(cub5_file, order=2, output=None):
    """Convert .cub5 to .vol with high-order curving via extract_curved_mesh.

    Returns dict with mesh info and CAD reference values.
    """
    # NGSolve FIRST
    from ngsolve import Mesh, Integrate, CF, BND

    # Setup Cubit
    _this_dir = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, _this_dir)
    from calc_common import setup_paths, setup_cubit
    setup_paths()
    cubit = setup_cubit(cub5_file)
    if cubit is None:
        return {"error": "Cubit not available"}

    # Export via C++ plugin (Path A)
    t0 = time.perf_counter()

    # Determine output path
    if output is None:
        base = os.path.splitext(cub5_file)[0]
        output = f"{base}_o{order}.vol"

    cubit.cmd(f'radia_export netgen "{output}" order {order} overwrite')
    t_curve = time.perf_counter() - t0
    t_save = t_curve

    # Verify
    mesh = Mesh(output)
    vol = Integrate(CF(1), mesh)
    area = Integrate(CF(1), mesh, BND)

    # CAD reference values
    cad_mats = {}
    for vid in cubit.get_entities("volume"):
        name = cubit.get_entity_name("volume", vid) or f"volume_{vid}"
        cad_mats[name] = cubit.volume(vid).volume()

    cad_bnds = {}
    for sid in cubit.parse_cubit_list("surface", "all"):
        cad_bnds[f"surface_{sid}"] = cubit.surface(sid).area()

    cad_edges = {}
    for cid in cubit.parse_cubit_list("curve", "all"):
        ps = list(cubit.parse_cubit_list("surface", f"in curve {cid}"))
        if len(ps) >= 2:
            cad_edges[f"curve_{cid}"] = cubit.curve(cid).length()

    cad_vol = sum(cad_mats.values())
    cad_area = sum(cad_bnds.values())
    vol_err = (vol - cad_vol) / cad_vol * 100 if cad_vol else 0
    area_err = (area - cad_area) / cad_area * 100 if cad_area else 0

    # Write companion JSON
    json_path = output + ".json"
    with open(json_path, "w") as f:
        json.dump({
            "materials": cad_mats,
            "boundaries": cad_bnds,
            "edges": cad_edges,
            "n_elements": mesh.ne,
            "order": order,
        }, f, indent=2)

    result = {
        "output": output,
        "order": order,
        "n_elements": mesh.ne,
        "ng_volume": vol,
        "ng_area": area,
        "cad_volume": cad_vol,
        "cad_area": cad_area,
        "vol_error_pct": vol_err,
        "area_error_pct": area_err,
        "time_curve": t_curve,
        "time_total": t_save,
    }
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Convert Cubit .cub5 to Netgen .vol (Path B)")
    parser.add_argument("--cub5", required=True, help="Input .cub5 file")
    parser.add_argument("--order", type=int, default=2, help="Curve order (2-5)")
    parser.add_argument("--output", help="Output .vol path (default: <cub5>_oN.vol)")

    args = parser.parse_args()
    if not os.path.exists(args.cub5):
        print(f"ERROR: {args.cub5} not found", file=sys.stderr)
        sys.exit(1)

    result = cub5_to_vol(args.cub5, args.order, args.output)

    if "error" in result:
        print(f"ERROR: {result['error']}", file=sys.stderr)
        sys.exit(1)

    print(f"Output: {result['output']}")
    print(f"Elements: {result['n_elements']}, Order: {result['order']}")
    print(f"Volume: {result['ng_volume']:.6e} (CAD: {result['cad_volume']:.6e}, "
          f"err: {result['vol_error_pct']:+.4e}%)")
    print(f"Area:   {result['ng_area']:.6e} (CAD: {result['cad_area']:.6e}, "
          f"err: {result['area_error_pct']:+.4e}%)")
    print(f"Time: {result['time_curve']:.2f}s (curve) / {result['time_total']:.2f}s (total)")

    # JSON on last line for programmatic use
    print(json.dumps(result))


if __name__ == "__main__":
    main()
