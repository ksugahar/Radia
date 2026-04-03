"""
Export Cubit mesh as Netgen .vol and/or curved .pkl.

Called as subprocess from Cubit Export Mesh menu:
    python calc_export_netgen.py --cub5 model.cub5 --order 3 --vol out.vol --pkl out.pkl

.vol: Netgen native format (linear only, curving NOT preserved)
.pkl: Python pickle (curving preserved via SetGeometry(None))

IMPORTANT: NGSolve must be imported BEFORE cubit.
Outputs JSON to stdout (last line).
"""

import argparse
import json
import os
import sys


def _setup():
    _this_dir = os.path.dirname(os.path.abspath(__file__))
    if _this_dir not in sys.path:
        sys.path.insert(0, _this_dir)
    from calc_common import setup_paths, setup_cubit
    return setup_paths, setup_cubit


def export_netgen(cub5_file, order, vol_path=None, pkl_path=None):
    setup_paths, setup_cubit = _setup()

    from ngsolve import Mesh as NGMesh, Integrate, CF, BND

    setup_paths()
    cubit = setup_cubit(cub5_file)
    if cubit is None:
        return {"error": "Cubit not available"}

    vol_ids = list(cubit.get_entities("volume"))
    if not vol_ids:
        return {"error": "No volumes found"}

    total_elems = sum(
        len(cubit.get_volume_tets(vid)) + len(cubit.get_volume_hexes(vid))
        for vid in vol_ids)
    if total_elems == 0:
        return {"error": "Volumes are not meshed"}

    from cubit_netgen_bridge import extract_curved_mesh

    # order=1: build with order=2 (C++ requires >=2), then reset to linear
    build_order = max(order, 2)
    ng_mesh = extract_curved_mesh(cubit, order=build_order)
    mesh = NGMesh(ng_mesh)
    if order == 1:
        mesh.Curve(1)

    result = {
        "order": order,
        "n_elements": total_elems,
        "volume": Integrate(CF(1), mesh),
        "area": Integrate(CF(1), mesh, BND),
    }

    # Save .vol (linear only — curving is lost)
    if vol_path:
        ng_mesh.Save(vol_path)
        result["vol_file"] = vol_path
        result["vol_note"] = "WARNING: .vol does not preserve curving"

    # Save .pkl (curving preserved)
    if pkl_path:
        import pickle
        ng_mesh.SetGeometry(None)
        with open(pkl_path, "wb") as f:
            pickle.dump(ng_mesh, f, protocol=2)
        result["pkl_file"] = pkl_path

    return result


def main():
    parser = argparse.ArgumentParser(description="Export Cubit mesh as Netgen .vol/.pkl")
    parser.add_argument("--cub5", required=True)
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--vol", default=None, help="Output .vol path")
    parser.add_argument("--pkl", default=None, help="Output .pkl path (curved)")

    args = parser.parse_args()
    result = export_netgen(args.cub5, args.order, args.vol, args.pkl)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
