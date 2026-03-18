"""
Volume calculator using NGSolve mesh integration + Cubit.

Called as subprocess from Cubit panel:
    python calc_volume.py --cub5 model.cub5 --order 3 [--step geometry.step]

IMPORTANT: NGSolve must be imported BEFORE cubit to avoid numpy DLL conflict.
"""

import argparse
import json
import os
import sys
import tempfile


def _setup_cubit():
    """Import cubit with path cleanup (NGSolve already imported by caller)."""
    cubit_path = os.environ.get(
        "CUBIT_PATH", r"C:\Program Files\Coreform Cubit 2025.3\bin"
    )
    if cubit_path not in sys.path:
        sys.path.append(cubit_path)

    # Block Cubit's site-packages (contains incompatible numpy for Python 3.10)
    cubit_site = os.path.join(cubit_path, "python3", "lib", "site-packages")
    if cubit_site in sys.path:
        sys.path.remove(cubit_site)

    import cubit
    cubit.init(["cubit", "-nojournal", "-batch"])

    # Clean up paths cubit.init() may have added
    for p in list(sys.path):
        if "Coreform Cubit" in p and "site-packages" in p:
            sys.path.remove(p)

    return cubit


def calculate_volume(cub5_file, order, step_file=None):
    """Calculate volume of Cubit model using NGSolve integration.

    Args:
        cub5_file: Path to .cub5 file (saved from Cubit GUI)
        order: Polynomial order for mesh.Curve()
        step_file: STEP file for geometry (auto-exported if omitted)
    """
    # 1. Import NGSolve FIRST
    from ngsolve import Mesh, Integrate, CF
    from netgen.occ import OCCGeometry

    # 2. Import cubit
    cubit = _setup_cubit()

    # 3. Load model
    cubit.cmd(f'open "{cub5_file}"')

    # 4. Get selected volumes (or all)
    selected = cubit.parse_cubit_list("volume", "selected")
    if not selected:
        selected = cubit.get_entities("volume")
    vol_ids = list(selected)

    if not vol_ids:
        return {"error": "No volumes found in model"}

    # 5. CAD volumes
    results = []
    for vid in vol_ids:
        v = cubit.volume(vid)
        results.append({
            "id": vid,
            "name": cubit.get_entity_name("volume", vid) or f"Volume {vid}",
            "cad_volume": v.volume(),
        })

    # 6. Check if meshed
    total_elems = sum(
        len(cubit.get_volume_tets(vid)) + len(cubit.get_volume_hexes(vid))
        for vid in vol_ids
    )
    if total_elems == 0:
        return {
            "volumes": results,
            "cad_total": sum(r["cad_volume"] for r in results),
            "error": "Volumes are not meshed. Mesh before calculating NGSolve volume.",
        }

    # 7. Register blocks
    cubit.cmd("delete block all")
    for i, vid in enumerate(vol_ids):
        cubit.cmd(f"block {i + 1} add volume {vid}")
    cubit.cmd(f"block {len(vol_ids) + 1} add tri all")

    # 8. Export STEP if not provided
    tmpdir = tempfile.mkdtemp(prefix="radia_vol_")
    if not step_file:
        step_file = os.path.join(tmpdir, "geometry.step")
        vol_list = " ".join(str(v) for v in vol_ids)
        cubit.cmd(f'export step "{step_file}" volume {vol_list} overwrite')

    # 9. Build Netgen mesh with geometry
    radia_src = os.path.join(os.path.dirname(__file__), "..")
    if os.path.abspath(radia_src) not in sys.path:
        sys.path.insert(0, os.path.abspath(radia_src))
    import cubit_mesh_export

    geo = OCCGeometry(step_file)
    ngmesh = cubit_mesh_export.export_netgen(cubit, geometry=geo)

    # 10. NGSolve volume integration
    mesh = Mesh(ngmesh)
    if order > 1:
        try:
            mesh.Curve(order)
        except Exception as e:
            # Still return order-1 result with warning
            vol_1 = Integrate(CF(1), mesh)
            for r in results:
                r["ngsolve_volume"] = None
            return {
                "volumes": results,
                "cad_total": sum(r["cad_volume"] for r in results),
                "ngsolve_total": vol_1,
                "order": 1,
                "warning": f"mesh.Curve({order}) failed: {e}. Showing order 1.",
            }

    total_vol = Integrate(CF(1), mesh)

    # Per-material volumes
    mats = mesh.GetMaterials()
    for i, mat in enumerate(mats):
        if i < len(results):
            try:
                vol = Integrate(CF(1), mesh, definedon=mesh.Materials(mat))
                results[i]["ngsolve_volume"] = vol
            except Exception:
                results[i]["ngsolve_volume"] = None

    cad_total = sum(r["cad_volume"] for r in results)

    return {
        "volumes": results,
        "cad_total": cad_total,
        "ngsolve_total": total_vol,
        "order": order,
    }


def main():
    parser = argparse.ArgumentParser(description="NGSolve volume calculator")
    parser.add_argument("--cub5", required=True, help="Cubit .cub5 model file")
    parser.add_argument("--order", type=int, default=1, help="Curve order (1-5)")
    parser.add_argument("--step", default=None, help="STEP file (auto-exported if omitted)")
    args = parser.parse_args()

    try:
        result = calculate_volume(args.cub5, args.order, args.step)
        print(json.dumps(result))
    except Exception as e:
        print(json.dumps({"error": str(e)}))
        sys.exit(1)


if __name__ == "__main__":
    main()
