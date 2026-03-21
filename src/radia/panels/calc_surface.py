"""
Surface area calculator using NGSolve mesh integration.

Called as subprocess from Cubit panel:
    python calc_surface.py --cub5 model.cub5 --order 1

Workflow:
  1. Open cub5, get CAD areas
  2. export_NGSolveCurvedMesh(cubit, order=N) -> NGSolve Mesh (curved via ACIS kernel)
  3. Integrate(CF(1), mesh, BND) for surface area

IMPORTANT: NGSolve must be imported BEFORE cubit.
Outputs JSON to stdout.
"""

import argparse
import json
import os
import sys


def _setup_cubit():
    """Import cubit with path cleanup (NGSolve already imported)."""
    radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    if os.path.abspath(radia_src) not in sys.path:
        sys.path.insert(0, os.path.abspath(radia_src))
    from install_panels import find_cubit_bin

    cubit_path = find_cubit_bin()
    if cubit_path and cubit_path not in sys.path:
        sys.path.append(cubit_path)

    for p in list(sys.path):
        if "site-packages" in p and ("cubit" in p.lower() or "Cubit" in p):
            sys.path.remove(p)

    import cubit
    import io, contextlib
    with contextlib.redirect_stderr(io.StringIO()):
        cubit.init(["cubit", "-nojournal", "-batch", "-noinit"])

    for p in list(sys.path):
        if "site-packages" in p and ("cubit" in p.lower() or "Cubit" in p):
            sys.path.remove(p)

    return cubit



def calculate_surface(cub5_file, order):
    """Calculate surface area using Cubit mesh + export_NGSolveCurvedMesh().

    Uses export_NGSolveCurvedMesh() which works directly with Cubit's ACIS kernel
    for high-order curving. No STEP files or OCC geometry needed.
    """
    from ngsolve import Integrate, CF, BND

    cubit = _setup_cubit()
    cubit.cmd(f'open "{cub5_file}"')

    vol_ids = list(cubit.get_entities("volume"))
    if not vol_ids:
        return {"error": "No volumes found in model"}

    # CAD surface areas
    results = []
    for vid in vol_ids:
        surfaces = cubit.get_relatives("volume", vid, "surface")
        cad_area = sum(cubit.surface(sid).area() for sid in surfaces)
        results.append({
            "id": vid,
            "name": cubit.get_entity_name("volume", vid) or f"Volume {vid}",
            "cad_area": cad_area,
        })
    cad_total = sum(r["cad_area"] for r in results)

    # Check if meshed
    total_elems = sum(
        len(cubit.get_volume_tets(vid)) + len(cubit.get_volume_hexes(vid))
        for vid in vol_ids
    )
    if total_elems == 0:
        return {"volumes": results, "cad_total": cad_total,
                "error": "Volumes are not meshed."}

    # Export to NGSolve mesh with curving (ACIS kernel, no STEP needed)
    radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    if os.path.abspath(radia_src) not in sys.path:
        sys.path.insert(0, os.path.abspath(radia_src))
    import cubit_mesh_export

    try:
        mesh = cubit_mesh_export.export_NGSolveCurvedMesh(cubit, order=order)
    except Exception as e:
        return {"volumes": results, "cad_total": cad_total,
                "error": f"export_NGSolveCurvedMesh(order={order}) failed: {e}"}

    # Integrate
    total_area = Integrate(CF(1), mesh, VOL_or_BND=BND)
    n_bnd_elements = sum(1 for _ in mesh.Elements(BND))

    bnd_areas = []
    for bnd in mesh.GetBoundaries():
        try:
            area = Integrate(CF(1), mesh, definedon=mesh.Boundaries(bnd))
            bnd_areas.append({"boundary": bnd, "area": area})
        except Exception:
            pass

    if len(results) == 1:
        results[0]["ngsolve_area"] = total_area
        results[0]["n_bnd_elements"] = n_bnd_elements

    return {
        "volumes": results,
        "cad_total": cad_total,
        "ngsolve_total": total_area,
        "n_bnd_elements": n_bnd_elements,
        "boundaries": bnd_areas,
        "order": order,
    }


def main():
    parser = argparse.ArgumentParser(description="NGSolve surface area calculator")
    parser.add_argument("--cub5", required=True, help="Cubit .cub5 model file")
    parser.add_argument("--order", type=int, default=1, help="Curve order (1-5)")
    args = parser.parse_args()

    try:
        result = calculate_surface(args.cub5, args.order)
        print(json.dumps(result))
    except Exception as e:
        print(json.dumps({"error": str(e)}))
        sys.exit(1)


if __name__ == "__main__":
    main()
