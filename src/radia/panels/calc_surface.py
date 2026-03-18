"""
Surface area calculator using NGSolve mesh integration.

Called as subprocess from Cubit panel:
    python calc_surface.py --cub5 model.cub5 --order 1

Workflow (STEP reimport for all orders):
  1. Open cub5, get CAD areas, export STEP
  2. Reset + reimport STEP heal (ACIS -> OCC seam fix)
  3. Remesh with same size
  4. export_netgen(cubit, geometry=OCCGeometry(step)) + SetGeomInfo
  5. mesh.Curve(order)
  6. Integrate(CF(1), mesh, BND) for surface area

IMPORTANT: NGSolve must be imported BEFORE cubit.
Outputs JSON to stdout.
"""

import argparse
import json
import os
import sys
import tempfile


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


def _get_mesh_size(cubit, vol_ids):
    """Estimate mesh size from existing mesh."""
    for vid in vol_ids:
        ne = len(cubit.get_volume_tets(vid)) + len(cubit.get_volume_hexes(vid))
        if ne > 0:
            bb = cubit.get_bounding_box("volume", vid)
            diag = bb[9] if len(bb) > 9 else 1.0
            return diag / max(ne ** (1.0 / 3.0), 1.0)
    return 0.1


def calculate_surface(cub5_file, order):
    """Calculate surface area using Cubit mesh + STEP reimport workflow.

    All orders use the same Cubit mesh with proper geometry mapping:
      STEP reimport -> export_netgen(geometry=OCC) -> Curve(order) -> Integrate
    """
    from ngsolve import Mesh, Integrate, CF, BND
    from netgen.occ import OCCGeometry

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

    # Get mesh size before reimport
    mesh_size = _get_mesh_size(cubit, vol_ids)

    # --- STEP reimport workflow (MCP knowledge) ---
    tmpdir = tempfile.mkdtemp(prefix="radia_surf_")
    step_file = os.path.join(tmpdir, "geometry.step").replace("\\", "/")
    vol_list = " ".join(str(v) for v in vol_ids)
    cubit.cmd(f'export step "{step_file}" volume {vol_list} overwrite')

    # Reimport STEP (ACIS -> OCC seam compatibility)
    cubit.cmd("reset")
    cubit.cmd(f'import step "{step_file}" heal')

    # Remesh with same size
    cubit.cmd("volume all scheme tetmesh")
    cubit.cmd(f"volume all size {mesh_size}")
    cubit.cmd("mesh volume all")

    # Register blocks
    new_vol_ids = list(cubit.get_entities("volume"))
    cubit.cmd("delete block all")
    for i, vid in enumerate(new_vol_ids):
        cubit.cmd(f"block {i + 1} add volume {vid}")
    cubit.cmd(f"block {len(new_vol_ids) + 1} add tri all")
    cubit.cmd(f'block {len(new_vol_ids) + 1} name "boundary"')

    # Export to Netgen with OCC geometry
    radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    if os.path.abspath(radia_src) not in sys.path:
        sys.path.insert(0, os.path.abspath(radia_src))
    import cubit_mesh_export

    geo = OCCGeometry(step_file)
    ngmesh = cubit_mesh_export.export_NetgenMesh(cubit, geometry=geo)

    # TODO: SetGeomInfo for curved surfaces (cylinder, torus, sphere)
    # cubit_mesh_export.set_torus_geominfo(ngmesh, ...)

    mesh = Mesh(ngmesh)

    # Curve (only effective with geometry + SetGeomInfo)
    if order > 1:
        try:
            mesh.Curve(order)
        except Exception as e:
            return {"volumes": results, "cad_total": cad_total,
                    "error": f"mesh.Curve({order}) failed: {e}"}

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
