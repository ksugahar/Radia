"""
Surface area calculator using NGSolve mesh integration.

Called as subprocess from Cubit panel:
    python calc_surface.py --cub5 model.cub5 --order 1

Workflow:
  1. Import NGSolve FIRST (before cubit)
  2. Open cub5, export surface mesh
  3. Integrate(CF(1), mesh, BND) for surface area

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


def calculate_surface(cub5_file, order):
    """Calculate surface area using NGSolve Integrate(CF(1), mesh, BND).

    For order 1: uses Cubit mesh directly.
    For order >1: uses Netgen OCC mesh with Curve(order).
    """
    from ngsolve import Mesh, Integrate, CF, BND

    cubit = _setup_cubit()
    cubit.cmd(f'open "{cub5_file}"')

    vol_ids = list(cubit.get_entities("volume"))
    if not vol_ids:
        return {"error": "No volumes found in model"}

    # CAD surface areas
    results = []
    for vid in vol_ids:
        v = cubit.volume(vid)
        surfaces = cubit.get_relatives("volume", vid, "surface")
        cad_area = sum(cubit.surface(sid).area() for sid in surfaces)
        results.append({
            "id": vid,
            "name": cubit.get_entity_name("volume", vid) or f"Volume {vid}",
            "cad_area": cad_area,
            "n_surfaces": len(surfaces),
        })
    cad_total = sum(r["cad_area"] for r in results)

    # Export STEP
    tmpdir = tempfile.mkdtemp(prefix="radia_surf_")
    step_file = os.path.join(tmpdir, "geometry.step").replace("\\", "/")
    vol_list = " ".join(str(v) for v in vol_ids)
    cubit.cmd(f'export step "{step_file}" volume {vol_list} overwrite')

    if order == 1:
        # Order 1: Use Cubit mesh
        total_elems = sum(
            len(cubit.get_volume_tets(vid)) + len(cubit.get_volume_hexes(vid))
            for vid in vol_ids
        )
        if total_elems == 0:
            return {
                "volumes": results,
                "cad_total": cad_total,
                "error": "Volumes are not meshed.",
            }

        cubit.cmd("delete block all")
        for i, vid in enumerate(vol_ids):
            cubit.cmd(f"block {i + 1} add volume {vid}")
        cubit.cmd(f"block {len(vol_ids) + 1} add tri all")

        radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
        if os.path.abspath(radia_src) not in sys.path:
            sys.path.insert(0, os.path.abspath(radia_src))
        import cubit_mesh_export

        ngmesh = cubit_mesh_export.export_NetgenMesh(cubit)
        mesh = Mesh(ngmesh)
        total_area = Integrate(CF(1), mesh, VOL_or_BND=BND)

        # Per-boundary areas
        bnd_areas = []
        for bnd in mesh.GetBoundaries():
            try:
                area = Integrate(CF(1), mesh, definedon=mesh.Boundaries(bnd))
                bnd_areas.append({"boundary": bnd, "area": area})
            except Exception:
                pass

        # Surface element info
        n_bnd_elements = sum(1 for el in mesh.Elements(BND))

    else:
        # Order > 1: Netgen OCC mesh with Curve
        from netgen.occ import OCCGeometry

        # Estimate mesh size from Cubit mesh
        for vid in vol_ids:
            ne = len(cubit.get_volume_tets(vid)) + len(cubit.get_volume_hexes(vid))
            if ne > 0:
                bb = cubit.get_bounding_box("volume", vid)
                diag = bb[9] if len(bb) > 9 else 1.0
                mesh_size = diag / max(ne ** (1.0 / 3.0), 1.0)
                break
        else:
            mesh_size = 0.1

        geo = OCCGeometry(step_file)
        try:
            ngmesh = geo.GenerateMesh(maxh=mesh_size)
        except Exception as e:
            return {"volumes": results, "cad_total": cad_total,
                    "error": f"Netgen OCC meshing failed: {e}"}

        mesh = Mesh(ngmesh)
        try:
            mesh.Curve(order)
        except Exception as e:
            return {"volumes": results, "cad_total": cad_total,
                    "error": f"mesh.Curve({order}) failed: {e}"}

        total_area = Integrate(CF(1), mesh, VOL_or_BND=BND)
        bnd_areas = []
        for bnd in mesh.GetBoundaries():
            try:
                area = Integrate(CF(1), mesh, definedon=mesh.Boundaries(bnd))
                bnd_areas.append({"boundary": bnd, "area": area})
            except Exception:
                pass
        n_bnd_elements = sum(1 for el in mesh.Elements(BND))

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
