"""
Volume calculator using NGSolve mesh integration + Cubit.

Called as subprocess from Cubit panel:
    python calc_volume.py --cub5 model.cub5 --order 3

Workflow:
  Order 1: Cubit mesh -> export_NetgenMesh -> Integrate (uses Cubit mesh as-is)
  Order >1: STEP export -> Netgen OCC GenerateMesh -> Curve(order) -> Integrate
            (Netgen OCC mesh has proper geometry mapping for high-order curving)

IMPORTANT: NGSolve must be imported BEFORE cubit to avoid numpy DLL conflict.
Outputs JSON to stdout (last line).
"""

import argparse
import json
import os
import sys
import tempfile


def _setup_cubit():
    """Import cubit with path cleanup (NGSolve already imported).

    Uses install_panels.find_cubit_bin() for cross-platform detection.
    Removes Cubit's bundled site-packages from sys.path to avoid
    numpy DLL conflicts between Cubit's Python and system Python.
    """
    # Find Cubit bin directory (cross-platform)
    radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    if os.path.abspath(radia_src) not in sys.path:
        sys.path.insert(0, os.path.abspath(radia_src))
    from install_panels import find_cubit_bin

    cubit_path = find_cubit_bin()
    if cubit_path and cubit_path not in sys.path:
        sys.path.append(cubit_path)

    # Block Cubit's bundled site-packages (may contain incompatible numpy)
    _remove_cubit_site_packages()

    import cubit
    # Suppress cubit.init() banner on stderr
    import io, contextlib
    with contextlib.redirect_stderr(io.StringIO()):
        cubit.init(["cubit", "-nojournal", "-batch"])

    # Clean up again (cubit.init may re-add paths)
    _remove_cubit_site_packages()

    return cubit


def _remove_cubit_site_packages():
    """Remove Cubit's bundled site-packages from sys.path."""
    for p in list(sys.path):
        # Match any Cubit installation's site-packages
        if "site-packages" in p and ("cubit" in p.lower() or "Cubit" in p):
            sys.path.remove(p)


def _get_mesh_size(cubit, vol_ids):
    """Estimate mesh size from existing mesh."""
    for vid in vol_ids:
        tets = cubit.get_volume_tets(vid)
        hexes = cubit.get_volume_hexes(vid)
        ne = len(tets) + len(hexes)
        if ne > 0:
            bb = cubit.get_bounding_box("volume", vid)
            diag = bb[9] if len(bb) > 9 else 1.0
            return diag / max(ne ** (1.0 / 3.0), 1.0)
    return 0.1


def calculate_volume(cub5_file, order):
    """Calculate volume using NGSolve integration.

    For order 1: uses Cubit mesh directly (fast, shows linear mesh accuracy).
    For order >1: uses Netgen OCC mesh (proper geometry mapping for Curve()).
    """
    # 1. Import NGSolve FIRST
    from ngsolve import Mesh, Integrate, CF
    from netgen.occ import OCCGeometry

    # 2. Import cubit
    cubit = _setup_cubit()

    # 3. Load model
    cubit.cmd(f'open "{cub5_file}"')

    # 4. Get volume IDs (all volumes in the loaded model)
    vol_ids = list(cubit.get_entities("volume"))

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
    cad_total = sum(r["cad_volume"] for r in results)

    # 6. Export STEP (needed for both order 1 and order > 1)
    tmpdir = tempfile.mkdtemp(prefix="radia_vol_")
    step_file = os.path.join(tmpdir, "geometry.step").replace("\\", "/")
    vol_list = " ".join(str(v) for v in vol_ids)
    cubit.cmd(f'export step "{step_file}" volume {vol_list} overwrite')

    if order == 1:
        # Order 1: Use Cubit mesh (shows linear mesh accuracy)
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

        # Register blocks
        cubit.cmd("delete block all")
        for i, vid in enumerate(vol_ids):
            cubit.cmd(f"block {i + 1} add volume {vid}")
        cubit.cmd(f"block {len(vol_ids) + 1} add tri all")

        # Export Cubit mesh to Netgen
        radia_src = os.path.join(os.path.dirname(__file__), "..")
        if os.path.abspath(radia_src) not in sys.path:
            sys.path.insert(0, os.path.abspath(radia_src))
        import cubit_mesh_export

        ngmesh = cubit_mesh_export.export_NetgenMesh(cubit)
        mesh = Mesh(ngmesh)
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

    else:
        # Order > 1: Use Netgen OCC mesh (proper geometry for Curve)
        # Netgen's OCC mesher creates mesh with correct UV mapping,
        # so mesh.Curve(order) works correctly.
        mesh_size = _get_mesh_size(cubit, vol_ids)
        geo = OCCGeometry(step_file)

        try:
            ngmesh = geo.GenerateMesh(maxh=mesh_size)
        except Exception as e:
            return {
                "volumes": results,
                "cad_total": cad_total,
                "error": f"Netgen OCC meshing failed: {e}",
            }

        mesh = Mesh(ngmesh)
        try:
            mesh.Curve(order)
        except Exception as e:
            return {
                "volumes": results,
                "cad_total": cad_total,
                "error": f"mesh.Curve({order}) failed: {e}",
            }

        total_vol = Integrate(CF(1), mesh)

        # For OCC mesh, per-material mapping may differ from Cubit blocks
        # Just report total for now
        if len(results) == 1:
            results[0]["ngsolve_volume"] = total_vol

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
    args = parser.parse_args()

    try:
        result = calculate_volume(args.cub5, args.order)
        print(json.dumps(result))
    except Exception as e:
        print(json.dumps({"error": str(e)}))
        sys.exit(1)


if __name__ == "__main__":
    main()
