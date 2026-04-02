"""
Volume calculator using NGSolve mesh integration.

Called as subprocess from Radia app or Cubit panel:
    python calc_volume.py --cub5 model.cub5 --order 3
    python calc_volume.py --step model.step --order 3 --maxh 0.01

Two mesh paths:
  .cub5: Cubit ACIS kernel -> extract_curved_mesh(order=N) -> NGSolve
  .step: OCC geometry -> Netgen mesh -> mesh.Curve(order) -> NGSolve

IMPORTANT: NGSolve must be imported BEFORE cubit to avoid numpy DLL conflict.
Outputs JSON to stdout (last line).
"""

import argparse
import json
import os
import sys


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
    # Use -noinit to prevent .cubit startup file from being played
    # (it loads register_toolbar.py which imports Qt, unavailable here)
    import io, contextlib
    with contextlib.redirect_stderr(io.StringIO()):
        cubit.init(["cubit", "-nojournal", "-batch", "-noinit"])

    # Clean up again (cubit.init may re-add paths)
    _remove_cubit_site_packages()

    return cubit


def _remove_cubit_site_packages():
    """Remove Cubit's bundled site-packages from sys.path."""
    for p in list(sys.path):
        # Match any Cubit installation's site-packages
        if "site-packages" in p and ("cubit" in p.lower() or "Cubit" in p):
            sys.path.remove(p)



def calculate_volume(cub5_file, order, output_pkl=None):
    """Calculate volume using NGSolve integration.

    Uses extract_curved_mesh() with Cubit's ACIS kernel for high-order curving.

    Args:
        cub5_file: Path to Cubit .cub5 model file.
        order: Polynomial order for curving (1-5).
        output_pkl: If set, save curved netgen.Mesh to this pickle file.
                    Uses SetGeometry(None) to detach ACIS callbacks before pickle.
    """
    # 1. Import NGSolve FIRST
    from ngsolve import Mesh as NGMesh, Integrate, CF

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

    # 6. Check if meshed
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

    # 7. Export Cubit mesh with curving (ACIS kernel, no STEP needed)
    from cubit_netgen_bridge import extract_curved_mesh

    try:
        ng_mesh = extract_curved_mesh(cubit, order=order)
    except Exception as e:
        return {
            "volumes": results,
            "cad_total": cad_total,
            "error": f"extract_curved_mesh(order={order}) failed: {e}",
        }

    # 8. Save curved mesh to pickle if requested
    if output_pkl:
        import pickle
        ng_mesh.SetGeometry(None)  # detach ACIS callbacks for serialization
        with open(output_pkl, "wb") as f:
            pickle.dump(ng_mesh, f, protocol=2)

    mesh = NGMesh(ng_mesh)
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

    result = {
        "volumes": results,
        "cad_total": cad_total,
        "ngsolve_total": total_vol,
        "order": order,
    }
    if output_pkl:
        result["output_pkl"] = output_pkl
    return result


def calculate_volume_step(step_file, order, maxh=0.01, output_pkl=None):
    """Calculate volume from STEP file using OCC + Netgen.

    No Cubit needed. Uses NGSolve OCC for geometry, Netgen for meshing,
    mesh.Curve(order) for high-order curving.
    """
    from ngsolve import Mesh as NGMesh, Integrate, CF
    from netgen.occ import OCCGeometry

    geo = OCCGeometry(step_file)
    ng = geo.GenerateMesh(maxh=maxh)

    # CAD volume from OCC
    shape = geo.shape
    cad_total = shape.mass

    results = []
    for i, solid in enumerate(shape.solids):
        results.append({
            "id": i + 1,
            "name": solid.name or f"Solid {i+1}",
            "cad_volume": solid.mass,
        })

    if order >= 2:
        ng.Curve(order)

    # Save curved mesh to pickle if requested
    if output_pkl:
        import pickle
        with open(output_pkl, "wb") as f:
            pickle.dump(ng, f, protocol=2)

    mesh = NGMesh(ng)
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

    result = {
        "volumes": results,
        "cad_total": cad_total,
        "ngsolve_total": total_vol,
        "order": order,
        "maxh": maxh,
    }
    if output_pkl:
        result["output_pkl"] = output_pkl
    return result


def main():
    # Import calc_common for shared main wrapper
    _panels_dir = os.path.dirname(os.path.abspath(__file__))
    if _panels_dir not in sys.path:
        sys.path.insert(0, _panels_dir)
    from calc_common import calc_main

    parser = argparse.ArgumentParser(description="NGSolve volume calculator")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--cub5", help="Cubit .cub5 model file")
    group.add_argument("--step", help="STEP geometry file (.step/.stp)")
    parser.add_argument("--order", type=int, default=2, help="Curve order (1-5)")
    parser.add_argument("--maxh", type=float, default=0.01,
                        help="Max mesh size (STEP path only)")
    parser.add_argument("--output-pkl", default=None,
                        help="Save curved netgen.Mesh to pickle file")

    def run(args):
        if args.step:
            return calculate_volume_step(args.step, args.order, args.maxh,
                                         args.output_pkl)
        else:
            return calculate_volume(args.cub5, args.order, args.output_pkl)

    calc_main(run, parser)


if __name__ == "__main__":
    main()
