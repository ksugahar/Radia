"""
Volume calculator using NGSolve mesh integration.

Called as subprocess from Radia app or Cubit panel:
    python calc_volume.py --vol model.vol
    python calc_volume.py --cub5 model.cub5 --order 3
    python calc_volume.py --step model.step --order 3 --maxh 0.01

Three mesh paths:
  .vol:  Pre-exported Netgen mesh (preferred, no Cubit needed)
  .cub5: Cubit ACIS kernel -> extract_curved_mesh(order=N) -> NGSolve
  .step: OCC geometry -> Netgen mesh -> mesh.Curve(order) -> NGSolve

IMPORTANT: NGSolve must be imported BEFORE cubit to avoid numpy DLL conflict.
Outputs JSON to stdout (last line).
"""

import argparse
import json
import os
import sys


def calculate_volume_vol(vol_file):
    """Calculate volume from pre-exported .vol file.

    No Cubit needed. Labels (materials, boundaries) are read from .vol.
    """
    from ngsolve import Mesh, Integrate, CF, BND

    mesh = Mesh(vol_file)
    total_vol = Integrate(CF(1), mesh)
    total_area = Integrate(CF(1), mesh, BND)

    results = []
    for mat in mesh.GetMaterials():
        vol = Integrate(CF(1), mesh, definedon=mesh.Materials(mat))
        results.append({"name": mat, "ngsolve_volume": vol})

    bnd_results = []
    checked = set()
    for bnd in mesh.GetBoundaries():
        if bnd in checked:
            continue
        checked.add(bnd)
        area = Integrate(CF(1), mesh, BND, definedon=mesh.Boundaries(bnd))
        bnd_results.append({"name": bnd, "ngsolve_area": area})

    return {
        "volumes": results,
        "boundaries": bnd_results,
        "ngsolve_total": total_vol,
        "ngsolve_area": total_area,
        "source": vol_file,
    }


def calculate_volume(cub5_file, order):
    """Calculate volume using Cubit ACIS kernel + extract_curved_mesh.

    Legacy path — prefer .vol export + calculate_volume_vol().
    """
    _panels_dir = os.path.dirname(os.path.abspath(__file__))
    if _panels_dir not in sys.path:
        sys.path.insert(0, _panels_dir)
    from calc_common import setup_paths, setup_cubit

    from ngsolve import Mesh as NGMesh, Integrate, CF, BND

    setup_paths()
    cubit = setup_cubit(cub5_file)
    if cubit is None:
        return {"error": "Cubit not available"}

    vol_ids = list(cubit.get_entities("volume"))
    if not vol_ids:
        return {"error": "No volumes found in model"}

    # CAD volumes
    results = []
    for vid in vol_ids:
        v = cubit.volume(vid)
        results.append({
            "id": vid,
            "name": cubit.get_entity_name("volume", vid) or f"Volume {vid}",
            "cad_volume": v.volume(),
        })
    cad_total = sum(r["cad_volume"] for r in results)

    total_elems = sum(
        len(cubit.get_volume_tets(vid)) + len(cubit.get_volume_hexes(vid))
        for vid in vol_ids)
    if total_elems == 0:
        return {
            "volumes": results, "cad_total": cad_total,
            "error": "Volumes are not meshed.",
        }

    from cubit_netgen_bridge import extract_curved_mesh

    try:
        build_order = max(order, 2)
        ng_mesh = extract_curved_mesh(cubit, order=build_order)
    except Exception as e:
        return {
            "volumes": results, "cad_total": cad_total,
            "error": f"extract_curved_mesh(order={order}) failed: {e}",
        }

    mesh = NGMesh(ng_mesh)
    if order == 1:
        mesh.Curve(1)

    total_vol = Integrate(CF(1), mesh)
    total_area = Integrate(CF(1), mesh, BND)

    mats = mesh.GetMaterials()
    for i, mat in enumerate(mats):
        if i < len(results):
            try:
                results[i]["ngsolve_volume"] = Integrate(
                    CF(1), mesh, definedon=mesh.Materials(mat))
            except Exception:
                results[i]["ngsolve_volume"] = None

    return {
        "volumes": results,
        "cad_total": cad_total,
        "ngsolve_total": total_vol,
        "ngsolve_area": total_area,
        "order": order,
    }


def calculate_volume_step(step_file, order, maxh=0.01):
    """Calculate volume from STEP file using OCC + Netgen."""
    from ngsolve import Mesh as NGMesh, Integrate, CF, BND
    from netgen.occ import OCCGeometry

    geo = OCCGeometry(step_file)
    ng = geo.GenerateMesh(maxh=maxh)

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

    mesh = NGMesh(ng)
    total_vol = Integrate(CF(1), mesh)
    total_area = Integrate(CF(1), mesh, BND)

    mats = mesh.GetMaterials()
    for i, mat in enumerate(mats):
        if i < len(results):
            try:
                results[i]["ngsolve_volume"] = Integrate(
                    CF(1), mesh, definedon=mesh.Materials(mat))
            except Exception:
                results[i]["ngsolve_volume"] = None

    return {
        "volumes": results,
        "cad_total": cad_total,
        "ngsolve_total": total_vol,
        "ngsolve_area": total_area,
        "order": order,
        "maxh": maxh,
    }


def main():
    _panels_dir = os.path.dirname(os.path.abspath(__file__))
    if _panels_dir not in sys.path:
        sys.path.insert(0, _panels_dir)
    from calc_common import calc_main

    parser = argparse.ArgumentParser(description="NGSolve volume calculator")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--vol", help="Netgen .vol mesh file (preferred)")
    group.add_argument("--cub5", help="Cubit .cub5 model file (legacy)")
    group.add_argument("--step", help="STEP geometry file (.step/.stp)")
    parser.add_argument("--order", type=int, default=2, help="Curve order (1-5)")
    parser.add_argument("--maxh", type=float, default=0.01,
                        help="Max mesh size (STEP path only)")

    def run(args):
        if args.vol:
            return calculate_volume_vol(args.vol)
        elif args.step:
            return calculate_volume_step(args.step, args.order, args.maxh)
        else:
            return calculate_volume(args.cub5, args.order)

    calc_main(run, parser)


if __name__ == "__main__":
    main()
