"""
Export Cubit mesh as Netgen .vol with curving and labels.

Called as subprocess from Cubit Export Mesh menu:
    python calc_export_netgen.py --cub5 model.cub5 --order 3 --vol out.vol

Performs pre-export consistency check:
  - Per-material volume: NGSolve vs Cubit CAD
  - Per-boundary area: NGSolve vs Cubit CAD
  - Warnings if error exceeds threshold

.vol preserves: curving (post26+), material names, boundary names.

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


def export_netgen(cub5_file, order, vol_path=None):
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

    # --- CAD reference values (from Cubit ACIS kernel) ---
    cad_materials = []
    for vid in vol_ids:
        name = cubit.get_entity_name("volume", vid)
        if not name or name.startswith("Volume"):
            # Check block name
            for bid in cubit.get_block_id_list():
                bvols = list(cubit.parse_cubit_list("volume", f"in block {bid}"))
                if vid in bvols:
                    etype = cubit.get_block_element_type(bid)
                    if etype in ("TETRA", "TETRA4", "TET", "TET4",
                                 "HEX", "HEX8", "WEDGE", "WEDGE6",
                                 "PYRAMID", "PYRAMID5"):
                        name = cubit.get_exodus_entity_name("block", bid)
                        break
        if not name or name.startswith("Volume"):
            name = f"volume_{vid}"
        cad_materials.append({
            "name": name,
            "cad_volume": cubit.volume(vid).volume(),
        })

    surface_ids = list(cubit.parse_cubit_list("surface", "all"))
    cad_boundaries = []
    # Build surface -> label mapping (same priority as _set_labels)
    surf_labels = {}
    try:
        for ssid in cubit.get_sideset_id_list():
            ssname = cubit.get_exodus_entity_name("sideset", ssid)
            for sid in cubit.parse_cubit_list("surface", f"in sideset {ssid}"):
                if sid in surface_ids:
                    surf_labels[sid] = ssname
    except Exception:
        pass
    for bid in cubit.get_block_id_list():
        btris = set(cubit.get_block_tris(bid))
        try:
            bfaces = set(cubit.get_block_faces(bid))
        except Exception:
            bfaces = set()
        if not btris and not bfaces:
            continue
        bname = cubit.get_exodus_entity_name("block", bid)
        for sid in surface_ids:
            if sid in surf_labels:
                continue
            stris = set(cubit.parse_cubit_list("tri", f"in surface {sid}"))
            sfaces = set(cubit.parse_cubit_list("face", f"in surface {sid}"))
            if (stris and stris & btris) or (sfaces and sfaces & bfaces):
                surf_labels[sid] = bname
    for sid in surface_ids:
        if sid not in surf_labels:
            ename = cubit.get_entity_name("surface", sid)
            if ename and not ename.startswith("Surface"):
                surf_labels[sid] = ename
            else:
                surf_labels[sid] = f"surface_{sid}"
        cad_boundaries.append({
            "name": surf_labels[sid],
            "cad_area": cubit.surface(sid).area(),
        })

    # --- CAD reference values: curves (edges / BBND) ---
    curve_ids = list(cubit.parse_cubit_list("curve", "all"))
    cad_edges = []
    # Build curve -> label mapping (same priority as _set_labels edge logic)
    curve_labels = {}
    try:
        for ssid in cubit.get_sideset_id_list():
            ssname = cubit.get_exodus_entity_name("sideset", ssid)
            for cid in cubit.parse_cubit_list("curve", f"in sideset {ssid}"):
                if cid not in curve_labels:
                    curve_labels[cid] = ssname
    except Exception:
        pass
    # Only include curves that have >= 2 parent surfaces (boundary curves)
    for cid in curve_ids:
        parent_surfs = list(cubit.parse_cubit_list("surface", f"in curve {cid}"))
        if len(parent_surfs) < 2:
            continue
        if cid not in curve_labels:
            ename = cubit.get_entity_name("curve", cid)
            if ename and not ename.startswith("Curve"):
                curve_labels[cid] = ename
            else:
                curve_labels[cid] = f"curve_{cid}"
        cad_edges.append({
            "name": curve_labels[cid],
            "cad_length": cubit.curve(cid).length(),
        })

    # --- Extract curved mesh with labels ---
    from cubit_mesh_export import extract_curved_mesh

    build_order = max(order, 2)
    ng_mesh = extract_curved_mesh(cubit, order=build_order)
    mesh = NGMesh(ng_mesh)
    if order == 1:
        mesh.Curve(1)

    # --- Consistency check: per-material volume ---
    # Use NGSolve label as display name (overrides CAD-side fallback)
    warnings = []
    materials = list(dict.fromkeys(mesh.GetMaterials()))  # unique, order preserved
    for i, mat in enumerate(materials):
        ng_vol = Integrate(CF(1), mesh, definedon=mesh.Materials(mat))
        if i < len(cad_materials):
            cad_materials[i]["name"] = mat  # use NGSolve label
            cad_materials[i]["ng_volume"] = ng_vol
            cad_v = cad_materials[i]["cad_volume"]
            if cad_v > 0:
                err = (ng_vol - cad_v) / cad_v * 100
                cad_materials[i]["error_pct"] = err
                if abs(err) > 1.0:
                    warnings.append(
                        f"Material \"{mat}\": V error {err:+.2e}% (>{1.0}%)")

    # --- Consistency check: per-boundary area ---
    # Match by index (surface order is consistent between Cubit and _set_labels)
    boundaries = list(dict.fromkeys(mesh.GetBoundaries()))  # unique
    for i, bnd in enumerate(boundaries):
        ng_area = Integrate(CF(1), mesh, BND, definedon=mesh.Boundaries(bnd))
        if i < len(cad_boundaries):
            cad_boundaries[i]["name"] = bnd  # use NGSolve label
            cad_boundaries[i]["ng_area"] = ng_area
            cad_a = cad_boundaries[i]["cad_area"]
            if cad_a > 0:
                err = (ng_area - cad_a) / cad_a * 100
                cad_boundaries[i]["error_pct"] = err
                if abs(err) > 1.0:
                    warnings.append(
                        f"Boundary \"{bnd}\": A error {err:+.2e}% (>{1.0}%)")

    # --- Consistency check: per-edge length (BBND) ---
    # BBND requires .vol round-trip for per-label check.
    # Here we do total BBND length check (reliable without round-trip).
    try:
        from ngsolve import BBND
        ng_length_total = Integrate(CF(1), mesh, BBND)
        cad_length_total = sum(e["cad_length"] for e in cad_edges)
        for ce in cad_edges:
            # Per-edge: best-effort (may not match due to Netgen topology)
            try:
                ng_len = Integrate(CF(1), mesh, BBND,
                                   definedon=mesh.BBoundaries(ce["name"]))
                ce["ng_length"] = ng_len
                if ce["cad_length"] > 0:
                    err = (ng_len - ce["cad_length"]) / ce["cad_length"] * 100
                    ce["error_pct"] = err
            except Exception:
                pass
        if cad_length_total > 0:
            total_err = (ng_length_total - cad_length_total) / cad_length_total * 100
            if abs(total_err) > 1.0:
                warnings.append(
                    f"Total edge length: {total_err:+.2e}% (>{1.0}%)")
    except Exception:
        pass

    # --- Save .vol + companion JSON ---
    result = {
        "order": order,
        "n_elements": total_elems,
        "materials": cad_materials,
        "boundaries": cad_boundaries,
        "edges": cad_edges,
        "warnings": warnings,
    }

    if vol_path:
        ng_mesh.Save(vol_path)
        result["vol_file"] = vol_path

        # Write companion JSON with CAD reference values
        json_path = vol_path + ".json"
        cad_ref = {
            "vol_file": os.path.basename(vol_path),
            "order": order,
            "materials": {m["name"]: m["cad_volume"]
                          for m in cad_materials},
            "boundaries": {b["name"]: b["cad_area"]
                           for b in cad_boundaries},
            "edges": {e["name"]: e["cad_length"]
                      for e in cad_edges},
        }
        with open(json_path, "w") as f:
            json.dump(cad_ref, f, indent=2)
        result["json_file"] = json_path

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Export Cubit mesh as Netgen .vol with consistency check")
    parser.add_argument("--cub5", required=True)
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--vol", default=None, help="Output .vol path")

    args = parser.parse_args()
    result = export_netgen(args.cub5, args.order, args.vol)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
