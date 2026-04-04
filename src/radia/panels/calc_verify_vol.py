"""
Verify exported Netgen .vol mesh against CAD reference JSON.

Called as subprocess from Cubit Export Mesh GUI (after C++ exports .vol + .json):
    python calc_export_netgen.py --vol model.vol

Reads .vol (NGSolve), reads companion .vol.json (CAD values from C++),
computes NGSolve Integrate per label, updates JSON with ng_* values and errors.

Outputs updated JSON to stdout (last line).
Does NOT require Cubit — only NGSolve.
"""

import argparse
import json
import os
import sys


def verify_vol(vol_path):
    """Read .vol + companion JSON, compute NGSolve values, return result dict."""
    from ngsolve import Mesh, Integrate, CF, BND

    json_path = vol_path + ".json"
    if not os.path.exists(json_path):
        return {"error": f"Companion JSON not found: {json_path}"}

    with open(json_path, "r") as f:
        cad_ref = json.load(f)

    mesh = Mesh(vol_path)
    warnings = []

    # --- Per-material volume ---
    cad_mats = cad_ref.get("materials", {})
    materials = []
    for mat in dict.fromkeys(mesh.GetMaterials()):
        ng_vol = Integrate(CF(1), mesh, definedon=mesh.Materials(mat))
        entry = {"name": mat, "ng_volume": ng_vol}
        if mat in cad_mats:
            cad_v = cad_mats[mat]
            entry["cad_volume"] = cad_v
            if cad_v > 0:
                err = (ng_vol - cad_v) / cad_v * 100
                entry["error_pct"] = err
                if abs(err) > 1.0:
                    warnings.append(
                        f"Material \"{mat}\": V error {err:+.2e}%")
        materials.append(entry)

    # --- Per-boundary area ---
    cad_bnds = cad_ref.get("boundaries", {})
    boundaries = []
    for bnd in dict.fromkeys(mesh.GetBoundaries()):
        ng_area = Integrate(CF(1), mesh, BND,
                            definedon=mesh.Boundaries(bnd))
        entry = {"name": bnd, "ng_area": ng_area}
        if bnd in cad_bnds:
            cad_a = cad_bnds[bnd]
            entry["cad_area"] = cad_a
            if cad_a > 0:
                err = (ng_area - cad_a) / cad_a * 100
                entry["error_pct"] = err
                if abs(err) > 1.0:
                    warnings.append(
                        f"Boundary \"{bnd}\": A error {err:+.2e}%")
        boundaries.append(entry)

    # --- Per-edge length (BBND) ---
    cad_edges = cad_ref.get("edges", {})
    edges = []
    try:
        from ngsolve import BBND
        for ename in dict.fromkeys(mesh.GetBBoundaries()):
            ng_len = Integrate(CF(1), mesh, BBND,
                               definedon=mesh.BBoundaries(ename))
            entry = {"name": ename, "ng_length": ng_len}
            if ename in cad_edges:
                cad_l = cad_edges[ename]
                entry["cad_length"] = cad_l
                if cad_l > 0:
                    err = (ng_len - cad_l) / cad_l * 100
                    entry["error_pct"] = err
            edges.append(entry)

        # Total length check
        cad_length_total = sum(cad_edges.values())
        if cad_length_total > 0:
            ng_length_total = Integrate(CF(1), mesh, BBND)
            total_err = (ng_length_total - cad_length_total) / cad_length_total * 100
            if abs(total_err) > 1.0:
                warnings.append(f"Total edge length: {total_err:+.2e}%")
    except Exception:
        pass

    # Also include CAD-only edges (in JSON but not in mesh BBoundaries)
    ng_edge_names = {e["name"] for e in edges}
    for ename, cad_l in cad_edges.items():
        if ename not in ng_edge_names:
            edges.append({"name": ename, "cad_length": cad_l})

    result = {
        "vol_file": vol_path,
        "order": cad_ref.get("order", 0),
        "n_elements": mesh.ne,
        "materials": materials,
        "boundaries": boundaries,
        "edges": edges,
        "warnings": warnings,
    }

    # Update companion JSON with NGSolve values
    cad_ref["ng_materials"] = {m["name"]: m["ng_volume"] for m in materials}
    cad_ref["ng_boundaries"] = {b["name"]: b["ng_area"] for b in boundaries}
    cad_ref["ng_edges"] = {e["name"]: e.get("ng_length", 0) for e in edges
                           if "ng_length" in e}
    cad_ref["warnings"] = warnings
    with open(json_path, "w") as f:
        json.dump(cad_ref, f, indent=2)

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Verify .vol mesh against CAD reference JSON")
    parser.add_argument("--vol", required=True, help="Path to .vol mesh file")

    args = parser.parse_args()
    result = verify_vol(args.vol)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
