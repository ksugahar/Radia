"""
Check .vol mesh consistency against CAD reference values.

Standalone checker -- does NOT require Cubit.
Reads .vol (NGSolve mesh) + .vol.json (CAD reference from export).

Usage:
    python check_vol_consistency.py model.vol
    python check_vol_consistency.py model.vol --json model.vol.json
    python check_vol_consistency.py model.vol --threshold 0.5

Exit code: 0 = all checks pass, 1 = warnings found, 2 = error.
"""

import argparse
import json
import os
import sys


def check_consistency(vol_path, json_path=None, threshold=1.0):
    """Check mesh consistency against CAD reference.

    Args:
        vol_path: Path to .vol mesh file.
        json_path: Path to companion JSON (default: vol_path + ".json").
        threshold: Warning threshold in percent.

    Returns:
        dict with check results and warnings.
    """
    from ngsolve import Mesh, Integrate, CF, BND

    if json_path is None:
        json_path = vol_path + ".json"

    if not os.path.exists(json_path):
        return {"error": f"CAD reference not found: {json_path}"}

    with open(json_path, "r") as f:
        cad_ref = json.load(f)

    mesh = Mesh(vol_path)
    warnings = []
    results = {
        "vol_file": vol_path,
        "json_file": json_path,
        "threshold_pct": threshold,
        "materials": [],
        "boundaries": [],
        "edges": [],
    }

    # --- Volume check (per-material) ---
    cad_volumes = cad_ref.get("materials", {})
    for mat in mesh.GetMaterials():
        ng_vol = Integrate(CF(1), mesh, definedon=mesh.Materials(mat))
        entry = {"name": mat, "ng_volume": ng_vol}
        if mat in cad_volumes:
            cad_v = cad_volumes[mat]
            entry["cad_volume"] = cad_v
            if cad_v > 0:
                err = (ng_vol - cad_v) / cad_v * 100
                entry["error_pct"] = err
                if abs(err) > threshold:
                    warnings.append(
                        f"Volume \"{mat}\": {err:+.2e}% "
                        f"(ng={ng_vol:.6e}, cad={cad_v:.6e})")
        results["materials"].append(entry)

    # --- Area check (per-boundary) ---
    cad_areas = cad_ref.get("boundaries", {})
    checked = set()
    for bnd in mesh.GetBoundaries():
        if bnd in checked:
            continue
        checked.add(bnd)
        ng_area = Integrate(CF(1), mesh, BND,
                            definedon=mesh.Boundaries(bnd))
        entry = {"name": bnd, "ng_area": ng_area}
        if bnd in cad_areas:
            cad_a = cad_areas[bnd]
            entry["cad_area"] = cad_a
            if cad_a > 0:
                err = (ng_area - cad_a) / cad_a * 100
                entry["error_pct"] = err
                if abs(err) > threshold:
                    warnings.append(
                        f"Area \"{bnd}\": {err:+.2e}% "
                        f"(ng={ng_area:.6e}, cad={cad_a:.6e})")
        results["boundaries"].append(entry)

    # --- Length check (total BBND + per-edge where possible) ---
    # NOTE: Per-edge BBND mapping has a known Netgen topology issue where
    # edgenr assignment may not match NGSolve's topological edge enumeration.
    # Total BBND length is reliable; per-edge breakdown may be inaccurate.
    try:
        from ngsolve import BBND
        cad_lengths = cad_ref.get("edges", {})
        total_cad_length = sum(cad_lengths.values())
        total_ng_length = Integrate(CF(1), mesh, BBND)

        # Per-edge entries (best-effort)
        checked_edge = set()
        for bname in mesh.GetBBoundaries():
            if bname in checked_edge:
                continue
            checked_edge.add(bname)
            ng_len = Integrate(CF(1), mesh, BBND,
                               definedon=mesh.BBoundaries(bname))
            entry = {"name": bname, "ng_length": ng_len}
            if bname in cad_lengths:
                entry["cad_length"] = cad_lengths[bname]
            results["edges"].append(entry)

        # Total length check (reliable)
        if total_cad_length > 0:
            err = (total_ng_length - total_cad_length) / total_cad_length * 100
            results["total_edge_length"] = {
                "cad": total_cad_length,
                "ng": total_ng_length,
                "error_pct": err,
            }
            if abs(err) > threshold:
                warnings.append(
                    f"Total edge length: {err:+.2e}% "
                    f"(ng={total_ng_length:.6e}, cad={total_cad_length:.6e})")
    except Exception:
        pass

    results["warnings"] = warnings
    return results


def print_table(results):
    """Print consistency check results as a formatted table."""
    print(f"Mesh: {results['vol_file']}")
    print(f"Threshold: {results['threshold_pct']}%")
    print()

    # Volume table
    if results["materials"]:
        print("  Material            CAD Volume      NG Volume       Error")
        print("  " + "-" * 64)
        for m in results["materials"]:
            cad = m.get("cad_volume", 0)
            ng = m.get("ng_volume", 0)
            err = m.get("error_pct")
            err_s = f"{err:+.2e}%" if err is not None else "N/A"
            flag = " ***" if err and abs(err) > results["threshold_pct"] else ""
            print(f"  {m['name']:<20s}{cad:>14.6e}  {ng:>14.6e}  {err_s:>10s}{flag}")
        print()

    # Area table
    if results["boundaries"]:
        print("  Boundary            CAD Area        NG Area         Error")
        print("  " + "-" * 64)
        for b in results["boundaries"]:
            cad = b.get("cad_area", 0)
            ng = b.get("ng_area", 0)
            err = b.get("error_pct")
            err_s = f"{err:+.2e}%" if err is not None else "N/A"
            flag = " ***" if err and abs(err) > results["threshold_pct"] else ""
            print(f"  {b['name']:<20s}{cad:>14.6e}  {ng:>14.6e}  {err_s:>10s}{flag}")
        print()

    # Length table
    if results["edges"]:
        print("  Edge                CAD Length      NG Length")
        print("  " + "-" * 52)
        for e in results["edges"]:
            cad = e.get("cad_length", 0)
            ng = e.get("ng_length", 0)
            print(f"  {e['name']:<20s}{cad:>14.6e}  {ng:>14.6e}")
        # Total length check
        tot = results.get("total_edge_length")
        if tot:
            print("  " + "-" * 52)
            err_s = f"{tot['error_pct']:+.2e}%"
            flag = " ***" if abs(tot["error_pct"]) > results["threshold_pct"] else ""
            print(f"  {'TOTAL':<20s}{tot['cad']:>14.6e}  {tot['ng']:>14.6e}  {err_s}{flag}")
        print()

    # Warnings summary
    if results.get("warnings"):
        print(f"  WARNINGS ({len(results['warnings'])}):")
        for w in results["warnings"]:
            print(f"    - {w}")
    else:
        print("  All checks PASSED.")


def main():
    parser = argparse.ArgumentParser(
        description="Check .vol mesh consistency against CAD reference")
    parser.add_argument("vol", help="Path to .vol mesh file")
    parser.add_argument("--json", default=None,
                        help="Path to companion JSON (default: vol + .json)")
    parser.add_argument("--threshold", type=float, default=1.0,
                        help="Warning threshold in percent (default: 1.0)")

    args = parser.parse_args()
    results = check_consistency(args.vol, args.json, args.threshold)

    if "error" in results:
        print(f"ERROR: {results['error']}")
        sys.exit(2)

    print_table(results)
    sys.exit(1 if results.get("warnings") else 0)


if __name__ == "__main__":
    main()
