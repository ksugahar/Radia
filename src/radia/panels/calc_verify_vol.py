"""
Verify exported Netgen .vol mesh against CAD reference JSON.

Called as subprocess from Cubit Export Mesh GUI (after C++ exports .vol + .json):
    python calc_verify_vol.py --vol model.vol

Delegates all checks to ``cubit_mesh_export.check``, then preserves the Cubit
toolbar's historical sidecar update with measured ``ng_*`` values.

Outputs the versioned check report to stdout (last line). Does not require
Cubit; it requires cubit-mesh-export and NGSolve in the external Python.
"""

import argparse
import json


def verify_vol(vol_path):
    """Run the canonical checker and preserve the toolbar sidecar update."""
    try:
        from cubit_mesh_export.check import check_consistency

        json_path = vol_path + ".json"
        result = check_consistency(vol_path, json_path=json_path)
    except Exception as exc:
        return {"error": str(exc)}

    # The Cubit export dialog historically augments the CAD sidecar with the
    # measured NGSolve values. Keep that UI contract while all validation logic
    # remains in cubit_mesh_export.check.
    with open(json_path, "r", encoding="utf-8-sig") as stream:
        cad_ref = json.load(stream)
    cad_ref["ng_materials"] = {
        row["name"]: row["ng_volume"] for row in result["materials"]
    }
    cad_ref["ng_boundaries"] = {
        row["name"]: row["ng_area"] for row in result["boundaries"]
    }
    cad_ref["ng_edges"] = {
        row["name"]: row["ng_length"]
        for row in result["edges"] if "ng_length" in row
    }
    cad_ref["warnings"] = result["warnings"]
    with open(json_path, "w", encoding="utf-8") as stream:
        json.dump(cad_ref, stream, indent=2)

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
