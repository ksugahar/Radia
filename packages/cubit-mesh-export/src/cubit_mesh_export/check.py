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

import numpy as np


def _labels(value):
    if value is None:
        return set()
    if isinstance(value, str):
        return {part.strip() for part in value.split(",") if part.strip()}
    return {str(part) for part in value}


def _node_number(node):
    number = getattr(node, "nr", None)
    return int(node) if number is None else int(number)


def _element_number(element):
    number = getattr(element, "nr", None)
    return int(getattr(element, "index")) if number is None else int(number)


def _classify_conductor_face(volume_records, conductive, air):
    conductive_records = [row for row in volume_records if row[1] in conductive]
    if not conductive_records:
        return "nonconductive"
    if len(conductive_records) >= 2:
        materials = {row[1] for row in conductive_records}
        return "conductor-conductor" if len(materials) == 1 else "conductive-interface"

    neighbors = [row[1] for row in volume_records if row[1] not in conductive]
    if any(material in air for material in neighbors):
        return "conductor-air"
    if neighbors:
        return "conductor-insulator"
    return "conductor-exterior"


def _face_adjacency_quality(mesh, conductive_materials, air_materials,
                            sibc_boundaries, require_all_sibc_labeled):
    import ngsolve as ng

    conductive = _labels(conductive_materials)
    air = _labels(air_materials)
    sibc = _labels(sibc_boundaries)
    if not conductive:
        return {
            "enabled": False,
            "conductive_materials": [],
            "air_materials": sorted(air),
        }, []

    face_to_volume = {}
    for element in mesh.Elements(ng.VOL):
        record = (_element_number(element), str(element.mat))
        for face in element.faces:
            face_to_volume.setdefault(_node_number(face), []).append(record)

    face_to_boundary = {}
    for element in mesh.Elements(ng.BND):
        record = (_element_number(element), str(element.mat))
        for face in element.faces:
            face_to_boundary.setdefault(_node_number(face), []).append(record)

    role_counts = {}
    marked_non_sibc = []
    unlabeled_sibc = []
    for face_nr in sorted(set(face_to_volume) | set(face_to_boundary)):
        volume_records = tuple(face_to_volume.get(face_nr, ()))
        boundary_records = tuple(face_to_boundary.get(face_nr, ()))
        role = _classify_conductor_face(volume_records, conductive, air)
        role_counts[role] = role_counts.get(role, 0) + 1
        boundary_labels = {row[1] for row in boundary_records}
        is_sibc_role = role in {"conductor-air", "conductor-exterior"}
        marked = bool(boundary_labels & sibc)
        if marked and not is_sibc_role:
            marked_non_sibc.append({
                "face_nr": face_nr,
                "role": role,
                "boundary_labels": sorted(boundary_labels),
            })
        if require_all_sibc_labeled and is_sibc_role and not marked:
            unlabeled_sibc.append({
                "face_nr": face_nr,
                "role": role,
                "boundary_labels": sorted(boundary_labels),
            })

    warnings = []
    if marked_non_sibc:
        warnings.append(
            f"{len(marked_non_sibc)} SIBC-labelled faces are not conductor-air/exterior faces")
    if unlabeled_sibc:
        warnings.append(
            f"{len(unlabeled_sibc)} conductor-air/exterior faces lack an SIBC label")

    return {
        "enabled": True,
        "conductive_materials": sorted(conductive),
        "air_materials": sorted(air),
        "sibc_boundaries": sorted(sibc),
        "face_count": sum(role_counts.values()),
        "face_role_counts": dict(sorted(role_counts.items())),
        "sibc_candidate_face_count": (
            role_counts.get("conductor-air", 0)
            + role_counts.get("conductor-exterior", 0)
        ),
        "loop_bridge_face_count": (
            role_counts.get("conductor-conductor", 0)
            + role_counts.get("conductive-interface", 0)
        ),
        "marked_non_sibc_faces": marked_non_sibc[:50],
        "unlabeled_sibc_faces": unlabeled_sibc[:50],
    }, warnings


def check_mesh_quality(mesh_or_path, *, min_scaled_jacobian=1.0e-6,
                       jacobian_tolerance=0.0, integration_order=None,
                       min_curve_order=1, require_tetrahedra=False,
                       required_materials=(), required_boundaries=(),
                       conductive_materials=(), air_materials=("air", "vacuum"),
                       sibc_boundaries=(), require_all_sibc_labeled=False):
    """Check curved-element mappings, labels, and conductor-face adjacency.

    Unlike a corner-only tetrahedron metric, this gate samples the actual
    NGSolve element transformation at a high-order integration rule.  It can
    therefore detect an inverted or nearly collapsed mapping inside a curved
    element even when the linear corner tetrahedron is valid.
    """
    import ngsolve as ng

    mesh = ng.Mesh(str(mesh_or_path)) if isinstance(mesh_or_path, (str, os.PathLike)) else mesh_or_path
    curve_order = int(mesh.GetCurveOrder())
    sample_order = (
        max(2, 2 * curve_order + 2)
        if integration_order is None else int(integration_order)
    )
    if sample_order < 1:
        raise ValueError("integration_order must be positive")
    if min_scaled_jacobian < 0.0 or min_scaled_jacobian > 1.0:
        raise ValueError("min_scaled_jacobian must lie in [0, 1]")

    min_det = float("inf")
    max_det = float("-inf")
    min_scaled = float("inf")
    sample_count = 0
    invalid_count = 0
    low_scaled_count = 0
    invalid_elements = set()
    low_scaled_elements = set()
    volume_element_count = 0
    tetrahedron_count = 0

    for element in mesh.Elements(ng.VOL):
        volume_element_count += 1
        if element.type == ng.ET.TET:
            tetrahedron_count += 1
        trafo = mesh.GetTrafo(element)
        rule = ng.IntegrationRule(element.type, sample_order)
        for point in rule:
            jacobian = np.asarray(trafo(point).jacobi, dtype=float)
            determinant = float(np.linalg.det(jacobian))
            column_norm_product = float(np.prod(np.linalg.norm(jacobian, axis=0)))
            scaled = (
                determinant / column_norm_product
                if column_norm_product > 0.0 else float("-inf")
            )
            sample_count += 1
            min_det = min(min_det, determinant)
            max_det = max(max_det, determinant)
            min_scaled = min(min_scaled, scaled)
            if not np.isfinite(determinant) or determinant <= jacobian_tolerance:
                invalid_count += 1
                invalid_elements.add(_element_number(element))
            if not np.isfinite(scaled) or scaled < min_scaled_jacobian:
                low_scaled_count += 1
                low_scaled_elements.add(_element_number(element))

    materials = sorted(set(mesh.GetMaterials()))
    boundaries = sorted(set(mesh.GetBoundaries()))
    missing_materials = sorted(_labels(required_materials) - set(materials))
    missing_boundaries = sorted(_labels(required_boundaries) - set(boundaries))
    adjacency, adjacency_warnings = _face_adjacency_quality(
        mesh,
        conductive_materials,
        air_materials,
        sibc_boundaries,
        require_all_sibc_labeled,
    )

    warnings = list(adjacency_warnings)
    if curve_order < int(min_curve_order):
        warnings.append(f"curve order {curve_order} is below required order {min_curve_order}")
    if invalid_count:
        warnings.append(f"{invalid_count} mapping samples have a non-positive Jacobian")
    if low_scaled_count:
        warnings.append(
            f"{low_scaled_count} mapping samples are below scaled Jacobian {min_scaled_jacobian:g}")
    if require_tetrahedra and tetrahedron_count != volume_element_count:
        warnings.append(
            f"{volume_element_count - tetrahedron_count} volume elements are not tetrahedra")
    if missing_materials:
        warnings.append("missing materials: " + ", ".join(missing_materials))
    if missing_boundaries:
        warnings.append("missing boundaries: " + ", ".join(missing_boundaries))

    return {
        "passed": not warnings,
        "curve_order": curve_order,
        "integration_order": sample_order,
        "volume_element_count": volume_element_count,
        "tetrahedron_count": tetrahedron_count,
        "mapping_sample_count": sample_count,
        "minimum_jacobian": None if sample_count == 0 else min_det,
        "maximum_jacobian": None if sample_count == 0 else max_det,
        "minimum_scaled_jacobian": None if sample_count == 0 else min_scaled,
        "invalid_jacobian_sample_count": invalid_count,
        "low_scaled_jacobian_sample_count": low_scaled_count,
        "invalid_jacobian_elements": sorted(invalid_elements)[:50],
        "low_scaled_jacobian_elements": sorted(low_scaled_elements)[:50],
        "materials": materials,
        "boundaries": boundaries,
        "missing_materials": missing_materials,
        "missing_boundaries": missing_boundaries,
        "adjacency": adjacency,
        "warnings": warnings,
    }


def check_consistency(vol_path, json_path=None, threshold=1.0, **quality_options):
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

    results["quality"] = check_mesh_quality(mesh, **quality_options)
    warnings.extend(results["quality"]["warnings"])

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

    quality = results.get("quality")
    if quality:
        print("  Curved mapping quality")
        print("  " + "-" * 64)
        print(
            f"  curve order={quality['curve_order']}, "
            f"elements={quality['volume_element_count']}, "
            f"samples={quality['mapping_sample_count']}")
        print(
            f"  min det(J)={quality['minimum_jacobian']:.6e}, "
            f"min scaled J={quality['minimum_scaled_jacobian']:.6e}")
        adjacency = quality.get("adjacency", {})
        if adjacency.get("enabled"):
            print(f"  face roles={adjacency['face_role_counts']}")
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
    parser.add_argument("--min-scaled-jacobian", type=float, default=1.0e-6,
                        help="Minimum sampled scaled Jacobian (default: 1e-6)")
    parser.add_argument("--integration-order", type=int, default=None,
                        help="Curved-map sampling order (default: 2*curve_order+2)")
    parser.add_argument("--min-curve-order", type=int, default=1)
    parser.add_argument("--require-tet", action="store_true")
    parser.add_argument("--required-materials", default="")
    parser.add_argument("--required-boundaries", default="")
    parser.add_argument("--conductors", default="",
                        help="Comma-separated conductive material labels")
    parser.add_argument("--air-materials", default="air,vacuum")
    parser.add_argument("--sibc-boundaries", default="",
                        help="Comma-separated SIBC boundary labels")
    parser.add_argument("--require-all-sibc-labeled", action="store_true")

    args = parser.parse_args()
    results = check_consistency(
        args.vol,
        args.json,
        args.threshold,
        min_scaled_jacobian=args.min_scaled_jacobian,
        integration_order=args.integration_order,
        min_curve_order=args.min_curve_order,
        require_tetrahedra=args.require_tet,
        required_materials=args.required_materials,
        required_boundaries=args.required_boundaries,
        conductive_materials=args.conductors,
        air_materials=args.air_materials,
        sibc_boundaries=args.sibc_boundaries,
        require_all_sibc_labeled=args.require_all_sibc_labeled,
    )

    if "error" in results:
        print(f"ERROR: {results['error']}")
        sys.exit(2)

    print_table(results)
    sys.exit(1 if results.get("warnings") else 0)


if __name__ == "__main__":
    main()
