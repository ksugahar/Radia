"""Validate Cubit-exported BDF files with the independent pyNastran reader.

This is a validation-class test, not a runtime dependency of
``cubit-mesh-export``.  Generate the corpus first with
``test_export_combinations.py``, then run this script in an environment that
contains pyNastran::

    python -m pip install "pyNastran==1.4.1"
    python validation_test/cubit/validate_nastran_with_pynastran.py

The exporter deliberately writes a mesh/property interchange deck without
inventing physical material constants.  Missing MAT cards are therefore
reported as a downstream completion requirement, not treated as a parse
failure.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

REPORT_SCHEMA = "cubit-mesh-export.nastran-interchange.v1"


def _material_id(prop):
    """Return the material ID referenced by a supported property card."""
    if prop.type == "PSOLID":
        return prop.mid
    if prop.type == "PSHELL":
        return prop.mid1
    return None


def validate_file(path: Path, bdf_class) -> dict:
    """Parse one BDF and check mesh-level referential integrity."""
    row = {"file": path.name, "passed": False}
    try:
        model = bdf_class(debug=False, log=None)
        model.read_bdf(str(path), xref=False, punch=False)

        missing_node_refs = sorted({
            node_id
            for element in model.elements.values()
            for node_id in element.node_ids
            if node_id is not None and node_id not in model.nodes
        })
        missing_property_refs = sorted({
            element.Pid()
            for element in model.elements.values()
            if element.Pid() not in model.properties
        })
        nonfinite_nodes = sorted(
            node_id
            for node_id, node in model.nodes.items()
            if not all(math.isfinite(float(value)) for value in node.xyz)
        )
        referenced_materials = sorted({
            material_id
            for prop in model.properties.values()
            if (material_id := _material_id(prop)) is not None
        })
        missing_materials = [
            material_id for material_id in referenced_materials
            if material_id not in model.materials
        ]

        errors = []
        if not model.nodes:
            errors.append("no GRID cards")
        if not model.elements:
            errors.append("no supported element cards")
        if missing_node_refs:
            errors.append(f"missing GRID references: {missing_node_refs}")
        if missing_property_refs:
            errors.append(f"missing property references: {missing_property_refs}")
        if nonfinite_nodes:
            errors.append(f"non-finite GRID coordinates: {nonfinite_nodes}")

        row.update({
            "passed": not errors,
            "nodes": len(model.nodes),
            "elements": dict(sorted(Counter(
                element.type for element in model.elements.values()
            ).items())),
            "properties": dict(sorted(Counter(
                prop.type for prop in model.properties.values()
            ).items())),
            "sets": {
                str(set_id): len(mesh_set.ids)
                for set_id, mesh_set in sorted(model.sets.items())
            },
            "missing_material_ids": missing_materials,
            "interchange_ready": not errors,
            "analysis_deck_complete": not missing_materials,
            "errors": errors,
        })
    except Exception as exc:  # noqa: BLE001 - retain all parser failures
        row.update({
            "interchange_ready": False,
            "analysis_deck_complete": False,
            "errors": [f"{type(exc).__name__}: {exc}"],
        })
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    default_dir = Path(__file__).resolve().parent / "export_test_output"
    parser.add_argument("input_dir", nargs="?", type=Path, default=default_dir)
    parser.add_argument("--output-json", type=Path)
    args = parser.parse_args()

    try:
        import pyNastran
        from pyNastran.bdf.bdf import BDF
    except ImportError as exc:
        print(
            "pyNastran is required for this independent validation; install "
            'it with: python -m pip install "pyNastran==1.4.1"',
            file=sys.stderr,
        )
        raise SystemExit(2) from exc

    files = sorted(args.input_dir.glob("*.bdf"))
    if not files:
        print(f"No BDF files found in {args.input_dir}", file=sys.stderr)
        return 2

    rows = [validate_file(path, BDF) for path in files]
    passed = sum(row["passed"] for row in rows)
    report = {
        "schema": REPORT_SCHEMA,
        "generated_utc": datetime.now(UTC).isoformat(),
        "reader": {"name": "pyNastran", "version": pyNastran.__version__},
        "input_dir": str(args.input_dir.resolve()),
        "summary": {
            "passed": passed,
            "failed": len(rows) - passed,
            "total": len(rows),
            "interchange_ready": sum(row["interchange_ready"] for row in rows),
            "analysis_deck_complete": sum(
                row["analysis_deck_complete"] for row in rows
            ),
            "mesh_interchange_only": sum(
                not row["analysis_deck_complete"] for row in rows
            ),
        },
        "files": rows,
    }

    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(rendered + "\n", encoding="utf-8")
    return 0 if passed == len(rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
