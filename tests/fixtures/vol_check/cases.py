"""Deterministic check-vol cases: committed fixtures plus small text mutations.

Shared by the check-vol contract test and the script that records its
baseline ``expected.json`` from radia.vol_check.
Every case is (vol_bytes, vol_suffix, contract_dict_or_None).
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
HEX = ROOT / "tests/fixtures/hdiv/hex_pair.vol"
HEX_LABELS = ROOT / "tests/fixtures/hdiv/hex_pair.labels.json"


def _hex():
    return HEX.read_text(encoding="utf-8")


def _replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f"mutation anchor not unique: {old!r}")
    return text.replace(old, new)


def _surface_rows(text, transform):
    head, rest = text.split("surfaceelements\n", 1)
    count, rest = rest.split("\n", 1)
    rows, tail = rest.split("\n\n", 1)
    rows = transform(rows.split("\n"))
    return f"{head}surfaceelements\n{len(rows)}\n" + "\n".join(rows) + "\n\n" + tail


def _contract(required_materials=("body",), required_boundaries=("exterior",),
              allowed_materials=("body",), allowed_boundaries=("exterior",)):
    return {"schema": "radia.vol-label-contract.v1",
            "required": {"materials": list(required_materials), "boundaries": list(required_boundaries)},
            "allowed": {"materials": list(allowed_materials), "boundaries": list(allowed_boundaries)}}


def cases():
    """Return {name: (vol_bytes, suffix, contract or None)}."""
    hex_text = _hex()
    hex_contract = json.loads(HEX_LABELS.read_text(encoding="utf-8"))
    out = {}

    def add(name, text_or_bytes, contract, suffix=".vol"):
        data = text_or_bytes if isinstance(text_or_bytes, bytes) else text_or_bytes.encode("utf-8")
        out[name] = (data, suffix, contract)

    # Committed fixtures, as stored.
    add("hex_pair", HEX.read_bytes(), hex_contract)
    tph = ROOT / "tests/fixtures/hdiv/tet_pyramid_hex"
    add("tet_pyramid_hex", (tph.with_suffix(".vol")).read_bytes(),
        json.loads((ROOT / "tests/fixtures/hdiv/tet_pyramid_hex.labels.json").read_text(encoding="utf-8")))
    prolate = ROOT / "tests/feec/data/prolate_named"
    add("prolate_named", (prolate / "prolate_named.vol").read_bytes(),
        json.loads((prolate / "prolate-label-contract.json").read_text(encoding="utf-8")))
    add("affine_field_tetra", (ROOT / "tests/fixtures/beam/affine_field_tetra.vol").read_bytes(), None)
    add("mixed_cylinder_o2", (ROOT / "tests/fixtures/cubit/mixed_cylinder_o2.vol.gz").read_bytes(), None,
        suffix=".vol.gz")

    # Label contract violations on a valid mesh.
    add("labels_missing_required_material", hex_text, _contract(required_materials=("body", "coil"),
                                                                allowed_materials=("body", "coil")))
    add("labels_unexpected_boundary", hex_text, _contract(required_boundaries=(), allowed_boundaries=("outer",)))
    add("labels_unnamed_material", _replace_once(hex_text, "materials\n1\n1 body\n", "materials\n0\n"),
        hex_contract)

    # Topology violations.
    add("topology_duplicate_surface_row", _surface_rows(hex_text, lambda r: r + [r[0]]), hex_contract)
    add("topology_unknown_surface_domain",
        _surface_rows(hex_text, lambda r: [r[0].replace(" 2 1 1 0 4", " 2 1 3 0 4", 1)] + r[1:]), hex_contract)
    add("topology_zero_zero_surface",
        _surface_rows(hex_text, lambda r: [r[0].replace(" 2 1 1 0 4", " 2 1 0 0 4", 1)] + r[1:]), hex_contract)
    add("topology_missing_exterior_face", _surface_rows(hex_text, lambda r: r[1:]), hex_contract)
    add("topology_inverted_hex",
        _replace_once(hex_text, "1 8 1 2 3 4 5 6 7 8\n", "1 8 5 6 7 8 1 2 3 4\n"), hex_contract)
    add("topology_unused_volume_domain",
        _replace_once(hex_text, "1 8 5 6 7 8 9 10 11 12\n", "2 8 5 6 7 8 9 10 11 12\n"), hex_contract)
    return out


def write_case(directory, name, case):
    """Write one case as mesh.vol[.gz] (+ labels.json); return (mesh, contract or None)."""
    data, suffix, contract = case
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    mesh = directory / ("mesh" + suffix)
    mesh.write_bytes(data)
    path = None
    if contract is not None:
        path = directory / "labels.json"
        path.write_text(json.dumps(contract), encoding="utf-8")
    return mesh, path


VOLATILE = ("generated_at_utc", "versions", "vol_file", "json_file")


def normalize(report):
    """Drop run-specific fields; keep every verdict, count and measured value."""
    return {k: v for k, v in report.items() if k not in VOLATILE}
