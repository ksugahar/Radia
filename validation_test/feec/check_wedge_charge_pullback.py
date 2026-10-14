"""Check saved WEDGE pullback lane samples, without importing either solver.

Consumes the native/radia lane JSON from the experimental
validate_hdiv_vim_wedge_charge_pullback_native driver. This only checks sampled
agreement: no basis-rank proof, runtime authentication, Gram root-cause claim,
or campaign acceptance follows. The native driver remains separately reviewed.

Usage: python check_wedge_charge_pullback.py --native-json native.json
       --radia-json radia.json --output checked.json
Exit 0 means sampled agreement only; exit 1 means mismatch or invalid input.
Input hashes bind the report to bytes read, not to an authenticated runtime.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

CASES = {"two_structured_wedges": 2, "explicit_full_four_wedges": 4}
VOLUME_POINTS = ((.17, .23, .19), (.41, .12, .67), (.08, .61, .83))
TRI_POINTS = ((.17, .23), (.41, .12), (.08, .61))
QUAD_POINTS = ((.17, .23), (.41, .67), (.83, .19))


def _number(value):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError("expected a finite real number")
    return value


def _positive_int(value):
    if type(value) is not int or value <= 0:
        raise ValueError("expected a positive integer")
    return value


def _vector(value, size):
    if not isinstance(value, (list, tuple)) or len(value) != size:
        raise ValueError(f"expected {size} coordinates")
    return tuple(_number(v) for v in value)


def _rows(case, family):
    count = _positive_int(case["elements" if family == "volume" else "boundary_faces"])
    rows = case[family]
    if not isinstance(rows, list) or len(rows) != count * 3:
        raise ValueError("incomplete sample coverage")
    indexed, face_types = {}, {}
    for row in rows:
        host = row["host"]
        if type(host) is not int or not 0 <= host < count:
            raise ValueError("invalid sample host")
        if family == "boundary":
            triangle = row["triangle"]
            if type(triangle) is not bool:
                raise ValueError("triangle must be boolean")
            if host in face_types and face_types[host] != triangle:
                raise ValueError("inconsistent face type")
            face_types[host] = triangle
            points = TRI_POINTS if triangle else QUAD_POINTS
        else:
            points = VOLUME_POINTS
        point = _vector(row["point"], len(points[0]))
        key = (host, point)
        if point not in points or key in indexed:
            raise ValueError("unexpected or duplicate sample")
        indexed[key] = (_vector(row["physical_point"], 3),
                        _number(row["reference_charge"]))
    return indexed, face_types


def compare(native, radia):
    """Return fail-closed sample diagnostics; never grant numerical acceptance."""
    result = {
        "schema": "radia.validation.wedge-pullback-samples.v1",
        "status": "invalid_input", "sample_checks_passed": False,
        "acceptance_evidence": False, "campaign_level": False,
        "runtime_provenance_verified": False, "basis_rank_verified": False,
        "gram_root_cause_established": False,
    }
    try:
        if native["lane"] != "ngsolve-native" or radia["lane"] != "radia-charge-pullback":
            raise ValueError("incorrect lane identity")
        for field in ("python", "ngsolve"):
            if not isinstance(native[field], str) or not native[field].strip():
                raise ValueError(f"missing {field} version")
            if native[field] != radia[field]:
                raise ValueError(f"different {field} versions")
        if set(native["cases"]) != set(CASES) or set(radia["cases"]) != set(CASES):
            raise ValueError("required case set mismatch")
        metrics = {"map_absolute": 0.0, "charge_absolute": 0.0, "charge_relative": 0.0}
        for name, elements in CASES.items():
            ref, actual = native["cases"][name], radia["cases"][name]
            for field in ("elements", "boundary_faces", "ndof"):
                if _positive_int(ref[field]) != _positive_int(actual[field]):
                    raise ValueError(f"{name}: {field} mismatch")
            if ref["elements"] != elements:
                raise ValueError("unexpected fixture element count")
            _positive_int(actual["charge_modes"])
            for family in ("volume", "boundary"):
                refs, ref_types = _rows(ref, family)
                acts, act_types = _rows(actual, family)
                if refs.keys() != acts.keys() or ref_types != act_types:
                    raise ValueError("sample identity mismatch")
                for key, (position, charge) in refs.items():
                    actual_position, actual_charge = acts[key]
                    absolute = abs(actual_charge - charge)
                    measured = {
                        "map_absolute": math.dist(position, actual_position),
                        "charge_absolute": absolute,
                        "charge_relative": absolute / max(abs(charge), abs(actual_charge), 1e-14),
                    }
                    for metric, value in measured.items():
                        metrics[metric] = max(metrics[metric], _number(value))
        checks = {"map_absolute": metrics["map_absolute"] < 1e-13,
                  "charge_absolute": metrics["charge_absolute"] < 1e-10,
                  "charge_relative": metrics["charge_relative"] < 1e-10}
        passed = all(checks.values())
        result.update(status="sampled_agreement" if passed else "sampled_mismatch",
                      sample_checks_passed=passed, metrics=metrics, checks=checks)
    except (KeyError, TypeError, ValueError, OverflowError) as error:
        result["error"] = str(error)
    return result


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-json", required=True, type=Path)
    parser.add_argument("--radia-json", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    input_paths = (args.native_json, args.radia_json)
    if any(args.output.resolve() == path.resolve() or (
            args.output.exists() and path.exists() and args.output.samefile(path))
           for path in input_paths):
        parser.error("output must not overwrite input evidence")
    try:
        inputs = [path.read_bytes() for path in (args.native_json, args.radia_json)]
        payloads = [json.loads(data, object_pairs_hook=_unique_object) for data in inputs]
        result = compare(*payloads)
        result["input_sha256"] = dict(zip(("native", "radia"),
                                         (hashlib.sha256(data).hexdigest() for data in inputs)))
    except (OSError, ValueError) as error:
        result = compare({}, {})
        result["error"] = str(error)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return 0 if result["sample_checks_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
