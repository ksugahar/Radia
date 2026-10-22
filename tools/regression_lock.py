"""Compare lightweight outputs with reviewed historical baselines, not physics oracles.

Baselines record source SHA, conditions, per-quantity tolerances and runtime_s.
Runtime is descriptive only; no runtime comparison or 30-second assertion occurs.
Regeneration writes a separate proposal and reason/diff, never an approved file.
Create the first baseline by hand with user-owned tolerances. Complex quantities
must be split into real/imaginary quantities. Conditions compare JSON types too.
Arrays passed from NumPy must use .tolist(); no NumPy dependency is required.
Proposal hashes normalize CRLF/CR to LF; original approved bytes are never edited.
"""
from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import math
from numbers import Real
import platform
import subprocess
from pathlib import Path
import re


SOURCE_ROOT = Path(__file__).resolve().parents[1]


def _json(text):
    def pairs(items):
        obj = {}
        for key, value in items:
            if key in obj:
                raise ValueError(f"duplicate JSON key: {key}")
            obj[key] = value
        return obj
    def invalid(value):
        raise ValueError(f"non-finite JSON value: {value}")
    return json.loads(text, object_pairs_hook=pairs, parse_constant=invalid)


def _number(value, label, *, nonnegative=False):
    try:
        valid = isinstance(value, Real) and not isinstance(value, bool) and math.isfinite(value)
    except OverflowError:
        valid = False
    if not valid or nonnegative and value < 0:
        raise ValueError(f"{label}: actual={value!r}; expected finite {'nonnegative ' if nonnegative else ''}real number (NumPy arrays: use .tolist(); complex: split real/imaginary)")


def _values(value, label):
    if isinstance(value, list):
        if not value:
            raise ValueError(f"{label}: empty quantity array")
        for index, item in enumerate(value):
            _values(item, f"{label}[{index}]")
    else:
        _number(value, label)


def _conditions(value):
    if not isinstance(value, dict) or not value:
        raise ValueError("conditions: expected a nonempty JSON object")
    # Reject NaN, Infinity, unsupported objects and non-string object keys.
    def check(item):
        if isinstance(item, dict):
            if any(not isinstance(key, str) for key in item):
                raise ValueError("conditions: keys must be strings")
            for child in item.values():
                check(child)
        elif isinstance(item, list):
            for child in item:
                check(child)
        elif item is not None and type(item) not in (str, bool):
            _number(item, "conditions")
    check(value)
    return json.dumps(value, sort_keys=True, allow_nan=False)


def _validate(baseline):
    keys = {"schema", "source_sha", "conditions", "runtime_s", "quantities"}
    if not isinstance(baseline, dict) or set(baseline) != keys:
        raise ValueError(f"baseline: expected keys={sorted(keys)}; actual={list(baseline) if isinstance(baseline, dict) else baseline!r}")
    if baseline["schema"] != "radia.regression-lock.v1":
        raise ValueError("baseline.schema: expected radia.regression-lock.v1")
    if not isinstance(baseline["source_sha"], str) or not re.fullmatch(r"[0-9a-f]{40}", baseline["source_sha"]):
        raise ValueError("baseline.source_sha: expected full 40-character source SHA")
    _conditions(baseline["conditions"])
    _number(baseline["runtime_s"], "baseline.runtime_s", nonnegative=True)
    quantities = baseline["quantities"]
    if not isinstance(quantities, dict) or not quantities:
        raise ValueError("baseline.quantities: expected nonempty object")
    for name, spec in quantities.items():
        if not isinstance(name, str) or not name or name == "runtime_s":
            raise ValueError("quantity names must be nonempty strings; runtime_s is metadata only")
        if not isinstance(spec, dict) or set(spec) != {"expected", "atol", "rtol", "unit"}:
            raise ValueError(f"quantities.{name}: expected expected/atol/rtol/unit")
        _values(spec["expected"], f"quantities.{name}.expected")
        _number(spec["atol"], f"quantities.{name}.atol", nonnegative=True)
        _number(spec["rtol"], f"quantities.{name}.rtol", nonnegative=True)
        if not isinstance(spec["unit"], str) or not spec["unit"].strip():
            raise ValueError(f"quantities.{name}.unit: expected declared unit")
    return baseline


def load_baseline(path):
    """Read and validate stored metadata without generating any reference values."""
    return _validate(_json(Path(path).read_text(encoding="utf-8")))


def _compare(actual, expected, atol, rtol, label):
    if isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            raise AssertionError(f"{label}: actual=type {type(actual).__name__}, length {len(actual) if hasattr(actual, '__len__') else 'scalar'}; expected=list length {len(expected)}; array shape differs; NumPy arrays: use .tolist()")
        for index, (got, wanted) in enumerate(zip(actual, expected)):
            _compare(got, wanted, atol, rtol, f"{label}[{index}]")
    else:
        _number(actual, label)
        bound = atol + rtol * abs(expected)
        if not math.isfinite(bound):
            raise ValueError(f"{label}: tolerance budget overflow")
        if abs(actual - expected) > bound:
            raise AssertionError(f"{label}: actual={actual!r}; expected={expected!r}; atol={atol!r}; rtol={rtol!r}; allowed={bound!r}")


def assert_regression(actual, baseline, *, conditions):
    """Compare exact quantity sets; bound = atol + rtol * abs(expected).

    This is a prior-output regression lock, not an independent validation.
    Call the same production code path used by the associated validation.
    Evidence of that path coverage is owned by the caller/lifecycle review.
    """
    _validate(baseline)
    if _conditions(conditions) != _conditions(baseline["conditions"]):
        raise AssertionError(f"conditions: actual={conditions!r}; expected={baseline['conditions']!r}")
    if not isinstance(actual, dict) or set(actual) != set(baseline["quantities"]):
        raise AssertionError(f"quantity names: actual={list(actual) if isinstance(actual, dict) else actual!r}; expected={list(baseline['quantities'])!r}")
    for name, spec in baseline["quantities"].items():
        _compare(actual[name], spec["expected"], spec["atol"], spec["rtol"], name)


def _source_identity():
    """Capture the repository containing this generator, never a caller assertion."""
    def git(*args):
        return subprocess.check_output(["git", "-c", f"safe.directory={SOURCE_ROOT}", *args],
            cwd=SOURCE_ROOT, text=True, encoding="utf-8").strip()
    head = git("rev-parse", "HEAD")
    dirty = git("status", "--porcelain", "--untracked-files=all")
    return {"head": head, "dirty": bool(dirty), "status_porcelain": dirty,
            "python_version": platform.python_version(), "python_implementation": platform.python_implementation()}


def _lf(data):
    return data.replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def _proposal_identity(source_sha, reason):
    if not isinstance(reason, str) or not reason.strip():
        raise ValueError("regeneration reason is required")
    identity = _source_identity()
    if identity["dirty"]:
        raise ValueError(f"proposal requires clean generator repository; actual dirty={identity['status_porcelain']!r}")
    if source_sha != identity["head"]:
        raise ValueError(f"proposal source_sha: actual={source_sha!r}; expected git HEAD={identity['head']!r}")
    return identity


def propose_regeneration(baseline_path, *, actual, source_sha, conditions, runtime_s, reason, out):
    """Write a new proposal directory. Never replace approved baseline bytes.

    Existing tolerances and units are preserved. No accept/apply command exists.
    The caller must rerun the retained validation before requesting re-baseline.
    """
    identity = _proposal_identity(source_sha, reason)
    original = Path(baseline_path).read_bytes()
    baseline = _validate(_json(original))
    proposed = _json(json.dumps(baseline))
    if not isinstance(actual, dict) or set(actual) != set(proposed["quantities"]):
        raise ValueError("proposal actual quantity names must match the approved baseline")
    for name in proposed["quantities"]:
        proposed["quantities"][name]["expected"] = actual[name]
    proposed.update(source_sha=source_sha, conditions=conditions, runtime_s=runtime_s)
    _validate(proposed)
    return _write_proposal(original, proposed, source_sha=source_sha, conditions=conditions,
        runtime_s=runtime_s, reason=reason, out=out, identity=identity,
        changed=json.dumps(proposed, sort_keys=True) != json.dumps(baseline, sort_keys=True),
        quantity_changes=[name for name in actual if json.dumps(actual[name], sort_keys=True)
                          != json.dumps(baseline["quantities"][name]["expected"], sort_keys=True)])


def _write_proposal(original, proposed, *, source_sha, conditions, runtime_s, reason, out,
                    identity, changed, quantity_changes):
    """The single proposal writer: provenance, hashes, exclusive output and diff."""
    encoded = ((json.dumps(proposed, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")
               if changed else original)
    directory = Path(out)
    directory.mkdir(parents=True, exist_ok=False)
    (directory / "baseline.json").write_bytes(encoded)
    report = {"reason": reason, "approval": "proposal-only; review and user approval required",
              "approved_sha256": hashlib.sha256(_lf(original)).hexdigest(),
              "proposed_sha256": hashlib.sha256(_lf(encoded)).hexdigest(),
              "source_sha": source_sha, "generation": identity,
              "hash_convention": "UTF-8 bytes with CRLF/CR normalized to LF",
              "runtime_s": runtime_s, "conditions": conditions, "changed": changed,
              "byte_changed": encoded != original, "quantity_changes": quantity_changes,
              "diff": list(difflib.unified_diff(original.decode("utf-8").splitlines(),
                  encoded.decode("utf-8").splitlines(), fromfile="approved/baseline.json",
                  tofile="proposed/baseline.json", lineterm=""))}
    (directory / "proposal.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return report


def flatten_json_record(record):
    """Separate dotted numerical leaves from exact JSON structure/literals.

    '~' and '.' in object keys are escaped as '~0' and '~1'. Container metadata
    distinguishes object keys from array indices, empty containers and literals.
    No field values, tuple-cleaning policy or units are transformed here.
    """
    values, structure = {}, {}
    def visit(item, path):
        if isinstance(item, dict):
            if any(not isinstance(key, str) for key in item):
                raise ValueError("record object keys must be strings")
            structure[path] = {"object": sorted(item)}
            for key, child in item.items():
                visit(child, path + "." + key.replace("~", "~0").replace(".", "~1"))
        elif isinstance(item, list):
            structure[path] = {"array_length": len(item)}
            for index, child in enumerate(item):
                visit(child, path + "." + str(index))
        elif type(item) in (int, float):
            _number(item, path)
            values[path] = item
        elif item is None or type(item) in (str, bool):
            structure[path] = {"literal": item}
        else:
            raise ValueError(f"{path}: expected JSON value, got {type(item).__name__}; clean tuples/arrays before flattening")
    visit(record, "root")
    return values, structure


def propose_json_regeneration(baseline_path, *, actual, source_sha, conditions, runtime_s, reason, out):
    """Propose an existing nested JSON record through the common writer.

    This is an adapter for unchanged legacy fixture schemas, not an independent
    acceptance route. Semantic equality preserves approved original bytes,
    including signed-zero spelling. No approved file is modified.
    """
    identity = _proposal_identity(source_sha, reason)
    _conditions(conditions)
    _number(runtime_s, "runtime_s", nonnegative=True)
    original = Path(baseline_path).read_bytes()
    before = _json(original)
    before_values, before_structure = flatten_json_record(before)
    after_values, after_structure = flatten_json_record(actual)
    changes = [key for key in sorted(before_values.keys() | after_values.keys())
               if key not in before_values or key not in after_values
               or before_values[key] != after_values[key]]
    return _write_proposal(original, actual, source_sha=source_sha, conditions=conditions,
        runtime_s=runtime_s, reason=reason, out=out, identity=identity,
        changed=actual != before or before_structure != after_structure, quantity_changes=changes)


def main(argv=None):
    parser = argparse.ArgumentParser(__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    proposal = sub.add_parser("propose", help="write a proposal; never apply a baseline")
    proposal.add_argument("--baseline", type=Path, required=True)
    proposal.add_argument("--actual", type=Path, required=True, help="JSON quantity-name to output-value object")
    proposal.add_argument("--conditions", type=Path, required=True)
    proposal.add_argument("--source-sha", required=True, help="Must equal the clean generator repository HEAD")
    proposal.add_argument("--runtime-s", type=float, required=True)
    proposal.add_argument("--reason", required=True)
    proposal.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    propose_regeneration(args.baseline, actual=_json(args.actual.read_text(encoding="utf-8")),
        source_sha=args.source_sha, conditions=_json(args.conditions.read_text(encoding="utf-8")),
        runtime_s=args.runtime_s, reason=args.reason, out=args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
