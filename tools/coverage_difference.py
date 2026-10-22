"""Manual evidence comparison before reducing validation to a regression lock.

Consumes explicit exports, never runs commands or imports coverage/solver modules.
Matching coverage is not a physics oracle or permission to retire a validation.
Native entry-point hits do not establish internal native/JIT branch coverage.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path, PurePosixPath
import re


def _require(ok, label, actual, expected):
    if not ok:
        raise ValueError(f"{label}: actual={actual!r}; expected={expected}")


def _hash(value, size, label):
    _require(isinstance(value, str) and re.fullmatch(f"[0-9a-f]{{{size}}}", value),
             label, value, f"{size} lowercase hex characters")


def _path(value):
    _require(isinstance(value, str) and value and "\\" not in value and ":" not in value
             and not PurePosixPath(value).is_absolute()
             and all(part not in {"", ".", ".."} for part in value.split("/")),
             "production path", value, "normalized repository-relative path")


def _integers(values, label, *, arcs=False):
    _require(isinstance(values, list), label, values, "array")
    entries = []
    for item in values:
        if arcs:
            _require(isinstance(item, list) and len(item) == 2
                     and all(type(n) is int and n != 0 for n in item),
                     label, item, "two nonzero line numbers (negative entry/exit allowed)")
            entries.append(tuple(item))
        else:
            _require(type(item) is int and item > 0, label, item, "positive integer line")
            entries.append(item)
    _require(len(entries) == len(set(entries)), label, values, "unique entries")
    return set(entries)


def load_export(path):
    """Read an explicit export; no git ancestry, source execution or path rewriting."""
    def pairs(items):
        result = {}
        for key, value in items:
            _require(key not in result, "JSON key", key, "unique key")
            result[key] = value
        return result
    def invalid(value):
        raise ValueError(f"nonfinite JSON token: {value}")
    return json.loads(Path(path).read_text(encoding="utf-8"),
                      object_pairs_hook=pairs, parse_constant=invalid)


def _validate(data, role):
    keys = {"schema", "role", "identity", "run", "files", "native"}
    _require(isinstance(data, dict) and set(data) == keys, role, data, sorted(keys))
    _require(data["schema"] == "radia.coverage-evidence.v1" and data["role"] == role,
             "schema/role", (data["schema"], data["role"]), f"radia.coverage-evidence.v1/{role}")
    identity = data["identity"]
    _require(isinstance(identity, dict) and set(identity) == {
        "source_sha", "artifacts", "backend", "route", "instrumentation", "production_scope"},
        "identity", identity, "source/artifacts/backend/route/instrumentation/production_scope")
    _hash(identity["source_sha"], 40, "source_sha")
    _require(identity["backend"] in {"cpu", "gpu"}, "backend", identity["backend"], "cpu or gpu")
    _require(isinstance(identity["artifacts"], dict) and identity["artifacts"],
             "artifacts", identity["artifacts"], "nonempty artifact-name to SHA256 map")
    for name, digest in identity["artifacts"].items():
        _path(name)
        _hash(digest, 64, "artifact hash")
    route = identity["route"]
    route_keys = {"entry_point", "formulation", "element", "order", "nonlinear", "boundary", "constitutive"}
    _require(isinstance(route, dict) and set(route) == route_keys, "route", route, sorted(route_keys))
    _require(type(route["order"]) is int and route["order"] > 0 and type(route["nonlinear"]) is bool
             and all(isinstance(route[k], str) and route[k].strip()
                     for k in route_keys - {"order", "nonlinear"}), "route values", route, "typed declared route")
    scope = identity["production_scope"]
    _require(isinstance(scope, list) and scope and len(scope) == len(set(scope)),
             "production_scope", scope, "nonempty unique path array")
    for path in scope:
        _path(path)
    ins = identity["instrumentation"]
    _require(isinstance(ins, dict) and set(ins) == {"tool", "version", "branch", "exclusions", "subprocess"},
             "instrumentation", ins, "tool/version/branch/exclusions/subprocess")
    _require(all(isinstance(ins[k], str) and ins[k] for k in ("tool", "version"))
             and type(ins["branch"]) is bool and isinstance(ins["exclusions"], list)
             and all(isinstance(x, str) for x in ins["exclusions"]), "instrumentation values", ins, "typed metadata")
    child = ins["subprocess"]
    _require(isinstance(child, dict) and set(child) == {"required", "sentinel_verified"}
             and all(type(v) is bool for v in child.values()), "subprocess", child, "required/sentinel_verified booleans")
    run = data["run"]
    _require(isinstance(run, dict) and set(run) == {"command", "conditions", "runtime_s", "gate"},
             "run", run, "command/conditions/runtime_s/gate")
    _require(isinstance(run["command"], list) and run["command"]
             and all(isinstance(x, str) and x for x in run["command"]), "command", run["command"], "argv array")
    _require(isinstance(run["conditions"], dict) and run["conditions"], "conditions", run["conditions"], "nonempty JSON object")
    _require(type(run["runtime_s"]) in (int, float) and math.isfinite(run["runtime_s"])
             and run["runtime_s"] >= 0, "runtime_s", run["runtime_s"], "nonnegative finite metadata")
    gate = run["gate"]
    _require(isinstance(gate, dict) and set(gate) == {"status", "evidence_sha256"}
             and gate["status"] in {"passed", "failed", "partial", "timed-out"}, "gate", gate, "status/evidence_sha256")
    _hash(gate["evidence_sha256"], 64, "gate evidence")
    _require(isinstance(data["files"], dict) and set(data["files"]) <= set(scope),
             "files", data["files"], "only declared production scope")
    for path, entry in data["files"].items():
        _require(isinstance(entry, dict) and set(entry) == {"sha256", "lines", "arcs"}, "file", entry, "sha256/lines/arcs")
        _hash(entry["sha256"], 64, path)
        _integers(entry["lines"], path)
        _integers(entry["arcs"], path, arcs=True)
    _require(isinstance(data["native"], dict), "native", data["native"], "named entry-point map")
    for name, entry in data["native"].items():
        _require(isinstance(name, str) and name and isinstance(entry, dict)
                 and set(entry) == {"calls", "artifact", "configuration"}, "native entry", entry, "calls/artifact/configuration")
        _require(type(entry["calls"]) is int and entry["calls"] >= 0
                 and entry["artifact"] in identity["artifacts"] and isinstance(entry["configuration"], dict),
                 "native values", entry, "nonnegative count, declared artifact, configuration object")
    # Also reject non-JSON inputs/NaN nested in caller-provided dictionaries.
    json.dumps(data, allow_nan=False)


def compare_exports(heavy, light):
    """Return explicit differences/unknowns, never a downgrade approval or CI gate."""
    _validate(heavy, "heavy")
    _validate(light, "light")
    hi, li = heavy["identity"], light["identity"]
    for key in ("source_sha", "artifacts", "backend", "route", "production_scope"):
        _require(hi[key] == li[key], f"pair identity.{key}", li[key], hi[key])
    for key in ("tool", "version", "exclusions"):
        _require(hi["instrumentation"][key] == li["instrumentation"][key],
                 f"pair instrumentation.{key}", li["instrumentation"][key], hi["instrumentation"][key])
    unknown = []
    for role, data in (("heavy", heavy), ("light", light)):
        ins = data["identity"]["instrumentation"]
        if not ins["branch"]:
            unknown.append(f"{role}: branch tracing absent")
        if ins["subprocess"]["required"] and not ins["subprocess"]["sentinel_verified"]:
            unknown.append(f"{role}: subprocess sentinel unverified")
        if data["run"]["gate"]["status"] != "passed":
            unknown.append(f"{role}: original numerical gate {data['run']['gate']['status']}")
        if not any(entry["lines"] for entry in data["files"].values()):
            unknown.append(f"{role}: no executed production lines")
        if data["identity"]["backend"] == "gpu" and not data["native"]:
            unknown.append(f"{role}: GPU entry-point evidence absent")
    differences = {}
    for path in hi["production_scope"]:
        if path not in heavy["files"] or path not in light["files"]:
            unknown.append(f"{path}: missing file capture")
            continue
        h, l = heavy["files"][path], light["files"][path]
        _require(h["sha256"] == l["sha256"], path, l["sha256"], h["sha256"])
        differences[path] = {}
        for kind in ("lines", "arcs"):
            hs, ls = _integers(h[kind], path, arcs=kind == "arcs"), _integers(l[kind], path, arcs=kind == "arcs")
            differences[path][kind] = {"both": sorted(hs & ls), "heavy_only": sorted(hs - ls), "light_only": sorted(ls - hs)}
    routes = []
    for name in sorted(heavy["native"].keys() | light["native"].keys()):
        h, l = heavy["native"].get(name), light["native"].get(name)
        if h is None or l is None:
            state = "heavy-only" if h else "light-only"
        elif h["artifact"] != l["artifact"] or h["configuration"] != l["configuration"]:
            state = "unknown"
            unknown.append(f"{name}: native entry-point configuration mismatch")
        else:
            state = "covered-both" if h["calls"] and l["calls"] else "unknown"
            if state == "unknown":
                unknown.append(f"{name}: no invocation observed in both runs")
        routes.append({"entry_point": name, "state": state, "internal_coverage": "UNKNOWN"})
    return {"schema": "radia.coverage-difference.v1", "identity": hi,
            "runs": {"heavy": heavy["run"], "light": light["run"]},
            "files": differences, "route_matrix": routes, "unknowns": unknown,
            "recommendation": "insufficient-evidence" if unknown else "owner-review-required",
            "automatic_downgrade": False,
            "limitations": ["Input provenance is declared, not independently attested.",
                            "No enclosing-function discovery or native internal branch measurement.",
                            "Coverage does not replace an independent numerical oracle."]}


def main(argv=None):
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--heavy", type=Path, required=True)
    parser.add_argument("--light", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    report = compare_exports(load_export(args.heavy), load_export(args.light))
    # Exclusive creation preserves previous reports and original evidence.
    with args.out.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
