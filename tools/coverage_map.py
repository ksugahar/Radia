"""Read immutable component/detector projections into a private incidence map.

Manual evidence reader, not a CI gate, adoption registry or independence oracle.
Consumes existing A4 exports; never runs the detector or modifies its source.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import coverage_difference as a4


def require(condition, message):
    if not condition:
        raise ValueError(message)


def keys(value, expected):
    require(isinstance(value, dict) and set(value) == set(expected),
            f"expected fields {sorted(expected)}; actual={value!r}")


def plain(path):
    require(not any(p.is_symlink() or p.is_junction() for p in [path, *path.parents]),
            f"links/junctions forbidden: {path}")
    return path


def decode(raw):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, f"duplicate JSON key: {key}")
            result[key] = value
        return result
    value = json.loads(raw.decode("utf-8"), object_pairs_hook=pairs)
    json.dumps(value, allow_nan=False)
    return value


def referenced(base, ref, inputs):
    keys(ref, {"path", "sha256"})
    a4._path(ref["path"])
    a4._hash(ref["sha256"], 64, "reference SHA256")
    path = plain(base / ref["path"])
    raw = path.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == ref["sha256"], f"hash mismatch: {path}")
    inputs.append(dict(ref))
    return raw


def text(value):
    require(isinstance(value, str) and value.strip(), "nonempty identifier/text required")
    return value


def components(data):
    keys(data, {"schema", "source_sha", "components"})
    require(data["schema"] == "radia.coverage-map.inventory.v1", "inventory schema")
    require(isinstance(data["components"], list) and data["components"], "empty inventory")
    seen = set()
    for item in data["components"]:
        kind = item.get("kind")
        common = {"id", "kind", "expected_detector_ids"}
        if kind in {"python_module", "python_line", "python_branch", "python_function"}:
            common |= {"path", "sha256"}
            a4._path(item["path"])
            a4._hash(item["sha256"], 64, "component file SHA256")
            if kind == "python_line":
                common.add("line")
                a4._integers([item["line"]], "component line")
            elif kind == "python_branch":
                common.add("arc")
                a4._integers([item["arc"]], "component arc", arcs=True)
            elif kind == "python_function":
                common.add("qualified_name")
                text(item["qualified_name"])
        elif kind in {"native_entry", "native_internal"}:
            common |= {"entry", "artifact", "configuration"}
            text(item["entry"])
            a4._path(item["artifact"])
            require(isinstance(item["configuration"], dict), "native configuration object")
        else:
            raise ValueError(f"unknown component kind: {kind!r}")
        keys(item, common)
        require(text(item["id"]) not in seen, "duplicate component ID")
        seen.add(item["id"])
        expected = item["expected_detector_ids"]
        require(isinstance(expected, list) and expected and
                all(isinstance(v, str) and v.strip() for v in expected) and
                len(set(expected)) == len(expected), "nonempty unique expected detectors")
    return data["components"]


def captures(base, detector, request, inputs):
    keys(detector, {"id", "nodeid", "test_source", "metadata", "captures"})
    text(detector["id"])
    text(detector["nodeid"])
    require(isinstance(detector["metadata"], dict), "detector metadata object")
    referenced(base, detector["test_source"], inputs)
    require(isinstance(detector["captures"], list), "captures array")
    result = []
    for pair in detector["captures"]:
        keys(pair, {"export", "attribution"})
        export = decode(referenced(base, pair["export"], inputs))
        a4._validate(export, export.get("role"))  # Reuse A4 v1 shape validation unchanged.
        require(export["role"] in {"heavy", "light"}, "A4 heavy/light role required")
        identity = export["identity"]
        require(identity["source_sha"] == request["source_sha"], "capture source mismatch")
        require(all(request["artifacts"].get(k) == v for k, v in
                    identity["artifacts"].items()), "capture artifact mismatch")
        binding = decode(referenced(base, pair["attribution"], inputs))
        keys(binding, {"schema", "source_sha", "export_sha256", "test_source_sha256",
                       "detector_ids", "nodeids", "reviewed", "capture_complete"})
        require(binding["schema"] == "radia.coverage-map.attribution.v1", "attribution schema")
        require(binding["source_sha"] == request["source_sha"] and
                binding["export_sha256"] == pair["export"]["sha256"] and
                binding["test_source_sha256"] == detector["test_source"]["sha256"],
                "attribution identity mismatch")
        ids, nodes = binding["detector_ids"], binding["nodeids"]
        require(isinstance(ids, list) and ids and len(ids) == len(set(ids)) and
                all(isinstance(v, str) and v for v in ids), "attribution detector IDs")
        require(isinstance(nodes, list) and len(nodes) == len(ids) and
                all(isinstance(v, str) and v for v in nodes), "attribution nodeids")
        require(detector["id"] in ids and nodes[ids.index(detector["id"])] == detector["nodeid"],
                "attribution does not bind this detector/nodeid")
        require(type(binding["reviewed"]) is bool and type(binding["capture_complete"]) is bool,
                "attribution flags must be booleans")
        result.append((export, binding))
    return result


def observation(component, export, binding, request):
    if len(binding["detector_ids"]) != 1 or not binding["reviewed"]:
        return False, False, "aggregate or unreviewed attribution"
    kind = component["kind"]
    if kind in {"python_function", "native_internal"}:
        return False, False, "function entry/native internal instrumentation unavailable"
    ins = export["identity"]["instrumentation"]
    complete = (binding["capture_complete"] and export["run"]["gate"]["status"] == "passed"
                and not ins["exclusions"] and not
                (ins["subprocess"]["required"] and not ins["subprocess"]["sentinel_verified"]))
    if kind == "native_entry":
        native = export["native"].get(component["entry"])
        if native is None:
            return False, False, "native entry capture missing"
        require(component["artifact"] in request["artifacts"] and
                native["artifact"] == component["artifact"] and
                native["configuration"] == component["configuration"], "native identity mismatch")
        return native["calls"] > 0, complete, "incomplete native capture" if not complete else ""
    file = export["files"].get(component["path"])
    if file is None:
        return False, False, "file capture missing/out of scope"
    require(file["sha256"] == component["sha256"], "component/capture file hash mismatch")
    if kind == "python_line":
        hit = component["line"] in file["lines"]
    elif kind == "python_branch":
        complete = complete and ins["branch"]
        hit = ins["branch"] and component["arc"] in file["arcs"]
    else:
        hit = bool(file["lines"])
    return hit, complete, "capture/gate/tracing incomplete" if not complete else ""


def build_map(path):
    path = plain(Path(path).absolute())
    raw = path.read_bytes()
    request = decode(raw)
    keys(request, {"schema", "source_sha", "artifacts", "inventory", "detectors"})
    require(request["schema"] == "radia.coverage-map.request.v1", "request schema")
    a4._hash(request["source_sha"], 40, "source SHA")
    require(isinstance(request["artifacts"], dict) and request["artifacts"], "artifact identities required")
    for name, digest in request["artifacts"].items():
        a4._path(name)
        a4._hash(digest, 64, "artifact SHA256")
    inputs = [{"path": path.name, "sha256": hashlib.sha256(raw).hexdigest()}]
    inventory = decode(referenced(path.parent, request["inventory"], inputs))
    require(inventory["source_sha"] == request["source_sha"], "inventory source mismatch")
    universe = components(inventory)
    for component in universe:
        if component["kind"].startswith("native_"):
            require(component["artifact"] in request["artifacts"], "inventory artifact undeclared")
    require(isinstance(request["detectors"], list), "detectors array")
    detectors, data, identities, capture_owners = {}, {}, set(), {}
    for detector in request["detectors"]:
        require(detector["id"] not in detectors, "duplicate detector declaration; use captures for repeats")
        detectors[detector["id"]] = detector
        data[detector["id"]] = captures(path.parent, detector, request, inputs)
        identity = (detector["test_source"]["sha256"], detector["nodeid"])
        require(identity not in identities, "same test source/nodeid under multiple detector IDs")
        identities.add(identity)
        for _, binding in data[detector["id"]]:
            if len(binding["detector_ids"]) == 1 and binding["reviewed"]:
                digest = binding["export_sha256"]
                require(capture_owners.get(digest, detector["id"]) == detector["id"],
                        "one capture cannot attest two single-detector identities")
                capture_owners[digest] = detector["id"]
    rows, incidence = [], {key: [] for key in detectors}
    for component in universe:
        hits, gaps = set(), {}
        for detector_id in component["expected_detector_ids"]:
            seen, complete, reasons = False, False, []
            for export, binding in data.get(detector_id, []):
                hit, done, reason = observation(component, export, binding, request)
                seen, complete = seen or hit, complete or done
                if reason:
                    reasons.append(reason)
            if seen:
                hits.add(detector_id)
                incidence[detector_id].append(component["id"])
            if not complete:
                gaps[detector_id] = sorted(set(reasons)) or ["detector/capture missing"]
        count = len(hits)
        category = "2+" if count >= 2 else str(count)
        rows.append({"component": component, "observed_test_count": count,
                     "observed_test_ids": sorted(hits), "observed_category": category,
                     "category": "UNKNOWN" if gaps else category, "complete": not gaps,
                     "gaps": gaps, "mutation_detection": "UNKNOWN", "change_activity": "UNKNOWN"})
    return {"schema": "radia.coverage-map.v1", "source_sha": request["source_sha"],
            "artifacts": request["artifacts"], "inputs": inputs, "components": rows,
            "detectors": list(detectors.values()), "detector_components": incidence,
            "qualified": False, "automatic_retirement": False,
            "limitations": ["Reviewed input declarations are not independently attested by this reader.",
                            "Distinct detector identities do not prove statistical independence.",
                            "Incidence is not chronological path order, accuracy or adoption.",
                            "Mutation, function entry, native internals and activity adapters are deferred."]}


def private_output(path):
    require(path.is_absolute(), "output must be an absolute private path")
    plain(path)
    require(path.parent.is_dir() and not path.exists(), "output parent must exist; output must be new")
    require(not any((p / ".git").exists() for p in path.parents), "output inside Git worktree forbidden")


def main(argv=None):
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    private_output(args.out)
    report = build_map(args.request)
    with args.out.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
