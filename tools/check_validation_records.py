"""Read-only historical validation-record contracts; never rerun numerical work."""
from __future__ import annotations

import argparse
from datetime import date
import hashlib
import ipaddress
import json
import math
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys

from ci_distribution_contract import matches
from validation_adoption import validate_adoption, validate_champion

SCHEMA = Path(__file__).with_name("validation_record.schema.json")
PRIVATE = re.compile(r"(?:\b[A-Za-z]:[\\/](?!/)|\\\\[A-Za-z0-9_.-]+\\|/(?:Users|home|temp|tmp)/|\bcommercial[ _-]+(?:comparison|reference)\b)", re.I)


def read_json(data):
    def invalid(value):
        raise ValueError(f"non-finite JSON value: {value}")
    def finite_float(value):
        result = float(value)
        if not math.isfinite(result):
            invalid(value)
        return result
    def unique_keys(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate JSON key: {key}")
            result[key] = value
        return result
    return json.loads(data, parse_constant=invalid, parse_float=finite_float, object_pairs_hook=unique_keys)


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], text=True, encoding="utf-8").strip()


def local_path(base, value):
    path = PurePosixPath(value)
    if not value or "\\" in value or ":" in value or path.is_absolute() or ".." in path.parts:
        raise ValueError(f"unsafe relative record path: {value!r}")
    resolved = (base / value).resolve()
    if not resolved.is_relative_to(base.resolve()):
        raise ValueError(f"record path escapes its base: {value!r}")
    return resolved


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def covered_hashes(root, patterns, *, head="HEAD"):
    """Hash exact bytes; normalize CRLF/CR to LF only for UTF-8 text."""
    tracked = git(root, "ls-tree", "-r", "--name-only", head).splitlines()
    paths = set()
    for pattern in patterns:
        local_path(root, pattern)
        paths.update(p for p in tracked if matches(p, pattern))
    result = {}
    for name in sorted(paths):
        path = local_path(root, name)
        if not path.is_file():
            continue
        data = path.read_bytes()
        try:
            data = data.decode("utf-8").replace("\r\n", "\n").replace("\r", "\n").encode("utf-8")
        except UnicodeDecodeError:
            pass
        result[name] = sha256(data)
    return result


def private_text(value):
    # Inspect decoded values, not JSON escaping. URLs and dotted versions are not paths.
    if isinstance(value, dict):
        return any(private_text(k) or private_text(v) for k, v in value.items())
    if isinstance(value, list):
        return any(private_text(v) for v in value)
    if not isinstance(value, str):
        return False
    if re.fullmatch(r"[0-9]+(?:\.[0-9]+){3}", value):
        try:
            ipaddress.ip_address(value)
            return True
        except ValueError:
            pass
    return bool(PRIVATE.search(value))


def check_schema_definition(schema):
    """Fail closed on keywords outside our fixed record-contract subset.

    This is not a general JSON Schema implementation. The repository owns this
    single schema and the explicitly supported subset below; no external refs.
    """
    allowed = {"$schema", "$id", "title", "type", "properties", "required",
               "additionalProperties", "items", "minItems", "uniqueItems",
               "minLength", "minimum", "maximum", "pattern", "const", "enum",
               "anyOf", "oneOf", "allOf", "not", "if", "then", "else"}
    if not isinstance(schema, dict) or set(schema) - allowed:
        raise ValueError("record schema definition uses unsupported keywords")
    children = list(schema.get("properties", {}).values())
    for name in ("additionalProperties", "items", "not", "if", "then", "else"):
        if isinstance(schema.get(name), dict):
            children.append(schema[name])
    for name in ("anyOf", "oneOf", "allOf"):
        children.extend(schema.get(name, []))
    for child in children:
        check_schema_definition(child)


def check_shape(value, schema, path="$"):
    """Validate only the repository's fixed contract using the standard library."""
    def fail(message):
        raise ValueError(f"record schema at {path}: {message}; actual={value!r}; expected={ {k: v for k, v in schema.items() if k not in ('properties', 'items')}!r}")
    def accepts(branch):
        try:
            check_shape(value, branch, path)
            return True
        except ValueError:
            return False
    def equal(a, b):
        if isinstance(a, bool) or isinstance(b, bool):
            return type(a) is type(b) and a == b
        return a == b
    if "const" in schema and not equal(value, schema["const"]):
        fail("constant mismatch")
    if "enum" in schema and not any(equal(value, choice) for choice in schema["enum"]):
        fail(f"actual={value!r}; expected one of {schema['enum']!r}")
    for keyword in ("anyOf", "oneOf", "allOf"):
        if keyword in schema:
            count = sum(accepts(branch) for branch in schema[keyword])
            if not (count > 0 if keyword == "anyOf" else count == 1 if keyword == "oneOf" else count == len(schema[keyword])):
                fail(keyword + " mismatch")
    if "not" in schema and accepts(schema["not"]):
        fail("forbidden alternative")
    if "if" in schema:
        branch = "then" if accepts(schema["if"]) else "else"
        if branch in schema:
            check_shape(value, schema[branch], path)
    kinds = schema.get("type", [])
    kinds = [kinds] if isinstance(kinds, str) else kinds
    validators = {"object": lambda: isinstance(value, dict), "array": lambda: isinstance(value, list),
                  "string": lambda: isinstance(value, str), "boolean": lambda: type(value) is bool,
                  "integer": lambda: type(value) is int, "null": lambda: value is None,
                  "number": lambda: type(value) is int or type(value) is float and math.isfinite(value)}
    if any(kind not in validators for kind in kinds):
        fail("unsupported type definition")
    if kinds and not any(validators[kind]() for kind in kinds):
        fail(f"actual={value!r} ({type(value).__name__}); expected type {kinds!r}")
    if isinstance(value, dict):
        if set(schema.get("required", [])) - set(value):
            fail(f"actual keys={sorted(value)!r}; expected required fields={schema['required']!r}; missing={sorted(set(schema['required']) - set(value))!r}")
        properties = schema.get("properties", {})
        additional = schema.get("additionalProperties", True)
        for key, item in value.items():
            if key in properties:
                check_shape(item, properties[key], path + "." + key)
            elif additional is False:
                fail("unknown field: " + key)
            elif isinstance(additional, dict):
                check_shape(item, additional, path + "." + key)
    if isinstance(value, list):
        if len(value) < schema.get("minItems", 0):
            fail("array too short")
        if schema.get("uniqueItems") and len({json.dumps(x, sort_keys=True) for x in value}) != len(value):
            fail("duplicate array entries")
        for i, item in enumerate(value):
            if "items" in schema:
                check_shape(item, schema["items"], f"{path}[{i}]")
    if isinstance(value, str):
        if len(value) < schema.get("minLength", 0):
            fail("string too short")
        if "pattern" in schema and re.fullmatch(schema["pattern"], value) is None:
            fail("pattern mismatch")
    if type(value) in (int, float):
        if type(value) is float and not math.isfinite(value):
            fail("non-finite number")
        if "minimum" in schema and value < schema["minimum"] or "maximum" in schema and value > schema["maximum"]:
            fail("number outside bounds")


def _validate_record(root, record_path, *, head="HEAD"):
    """Validate source identity, hashes and staleness without importing its entry."""
    root = root.resolve()
    record_path = record_path.resolve()
    if not record_path.is_relative_to(root):
        raise ValueError("record is outside repository")
    record = read_json(record_path.read_text(encoding="utf-8"))
    schema = read_json(SCHEMA.read_text(encoding="utf-8"))
    check_schema_definition(schema)
    check_shape(record, schema)
    date.fromisoformat(record["date"])  # Reject impossible calendar dates.
    if not record["claim"].strip() or not record["milestone"].strip():
        raise ValueError("milestone and claim must contain text")
    date.fromisoformat(record["id"].split("/")[-1][:8])
    candidates = [root / "validation_records", *sorted((root / "packages").glob("*/validation_records"))]
    base = next((p for p in candidates if record_path.is_relative_to(p)), None)
    if base is None or record["id"] != record_path.parent.relative_to(base).as_posix():
        raise ValueError("record id must match its records-root-relative directory")
    # Relations are documentary, not automatic performance ranking or promotion.
    def predecessors(item, stack):
        ids = [row["id"] for row in item["supersedes"]]
        if len(set(ids)) != len(ids):
            raise ValueError(f"duplicate supersedes ids: actual={ids!r}; expected=unique")
        for identity in ids:
            if identity in stack:
                raise ValueError(f"supersedes cycle: actual={stack + [identity]!r}; expected=acyclic")
            target = local_path(base, identity) / "record.json"
            other = read_json(target.read_text(encoding="utf-8"))
            check_shape(other, schema)
            if other["id"] != identity:
                raise ValueError(f"supersedes target id: actual={other['id']!r}; expected={identity!r}")
            if other["superseded_by"] != item["id"]:
                raise ValueError(f"supersedes backlink: actual={other['superseded_by']!r}; expected={item['id']!r}")
            predecessors(other, stack + [identity])
    predecessors(record, [record["id"]])
    if private_text(record):
        raise ValueError("record contains a private path or comparison marker")
    source = record["source"]["repo_sha"]
    provenance = record["source"]["provenance"]
    if provenance == "exact" and source is None:
        raise ValueError("source.provenance=exact requires a full source.repo_sha; actual=null")
    if provenance == "pre-rewrite":
        if source is not None:
            raise ValueError("source.provenance=pre-rewrite requires source.repo_sha=null; actual=non-null identity (omitted)")
        if record["status"] != "stale":
            raise ValueError(f"pre-rewrite evidence requires status=stale; actual={record['status']!r}")
        if not record["stale_reason"].startswith("pre-rewrite:") or not record["stale_reason"][len("pre-rewrite:"):].strip():
            raise ValueError("pre-rewrite stale_reason requires 'pre-rewrite:' followed by a nonempty explanation")
    current = covered_hashes(root, record["covered_paths"], head=head)
    expected = record["covered_sha256"]
    if not expected:
        raise ValueError("covered_sha256 must contain at least one source file")
    for pattern in record["covered_paths"]:
        if not any(matches(name, pattern) for name in current.keys() | expected.keys()):
            raise ValueError(f"covered_paths pattern matches no tracked or recorded file: {pattern}")
    for name in expected:
        local_path(root, name)
        if not any(matches(name, pattern) for pattern in record["covered_paths"]):
            raise ValueError(f"covered_sha256 path outside covered_paths: {name}")
    touched = {name for name in current.keys() | expected.keys() if current.get(name) != expected.get(name)}
    if record["status"] == "current" and touched:
        details = {name: {"actual": current.get(name), "expected": expected.get(name)} for name in sorted(touched)}
        raise ValueError(f"record claims current after covered content changed: {details}")
    covered = current
    source_only = record["source"]["source_only"]
    if not source_only and not (record["source"]["native_sha256"] or record["source"]["wheel_sha256"]):
        raise ValueError("source_only=false requires native_sha256 or wheel_sha256 artifact identity")
    if record["superseded_by"] is not None:
        target = local_path(base, record["superseded_by"]) / "record.json"
        if target == record_path or not target.is_file():
            raise ValueError(f"superseded_by target must exist and differ: {target}")
        successor = read_json(target.read_text(encoding="utf-8"))
        check_shape(successor, schema)
        if successor["id"] != record["superseded_by"] or not any(row["id"] == record["id"] for row in successor["supersedes"]):
            raise ValueError(f"superseded_by forward link: actual={successor['supersedes']!r}; expected={record['id']!r}")

    entry = record["entry"]
    script_path = entry["script"]
    data = local_path(record_path.parent, script_path).read_bytes()
    if sha256(data) != entry["script_sha256"]:
        raise ValueError(f"entry script sha256 mismatch: actual={sha256(data)} expected={entry['script_sha256']}")
    result = record["result"]
    result_path = local_path(record_path.parent, result["path"])
    payload = result_path.read_bytes()
    if sha256(payload) != result["sha256"]:
        raise ValueError(f"result sha256 mismatch: actual={sha256(payload)} expected={result['sha256']}")
    parsed_result = read_json(payload)
    if private_text(parsed_result):
        raise ValueError("result contains a private path, address or comparison marker")
    if private_text(data.decode("utf-8-sig")):
        raise ValueError("entry contains a private path, address or comparison marker")
    adoption = validate_adoption(root, record_path.parent, record, parsed_result,
                                 local_path=local_path, sha256=sha256)
    champion = validate_champion(root, base, record_path.parent, record, parsed_result,
                                 local_path=local_path, sha256=sha256, read_json=read_json,
                                 load_record=lambda repo, path: _validate_record(repo, path, head=head))
    return {"path": record_path.relative_to(root).as_posix(), "status": record["status"],
            "source": source, "provenance": provenance, "covered_files": len(covered), "changed_covered_files": sorted(touched),
            "numerical_work_executed": False, **adoption, **champion}


def validate_record(root, record_path, *, head="HEAD"):
    try:
        return _validate_record(root, record_path, head=head)
    except (ValueError, OSError, OverflowError, subprocess.CalledProcessError) as exc:
        raise ValueError(f"record {record_path}: {exc}") from exc


def validate_tree(root, *, head="HEAD"):
    """Root and independently owned package record directories share one schema."""
    bases = [root / "validation_records", *sorted((root / "packages").glob("*/validation_records"))]
    results = []
    for base in bases:
        if not base.exists():
            continue
        for payload in base.rglob("*.json"):
            if payload.name == "record.json":
                continue
            if not any((parent / "record.json").is_file() for parent in payload.parents if parent == base or parent.is_relative_to(base)):
                raise ValueError(f"orphan record JSON has no record.json wrapper: {payload.relative_to(root)}")
        for record in sorted(base.rglob("record.json")):
            results.append(validate_record(root, record, head=head))
    return {"passed": True, "records": results, "numerical_work_executed": False}


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--head", default="HEAD")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    try:
        report = validate_tree(args.root, head=args.head)
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        report = {"passed": False, "error": str(exc), "numerical_work_executed": False}
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
