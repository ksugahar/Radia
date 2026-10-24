"""Read-only ownership/trigger contracts; not a CI runtime impact selector.

Ordered rules replace parent ownership, never union it. Workflow projections
describe current behavior; migration exceptions describe debt, not ownership.
Existing sparsesolv and MCP selectors remain the execution authority in PR1.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import json
from pathlib import Path
import re
import subprocess

SCHEMA = "radia.ci-distributions.v1"
DISTRIBUTIONS = ("radia", "radia-simulink", "radia-mcp", "radia-optuna",
                 "cubit-mesh-export")
ROOT_CONFIGS = ("CMakeLists.txt", "Build.ps1", "pyproject.toml", "MANIFEST.in",
                "tests/test_tier_manifest.json")
KINDS = {"shared_abi", "shared_build", "integration"}


@lru_cache(maxsize=None)
def _glob(pattern):
    """Anchored Actions path glob subset. Reject unsupported syntax explicitly."""
    if not pattern or pattern.startswith(("/", "!")) or "\\" in pattern or ".." in pattern.split("/"):
        raise ValueError(f"invalid relative pattern: {pattern}")
    result = ""; i = 0
    while i < len(pattern):
        char = pattern[i]
        if char == "*":
            if pattern[i:i+2] == "**":
                if pattern[i:i+3] == "**/":
                    result += "(?:.*/)?"; i += 3
                else:
                    result += ".*"; i += 2
            else:
                result += "[^/]*"; i += 1
        elif char == "?":
            result += "[^/]"; i += 1
        elif char in "[]{}":
            raise ValueError(f"unsupported glob syntax: {pattern}")
        else:
            result += re.escape(char); i += 1
    return re.compile("^" + result + "$")


def matches(path, pattern):
    return bool(_glob(pattern).fullmatch(path))


def owner_rule(manifest, path):
    selected = None
    for rule in manifest["rules"]:
        if matches(path, rule["pattern"]):
            selected = rule
    if selected is None:
        raise ValueError(f"unowned path: {path}")
    return selected


def selected(path, projection):
    if "paths" in projection:
        included = False
        for pattern in projection["paths"]:
            if matches(path, pattern.lstrip("!")):
                included = not pattern.startswith("!")
        return included
    return not any(matches(path, p) for p in projection.get("paths-ignore", []))


def tracked(root):
    return subprocess.check_output(["git", "-C", str(root), "ls-files", "-z"],
                                   encoding="utf-8").rstrip("\0").split("\0")


def load(root):
    return json.loads((root / "tools/ci_distribution_manifest.json").read_text(encoding="utf-8"))


def validate_schema(manifest):
    if manifest.get("schema") != SCHEMA or manifest.get("distributions") != list(DISTRIBUTIONS):
        raise ValueError("unknown schema or distributions")
    if manifest.get("infrastructure_owner") != "repo":
        raise ValueError("infrastructure owner must be repo")
    if manifest.get("rule_semantics") != "ordered_last_match_replaces":
        raise ValueError("unknown rule precedence")
    patterns = set()
    for rule in manifest["rules"]:
        pattern = rule["pattern"]; _glob(pattern)
        if pattern in ("*", "**", "**/*"):
            raise ValueError("root catch-all hides unclassified files")
        if pattern in patterns:
            raise ValueError(f"duplicate rule: {pattern}")
        patterns.add(pattern)
        owners = rule["owners"]
        if not owners or len(set(owners)) != len(owners) or not set(owners) <= set(DISTRIBUTIONS) | {"repo"}:
            raise ValueError(f"invalid owners: {pattern}")
        if len(owners) > 1:
            reason = rule.get("reason", {})
            if reason.get("kind") not in KINDS or not reason.get("detail", "").strip():
                raise ValueError(f"shared reason required: {pattern}")
            if rule.get("selection") not in ("all_owners", "content_diff", "content_diff:tools/sparsesolv_ci_impact.py"):
                raise ValueError(f"shared selection required: {pattern}")
        elif "reason" in rule or "selection" in rule:
            raise ValueError(f"unshared path carries shared metadata: {pattern}")
    for path in ROOT_CONFIGS:
        rule = owner_rule(manifest, path)
        if set(rule["owners"]) != set(DISTRIBUTIONS) or rule.get("selection") != "content_diff:tools/sparsesolv_ci_impact.py":
            raise ValueError(f"root content-diff contract missing: {path}")
    routing = manifest["root_content_diff"]
    if routing.get("router") != "tools/sparsesolv_ci_impact.py" or routing.get("unknown") != "all_owners":
        raise ValueError("unknown root router/fallback")
    attribution = routing["attribution"]
    if set(attribution) != set(ROOT_CONFIGS):
        raise ValueError("incomplete root attribution")
    maps = [attribution["Build.ps1"], attribution["tests/test_tier_manifest.json"]]
    maps += [attribution["CMakeLists.txt"][key] for key in ("targets", "options")]
    maps += [attribution["pyproject.toml"][key] for key in ("extras", "packages")]
    lists = [attribution["pyproject.toml"]["build_system"]]
    lists += [owners for mapping in maps for owners in mapping.values()]
    for owners in lists:
        if not isinstance(owners, list) or not owners or not set(owners) <= set(DISTRIBUTIONS) | {"repo"}:
            raise ValueError("invalid root attribution owners")
    for name, lane in manifest["lanes"].items():
        if not lane["owners"] or not set(lane["owners"]) <= set(DISTRIBUTIONS) | {"repo"}:
            raise ValueError(f"invalid lane owners: {name}")
        if lane["workflow"] != f".github/workflows/{name}.yml":
            raise ValueError(f"invalid workflow: {name}")
        for event, projection in lane["events"].items():
            if event == "schedule":
                if not isinstance(projection, list) or not all(set(x) == {"cron"} for x in projection):
                    raise ValueError("invalid schedule")
                continue
            if "paths" in projection and "paths-ignore" in projection:
                raise ValueError("paths and paths-ignore cannot coexist")
            for key in ("paths", "paths-ignore"):
                patterns = projection.get(key, [])
                if not isinstance(patterns, list):
                    raise ValueError("ordered path array required")
                for pattern in patterns:
                    _glob(pattern.lstrip("!"))
    ids = set()
    for exception in manifest["migration_exceptions"]:
        if exception["id"] in ids:
            raise ValueError("duplicate exception id")
        ids.add(exception["id"])
        if exception["lane"] not in manifest["lanes"] or exception["owner"] not in (*DISTRIBUTIONS, "repo"):
            raise ValueError("unknown exception lane/owner")
        if exception["kind"] not in {"foreign_trigger", "markdown_native", "content_diff_pending", "missing_trigger"}:
            raise ValueError("unknown exception kind")
        for key in ("id", "pattern", "current", "target", "follow_up"):
            if not isinstance(exception[key], str) or not exception[key].strip():
                raise ValueError(f"missing exception {key}")
        _glob(exception["pattern"])


def workflow_events(document):
    """Preserve release/manual payloads too, separately from branch path events."""
    events = document["on"]
    return {name: value or {} for name, value in events.items()}


def numerical_environment_input(path):
    """Cross-distribution environment edits require full CPU numerical validation.

    This safety gate does not transfer distribution or release ownership.
    """
    name = Path(path).name.lower()
    return (name in {'pyproject.toml', 'build.ps1', 'cmakelists.txt', 'test_tier_manifest.json'}
            or name.startswith(('requirements', 'constraints'))
            or path.startswith(('.github/workflows/', '.github/actions/', 'cmake/'))
            or path.endswith(('.cmake', '.toolchain')))


def issues(manifest, files):
    """Target-contract violations as (lane, kind, path), with no runtime routing."""
    result = set()
    for name, lane in manifest["lanes"].items():
        branch_events = [v for k, v in lane["events"].items() if k in ("push", "pull_request")
                         and ("branches" in v or k == "pull_request")]
        if not branch_events:
            continue  # Tag-only/manual/release completion is a separate contract.
        for path in files:
            owners = owner_rule(manifest, path)["owners"]
            triggered = any(selected(path, event) for event in branch_events)
            if triggered and lane["scope"] != "repository_contract" and not set(owners).intersection(lane["owners"]):
                if not (name == "sparsesolv" and numerical_environment_input(path)):
                    result.add((name, "foreign_trigger", path))
            if triggered and lane["native"] and path.endswith(".md") and path not in lane["markdown_inputs"]:
                result.add((name, "markdown_native", path))
            if triggered and path in ROOT_CONFIGS and lane["native"] and path not in lane["content_diff_active"]:
                result.add((name, "content_diff_pending", path))
        for path in lane["required_paths"]:
            # Exact paths, including markdown or build closure omissions.
            if path not in files:
                raise ValueError(f"required lane input vanished: {name}: {path}")
            if not all(selected(path, event) for event in branch_events):
                result.add((name, "missing_trigger", path))
    return result


def exception_coverage(manifest, violations):
    remaining = set(violations)
    for exception in manifest["migration_exceptions"]:
        covered = {v for v in violations if v[:2] == (exception["lane"], exception["kind"])
                   and matches(v[2], exception["pattern"])}
        if not covered:
            raise ValueError(f"resolved/unused exception: {exception['id']}")
        remaining -= covered
    if remaining:
        raise ValueError(f"undeclared CI ownership debt: {sorted(remaining)[:12]}")


def exception_digest(exception):
    return hashlib.sha256(json.dumps(exception, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def check_shrink_only(current, ceiling):
    """Removal is permitted; no new/reworded/broadened exception is implicit approval."""
    for exception in current:
        if ceiling.get(exception["id"]) != exception_digest(exception):
            raise ValueError(f"new or expanded exception: {exception['id']}")


def include_closure(root, seeds, prefix):
    """Resolve local quote includes (relative or include-root), not external SDKs."""
    pending = list(seeds); seen = set()
    while pending:
        path = pending.pop()
        if path in seen:
            continue
        source = root / path
        if not source.is_file():
            raise ValueError(f"include seed missing: {path}")
        seen.add(path)
        for include in re.findall(r'^\s*#\s*include\s*"([^"]+)"', source.read_text(encoding="utf-8"), re.M):
            choices = [source.parent / include, root / prefix / include]
            candidates = [p for p in choices if p.is_file()]
            if candidates:
                target = candidates[0].resolve().relative_to(root.resolve()).as_posix()
                if target.startswith(prefix + "/"):
                    pending.append(target)
    return sorted(seen)

