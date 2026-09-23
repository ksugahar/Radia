"""Read-only exact-artifact publication gate; never execute acceptance data."""
from __future__ import annotations

import argparse
import base64
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

SHA = re.compile(r"[0-9a-f]{40}\Z")
HASH = re.compile(r"[0-9a-f]{64}\Z")
STEPS = {"pip-check", "full6", "focused", "pip-freeze"}
CASES = {(kind, order) for kind in ("tet", "hex", "wedge") for order in (1, 2)}
TESTS = {
    "test_negative_jacobi_diagonal_raises", "test_positive_jacobi_diagonal_still_solves",
    "test_frozen_wedge_jacobi_reaches_true_residual",
    *(f"test_nonfinite_scalar_pcg_cannot_report_convergence[{where}-{value}]"
      for where in ("rhs", "x0") for value in ("nan", "inf")),
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def no_duplicates(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def strict_json(data):
    def invalid(value):
        raise ValueError(f"Nonfinite JSON value: {value}")
    result = json.loads(data, object_pairs_hook=no_duplicates, parse_constant=invalid)
    def visit(value):
        if isinstance(value, float):
            require(math.isfinite(value), "Nonfinite numeric value")
        elif isinstance(value, dict):
            for item in value.values():
                visit(item)
        elif isinstance(value, list):
            for item in value:
                visit(item)
    visit(result)
    return result


class GitHub:
    def __init__(self, repo):
        require(re.fullmatch(r"[\w.-]+/[\w.-]+", repo), "Invalid repository")
        self.repo = repo
        self.trees = {}

    def get(self, suffix):
        request = urllib.request.Request(
            f"https://api.github.com/repos/{self.repo}/{suffix}",
            headers={"Authorization": f"Bearer {os.environ['GH_TOKEN']}",
                     "Accept": "application/vnd.github+json"})
        with urllib.request.urlopen(request, timeout=60) as response:
            return strict_json(response.read())

    def file(self, path, commit):
        if commit not in self.trees:
            tree = self.get(f"git/trees/{commit}?recursive=1")
            require(tree.get("truncated") is False, "Cannot audit truncated acceptance tree")
            self.trees[commit] = {entry["path"]: entry for entry in tree["tree"]}
        entry = self.trees[commit].get(path, {})
        require(entry.get("type") == "blob" and entry.get("mode") == "100644",
                "Acceptance must be a tracked non-executable regular file")
        record = self.get(f"contents/{path}?ref={commit}")
        require(record.get("type") == "file" and record.get("encoding") == "base64",
                "Acceptance must be a tracked regular file")
        require(record.get("sha") == entry["sha"], "Acceptance blob identity mismatch")
        return base64.b64decode(record["content"], validate=False)


def verify_run(api, run_id, commit):
    require(SHA.fullmatch(commit), "Acceptance commit must be a full lowercase SHA")
    repo = api.get("")
    run = api.get(f"actions/runs/{run_id}")
    workflow = api.get("actions/workflows/build-test.yml")
    require(run.get("id") == run_id, "Run ID mismatch")
    require(run.get("repository", {}).get("id") == repo["id"], "Foreign repository")
    require(run.get("head_repository", {}).get("id") == repo["id"], "Foreign head repository")
    require(run.get("workflow_id") == workflow["id"], "Wrong build workflow")
    require(run.get("event") == "push" and run.get("status") == "completed"
            and run.get("conclusion") == "success", "Not a successful tag push build")
    require(SHA.fullmatch(run.get("head_sha", "")), "Invalid source SHA")
    reachable = api.get(f"compare/main...{commit}")
    require(reachable.get("status") in ("behind", "identical")
            and reachable.get("ahead_by") == 0, "Acceptance is not reachable from main")
    return run


def verify_context(api, context, run, version):
    tag = f"v{version}"
    require(context.get("schema") == "radia.ci-release-context.v1", "Unknown CI context")
    require(str(context.get("run_id")) == str(run["id"])
            and context.get("sha") == run["head_sha"], "Context identity mismatch")
    require(context.get("ref_type") == "tag" and context.get("ref") == f"refs/tags/{tag}",
            "CI did not build this exact release tag")
    require(tag in context.get("release_tags", []), "Missing immutable release tag")
    obj = api.get(f"git/ref/tags/{urllib.parse.quote(tag, safe='')}")["object"]
    for _ in range(10):
        if obj["type"] != "tag":
            break
        obj = api.get(f"git/tags/{obj['sha']}")["object"]
    require(obj["type"] == "commit" and obj["sha"] == run["head_sha"], "Remote tag moved")


def verify_host(acceptance, full, junit, identity, host):
    require(acceptance.get("passed") is True, "Host did not pass")
    for key, value in identity.items():
        require(acceptance.get(key) == value, f"Host identity mismatch: {key}")
    require(acceptance.get("sources_verified") == identity["sources_verified"], "Source inventory mismatch")
    steps = acceptance.get("steps", [])
    require(len(steps) == len(STEPS) and {s.get("name") for s in steps} == STEPS,
            "Missing or duplicate acceptance step")
    require(all(type(s.get("exit_code")) is int and s["exit_code"] == 0 for s in steps),
            "Acceptance command failed")
    require(full.get("completed") is True and full.get("version") == identity["version"],
            "Incomplete native smoke")
    require(full.get("drive") == 300000, "Wrong strong-field drive")
    rows = full.get("rows", [])
    require(len(rows) == 6 and {(r.get("kind"), r.get("order")) for r in rows} == CASES,
            "Missing or duplicate native case")
    for row in rows:
        require(type(row["order"]) is int, "Invalid element order")
        require(all(type(row[k]) in (int, float) and math.isfinite(row[k])
                    for k in ("linear_parity", "residual")), "Invalid numerical result")
        require(0 <= row["linear_parity"] < 1e-6 and 0 <= row["residual"] <= 2e-5,
                "Native numerical gate failed")
        stats = row["nonlinear_stats"]
        require(stats.get("nonlinear_final_relative_residual") == row["residual"],
                "Native residual summary disagrees with solver statistics")
        require(stats.get("nonlinear_converged_final_stage") is True
                and stats.get("nonlinear_tangent_assemblies", 0) > 0,
                "Nonlinear solve not exercised or not converged")
    require(b"<!DOCTYPE" not in junit.upper(), "DTD is forbidden")
    root = ET.fromstring(junit)
    cases = root.findall(".//testcase")
    require(len(cases) == 7 and {c.get("name") for c in cases} == TESTS,
            "Missing or duplicate focused case")
    require(not any(c.find(tag) is not None for c in cases for tag in ("failure", "error", "skipped")),
            "Focused case did not pass")
    suites = root.findall(".//testsuite")
    require(len(suites) == 1 and suites[0].get("tests") == "7"
            and all(suites[0].get(k) == "0" for k in ("failures", "errors", "skipped")),
            "Focused suite totals mismatch")
    require(suites[0].get("hostname", "").lower() == host, "Wrong acceptance host")
    require(all(math.isfinite(float(node.get("time", "nan")))
                and float(node.get("time", "nan")) >= 0 for node in suites + cases),
            "Nonfinite or invalid focused duration")


def verify_artifact(api, run, commit, wheel_dir, context_dir, expected_hash):
    require(HASH.fullmatch(expected_hash), "Wheel hash must be lowercase SHA-256")
    wheels = list(wheel_dir.rglob("*.whl"))
    require(len(wheels) == 1, "Exactly one wheel required")
    wheel = wheels[0]
    require(hashlib.sha256(wheel.read_bytes()).hexdigest() == expected_hash, "Wheel hash mismatch")
    require(re.fullmatch(r"radia-\d+\.\d+\.\d+-cp312-cp312-win_amd64.whl", wheel.name), "Unexpected wheel")
    version = wheel.name.split("-")[1]
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        require(len(names) == len(set(names)), "Duplicate wheel entry")
        native = hashlib.sha256(archive.read("radia/_radia_pybind.pyd")).hexdigest()
        require("radia/radia_motor_rom.dll" in names, "Missing motor C ABI")
        init = archive.read("radia/__init__.py").decode("utf-8")
        require(re.search(r'^__version__\s*=\s*"' + re.escape(version) + '"', init, re.M),
                "Wheel version mismatch")
        sources = sum(n.startswith("radia/") and n.endswith(".py") for n in names)
    context = strict_json((context_dir / "ci-release-context.json").read_bytes())
    verify_context(api, context, run, version)
    identity = dict(wheel_sha256=expected_hash, native_sha256=native,
                    source_commit=run["head_sha"], ci_run=run["id"],
                    sources_verified=sources, version=version)
    prefix = f"validation_test/esrf_three_engine/results/candidate_{run['head_sha'][:9]}"
    # The hosts come from the same tuple release-quad produces evidence for,
    # so this gate cannot quietly ask for fewer machines than policy requires.
    # It used to name two hosts of its own, one of which is not an acceptance
    # target at all, and would have published on two of the four.
    spec = importlib.util.spec_from_file_location(
        "radia_release_acceptance", Path(__file__).resolve().with_name("release_acceptance.py"))
    acceptance_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(acceptance_module)
    RELEASE_ACCEPTANCE_HOSTS = acceptance_module.RELEASE_ACCEPTANCE_HOSTS
    missing = []
    for host in RELEASE_ACCEPTANCE_HOSTS:
        try:
            acceptance = strict_json(api.file(f"{prefix}/{host}/acceptance.json", commit))
            full = strict_json(api.file(f"{prefix}/{host}/full6.json", commit))
            focused = api.file(f"{prefix}/{host}/focused.xml", commit)
        except Exception as exc:  # noqa: BLE001 -- report every missing host, then fail
            missing.append(f"{host} ({type(exc).__name__})")
            continue
        verify_host(acceptance, full, focused, identity, host)
    require(not missing,
            f"Acceptance evidence missing for {', '.join(missing)}; release-quad "
            f"requires all of {', '.join(RELEASE_ACCEPTANCE_HOSTS)} at {prefix}")
    return dict(identity, acceptance_commit=commit,
                hosts=list(RELEASE_ACCEPTANCE_HOSTS), passed=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", required=True)
    p.add_argument("--run-id", type=int, required=True)
    p.add_argument("--acceptance-commit", required=True)
    p.add_argument("--wheel-sha256", required=True)
    p.add_argument("--wheel-dir", type=Path, required=True)
    p.add_argument("--context-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    api = GitHub(args.repo)
    run = verify_run(api, args.run_id, args.acceptance_commit)
    report = verify_artifact(api, run, args.acceptance_commit, args.wheel_dir,
                             args.context_dir, args.wheel_sha256)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
