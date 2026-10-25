"""Read-only exact-artifact publication gate; never execute acceptance data."""
from __future__ import annotations

import argparse
import base64
import hashlib
import hmac
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

try:
    from release_acceptance import RELEASE_ACCEPTANCE_HOSTS
    from release_acceptance import PUBLIC_RELEASE_ACCEPTANCE_ROLES
    from release_acceptance import RELEASE_ACCEPTANCE_HOSTNAMES
    from release_acceptance import CURRENT_RELEASE_ACCEPTANCE_HOSTNAMES
    from release_acceptance import RELEASE_ACCEPTANCE_HOST_LABELS
except ModuleNotFoundError:  # Loaded by path in unit tests rather than run as a script.
    _acceptance_spec = importlib.util.spec_from_file_location(
        "radia_release_acceptance", Path(__file__).resolve().with_name("release_acceptance.py"))
    _acceptance_module = importlib.util.module_from_spec(_acceptance_spec)
    _acceptance_spec.loader.exec_module(_acceptance_module)
    RELEASE_ACCEPTANCE_HOSTNAMES = _acceptance_module.RELEASE_ACCEPTANCE_HOSTNAMES
    CURRENT_RELEASE_ACCEPTANCE_HOSTNAMES = _acceptance_module.CURRENT_RELEASE_ACCEPTANCE_HOSTNAMES
    RELEASE_ACCEPTANCE_HOST_LABELS = _acceptance_module.RELEASE_ACCEPTANCE_HOST_LABELS
    RELEASE_ACCEPTANCE_HOSTS = _acceptance_module.RELEASE_ACCEPTANCE_HOSTS
    PUBLIC_RELEASE_ACCEPTANCE_ROLES = _acceptance_module.PUBLIC_RELEASE_ACCEPTANCE_ROLES

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
STRENGTHENED_TESTS = TESTS | {
    "test_indefinite_operator_with_positive_jacobi_diagonal_raises_in_batched_pcg",
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
        url = f"https://api.github.com/repos/{self.repo}"
        if suffix:
            url += "/" + suffix
        request = urllib.request.Request(
            url,
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


def verify_acceptance_ci(api, commit):
    """Require successful main CI on the commit containing acceptance records."""
    for workflow in ("radia.yml", "policy-lint.yml"):
        query = f"actions/workflows/{workflow}/runs?head_sha={commit}&event=push&per_page=100"
        runs = api.get(query).get("workflow_runs", [])
        matching = [run for run in runs if run.get("head_sha") == commit
                    and run.get("head_branch") == "main"]
        require(matching, f"No main acceptance CI: {workflow}")
        latest = max(matching, key=lambda run: (run.get("run_number", 0), run.get("run_attempt", 0)))
        require(latest.get("status") == "completed" and latest.get("conclusion") == "success",
                f"Acceptance CI not successful: {workflow}")


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
    verify_acceptance_ci(api, commit)
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


def _json_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def _xml_digest(value):
    return hashlib.sha256(ET.canonicalize(value.decode()).encode()).hexdigest()


def verify_public_attestation(acceptance, full, junit, host):
    """Check the maintainer's private-machine verification and public byte binding.

    Like the original committed acceptance records, this is a maintainer
    attestation, not remote hardware authentication. A release-private HMAC key
    gives the same physical machine one opaque token across both target roles.
    Raw records and the key remain in private release storage.
    """
    proof = acceptance.get("private_runtime_verification", {})
    require(proof.get("schema") == "radia.private-runtime-verification.v1",
            "Missing private runtime verification")
    require(proof.get("execution_role") == PUBLIC_RELEASE_ACCEPTANCE_ROLES[host]
            and proof.get("current_machine_verified") is True,
            "Wrong or unverified execution role")
    for key in ("verification_batch", "machine_token"):
        require(isinstance(proof.get(key), str) and HASH.fullmatch(proof[key]),
                "Invalid private verification token")
    raw = proof.get("raw_records_sha256", {})
    require(set(raw) == {"acceptance", "full6", "focused"}
            and all(isinstance(v, str) and HASH.fullmatch(v) for v in raw.values()),
            "Missing raw evidence hashes")
    clean = {k: v for k, v in acceptance.items() if k != "private_runtime_verification"}
    require(proof.get("public_records_sha256") == dict(
        acceptance=_json_digest(clean), full6=_json_digest(full), focused=_xml_digest(junit)),
        "Public evidence differs from verified export")
    require(not any(k.lower() in {"hostname", "host", "host_name"}
                    for node in ET.fromstring(junit).iter() for k in node.attrib),
            "Public XML contains private host metadata")
    return proof


def prepare_public_host(acceptance, full, junit, identity, host, batch_key):
    """Verify raw records on the current fleet before removing identifying data."""
    require(len(batch_key) >= 32, "Private verification key is too short")
    verify_host(acceptance, full, junit, identity, host)
    tree = ET.fromstring(junit)
    suite = tree if tree.tag == "testsuite" else tree.find("testsuite")
    hostname = suite.get("hostname", "").lower()
    require(hostname == CURRENT_RELEASE_ACCEPTANCE_HOSTNAMES[host],
            "Current acceptance host required")
    a, f = json.loads(json.dumps(acceptance)), json.loads(json.dumps(full))
    require("private_runtime_verification" not in a, "Expected original raw acceptance")
    f["radia"] = "<installed-wheel>/radia/__init__.py"
    for node in tree.iter():
        for key in list(node.attrib):
            if key.lower() in {"hostname", "host", "host_name"}:
                del node.attrib[key]
    public_xml = ET.tostring(tree, encoding="utf-8", xml_declaration=True)
    a["private_runtime_verification"] = dict(
        schema="radia.private-runtime-verification.v1",
        execution_role=PUBLIC_RELEASE_ACCEPTANCE_ROLES[host], current_machine_verified=True,
        verification_batch=hashlib.sha256(batch_key).hexdigest(),
        machine_token=hmac.new(batch_key, hostname.encode(), hashlib.sha256).hexdigest(),
        raw_records_sha256=dict(acceptance=_json_digest(acceptance), full6=_json_digest(full),
                               focused=_xml_digest(junit)),
        public_records_sha256=dict(acceptance=_json_digest(a), full6=_json_digest(f),
                                  focused=_xml_digest(public_xml)))
    verify_host(a, f, public_xml, identity, host, public=True)
    return a, f, public_xml


def verify_private_batch(proofs):
    require(len(proofs) == len(RELEASE_ACCEPTANCE_HOSTS)
            and {p["execution_role"] for p in proofs} == set(PUBLIC_RELEASE_ACCEPTANCE_ROLES.values()),
            "All execution roles required")
    require(len({p["verification_batch"] for p in proofs}) == 1,
            "Acceptance records were not verified in one private batch")
    require(len({p["machine_token"] for p in proofs}) == len(RELEASE_ACCEPTANCE_HOSTS),
            "Distinct physical acceptance machines required")


def verify_host(acceptance, full, junit, identity, host, *, public=False):
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
    case_names = {c.get("name") for c in cases}
    version_parts = tuple(int(part) for part in identity["version"].split("."))
    allowed_suites = (STRENGTHENED_TESTS,) if version_parts >= (5, 1, 0) else (TESTS, STRENGTHENED_TESTS)
    require(len(cases) == len(case_names) and case_names in allowed_suites,
            "Missing or duplicate focused case")
    require(not any(c.find(tag) is not None for c in cases for tag in ("failure", "error", "skipped")),
            "Focused case did not pass")
    suites = root.findall(".//testsuite")
    require(len(suites) == 1 and suites[0].get("tests") == str(len(cases))
            and all(suites[0].get(k) == "0" for k in ("failures", "errors", "skipped")),
            "Focused suite totals mismatch")
    if public:
        verify_public_attestation(acceptance, full, junit, host)
    else:
        hostnames = RELEASE_ACCEPTANCE_HOSTNAMES.get(host, {host})
        require(suites[0].get("hostname", "").lower() in hostnames, "Wrong acceptance host")
    require(all(math.isfinite(float(node.get("time", "nan")))
                and float(node.get("time", "nan")) >= 0 for node in suites + cases),
            "Nonfinite or invalid focused duration")


def verify_application_acceptance(proof, identity):
    require(proof.get("schema") == "radia.release-quad.simulink-candidate.v1",
            "Unknown MATLAB/Simulink acceptance schema")
    require(proof.get("commit") == identity["source_commit"]
            and proof.get("version") == identity["version"], "MATLAB/Simulink source identity mismatch")
    for key in ("package_sha256", "mex_sha256"):
        require(isinstance(proof.get(key), str) and HASH.fullmatch(proof[key]),
                f"Missing MATLAB/Simulink artifact identity: {key}")
    targets = proof.get("targets", {})
    require(set(targets) == set(RELEASE_ACCEPTANCE_HOSTS), "MATLAB/Simulink requires both hosts")
    for host, row in targets.items():
        # Rows are keyed by host and carry the human label simulink-candidate prints.
        require(row.get("label") == RELEASE_ACCEPTANCE_HOST_LABELS[host]
                and row.get("status") == "passed",
                f"MATLAB/Simulink failed or missing: {host}")
        require(row.get("mex_sha256") == proof["mex_sha256"], f"MEX identity mismatch: {host}")
        require(bool(row.get("verified_at_utc")), f"Missing MATLAB/Simulink verification time: {host}")


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
    missing = []
    public = tuple(map(int, version.split("."))) >= (5, 3, 1)
    private_proofs = []
    for host in RELEASE_ACCEPTANCE_HOSTS:
        directory = PUBLIC_RELEASE_ACCEPTANCE_ROLES[host] if public else host
        try:
            acceptance = strict_json(api.file(f"{prefix}/{directory}/acceptance.json", commit))
            full = strict_json(api.file(f"{prefix}/{directory}/full6.json", commit))
            focused = api.file(f"{prefix}/{directory}/focused.xml", commit)
        except Exception as exc:  # noqa: BLE001 -- report every missing host, then fail
            missing.append(f"{host} ({type(exc).__name__})")
            continue
        if public:
            verify_host(acceptance, full, focused, identity, host, public=True)
            private_proofs.append(acceptance["private_runtime_verification"])
        else:
            verify_host(acceptance, full, focused, identity, host)
            hostname = ET.fromstring(focused).find(".//testsuite").get("hostname").lower()
            require(hostname == CURRENT_RELEASE_ACCEPTANCE_HOSTNAMES[host],
                    "Current acceptance host required")
    require(not missing,
            f"Acceptance evidence missing for {', '.join(missing)}; release-quad "
            f"requires all of {', '.join(RELEASE_ACCEPTANCE_HOSTS)} at {prefix}")
    if public:
        verify_private_batch(private_proofs)
    # MATLAB/Simulink acceptance belongs to the independent radia-simulink publisher.
    roles = ({"execution_roles": list(PUBLIC_RELEASE_ACCEPTANCE_ROLES.values())} if public
             else {"hosts": list(RELEASE_ACCEPTANCE_HOSTS)})
    return dict(identity, acceptance_commit=commit, **roles, passed=True)



def publication_identity(api, env):
    """Bind the publishing attempt, independently of the native producer run."""
    require(env.get("GITHUB_EVENT_NAME") == "workflow_dispatch", "Tag dispatch required")
    ref = env.get("GITHUB_REF", "")
    require(env.get("GITHUB_REF_TYPE") == "tag" and re.fullmatch(r"refs/tags/v\d+\.\d+\.\d+", ref), "Exact release tag ref required")
    sha = env.get("GITHUB_SHA", "")
    require(SHA.fullmatch(sha), "Invalid publishing source SHA")
    require(env.get("GITHUB_REPOSITORY") == api.repo, "Publishing repository mismatch")
    for key in ("GITHUB_RUN_ID", "GITHUB_RUN_ATTEMPT"):
        require(re.fullmatch(r"[1-9]\d*", env.get(key, "")), "Missing publishing run identity")
    obj = api.get("git/ref/tags/" + urllib.parse.quote(ref.removeprefix("refs/tags/"), safe=""))["object"]
    for _ in range(10):
        if obj["type"] != "tag":
            break
        obj = api.get(f"git/tags/{obj['sha']}")["object"]
    require(obj["type"] == "commit" and obj["sha"] == sha, "Publishing tag moved or deleted")
    return {key: env[key] for key in ("GITHUB_REPOSITORY", "GITHUB_REF", "GITHUB_SHA", "GITHUB_RUN_ID", "GITHUB_RUN_ATTEMPT")}


def native_bundle(wheel_dir, native_dir):
    """Only wheel-identical native files may be written to the binaries release."""
    wheels = list(wheel_dir.rglob("*.whl"))
    require(len(wheels) == 1, "Exactly one accepted wheel required")
    files = sorted(p for p in native_dir.rglob("*") if p.is_file())
    binaries = [p for p in files if p.suffix in {".pyd", ".dll", ".so", ".dylib"}]
    require(binaries and len({p.name for p in binaries}) == len(binaries), "Missing or duplicate native binaries")
    result = {}
    with zipfile.ZipFile(wheels[0]) as archive:
        for path in binaries:
            data = path.read_bytes()
            require(data == archive.read("radia/" + path.name), "Native artifact differs from accepted wheel: " + path.name)
            result[path.relative_to(native_dir).as_posix()] = hashlib.sha256(data).hexdigest()
    return result


def verify_publication(api, report, env, wheel_dir, native_dir, run_id, commit, expected_hash):
    require(report.get("passed") is True, "Promotion did not pass")
    require(report.get("publication") == publication_identity(api, env), "Publishing attempt identity changed; rerun all jobs and reapprove")
    require(report["source_commit"] == env["GITHUB_SHA"] and env["GITHUB_REF"] == "refs/tags/v" + report["version"], "Publishing tag/source mismatch")
    require(report["ci_run"] == run_id and report["acceptance_commit"] == commit and report["wheel_sha256"] == expected_hash, "Promotion inputs changed")
    wheels = list(wheel_dir.rglob("*.whl"))
    require(len(wheels) == 1 and hashlib.sha256(wheels[0].read_bytes()).hexdigest() == expected_hash, "Approved wheel changed")
    require(report["native_bundle"] == native_bundle(wheel_dir, native_dir), "Approved native bundle changed")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", required=True)
    p.add_argument("--run-id", type=int, required=True)
    p.add_argument("--acceptance-commit", required=True)
    p.add_argument("--wheel-sha256", required=True)
    p.add_argument("--wheel-dir", type=Path, required=True)
    p.add_argument("--context-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--native-dir", type=Path, required=True)
    p.add_argument("--after-approval", action="store_true")
    args = p.parse_args()
    api = GitHub(args.repo)
    publication = publication_identity(api, os.environ)
    if args.after_approval:
        report = strict_json(args.output.read_bytes())
        verify_publication(api, report, os.environ, args.wheel_dir, args.native_dir,
                           args.run_id, args.acceptance_commit, args.wheel_sha256)
        print("Approved publication identity and bytes verified")
        return
    run = verify_run(api, args.run_id, args.acceptance_commit)
    report = verify_artifact(api, run, args.acceptance_commit, args.wheel_dir,
                             args.context_dir, args.wheel_sha256)
    report["publication"] = publication
    report["native_bundle"] = native_bundle(args.wheel_dir, args.native_dir)
    verify_publication(api, report, os.environ, args.wheel_dir, args.native_dir,
                       args.run_id, args.acceptance_commit, args.wheel_sha256)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
