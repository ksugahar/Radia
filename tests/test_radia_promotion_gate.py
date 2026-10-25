"""Guards: publication rejects unapproved tags, substituted bytes and incomplete evidence.

Assumes: remote GitHub identities are mocked only in offline contracts; existing
host acceptance records retain their independently approved numerical criteria.
"""
import copy
import importlib.util
import hashlib
import json
import io
from pathlib import Path
import zipfile
import xml.etree.ElementTree as ET

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("promotion_gate", ROOT / "tools/verify_radia_promotion.py")
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)
EVIDENCE = ROOT / "validation_test/esrf_three_engine/results/candidate_ed58b5deb"


@pytest.fixture(autouse=True)
def configured_acceptance_roles(monkeypatch):
    names = {"validation": "validation-machine", "100": "build-machine"}
    monkeypatch.setattr(gate, "CURRENT_RELEASE_ACCEPTANCE_HOSTNAMES", names)
    monkeypatch.setattr(gate, "RELEASE_ACCEPTANCE_HOSTNAMES", {k: {v} for k, v in names.items()})


@pytest.mark.parametrize("suffix", ["", "actions/runs/123"])
def test_github_endpoint_has_no_empty_trailing_path(monkeypatch, suffix):
    seen = []
    def open_request(request, timeout):
        seen.append(request.full_url)
        return io.BytesIO(b'{"id": 123}')
    monkeypatch.setenv("GH_TOKEN", "test-token")
    monkeypatch.setattr(gate.urllib.request, "urlopen", open_request)
    assert gate.GitHub("owner/repo").get(suffix) == {"id": 123}
    assert seen == ["https://api.github.com/repos/owner/repo" + ("/"+suffix if suffix else "")]


def evidence(record="validation"):
    """Use public numeric records with synthetic host identity for gate tests."""
    a = gate.strict_json((EVIDENCE / record / "acceptance.json").read_bytes())
    f = gate.strict_json((EVIDENCE / record / "full6.json").read_bytes())
    x = (EVIDENCE / record / "focused.xml").read_bytes()
    tree = ET.fromstring(x)
    suite = tree if tree.tag == "testsuite" else tree.find("testsuite")
    suite.set("hostname", "validation-machine")  # Synthetic operational-role fixture only.
    x = ET.tostring(tree)
    identity = {k: a[k] for k in ("wheel_sha256", "native_sha256", "source_commit", "ci_run", "sources_verified", "version")}
    return a, f, x, identity


@pytest.mark.parametrize("record", ["validation", "compute"])
def test_committed_numeric_record_with_synthetic_host(record):
    gate.verify_host(*evidence(record), "validation")


@pytest.mark.parametrize("key", ["wheel_sha256", "native_sha256", "source_commit", "ci_run", "sources_verified", "version"])
def test_reject_substituted_identity(key):
    a, f, x, identity = evidence()
    a[key] = "wrong"
    with pytest.raises(ValueError):
        gate.verify_host(a, f, x, identity, "validation")


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "failed", "bool_exit"])
def test_reject_step_contract(mutation):
    a, f, x, identity = evidence()
    if mutation == "missing":
        a["steps"].pop()
    elif mutation == "duplicate":
        a["steps"][-1] = a["steps"][0]
    elif mutation == "bool_exit":
        a["steps"][0]["exit_code"] = False
    else:
        a["steps"][0]["exit_code"] = 1
    with pytest.raises(ValueError):
        gate.verify_host(a, f, x, identity, "validation")


@pytest.mark.parametrize("text", ['{"x":NaN}', '{"x":Infinity}', '{"x":1e999}', '{"x":1,"x":2}'])
def test_strict_json(text):
    with pytest.raises(ValueError):
        gate.strict_json(text)


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "residual", "not_converged", "no_tangent"])
def test_reject_native_failure(mutation):
    a, f, x, identity = evidence()
    if mutation == "missing":
        f["rows"].pop()
    elif mutation == "duplicate":
        f["rows"][-1] = f["rows"][0]
    elif mutation == "residual":
        f["rows"][0]["residual"] = 1
    elif mutation == "not_converged":
        f["rows"][0]["nonlinear_stats"]["nonlinear_converged_final_stage"] = False
    else:
        f["rows"][0]["nonlinear_stats"]["nonlinear_tangent_assemblies"] = 0
    with pytest.raises(ValueError):
        gate.verify_host(a, f, x, identity, "validation")


def test_host_and_focused_failure():
    a, f, x, identity = evidence()
    with pytest.raises(ValueError):
        gate.verify_host(a, f, x, identity, "hibino")
    with pytest.raises(ValueError):
        gate.verify_host(a, f, x.replace(b'failures="0"', b'failures="1"'), identity, "validation")


def test_strengthened_focused_suite_and_reject_unknown_case():
    a, f, x, identity = evidence()
    root = ET.fromstring(x)
    suite = root if root.tag == "testsuite" else root.find("testsuite")
    extra = ET.SubElement(suite, "testcase")
    extra.set("name", "test_indefinite_operator_with_positive_jacobi_diagonal_raises_in_batched_pcg")
    extra.set("time", "0.1")
    suite.set("tests", str(int(suite.get("tests")) + 1))
    strengthened = ET.tostring(root)
    gate.verify_host(a, f, strengthened, identity, "validation")
    extra.set("name", "test_unreviewed_extra_case")
    with pytest.raises(ValueError, match="Missing or duplicate focused case"):
        gate.verify_host(a, f, ET.tostring(root), identity, "validation")


@pytest.mark.parametrize("host,hostname", [("validation", "validation-machine"), ("100", "build-machine")])
def test_configured_os_hostnames(host, hostname):
    a, f, xml, identity = evidence()
    root = ET.fromstring(xml)
    root.find('.//testsuite').set('hostname', hostname)
    actual = ET.tostring(root)
    gate.verify_host(a, f, actual, identity, host)
    other = "100" if host == "validation" else "validation"
    with pytest.raises(ValueError, match='Wrong acceptance host'):
        gate.verify_host(a, f, actual, identity, other)


def private_current_evidence(host):
    a, f, xml, identity = evidence()
    a["version"] = f["version"] = identity["version"] = "5.3.1"
    tree = ET.fromstring(xml)
    suite = tree.find(".//testsuite")
    suite.set("hostname", gate.CURRENT_RELEASE_ACCEPTANCE_HOSTNAMES[host])
    names = {case.get("name") for case in suite.findall("testcase")}
    for name in gate.STRENGTHENED_TESTS - names:
        ET.SubElement(suite, "testcase", name=name, time="0.0")
    suite.set("tests", str(len(gate.STRENGTHENED_TESTS)))
    return a, f, ET.tostring(tree), identity


def test_public_export_preserves_numbers_and_passes_privacy_policy():
    spec = importlib.util.spec_from_file_location("public_lint", ROOT / "tools/policy_lint.py")
    lint = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(lint)
    proofs = []
    for host in gate.RELEASE_ACCEPTANCE_HOSTS:
        a, f, xml, identity = private_current_evidence(host)
        public_a, public_f, public_xml = gate.prepare_public_host(a, f, xml, identity, host, b"x" * 32)
        gate.verify_host(public_a, public_f, public_xml, identity, host, public=True)
        assert public_f["rows"] == f["rows"]
        assert {k: v for k, v in public_a.items() if k != "private_runtime_verification"} == a
        prefix = "validation_test/results/" + gate.PUBLIC_RELEASE_ACCEPTANCE_ROLES[host]
        for name, content in (("acceptance.json", json.dumps(public_a)),
                              ("full6.json", json.dumps(public_f)), ("focused.xml", public_xml.decode())):
            assert not lint.public_runtime_identity_hits(prefix + "/" + name, content)
        proofs.append(public_a["private_runtime_verification"])
    gate.verify_private_batch(proofs)
    assert proofs[0]["machine_token"] != proofs[1]["machine_token"]


def test_public_export_rejects_historical_machine_for_current_role():
    a, f, xml, identity = private_current_evidence("100")
    tree = ET.fromstring(xml)
    tree.find(".//testsuite").set("hostname", "other-machine")
    with pytest.raises(ValueError, match="Wrong acceptance host"):
        gate.prepare_public_host(a, f, ET.tostring(tree), identity, "100", b"x" * 32)


@pytest.mark.parametrize("changed", ["acceptance", "full6", "focused", "unverified", "missing"])
def test_public_gate_rejects_changed_or_unverified_export(changed):
    a, f, xml, identity = private_current_evidence("validation")
    a, f, xml = gate.prepare_public_host(a, f, xml, identity, "validation", b"x" * 32)
    if changed == "acceptance":
        a["steps"][0]["elapsed_s"] = a["steps"][0].get("elapsed_s", 0) + 1
    elif changed == "full6":
        f["rows"][0]["seconds"] += 1
    elif changed == "focused":
        tree = ET.fromstring(xml)
        tree.find(".//testcase").set("time", "123.0")
        xml = ET.tostring(tree)
    elif changed == "unverified":
        a["private_runtime_verification"]["current_machine_verified"] = False
    else:
        del a["private_runtime_verification"]
    with pytest.raises(ValueError):
        gate.verify_host(a, f, xml, identity, "validation", public=True)


@pytest.mark.parametrize("changed", ["machine_token", "verification_batch", "missing_role"])
def test_public_gate_requires_distinct_machines_in_one_complete_batch(changed):
    proofs = []
    for host in gate.RELEASE_ACCEPTANCE_HOSTS:
        a, f, xml, identity = private_current_evidence(host)
        a, _, _ = gate.prepare_public_host(a, f, xml, identity, host, b"x" * 32)
        proofs.append(a["private_runtime_verification"])
    if changed == "machine_token":
        proofs[1][changed] = proofs[0][changed]
    elif changed == "verification_batch":
        proofs[1][changed] = "a" * 64
    else:
        proofs.pop()
    with pytest.raises(ValueError):
        gate.verify_private_batch(proofs)


@pytest.mark.parametrize("key", ["residual", "linear_parity", "order"])
def test_bool_is_not_a_numerical_result(key):
    a, f, x, identity = evidence()
    f["rows"][0][key] = False if key != "order" else True
    with pytest.raises(ValueError):
        gate.verify_host(a, f, x, identity, "validation")


class FakeAPI:
    def __init__(self):
        self.records = {
            "": {"id": 42},
            "actions/workflows/build-test.yml": {"id": 8},
            "actions/runs/9": {"id": 9, "repository": {"id": 42}, "head_repository": {"id": 42}, "workflow_id": 8, "event": "push", "status": "completed", "conclusion": "success", "head_sha": "a" * 40},
            "compare/main..." + "b" * 40: {"status": "behind", "ahead_by": 0},
            "git/ref/tags/v5.0.0": {"object": {"type": "tag", "sha": "c" * 40}},
            "git/tags/" + "c" * 40: {"object": {"type": "commit", "sha": "a" * 40}},
        }

        for workflow in ("radia.yml", "policy-lint.yml"):
            self.records[f"actions/workflows/{workflow}/runs?head_sha={'b' * 40}&event=push&per_page=100"] = {
                "workflow_runs": [{"head_sha": "b" * 40, "head_branch": "main",
                                   "run_number": 1, "run_attempt": 1,
                                   "status": "completed", "conclusion": "success"}]}

    def get(self, path):
        return copy.deepcopy(self.records[path])


@pytest.mark.parametrize("key,value", [("event", "workflow_dispatch"), ("conclusion", "failure"), ("status", "in_progress"), ("workflow_id", 10), ("head_sha", "main"), ("repository", {"id": 1}), ("head_repository", {"id": 1})])
def test_run_api_rejection(key, value):
    api = FakeAPI()
    api.records["actions/runs/9"][key] = value
    with pytest.raises(ValueError):
        gate.verify_run(api, 9, "b" * 40)


def test_acceptance_commit_and_tag_identity():
    api = FakeAPI()
    run = gate.verify_run(api, 9, "b" * 40)
    ctx = dict(schema="radia.ci-release-context.v1", run_id=9, sha="a" * 40,
               ref_type="tag", ref="refs/tags/v5.0.0", release_tags=["v5.0.0"])
    gate.verify_context(api, ctx, run, "5.0.0")
    api.records["git/tags/" + "c" * 40]["object"]["sha"] = "d" * 40
    with pytest.raises(ValueError):
        gate.verify_context(api, ctx, run, "5.0.0")
    with pytest.raises(ValueError):
        gate.verify_run(api, 9, "main")
    api.records["compare/main..." + "b" * 40] = {"status": "ahead", "ahead_by": 1}
    with pytest.raises(ValueError):
        gate.verify_run(api, 9, "b" * 40)


def test_publication_has_no_automatic_bypass():
    workflow = (ROOT / ".github/workflows/release.yml").read_text(encoding="utf-8")
    before, publish = workflow.split("  publish-pypi:")
    assert "id-token: write" not in before
    assert "pypa/gh-action-pypi-publish" not in before
    assert "needs: verify-promotion" in publish
    assert "github.event_name == 'workflow_dispatch'" in publish
    assert "needs.verify-promotion.result == 'success'" in publish
    assert "environment: pypi" in publish
    assert "ref: ${{ github.sha }}" in publish
    assert "persist-credentials: false" in publish
    assert "--after-approval" in publish
    assert "Build" not in publish


def test_public_exact_wheel_end_to_end(tmp_path):
    wheel_dir = tmp_path / "dist"
    wheel_dir.mkdir()
    wheel = wheel_dir / "radia-5.3.1-cp312-cp312-win_amd64.whl"
    with zipfile.ZipFile(wheel, "w") as z:
        z.writestr("radia/__init__.py", '__version__ = "5.3.1"\n')
        z.writestr("radia/_radia_pybind.pyd", b"native")
        z.writestr("radia/radia_motor_rom.dll", b"motor")
    digest = hashlib.sha256(wheel.read_bytes()).hexdigest()
    api = FakeAPI()
    api.records["git/ref/tags/v5.3.1"] = api.records["git/ref/tags/v5.0.0"]
    run = gate.verify_run(api, 9, "b" * 40)
    ctx = dict(schema="radia.ci-release-context.v1", run_id=9, sha="a" * 40,
               ref_type="tag", ref="refs/tags/v5.3.1", release_tags=["v5.3.1"])
    (tmp_path / "ci-release-context.json").write_text(json.dumps(ctx))
    records = {}
    for host, role in gate.PUBLIC_RELEASE_ACCEPTANCE_ROLES.items():
        a, f, x, identity = private_current_evidence(host)
        identity.update(wheel_sha256=digest, native_sha256=hashlib.sha256(b"native").hexdigest(),
                        source_commit="a" * 40, ci_run=9, sources_verified=1)
        a.update(identity)
        a, f, x = gate.prepare_public_host(a, f, x, identity, host, b"x" * 32)
        records[role] = {"acceptance.json": json.dumps(a).encode(),
                         "full6.json": json.dumps(f).encode(), "focused.xml": x}

    def read_file(path, commit):
        assert commit == "b" * 40
        role, name = path.split("/")[-2:]
        return records[role][name]

    api.file = read_file
    report = gate.verify_artifact(api, run, "b" * 40, wheel_dir, tmp_path, digest)
    assert report["passed"] and report["execution_roles"] == ["validation", "build"]


def test_exact_wheel_end_to_end(tmp_path):
    wheel_dir = tmp_path / "dist"
    wheel_dir.mkdir()
    wheel = wheel_dir / "radia-5.0.0-cp312-cp312-win_amd64.whl"
    with zipfile.ZipFile(wheel, "w") as z:
        z.writestr("radia/__init__.py", '__version__ = "5.0.0"\n')
        z.writestr("radia/_radia_pybind.pyd", b"native")
        z.writestr("radia/radia_motor_rom.dll", b"motor")
    digest = hashlib.sha256(wheel.read_bytes()).hexdigest()
    api = FakeAPI()
    run = gate.verify_run(api, 9, "b" * 40)
    ctx = dict(schema="radia.ci-release-context.v1", run_id=9, sha="a" * 40,
               ref_type="tag", ref="refs/tags/v5.0.0", release_tags=["v5.0.0"])
    (tmp_path / "ci-release-context.json").write_text(json.dumps(ctx))

    def read_file(path, commit):
        assert commit == "b" * 40
        assert path.startswith("validation_test/esrf_three_engine/results/candidate_aaaaaaaaa/")
        host, name = path.split("/")[-2:]
        # Synthetic two-host API fixture, not a claim of machine acceptance.
        a, f, x, _ = evidence()
        tree = ET.fromstring(x)
        suite = tree if tree.tag == "testsuite" else tree.find("testsuite")
        suite.set("hostname", gate.CURRENT_RELEASE_ACCEPTANCE_HOSTNAMES[host])
        x = ET.tostring(tree)
        if name == "acceptance.json":
            a.update(wheel_sha256=digest, native_sha256=hashlib.sha256(b"native").hexdigest(),
                     source_commit="a" * 40, ci_run=9, sources_verified=1)
            return json.dumps(a).encode()
        return json.dumps(f).encode() if name == "full6.json" else x

    api.file = read_file
    report = gate.verify_artifact(api, run, "b" * 40, wheel_dir, tmp_path, digest)
    assert report["hosts"] == ["validation", "100"] and report["passed"]
    for missing_host in report["hosts"]:
        def incomplete(path, commit):
            if path.split("/")[-2] == missing_host:
                raise FileNotFoundError(path)
            return read_file(path, commit)
        api.file = incomplete
        with pytest.raises(ValueError, match="Acceptance evidence missing for " + missing_host):
            gate.verify_artifact(api, run, "b" * 40, wheel_dir, tmp_path, digest)
    api.file = read_file
    with pytest.raises(ValueError, match="hash mismatch"):
        gate.verify_artifact(api, run, "b" * 40, wheel_dir, tmp_path, "0" * 64)
    def same_machine(path, commit):
        data = read_file(path, commit)
        if path.endswith("focused.xml"):
            tree = ET.fromstring(data)
            tree.find(".//testsuite").set("hostname", "other-machine")
            return ET.tostring(tree)
        return data
    def old_fleet(path, commit):
        data = read_file(path, commit)
        if path.endswith("focused.xml"):
            tree = ET.fromstring(data)
            host = path.split("/")[-2]
            tree.find(".//testsuite").set("hostname", "validation-machine" if host == "validation" else "other-machine")
            return ET.tostring(tree)
        return data
    api.file = old_fleet
    with pytest.raises(ValueError, match="Wrong acceptance host"):
        gate.verify_artifact(api, run, "b" * 40, wheel_dir, tmp_path, digest)
    api.file = same_machine
    with pytest.raises(ValueError, match="Wrong acceptance host"):
        gate.verify_artifact(api, run, "b" * 40, wheel_dir, tmp_path, digest)
    api.file = read_file
    ctx["ref_type"] = "branch"
    (tmp_path / "ci-release-context.json").write_text(json.dumps(ctx))
    with pytest.raises(ValueError, match="exact release tag"):
        gate.verify_artifact(api, run, "b" * 40, wheel_dir, tmp_path, digest)


@pytest.mark.parametrize("mode", ["120000", "160000", "100755"])
def test_acceptance_not_symlink_submodule_or_executable(mode):
    api = gate.GitHub("ksugahar/Radia")
    api.get = lambda _: {"truncated": False, "tree": [{"path": "evidence.json", "type": "blob", "mode": mode}]}
    with pytest.raises(ValueError, match="regular file"):
        api.file("evidence.json", "b" * 40)


def test_new_releases_require_the_strengthened_focused_suite():
    a, f, xml, identity = evidence()
    for version in ("5.1.0", "5.1.1", "6.0.0"):
        a["version"] = f["version"] = identity["version"] = version
        with pytest.raises(ValueError, match="Missing or duplicate focused case"):
            gate.verify_host(a, f, xml, identity, "validation")
        root = ET.fromstring(xml)
        suite = root if root.tag == "testsuite" else root.find("testsuite")
        extra = ET.SubElement(suite, "testcase")
        extra.set("name", "test_indefinite_operator_with_positive_jacobi_diagonal_raises_in_batched_pcg")
        extra.set("time", "0.1")
        suite.set("tests", str(int(suite.get("tests")) + 1))
        gate.verify_host(a, f, ET.tostring(root), identity, "validation")


@pytest.mark.parametrize("status,conclusion", [("in_progress", None), ("completed", "failure"),
                                               ("completed", "cancelled"), ("completed", "skipped")])
def test_publication_waits_for_acceptance_commit_ci(status, conclusion):
    api = FakeAPI()
    key = f"actions/workflows/radia.yml/runs?head_sha={'b' * 40}&event=push&per_page=100"
    api.records[key]["workflow_runs"][0].update(status=status, conclusion=conclusion)
    with pytest.raises(ValueError, match="Acceptance CI not successful"):
        gate.verify_run(api, 9, "b" * 40)


def test_old_success_cannot_hide_a_failed_ci_rerun():
    api = FakeAPI()
    key = f"actions/workflows/radia.yml/runs?head_sha={'b' * 40}&event=push&per_page=100"
    failed = dict(api.records[key]["workflow_runs"][0], run_attempt=2, conclusion="failure")
    api.records[key]["workflow_runs"].append(failed)
    with pytest.raises(ValueError, match="Acceptance CI not successful"):
        gate.verify_run(api, 9, "b" * 40)


def application_proof():
    return {"schema": "radia.release-quad.simulink-candidate.v1", "version": "5.1.0",
            "commit": "a"*40, "package_sha256": "b"*64, "mex_sha256": "c"*64,
            "targets": {host: {"label": gate.RELEASE_ACCEPTANCE_HOST_LABELS[host],
                              "status": "passed", "mex_sha256": "c"*64,
                              "verified_at_utc": "2026-09-30T00:00:00Z"}
                        for host in gate.RELEASE_ACCEPTANCE_HOSTS}}


def test_application_acceptance_requires_two_exact_mex_results():
    identity = {"source_commit": "a"*40, "version": "5.1.0"}
    gate.verify_application_acceptance(application_proof(), identity)
    for mutation in ("missing_host", "failed", "wrong_mex", "wrong_source", "missing_time",
                     "wrong_label"):
        proof = application_proof()
        if mutation == "missing_host": del proof["targets"]["validation"]
        elif mutation == "wrong_label": proof["targets"]["100"]["label"] = "mdx1"
        elif mutation == "failed": proof["targets"]["100"]["status"] = "failed"
        elif mutation == "wrong_mex": proof["targets"]["validation"]["mex_sha256"] = "d"*64
        elif mutation == "wrong_source": proof["commit"] = "e"*40
        else: del proof["targets"]["validation"]["verified_at_utc"]
        with pytest.raises(ValueError):
            gate.verify_application_acceptance(proof, identity)


def test_cpu_promotion_does_not_request_matlab_acceptance(monkeypatch, tmp_path):
    wheel_dir = tmp_path / "wheel"
    wheel_dir.mkdir()
    wheel = wheel_dir / "radia-5.3.1-cp312-cp312-win_amd64.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr("radia/__init__.py", '__version__ = "5.3.1"\n')
        archive.writestr("radia/_radia_pybind.pyd", b"native")
        archive.writestr("radia/radia_motor_rom.dll", b"motor")
    (tmp_path / "ci-release-context.json").write_text("{}")
    monkeypatch.setattr(gate, "verify_context", lambda *args: None)
    monkeypatch.setattr(gate, "verify_host", lambda *args, **kwargs: None)

    class CPUAPI:
        def file(self, path, commit):
            assert path.rsplit("/", 1)[-1] in {"acceptance.json", "full6.json", "focused.xml"}
            role = path.split("/")[-2]
            assert role in gate.PUBLIC_RELEASE_ACCEPTANCE_ROLES.values()
            if path.endswith("/focused.xml"):
                return b'<testsuites><testsuite/></testsuites>'
            if path.endswith("/acceptance.json"):
                return json.dumps({"private_runtime_verification": dict(
                    execution_role=role, verification_batch="a" * 64,
                    machine_token=hashlib.sha256(role.encode()).hexdigest())}).encode()
            return b"{}"

    report = gate.verify_artifact(CPUAPI(), {"head_sha": "a" * 40, "id": 9},
                                 "b" * 40, wheel_dir, tmp_path,
                                 hashlib.sha256(wheel.read_bytes()).hexdigest())
    assert report["passed"]


def _publication_attempt():
    class API:
        repo = "owner/repo"
        def get(self, path):
            assert path == "git/ref/tags/v5.3.1"
            return {"object": {"type": "commit", "sha": "a" * 40}}
    env = dict(GITHUB_EVENT_NAME="workflow_dispatch", GITHUB_REF_TYPE="tag",
               GITHUB_REF="refs/tags/v5.3.1", GITHUB_SHA="a" * 40,
               GITHUB_REPOSITORY="owner/repo", GITHUB_RUN_ID="99", GITHUB_RUN_ATTEMPT="1")
    return API(), env


@pytest.mark.parametrize("key,value", [("GITHUB_EVENT_NAME", "workflow_run"),
    ("GITHUB_REF_TYPE", "branch"), ("GITHUB_REF", "refs/heads/main"),
    ("GITHUB_REF", "refs/tags/cubit-mesh-export-v2.2.7"),
    ("GITHUB_SHA", "b" * 40), ("GITHUB_REPOSITORY", "fork/repo"),
    ("GITHUB_RUN_ID", ""), ("GITHUB_RUN_ATTEMPT", "")])
def test_actual_publication_identity_rejects_non_tag_or_missing_identity(key, value):
    api, env = _publication_attempt()
    env[key] = value
    with pytest.raises(ValueError):
        gate.publication_identity(api, env)


@pytest.mark.parametrize("change", [None, "attempt", "wheel", "native", "source", "acceptance", "producer"])
def test_approved_publication_rechecks_attempt_and_exact_artifacts(tmp_path, change):
    api, env = _publication_attempt()
    wheel_dir, native_dir = tmp_path / "dist", tmp_path / "native"
    wheel_dir.mkdir(); native_dir.mkdir()
    wheel = wheel_dir / "radia-5.3.1-cp312-cp312-win_amd64.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr("radia/_radia_pybind.pyd", b"accepted native")
    native = native_dir / "_radia_pybind.pyd"
    native.write_bytes(b"accepted native")
    expected = hashlib.sha256(wheel.read_bytes()).hexdigest()
    report = dict(passed=True, publication=gate.publication_identity(api, env),
                  source_commit="a" * 40, version="5.3.1", ci_run=9,
                  acceptance_commit="b" * 40, wheel_sha256=expected,
                  native_bundle=gate.native_bundle(wheel_dir, native_dir))
    if change == "attempt": env["GITHUB_RUN_ATTEMPT"] = "2"
    if change == "wheel": wheel.write_bytes(b"substitution")
    if change == "native": native.write_bytes(b"substitution")
    if change == "source": report["source_commit"] = "c" * 40
    if change == "acceptance": report["acceptance_commit"] = "c" * 40
    if change == "producer": report["ci_run"] = 10
    if change:
        with pytest.raises((ValueError, zipfile.BadZipFile)):
            gate.verify_publication(api, report, env, wheel_dir, native_dir, 9, "b" * 40, expected)
    else:
        gate.verify_publication(api, report, env, wheel_dir, native_dir, 9, "b" * 40, expected)


def test_publication_cli_binds_before_and_after_approval(tmp_path, monkeypatch):
    import sys
    api, env = _publication_attempt()
    for key, value in env.items(): monkeypatch.setenv(key, value)
    monkeypatch.setattr(gate, "GitHub", lambda repo: api)
    monkeypatch.setattr(gate, "verify_run", lambda *args: {"head_sha": "a" * 40})
    wheel_dir, native_dir = tmp_path / "dist", tmp_path / "native"
    wheel_dir.mkdir(); native_dir.mkdir()
    wheel = wheel_dir / "radia-5.3.1-cp312-cp312-win_amd64.whl"
    with zipfile.ZipFile(wheel, "w") as archive: archive.writestr("radia/_radia_pybind.pyd", b"accepted")
    (native_dir / "_radia_pybind.pyd").write_bytes(b"accepted")
    digest = hashlib.sha256(wheel.read_bytes()).hexdigest()
    monkeypatch.setattr(gate, "verify_artifact", lambda *args: dict(passed=True, source_commit="a" * 40, version="5.3.1", ci_run=9, acceptance_commit="b" * 40, wheel_sha256=digest))
    args = ["gate", "--repo", api.repo, "--run-id", "9", "--acceptance-commit", "b" * 40, "--wheel-sha256", digest,
            "--wheel-dir", str(wheel_dir), "--context-dir", str(tmp_path), "--native-dir", str(native_dir), "--output", str(tmp_path / "promotion.json")]
    monkeypatch.setattr(sys, "argv", args)
    gate.main()
    monkeypatch.setattr(sys, "argv", args + ["--after-approval"])
    gate.main()
    monkeypatch.setenv("GITHUB_RUN_ATTEMPT", "2")
    with pytest.raises(ValueError, match="attempt identity changed"):
        gate.main()
