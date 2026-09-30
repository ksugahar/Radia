"""Publication gates reject incomplete, substituted and nonfinite evidence."""
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


def evidence(host="lab"):
    a = gate.strict_json((EVIDENCE / host / "acceptance.json").read_bytes())
    f = gate.strict_json((EVIDENCE / host / "full6.json").read_bytes())
    x = (EVIDENCE / host / "focused.xml").read_bytes()
    identity = {k: a[k] for k in ("wheel_sha256", "native_sha256", "source_commit", "ci_run", "sources_verified", "version")}
    return a, f, x, identity


@pytest.mark.parametrize("host", ["lab", "hibino"])
def test_real_committed_host_evidence(host):
    gate.verify_host(*evidence(host), host)


@pytest.mark.parametrize("key", ["wheel_sha256", "native_sha256", "source_commit", "ci_run", "sources_verified", "version"])
def test_reject_substituted_identity(key):
    a, f, x, identity = evidence()
    a[key] = "wrong"
    with pytest.raises(ValueError):
        gate.verify_host(a, f, x, identity, "lab")


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
        gate.verify_host(a, f, x, identity, "lab")


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
        gate.verify_host(a, f, x, identity, "lab")


def test_host_and_focused_failure():
    a, f, x, identity = evidence()
    with pytest.raises(ValueError):
        gate.verify_host(a, f, x, identity, "hibino")
    with pytest.raises(ValueError):
        gate.verify_host(a, f, x.replace(b'failures="0"', b'failures="1"'), identity, "lab")


def test_strengthened_focused_suite_and_reject_unknown_case():
    a, f, x, identity = evidence()
    root = ET.fromstring(x)
    suite = root if root.tag == "testsuite" else root.find("testsuite")
    extra = ET.SubElement(suite, "testcase")
    extra.set("name", "test_indefinite_operator_with_positive_jacobi_diagonal_raises_in_batched_pcg")
    extra.set("time", "0.1")
    suite.set("tests", str(int(suite.get("tests")) + 1))
    strengthened = ET.tostring(root)
    gate.verify_host(a, f, strengthened, identity, "lab")
    extra.set("name", "test_unreviewed_extra_case")
    with pytest.raises(ValueError, match="Missing or duplicate focused case"):
        gate.verify_host(a, f, ET.tostring(root), identity, "lab")


def test_100_role_accepts_actual_intel11_hostname_only_for_that_target():
    a, f, xml, identity = evidence()
    root = ET.fromstring(xml)
    root.find('.//testsuite').set('hostname', 'INTEL11')
    actual = ET.tostring(root)
    gate.verify_host(a, f, actual, identity, '100')
    for wrong_target in ('lab', 'mdx1', 'mdx2'):
        with pytest.raises(ValueError, match='Wrong acceptance host'):
            gate.verify_host(a, f, actual, identity, wrong_target)


@pytest.mark.parametrize("key", ["residual", "linear_parity", "order"])
def test_bool_is_not_a_numerical_result(key):
    a, f, x, identity = evidence()
    f["rows"][0][key] = False if key != "order" else True
    with pytest.raises(ValueError):
        gate.verify_host(a, f, x, identity, "lab")


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

        for workflow in ("radia-fast.yml", "policy-lint.yml"):
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
    workflow = (ROOT / ".github/workflows/release.yml").read_text()
    before, publish = workflow.split("  publish-pypi:")
    assert "id-token: write" not in before
    assert "pypa/gh-action-pypi-publish" not in before
    assert "needs: verify-promotion" in publish
    assert "github.event_name == 'workflow_dispatch'" in publish
    assert "needs.verify-promotion.result == 'success'" in publish
    assert "environment: pypi" in publish
    assert "actions/checkout" not in publish
    assert "Build" not in publish


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
        # Synthetic four-host API fixture, not a claim of machine acceptance.
        a, f, x, _ = evidence("lab")
        tree = ET.fromstring(x)
        suite = tree if tree.tag == "testsuite" else tree.find("testsuite")
        suite.set("hostname", host)
        x = ET.tostring(tree)
        if name == "acceptance.json":
            a.update(wheel_sha256=digest, native_sha256=hashlib.sha256(b"native").hexdigest(),
                     source_commit="a" * 40, ci_run=9, sources_verified=1)
            return json.dumps(a).encode()
        return json.dumps(f).encode() if name == "full6.json" else x

    api.file = read_file
    report = gate.verify_artifact(api, run, "b" * 40, wheel_dir, tmp_path, digest)
    assert report["hosts"] == ["lab", "100", "mdx1", "mdx2"] and report["passed"]
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
            gate.verify_host(a, f, xml, identity, "lab")
        root = ET.fromstring(xml)
        suite = root if root.tag == "testsuite" else root.find("testsuite")
        extra = ET.SubElement(suite, "testcase")
        extra.set("name", "test_indefinite_operator_with_positive_jacobi_diagonal_raises_in_batched_pcg")
        extra.set("time", "0.1")
        suite.set("tests", str(int(suite.get("tests")) + 1))
        gate.verify_host(a, f, ET.tostring(root), identity, "lab")


@pytest.mark.parametrize("status,conclusion", [("in_progress", None), ("completed", "failure"),
                                               ("completed", "cancelled"), ("completed", "skipped")])
def test_publication_waits_for_acceptance_commit_ci(status, conclusion):
    api = FakeAPI()
    key = f"actions/workflows/radia-fast.yml/runs?head_sha={'b' * 40}&event=push&per_page=100"
    api.records[key]["workflow_runs"][0].update(status=status, conclusion=conclusion)
    with pytest.raises(ValueError, match="Acceptance CI not successful"):
        gate.verify_run(api, 9, "b" * 40)


def test_old_success_cannot_hide_a_failed_ci_rerun():
    api = FakeAPI()
    key = f"actions/workflows/radia-fast.yml/runs?head_sha={'b' * 40}&event=push&per_page=100"
    failed = dict(api.records[key]["workflow_runs"][0], run_attempt=2, conclusion="failure")
    api.records[key]["workflow_runs"].append(failed)
    with pytest.raises(ValueError, match="Acceptance CI not successful"):
        gate.verify_run(api, 9, "b" * 40)


def application_proof():
    return {"schema": "radia.release-quad.simulink-candidate.v1", "version": "5.1.0",
            "commit": "a"*40, "package_sha256": "b"*64, "mex_sha256": "c"*64,
            "targets": {host: {"label": host, "status": "passed", "mex_sha256": "c"*64,
                              "verified_at_utc": "2026-09-30T00:00:00Z"}
                        for host in gate.RELEASE_ACCEPTANCE_HOSTS}}


def test_application_acceptance_requires_four_exact_mex_results():
    identity = {"source_commit": "a"*40, "version": "5.1.0"}
    gate.verify_application_acceptance(application_proof(), identity)
    for mutation in ("missing_host", "failed", "wrong_mex", "wrong_source", "missing_time"):
        proof = application_proof()
        if mutation == "missing_host": del proof["targets"]["mdx2"]
        elif mutation == "failed": proof["targets"]["100"]["status"] = "failed"
        elif mutation == "wrong_mex": proof["targets"]["lab"]["mex_sha256"] = "d"*64
        elif mutation == "wrong_source": proof["commit"] = "e"*40
        else: del proof["targets"]["mdx1"]["verified_at_utc"]
        with pytest.raises(ValueError):
            gate.verify_application_acceptance(proof, identity)
