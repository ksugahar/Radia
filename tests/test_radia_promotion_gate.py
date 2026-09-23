"""Publication gates reject incomplete, substituted and nonfinite evidence."""
import copy
import importlib.util
import hashlib
import json
from pathlib import Path
import zipfile
import xml.etree.ElementTree as ET

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("promotion_gate", ROOT / "tools/verify_radia_promotion.py")
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)
EVIDENCE = ROOT / "validation_test/esrf_three_engine/results/candidate_ed58b5deb"


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
