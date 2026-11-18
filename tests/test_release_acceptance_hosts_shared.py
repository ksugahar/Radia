"""The gate that publishes and the tool that produces evidence name one host list.

Release-quad requires validation role and 100号機 for the same release commit.
The promotion verifier used to carry its own two-host list -- one of them a
computation host that is not an acceptance target -- so a release with
evidence from two of the four machines would have passed the public gate.
Both tools now import the list from `tools/release_acceptance.py`; these
tests pin that the list is the policy's, that both tools use it, and that
neither has grown a private copy again.
"""
from __future__ import annotations

import pytest

import importlib.util
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"


def _load(name):
    if str(TOOLS) not in sys.path:
        sys.path.insert(0, str(TOOLS))
    spec = importlib.util.spec_from_file_location(name, TOOLS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_the_shared_list_is_the_policy_list():
    shared = _load("release_acceptance")
    assert shared.RELEASE_ACCEPTANCE_HOSTS == ("validation", "100")
    assert "hibino" not in shared.RELEASE_ACCEPTANCE_HOSTS, (
        "hibino is a computation host, not a release acceptance target")



def test_release_quad_produces_evidence_for_exactly_those_hosts():
    """The producer's target table is the shared tuple, key for key."""
    source = (TOOLS / "release_quad.py").read_text(encoding="utf-8")
    assert "from release_acceptance import RELEASE_ACCEPTANCE_HOSTS" in source
    assert "tuple(SOLVER_TARGETS) == RELEASE_ACCEPTANCE_HOSTS" in source
    keys = re.findall(r'^\s{4}"([a-z0-9]+)": \("', source, re.M)
    assert tuple(_load("release_quad").SOLVER_TARGETS) == ("validation", "100")


def test_the_promotion_gate_has_no_private_host_list():
    """No literal host tuple in the verifier: it must import the shared one."""
    source = (TOOLS / "verify_radia_promotion.py").read_text(encoding="utf-8")
    assert "from release_acceptance import RELEASE_ACCEPTANCE_HOSTS" in source
    assert "for host in RELEASE_ACCEPTANCE_HOSTS" in source
    # And it reports every missing host rather than stopping at the first.
    assert "Acceptance evidence missing for" in source


def test_retired_host_has_no_release_route():
    quad = _load("release_quad")
    assert not hasattr(quad, "SSH_MDX1")
    assert quad._deploy_mdx("md" + "x") == 2
    assert quad.cmd_phase8e(None) == 2
    assert "mdx2" not in quad.SIMULINK_TARGETS
    assert not hasattr(quad, "SSH_MDX2")


def test_release_all_is_read_only_candidate_verification(monkeypatch):
    quad = _load("release_quad")
    args = object()
    calls = []

    def forbidden(*args):
        pytest.fail("candidate handoff attempted a deployment")

    monkeypatch.setattr(quad, "cmd_phase8", forbidden)
    monkeypatch.setattr(quad, "cmd_phase9", forbidden)
    monkeypatch.setattr(quad, "cmd_done", lambda received: calls.append(received) or 4)
    assert quad.cmd_all(args) == 4
    assert calls == [args]


def test_done_has_no_compute_refresh_or_remote_cleanup_dependency():
    import inspect
    quad = _load("release_quad")
    source = inspect.getsource(quad.cmd_done)
    assert "cmd_phase8e(" not in source
    assert "cmd_temp_shadows(" not in source



@pytest.mark.parametrize("defect", [None, "source", "role", "artifact", "acceptance", "ci"])
def test_candidate_handoff_binds_all_four_receipts_to_exact_artifacts(tmp_path, monkeypatch, defect):
    from argparse import Namespace
    import hashlib
    import json

    quad = _load("release_quad")
    source = "a" * 40

    def evidence(name):
        path = tmp_path / name
        path.write_bytes(name.encode("ascii"))
        return {"path": name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}

    receipts = {}
    for name, (build_role, validation_role) in quad.RELEASE_LANES.items():
        artifact = evidence(name + ".whl")
        receipts[name] = {
            "build": {"role": build_role, "source_sha": source, "evidence": evidence(name + ".build"),
                      "artifacts": [artifact]},
            "validation": {"role": validation_role, "source_sha": source, "passed": True,
                           "evidence": evidence(name + ".validation"),
                           "artifact_sha256": {artifact["path"]: artifact["sha256"]}},
        }
    data = {"schema": 1, "source_sha": source, "lanes": receipts,
            "native_ci": {"source_sha": source, "passed": True, "evidence": evidence("ci.json")}}
    if defect == "source":
        receipts["matlab"]["validation"]["source_sha"] = "b" * 40
    elif defect == "role":
        receipts["windows-python-gpu"]["build"]["role"] = "cpu-build"
    elif defect == "artifact":
        (tmp_path / "linux-python-gpu.whl").write_bytes(b"changed")
    elif defect == "acceptance":
        del receipts["cpu"]
    elif defect == "ci":
        data["native_ci"]["passed"] = False

    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps(data), encoding="utf-8")
    monkeypatch.setattr(quad, "_verify_local_release_source",
                        lambda path, sha: 0 if path == str(tmp_path.resolve()) and sha == source else 4)
    args = Namespace(source_sha=source, manifest=manifest, release_source=tmp_path)
    assert quad.cmd_done(args) == (0 if defect is None else 4)
