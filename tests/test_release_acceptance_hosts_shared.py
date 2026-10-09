"""The gate that publishes and the tool that produces evidence name one host list.

Release-quad requires LAB and 100号機 for the same release commit.
The promotion verifier used to carry its own two-host list -- one of them a
computation host that is not an acceptance target -- so a release with
evidence from two of the four machines would have passed the public gate.
Both tools now import the list from `tools/release_acceptance.py`; these
tests pin that the list is the policy's, that both tools use it, and that
neither has grown a private copy again.
"""
from __future__ import annotations

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
    assert shared.RELEASE_ACCEPTANCE_HOSTS == ("lab", "100")
    assert "hibino" not in shared.RELEASE_ACCEPTANCE_HOSTS, (
        "hibino is a computation host, not a release acceptance target")


def test_the_policy_text_names_the_same_two_machines():
    policy = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    assert re.search(r"Release-quad requires LAB and 100号機", policy)
    assert "hibino remains a computation host" in policy


def test_release_quad_produces_evidence_for_exactly_those_hosts():
    """The producer's target table is the shared tuple, key for key."""
    source = (TOOLS / "release_quad.py").read_text(encoding="utf-8")
    assert "from release_acceptance import RELEASE_ACCEPTANCE_HOSTS" in source
    assert "tuple(SIMULINK_TARGETS) == RELEASE_ACCEPTANCE_HOSTS" in source
    keys = re.findall(r'^\s{4}"([a-z0-9]+)": \("', source, re.M)
    assert tuple(_load("release_quad").SIMULINK_TARGETS) == ("lab", "100")


def test_the_promotion_gate_has_no_private_host_list():
    """No literal host tuple in the verifier: it must import the shared one."""
    source = (TOOLS / "verify_radia_promotion.py").read_text(encoding="utf-8")
    assert "from release_acceptance import RELEASE_ACCEPTANCE_HOSTS" in source
    assert '("lab", "hibino")' not in source
    assert '["lab", "hibino"]' not in source
    assert "for host in RELEASE_ACCEPTANCE_HOSTS" in source
    # And it reports every missing host rather than stopping at the first.
    assert "Acceptance evidence missing for" in source


def test_retired_host_has_no_release_route():
    quad = _load("release_quad")
    assert quad.SSH_MDX1 == "mdx"
    assert "mdx2" not in quad.SIMULINK_TARGETS
    assert not hasattr(quad, "SSH_MDX2")


def test_release_all_never_refreshes_compute_hosts(monkeypatch):
    quad = _load("release_quad")
    calls = []
    def forbidden(*args):
        raise AssertionError("release must not contact compute hosts")
    monkeypatch.setattr(quad, "cmd_phase8e", forbidden)
    monkeypatch.setattr(quad, "cmd_temp_shadows", forbidden)
    monkeypatch.setattr(quad, "cmd_phase8", lambda args: calls.append(args.target) or 0)
    monkeypatch.setattr(quad, "cmd_phase9", lambda args: calls.append("phase9") or 0)
    assert quad.cmd_all(None) == 0
    assert calls == ["lab,100", "phase9"]


def test_optional_refresh_uses_only_the_mdx_alias(monkeypatch):
    quad = _load("release_quad")
    calls = []
    monkeypatch.setattr(quad, "_deploy_pypi", lambda host, label: calls.append((host, label)) or 0)
    assert quad.cmd_phase8e(None) == 0
    assert calls == [("mdx", "mdx1")]


def test_done_has_no_compute_refresh_or_remote_cleanup_dependency():
    import inspect
    quad = _load("release_quad")
    source = inspect.getsource(quad.cmd_done)
    assert "cmd_phase8e(" not in source
    assert "cmd_temp_shadows(" not in source
