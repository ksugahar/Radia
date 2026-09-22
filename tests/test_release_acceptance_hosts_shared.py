"""The gate that publishes and the tool that produces evidence name one host list.

Release-quad requires LAB, 100号機, mdx1 and mdx2 for the same release commit.
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
    assert shared.RELEASE_ACCEPTANCE_HOSTS == ("lab", "100", "mdx1", "mdx2")
    assert "hibino" not in shared.RELEASE_ACCEPTANCE_HOSTS, (
        "hibino is a computation host, not a release acceptance target")


def test_the_policy_text_names_the_same_four_machines():
    policy = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    assert re.search(r"Release-quad requires LAB, 100号機, mdx1, and mdx2", policy)
    assert "hibino remains a computation host" in policy


def test_release_quad_produces_evidence_for_exactly_those_hosts():
    """The producer's target table is the shared tuple, key for key."""
    source = (TOOLS / "release_quad.py").read_text(encoding="utf-8")
    assert "from release_acceptance import RELEASE_ACCEPTANCE_HOSTS" in source
    assert "tuple(SIMULINK_TARGETS) == RELEASE_ACCEPTANCE_HOSTS" in source
    keys = re.findall(r'^\s{4}"([a-z0-9]+)": \("', source, re.M)
    assert tuple(keys) == ("lab", "100", "mdx1", "mdx2")


def test_the_promotion_gate_has_no_private_host_list():
    """No literal host tuple in the verifier: it must import the shared one."""
    source = (TOOLS / "verify_radia_promotion.py").read_text(encoding="utf-8")
    assert "from release_acceptance import RELEASE_ACCEPTANCE_HOSTS" in source
    assert '("lab", "hibino")' not in source
    assert '["lab", "hibino"]' not in source
    assert "for host in RELEASE_ACCEPTANCE_HOSTS" in source
    # And it reports every missing host rather than stopping at the first.
    assert "Acceptance evidence missing for" in source
