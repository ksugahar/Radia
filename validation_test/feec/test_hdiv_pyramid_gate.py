"""TRIPWIRE for the HDiv-VIM tet/hex COUPLING gate: NGSolve H(div)-PYRAMID support.

A conforming mixed tet+hex mesh needs pyramid transition elements carrying an H(div) flux; that needs
NGSolve to implement HDiv on pyramids (Joachim/Schoberl committed to add it -- the lab WAITS, does not
reimplement).  As of NGSolve 6.2.2604 HDiv-pyramid is alloc-but-unimplemented (the first Assemble raises
"HDivHighOrderFESpace: Pyramid elements not implemented yet!").

This test is an executable tripwire, not a normal desired-state assertion:
NOT_IMPLEMENTED is an expected xfail today, while IMPLEMENTED / ALLOC_BUT_BROKEN / ERROR go RED loudly.
A red IMPLEMENTED here is GOOD NEWS: run the `ngsolve-hdiv-pyramid-check` skill and follow its playbook to
add the Radia pyramid charge-Gram mode (mirror the wedge port, commit af5ab64d) and enable mixed meshes.
Detection logic lives in tools/probe_hdiv_pyramid.py (single source of truth).
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")

_TOOLS = os.path.join(os.path.dirname(__file__), "..", "..", "tools")
sys.path.insert(0, os.path.abspath(_TOOLS))
from probe_hdiv_pyramid import probe   # noqa: E402


def test_ngsolve_hdiv_pyramid_tripwire():
    r = probe()
    # ERROR (probe broke) must NOT masquerade as the block still being present -- surface it loudly.
    if r["verdict"] == "ERROR":
        pytest.fail(f"pyramid probe ERROR (fix the probe): {r['detail']}")
    if r["verdict"] == "NOT_IMPLEMENTED":
        pytest.xfail("NGSolve HDiv-pyramid not implemented yet; HDiv-VIM mixed tet/hex remains blocked.")
    pytest.fail(
        f"NGSolve {r['ngsolve_version']} HDiv-pyramid verdict={r['verdict']}: {r['detail']} -- "
        "run the ngsolve-hdiv-pyramid-check playbook before treating this suite as green.")
