"""TRIPWIRE for the HDiv-VIM tet/hex COUPLING gate: NGSolve H(div)-PYRAMID support.

A conforming mixed tet+hex mesh needs pyramid transition elements carrying an H(div) flux; that needs
NGSolve to implement HDiv on pyramids (Joachim/Schoberl committed to add it -- the lab WAITS, does not
reimplement).  As of NGSolve 6.2.2604 HDiv-pyramid is alloc-but-unimplemented (the first Assemble raises
"HDivHighOrderFESpace: Pyramid elements not implemented yet!").

This test asserts the DESIRED end-state (HDiv-pyramid is functional) under xfail(strict=True): it is an
EXPECTED FAILURE today (suite stays green), and it XPASSes -> STRICT-FAIL (suite goes RED, loudly) the day
an NGSolve bump implements it.  A red here is GOOD NEWS: run the `ngsolve-hdiv-pyramid-check` skill and
follow its playbook to add the Radia pyramid charge-Gram mode (mirror the wedge port, commit af5ab64d) and
enable mixed meshes.  Detection logic lives in tools/probe_hdiv_pyramid.py (single source of truth).
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")

_TOOLS = os.path.join(os.path.dirname(__file__), "..", "..", "tools")
sys.path.insert(0, os.path.abspath(_TOOLS))
from probe_hdiv_pyramid import probe   # noqa: E402


@pytest.mark.xfail(reason="NGSolve HDiv-pyramid not implemented yet (6.2.2604 raises 'Pyramid elements not "
                          "implemented yet'); when this XPASSes the block has LIFTED -- run the "
                          "ngsolve-hdiv-pyramid-check skill and unblock HDiv-VIM mixed meshes.",
                   strict=True)
def test_ngsolve_hdiv_pyramid_is_implemented():
    r = probe()
    # ERROR (probe broke) must NOT masquerade as the block still being present -- surface it loudly.
    assert r["verdict"] != "ERROR", f"pyramid probe ERROR (fix the probe): {r['detail']}"
    assert r["verdict"] == "IMPLEMENTED", \
        f"NGSolve {r['ngsolve_version']} HDiv-pyramid verdict={r['verdict']}: {r['detail']}"
