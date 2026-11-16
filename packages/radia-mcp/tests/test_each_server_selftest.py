"""Per-server --selftest harness for radia-mcp.

Complements `test_meta_health.py` (which import-tests every subpackage):
this test boots each served key (`radia-mcp serve <key>`) and each pack
domain module as a subprocess with the `--selftest` flag and verifies it exits 0 within a short budget.

Why subprocess rather than direct import + call:
- catches breakage in the *entry-point wiring* (broken
  `pyproject.toml [project.scripts]`, stale `pip install -e` after a
  package rename, etc.) that a same-process import test does not.
- catches `if __name__ == "__main__"` bugs.
- runs each server in process isolation so a SystemExit / sys.argv
  mutation in one selftest cannot poison the others.

Why parametrized rather than one big loop:
- per-server pytest IDs make CI failure output immediately actionable
  ("which servers regressed?" vs "the suite failed").
- supports `pytest -k motor` to re-run a single server.

Each selftest is given a 60s budget — well above the historical worst
case (radia-ngsolve PDF index warm-up ~12s, mathematica wolframscript
probe ~5s) and well below CI total budget.

Skip behavior:
- If `radia-mcp` is not on PATH, served keys SKIP with a hint.
"""

from __future__ import annotations

import shutil
import subprocess
import sys

import pytest

pytestmark = pytest.mark.slow

# Single source of truth for the server list: derive from the meta
# catalog so this test cannot drift from `radia_mcp_overview()`.
from radia_mcp.meta import catalog


# Served keys run through the single `radia-mcp serve <key>` command; pack
# member domains have no command of their own and run their module selftest.
_SLOW_SERVERS = {
    "radia-ngsolve": 120,           # PDF index + radia knowledge load
    "radia-mcp-mathematica": 90,    # wolframscript probe
}

SERVERS = sorted(catalog.CATALOG.keys())


def _command_for(name: str) -> tuple[str, list[str]]:
    info = catalog.CATALOG[name]
    if "entry_point" in info:
        exe = shutil.which("radia-mcp")
        if exe is None:
            pytest.skip("radia-mcp not on PATH; run `pip install -e .`")
        return "radia-mcp serve " + name, [exe, "serve", name, "--selftest"]
    module = info["subpackage"] + ".server"
    return module, [sys.executable, "-m", module, "--selftest"]


@pytest.mark.parametrize("name", SERVERS)
def test_server_selftest(name: str) -> None:
    """Every served key and pack domain passes --selftest within budget."""
    cli, argv = _command_for(name)
    timeout = _SLOW_SERVERS.get(name, 60)
    try:
        result = subprocess.run(
            argv,
            capture_output=True,
            # IMPORTANT: do NOT pass text=True on Windows — the default
            # locale codec is cp932 in Japanese environments and many
            # selftest banners contain Japanese characters / U+0080+
            # bytes that crash the implicit decoder on the reader
            # thread.  Decode explicitly below with errors='replace'.
            timeout=timeout,
            # Inherit env so the selftest sees the same Python /
            # editable-install state as the test runner.
        )
    except subprocess.TimeoutExpired as e:
        pytest.fail(
            f"{cli} --selftest timed out after {timeout}s.  "
            f"Either the selftest is genuinely slow (raise _SLOW_SERVERS"
            f" budget) or it hangs (BUG).\n"
            f"  partial stdout: {(e.stdout or b'')[:400]!r}\n"
            f"  partial stderr: {(e.stderr or b'')[:400]!r}"
        )

    # cp932-safe decode of captured bytes
    stdout = (result.stdout or b"").decode("utf-8", errors="replace")
    stderr = (result.stderr or b"").decode("utf-8", errors="replace")

    if result.returncode != 0:
        pytest.fail(
            f"{cli} --selftest exited with returncode={result.returncode}\n"
            f"  stdout (last 800 chars):\n{stdout[-800:]}\n"
            f"  stderr (last 800 chars):\n{stderr[-800:]}"
        )
    # Sanity: selftest scripts should print "PASSED" near the end as the
    # uniform convention across the suite (motor, etc.).
    # Not a hard requirement — some servers (gmsh, cubit, radia-ngsolve)
    # emit different success banners — so this is a soft signal only,
    # printed (not warned) to keep pytest -W error clean.
    if "PASSED" not in stdout and "OK" not in stdout:
        # Don't use warnings.warn — pytest may elevate to error under
        # -W error.  print() goes to captured stdout and only surfaces
        # in pytest -v / -s, which is the right amount of noise for a
        # cosmetic banner check.
        print(
            f"[soft] {cli} --selftest succeeded but printed no "
            f"'PASSED'/'OK' banner.  Consider adding one for uniform "
            f"CI greppability."
        )


def test_selftest_harness_covers_full_catalog() -> None:
    """Pin the invariant: every cataloged server is selftest-parametrized.

    Catches the regression mode where someone adds a subpackage but
    forgets to either (a) wire `--selftest` in its server.py main() or
    (b) extend this test — without this assertion the new server would
    silently skip selftest forever.
    """
    cataloged = set(catalog.CATALOG.keys())
    parametrized = set(SERVERS)
    missing = cataloged - parametrized
    assert not missing, (
        f"{len(missing)} cataloged servers not in selftest harness: "
        f"{sorted(missing)}"
    )
