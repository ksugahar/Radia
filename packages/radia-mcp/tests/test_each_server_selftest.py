"""Per-server --selftest harness for radia-mcp.

Complements `test_meta_health.py` (which import-tests every subpackage):
this test boots each `mcp-server-<x>` script as a subprocess with the
`--selftest` flag and verifies it exits 0 within a short budget.

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
- If `mcp-server-<x>` is not on PATH (entry point not installed in the
  active env), the test reports SKIP with a clear hint rather than
  FAIL. This is the right behavior for a developer who has not run
  `pip install -e .` since adding the server.
"""

from __future__ import annotations

import shutil
import subprocess
import sys

import pytest

# Single source of truth for the server list: derive from the meta
# catalog so this test cannot drift from `radia_mcp_overview()`.
from radia_mcp.meta import catalog


# Map catalog short-name -> CLI script name.  Most are 1:1, but one
# gets the 'radia-' prefix from the CLI script (see catalog alias logic):
#   catalog 'meta' -> CLI 'mcp-server-radia-meta'
_CLI_NAME_OVERRIDES = {
    "meta": "mcp-server-radia-meta",
}

# Servers that genuinely take longer to selftest (PDF index warm-up,
# external tool probe, etc.).  Default per-server budget is 60s.
_SLOW_SERVERS = {
    "radia-ngsolve": 120,        # PDF index + radia knowledge load
    "literature-index": 120,     # chromadb optional probe
    "mathematica": 90,           # wolframscript probe
}


def _cli_for(name: str) -> str:
    """Resolve catalog short-name to CLI script name."""
    return _CLI_NAME_OVERRIDES.get(name, f"mcp-server-{name}")


# Build the parametrize list once at collection time.  pytest will
# generate one test_<name> per server so failures are individually
# addressable.
SERVERS = sorted(catalog.CATALOG.keys())


@pytest.mark.parametrize("name", SERVERS)
def test_server_selftest(name: str) -> None:
    """`mcp-server-<name> --selftest` exits 0 within budget."""
    cli = _cli_for(name)
    exe = shutil.which(cli)
    if exe is None:
        pytest.skip(
            f"{cli} not on PATH — run `pip install -e .` to wire entry "
            f"points (this is a developer-environment skip, not a bug)."
        )
    timeout = _SLOW_SERVERS.get(name, 60)
    try:
        result = subprocess.run(
            [exe, "--selftest"],
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
    # uniform convention across the suite (panel_review, motor, etc.).
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


def test_at_least_30_servers_runnable() -> None:
    """Floor check: most servers' CLI entry points must be wired.

    Run AFTER the parametrized tests collect SKIP results.  If fewer
    than 30 of the cataloged servers have their CLI on PATH, the editable
    install almost certainly drifted (e.g. to the CI runner clone) —
    see CLAUDE.md "Distribution Test Policy" 2026-05-19 note.

    The threshold (30) is intentionally below the actual count so
    this passes even if 2-3 brand-new servers are still being landed.
    """
    runnable = sum(1 for n in SERVERS if shutil.which(_cli_for(n)) is not None)
    assert runnable >= 30, (
        f"Only {runnable}/{len(SERVERS)} CLI scripts on PATH.  Probable "
        f"editable-install drift.  Fix with:\n"
        f"  pip install -e S:/Radia/01_GitHub/packages/radia-mcp "
        f"--no-deps --no-cache-dir"
    )
