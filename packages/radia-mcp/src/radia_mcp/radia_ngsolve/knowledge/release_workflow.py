"""
Triple-package release workflow for the Radia monorepo
(radia / cubit-mesh-export / radia-mcp).

The Radia monorepo ships three independent PyPI packages from one
git tree.  Releases are atomic across the three so users on 100号機
and mdx never see a mismatched (radia X, radia-mcp Y, cubit-mesh-
export Z) combination that would silently break the panel <-> MCP
<-> Cubit pipeline.

This module is the canonical AI-readable record of the workflow,
the gates that must pass at each phase, the failure modes that have
historically tripped releases, and the recovery patterns.

Read this when:
* The user asks for a release / version bump / PyPI publish.
* CI on a tag ref fails after `git push --tags`.
* You see "Release / Release radia-mcp / Release cubit-mesh-export"
  workflows skip in `gh run list` (means CI never went green).
* The user reports "is X.Y.Z on PyPI yet?" and propagation seems
  stuck.

The MCP server exposes this via release_workflow(topic=...). Topics:
overview, phases, preflight_gates, ci_failure_modes, recovery,
patch_bump_protocol, lab_lock_release, monorepo_lockstep,
ci_monitor_skill.
"""

RELEASE_WORKFLOW = """\
# Triple-package release workflow (radia + radia-mcp + cubit-mesh-export)

This document is the AI-readable canonical reference for the Radia
release flow.  It is the persistent equivalent of the Claude Code
``release-triple`` skill (`.claude/skills/release-triple/SKILL.md`).
Topics: overview, phases, preflight_gates, ci_failure_modes,
recovery, patch_bump_protocol, lab_lock_release, monorepo_lockstep,
ci_monitor_skill.

## ===
## overview — what gets released and why atomically
## ===

The Radia monorepo (`ksugahar/Radia`) ships three independent
packages to PyPI:

| Package           | Tag prefix              | What it ships                                     |
|-------------------|-------------------------|---------------------------------------------------|
| radia             | `v`                     | C++ core (.pyd) + Python (panels, MCP, BEM, PEEC) |
| cubit-mesh-export | `cubit-mesh-export-v`   | Cubit plugin .ccm/.pyd + check-vol CLI            |
| radia-mcp         | `radia-mcp-v`           | MCP servers (radia-ngsolve, cubit, build123d,     |
|                   |                         | gmsh, electromagnet, ih, peec, ...)               |

All three are released in lock-step from one composite git commit
because:
* radia-mcp imports radia at runtime for several tools, and a
  schema mismatch is silent until a tool crashes.
* cubit-mesh-export ships the Cubit C++ plugin binaries that
  radia's panels expect via APREPRO commands.  A 0.7.4 plugin with
  a 4.27.x panel that emits a flag the plugin doesn't know is a
  silent mismatch that only surfaces at user clicktime.
* End-to-end production target machines (100号機, mdx) install all
  three together with `pip install 'radia[cubit,gui]==X.Y.Z'`,
  which pulls a compatible cubit-mesh-export and a compatible
  radia-mcp.  Users never deal with version arithmetic.

## ===
## phases — the 9-phase pipeline
## ===

Phases are the INTERNAL CONTRACT of `tools/release_triple.py`.  The
skill's narrative description maps 1:1 to script subcommands; this
table is the AI-readable summary.

| Phase | Action | Mandatory? | Notes |
|------:|--------|-----------|-------|
| 0 | Clean rebuild of Cubit plugin (.ccm/.pyd) | If `src/cubit_plugin/` changed | ~2-5 min targeted, ~10 min full |
| 1 | Decide minor vs patch per package | always | git log per package since last tag |
| 2 | Bump 4 version files (radia: pyproject + __init__; radia-mcp: pyproject + __init__; cubit-mesh-export: pyproject) | always | strictly lock-step; mismatch = wheel install bug |
| **2.5** | **Pre-flight CI validation (4 gates) — ADDED 2026-05-03** | **always** | local equivalent of CI; saves 2-3 round-trips |
| 3 | Stage exactly the release files | always | NO `git add -A`; user has WIP |
| 4 | Composite commit (HEREDOC, all packages in title) | always | Co-Authored-By trailer required |
| 5 | Three (or two) annotated tags | always | only bump packages with changes |
| 6 | Push main + all tags | always | tag push triggers CI; CI workflow_run gates Release workflows |
| 7 | Monitor CI propagation to PyPI | always | use ci-monitor skill |
| 8 | Deploy 100号機 + mdx (PyPI install + cubit-plugin-install) | always | shared lab box must not be left to users |
| 9 | Cross-machine consistency probe (LAB / 100号機 / mdx hashes) | always | DRIFT row = a phase 2/3/CI bug |

## ===
## preflight_gates — Phase 2.5 4-gate pre-push validation (2026-05-03)
## ===

Added 2026-05-03 after v4.27.0 / v4.27.1 BOTH round-tripped through
tag CI failures, burning 3 patch numbers for one set of features.
Every gate below was a CI failure on at least one of those round-
trips; running the gates locally takes ~1 minute total.

### Gate 1: regenerate radia-mcp TOOLS.md and verify

```
python packages/radia-mcp/scripts/gen_tools_doc.py
git diff --stat packages/radia-mcp/docs/TOOLS.md
python -m pytest tests/mcp_server/test_tools_doc.py -xvs
```

If the regen produced "(import failed: ModuleNotFoundError)" lines,
edit `packages/radia-mcp/scripts/gen_tools_doc.py` SERVERS list to
remove the dead subpackage entry, regen, and verify again.  This was
the v4.27.1 -> v4.27.2 trip.

### Gate 2: collect-only sweep with the CI exclusion set

```
python -m pytest \\
  --ignore=tests/cubit \\
  --ignore=tests/panels \\
  --ignore=tests/test_far_field_accuracy.py \\
  --ignore=tests/test_rad_ngsolve_function.py \\
  --ignore=tests/test_tetrahedral_solver.py \\
  --ignore=tests/test_batch_evaluation.py \\
  --ignore=tests/test_curlA_equals_B.py \\
  --ignore=tests/test_mesh_import.py \\
  --ignore=tests/test_scalar_bie_sibc.py \\
  --collect-only -q
```

Expected: ``N tests collected``, **NO** ``error``, **NO** ``Interrupted``.
Any import error = a stale `from <removed_module>` somewhere in
tests/.  The v4.27.0 -> v4.27.1 trip was a stale
`tests/mcp_server/test_build123d_gmsh_elf_responses.py` importing
`radia_mcp.elf` after the elf subpackage was extracted.

### Gate 3: actually run the affected mcp_server tests

```
python -m pytest tests/mcp_server/ -x
```

Catches schema regressions in the MCP tool table that the doc
generator doesn't.

### Gate 4: workflow-runner --selftest matrix import check

If `--selftest` matrix in `.github/workflows/build-test.yml` was
edited (subpackage removed/added):

```
python -c "
import importlib
servers = ['radia_mcp.radia_ngsolve.server',
           'radia_mcp.cubit.server',
           'radia_mcp.build123d.server',
           'radia_mcp.gmsh.server']
for s in servers:
    try: importlib.import_module(s); print(f'OK   {s}')
    except Exception as e: print(f'FAIL {s}: {e}')"
```

Every server must print OK.  Subpackages no longer in the codebase
must NOT be in the .yml matrix.

**Only proceed to Phase 3 once all four gates above are green.**

## ===
## ci_failure_modes — known historical CI failures + cause + fix
## ===

| Symptom (CI log line)                                                   | Root cause                                                       | Local gate that catches it     | First seen   |
|-------------------------------------------------------------------------|------------------------------------------------------------------|--------------------------------|--------------|
| `ModuleNotFoundError: No module named 'radia_mcp.elf'`                  | stale test importing extracted subpackage                        | Gate 2 (collect-only sweep)    | v4.27.0      |
| `AssertionError: ... TOOLS.md is stale`                                 | TOOLS.md generator hardcodes a subpackage that no longer exists  | Gate 1 (regen + diff)          | v4.27.1      |
| `FAIL: examples/<dir>/ missing README.md`                               | Policy Lint workflow; new examples/ subdir without README        | manual: `ls examples/*/README.md` | (multiple)|
| `--selftest: ModuleNotFoundError: No module named 'radia_mcp.<sub>'`    | .yml matrix lists removed subpackage                             | Gate 4 (--selftest import)     | (multiple)   |
| `gh release download ... no assets match the file pattern`              | binaries-release upload race vs tag-CI start                     | (CI has 6-attempt retry as of v4.26)| v4.25.1, v4.27.0|
| `Basic tests failed` w/ collected 0 tests                               | pytest config conflict (pyproject.toml vs pytest.ini)            | Gate 2 (collect-only sweep)    | (rare)       |
| `_radia_pybind import failure: DLL load failed`                         | Cubit plugin .pyd / cubit-mesh-export .pyd built against wrong Python ABI | Phase 0 clean rebuild | release-triple Phase 0 |

## ===
## recovery — when CI on a tag fails AFTER push
## ===

CI on a tag ref running on a broken commit cannot be rescued by
pushing a fix to main.  The `Release` workflow is gated on
`workflow_run` of CI ON THE TAG REF, and the tag still points at
the broken commit.

**The only path forward is patch bump and re-tag**, per the skill's
"PyPI is immutable" policy.  Do NOT delete + recreate a pushed tag
unless the user explicitly authorizes it (and even then, only if
the version has not gone to PyPI -- which it hasn't, since CI
failed; PyPI never sees it).

The 2026-05-03 v4.27.x round-trip went:
  v4.27.0 (CI fail: stale elf imports)
    -> fix-forward commit on main (deletes test, cleans __init__,
       fixes .yml matrix)
    -> v4.27.1 (CI fail: TOOLS.md stale)
    -> fix-forward (regen + edit gen_tools_doc.py SERVERS list)
    -> v4.27.2 (CI green, PyPI propagated)

The Phase 2.5 pre-flight gates were added after this episode to
prevent the round-trip pattern.  Total cost was 3 patch numbers
and ~30 min of CI time per cycle that the gates would have spent
~1 min locally.

## ===
## patch_bump_protocol — exact steps for retry after CI failure
## ===

```
# 1. Fix the issue on main first (commit + push). Verify locally.
# 2. Bump 4 version files +1 in patch position (X.Y.Z -> X.Y.Z+1).
# 3. Compose a release commit on top with the bumps:

git add pyproject.toml src/radia/__init__.py \\
        packages/radia-mcp/pyproject.toml \\
        packages/radia-mcp/src/radia_mcp/__init__.py
git commit -m "Release vX.Y.Z+1 / radia-mcp-v0.M.N+1
<paragraph: '<previous bump> failed at <gate>; this is a re-tag>'>"

# 4. Tag at the new commit:
git tag -a vX.Y.Z+1           -m "..."
git tag -a radia-mcp-v0.M.N+1 -m "..."

# 5. Push main (already up-to-date) + new tags:
git push origin main
git push origin vX.Y.Z+1 radia-mcp-v0.M.N+1

# 6. Monitor with ci-monitor.
```

Old broken tags (vX.Y.Z, radia-mcp-v0.M.N) STAY in the repo as a
record of the failed attempt.  PyPI never saw them so there's no
immutability conflict.

## ===
## lab_lock_release — pre-release: stop processes that hold .pyd / .ccm files
## ===

Before any deploy that touches the Cubit plugin .ccm/.pyd or
the radia-mcp Scripts/mcp-server-*.exe, every machine in the deploy
target must release file locks:

```
pwsh -Command "
Get-Process -ErrorAction SilentlyContinue | Where-Object {
  \\$_.Name -like 'mcp-server*' -or \\$_.ProcessName -eq 'coreform_cubit'
} | ForEach-Object { Stop-Process -Id \\$_.Id -Force }
Start-Sleep -Seconds 2"
```

For 100号機 / mdx, prefix with `cat << 'PS' | ssh ... 'pwsh ...'`.
Do not wait for human users to close Cubit; on shared lab machines
the deploy responsibility is Claude's per the lab policy
(feedback_deploy_responsibility.md), and Cubit is re-launchable.

## ===
## monorepo_lockstep — version files that MUST stay in sync
## ===

| File | Field | Read by |
|------|-------|---------|
| `pyproject.toml` | `version = "X.Y.Z"` | setuptools, pip, PyPI metadata |
| `src/radia/__init__.py` | `__version__ = "X.Y.Z"` | runtime introspection, --version flags, tests |
| `packages/radia-mcp/pyproject.toml` | `version = "X.Y.Z"` | radia-mcp wheel metadata |
| `packages/radia-mcp/src/radia_mcp/__init__.py` | `__version__ = "X.Y.Z"` | runtime, --selftest --version |
| `packages/cubit-mesh-export/pyproject.toml` | `version = "X.Y.Z"` | cubit-mesh-export wheel metadata |
| `packages/cubit-mesh-export/src/cubit_mesh_export/__init__.py` | `__version__ = "X.Y.Z"` | runtime, COMPAT constants |

A mismatch between `pyproject.toml` and the matching `__init__.py`
causes:
* Cryptic wheel-install bugs (pip caches the wrong version).
* `tests/test_*_version_consistency.py` failures on CI.
* Users seeing `radia.__version__` not match `pip show radia`.

The `tools/release_triple.py preflight` subcommand checks this
lockstep and exits non-zero if violated.

## ===
## ci_monitor_skill — companion skill for Phase 7
## ===

The `.claude/skills/ci-monitor/SKILL.md` skill (added 2026-05-03)
watches GitHub Actions runs to completion and on failure auto-
fetches `gh run view <id> --log-failed` so the AI agent has the
diagnosis context immediately.

```
python .claude/skills/ci-monitor/monitor.py 25275890624 25275890613
# or, auto-discover:
python .claude/skills/ci-monitor/monitor.py --auto 3
```

Exit 0 = all green, exit 1 = at least one failure (with log tail
already printed for diagnosis).

This is the standard companion to release-triple Phase 7 -- after
`git push --tags` the agent should immediately call ci-monitor with
the 3 (or 2) new run IDs and not declare the release "done" until
the monitor exits 0.
"""


_TOPICS = (
    "overview",
    "phases",
    "preflight_gates",
    "ci_failure_modes",
    "recovery",
    "patch_bump_protocol",
    "lab_lock_release",
    "monorepo_lockstep",
    "ci_monitor_skill",
)


def get_release_workflow_documentation(topic: str = "") -> str:
    """Return the triple-package release workflow knowledge.

    Args:
        topic: empty for the full document, or one of the entries in
            ``_TOPICS`` for a single section.
    """
    if not topic:
        return RELEASE_WORKFLOW

    if topic not in _TOPICS:
        return (f"Unknown topic: {topic!r}. Available topics:\n"
                f"  {', '.join(_TOPICS)}\n\n"
                f"Pass empty string for the full document.")

    headers = []
    pos = 0
    while True:
        next_pos = RELEASE_WORKFLOW.find("\n## ", pos)
        if next_pos < 0:
            break
        line_end = RELEASE_WORKFLOW.find("\n", next_pos + 1)
        line = RELEASE_WORKFLOW[next_pos + 1:line_end]
        for t in _TOPICS:
            if line.startswith(f"## {t} "):
                headers.append((t, next_pos + 1))
                break
        pos = next_pos + 1

    req_starts = [off for kw, off in headers if kw == topic]
    if not req_starts:
        return f"Topic {topic!r} declared but not found in document."
    section_start = RELEASE_WORKFLOW.rfind("## ===", 0, req_starts[0])
    if section_start < 0:
        section_start = req_starts[0]

    section_end = len(RELEASE_WORKFLOW)
    last_req = req_starts[-1]
    for kw, off in headers:
        if kw != topic and off > last_req:
            delim = RELEASE_WORKFLOW.rfind("## ===", 0, off)
            section_end = delim if delim > 0 else off
            break

    return RELEASE_WORKFLOW[section_start:section_end].rstrip() + "\n"
