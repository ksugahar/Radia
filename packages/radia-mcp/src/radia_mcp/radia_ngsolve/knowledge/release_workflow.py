"""
Package release workflows for the Radia monorepo.

The Radia monorepo ships independently versioned solver, MCP and Simulink
distributions from one git tree. The numerical Radia solver waits for four
lane receipts of one SHA (cpu, linux-python-gpu, windows-python-gpu, matlab).
radia-mcp publishes independently and deploys a verified wheel to the
validation host and an editable source to the build host through its DUAL lane.
These release boundaries must not be collapsed into one coupled version or
deployment requirement.

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
overview, mcp_release, phases, simulink_candidate, optuna_candidate, preflight_gates, mcp_quality_review,
ci_failure_modes, recovery, patch_bump_protocol, lab_lock_release,
version_pairs, ci_monitor_skill.
"""

RELEASE_WORKFLOW = """\
# Package-scoped releases: Radia four-lane solver release and MCP release-dual

radia-mcp releases independently of Radia. Its release-dual deployment targets
are the validation host (verified wheel) and the build host (editable) only.
No retired host receives any deployment. Isolated CI is not deployment.
For MCP-only updates use update source -> reconnect -> check an affected live
tool, following packages/radia-mcp/docs/maintenance.md, not the solver release.
release-dual is the two-host scope, not a separate release_dual.py command.
The `mcp_release` topic is the authoritative MCP release contract. Other
sections describe the numerical solver, Simulink, or their named package lane;
they do not add MCP deployment targets or completion gates.

This document is the AI-readable canonical reference for the Radia
release flow.  Its canonical local orchestrator is `tools/release_quad.py`.
A Radia release waits for four lane receipts of one source SHA:

| Lane | Built on | Validated on |
|---|---|---|
| `cpu` | CPU-only build host (`cpu-build`) | Windows validation host |
| `linux-python-gpu` | supercomputer (`linux-gpu`), outside CI | supercomputer |
| `windows-python-gpu` | Windows GPU validation host (`windows-gpu`, runner `gpu-a100`) | same host |
| `matlab` | Windows GPU validation host (`windows-gpu`) | same host |

The build host is CPU-only: CPU builds, CPU CI (`100`, `ci-light`), early
debugging and CPU release-candidate preparation. It compiles no CUDA, installs
no NVIDIA/cuDSS packages and runs no GPU or MATLAB GPU work. All Windows GPU
builds, their CI and their release qualification run on the Windows GPU
validation host. `release_quad done --source-sha SHA --manifest FILE` verifies
the four receipts and their exact artifacts and writes the public lane summary
under validation_test/esrf_three_engine/results/candidate_<source-sha-first-9>/.
Dispatch the Release workflow on the exact tag with ci_run_id, wheel_sha256 and
the full acceptance_commit. Its read-only gate checks the GitHub
run/workflow/repository, tag context and peeled SHA, main reachability, the
four-lane summary, and wheel/native hashes; the published wheel is the
`windows-python-gpu` lane artifact. Publication uses the same artifact without
rebuilding. Changed bytes require new receipts. Linux wheels are not published.
Topics: overview, mcp_release, phases, simulink_candidate, optuna_candidate, preflight_gates, mcp_quality_review,
ci_failure_modes, recovery, patch_bump_protocol, lab_lock_release,
version_pairs, ci_monitor_skill.

## ===
## overview — package ownership and verification targets
## ===

The Radia monorepo (`ksugahar/Radia`) ships independent distributions
to PyPI. Each distribution follows its own declared verification targets:

| Package           | Tag prefix              | What it ships                                     |
|-------------------|-------------------------|---------------------------------------------------|
| radia             | `v`                     | C++ core (.pyd) + Python/pybind11 support, with   |
|                   |                         | the Python GPU payload from the Windows GPU lane  |
| radia-mcp         | `radia-mcp-v`           | MCP guidance and workflows, including mixed omega;|
|                   |                         | Cubit MCP is distributed by cubit-mesh-export     |
| radia-simulink | `radia-simulink-v` | Independently released mexcuda, MATLAB, LTspice and Simulink support |

cubit-mesh-export and radia-optuna release from their own repositories and
lanes. Packages are independently released. radia-mcp uses release-dual
(validation host wheel, build host editable). Compatibility is still checked
where an affected tool uses another package:
* radia-mcp imports radia at runtime for several tools, and a
  schema mismatch is silent until a tool crashes.
* cubit-mesh-export owns Cubit meshing, export, and its C++ plugin binaries.
  Radia consumes the checked `.vol` boundary; it does not couple solver release
  acceptance to the exporter version.
* The build host keeps the only, CPU-only editable development install in a
  dedicated venv. The validation host runs verified wheels only, never an
  editable install: `release_quad phase8` installs the published-lane wheel
  whose SHA-256 matches the verified candidate receipts, then a fresh process
  refuses editable registration and binds the installed RECORD hashes to that
  wheel, checks NGSolve/Netgen and GPU runtime pins against package metadata,
  and runs `pip check`. A matching Python version string alone is not
  acceptance. radia-mcp keeps its independent maintained-source policy.
Solver compatibility is checked in the owning numerical lane. It does not make
radia-mcp part of solver Phase 9 or require lock-step package versions.

## ===
## mcp_release — independent MCP publication and release-dual
## ===

radia-mcp publishes through its own reviewed main commit, `radia-mcp-v<VERSION>`
tag, package CI, verified wheel, and PyPI trusted-publishing job. Its wheel must
not contain the retired `radia_mcp.cubit` namespace or `mcp-server-cubit` entry
point. Mixed omega remains in radia-mcp.

After publication, install the verified wheel on the validation host and update
the maintained editable environment on the build host; verify registration plus
a fresh import on both. Do not replace newer reviewed source with an older public
wheel. The validation host must also confirm the affected live source and one
harmless affected tool. Existing build-host clients may remain
`next-launch-pending` until their normal restart; this does not block
release completion. Immediate all-user reconnection is a separate explicit
operation. Do not deploy radia-mcp to any other host, invoke the solver
release installer, restore an older tree, or couple its version to another package.

## ===
## phases — numerical solver release pipeline
## ===

These phases are the internal contract of `tools/release_quad.py`. The table is
the AI-readable summary; independent package lanes are excluded.

| Phase | Action | Mandatory? | Notes |
|------:|--------|-----------|-------|
| Pre | `sync-main`: fetch -> twin-aware rebase (`--empty=drop` skips commits whose patches landed on origin as rebased twins) -> impact preflight -> push | when NAS main diverged from origin/main | added 2026-08-07 after the recurring 40-min manual rebase-archaeology sessions; refuses on a dirty tree, leaves genuine conflicts in place with instructions |
| Pre | `evidence-motor --check`: verify the motor evidence SHA pins | whenever the pinned motor sources change | CPU MEX regeneration is retired; the command refuses to regenerate |
| Pre | `deployment-plan`: print the four lanes and the solver runtime role | anytime | read-only; never contacts a host |
| 1 | Decide the solver minor vs patch increment | always | inspect solver changes since the last Radia tag |
| 2 | Bump the Radia solver's own metadata/runtime version pair | always | independent distributions keep their own version pairs and release lanes |
| **2.5** | **Pre-flight CI validation (4 gates) — ADDED 2026-05-03** | **always** | local equivalent of CI; saves 2-3 round-trips |
| 3 | Stage exactly the release files | always | NO `git add -A`; user has WIP |
| 4 | Commit the reviewed solver release files | always | do not mix unrelated package bumps into the commit |
| 5 | Create the Radia solver tag | always | package tags are created only in their owning independent lane |
| 6 | Push approved main and the solver tag | after solver gates | Radia tag CI uploads `ci-release-context` but leaves PyPI on HOLD; other packages retain independent publication gates |
| 7 | Monitor the exact solver CI run | always | use ci-monitor skill |
| 8 | `phase8 --source-sha SHA --manifest FILE`: install and verify the exact hash-verified release wheel on the validation role | when the runtime must move | never an editable install; the CPU-only build host is never a target; this phase does not install radia-mcp or cubit-mesh-export |
| 9 | `phase9`: read-only verified-wheel runtime consistency probe | after phase 8 | refuses an editable install or a non-release version; package-independent MCP state is outside this probe |
| DoD | `done --source-sha SHA --manifest FILE`: verify the four lane receipts and exact artifacts, then write the public lane summary | always | read-only; does not deploy. Only exit 0 authorizes the release-owner handoff. The release source (the controller by default, or `--release-source`) must be the exact tracked-clean release SHA. A later validation runtime change deploys a newer verified release wheel; never restore an older tree or an editable source. |

## ===
## simulink_candidate — independent GPU distribution boundary
## ===

The `radia` CPU lane does not ship mexcuda, MATLAB, LTspice or Simulink
artifacts. The retired `release_quad simulink-candidate` route has been deleted.
Use the independent `packages/radia-simulink` build, installed-wheel acceptance
and `radia-simulink-release.yml` publication workflow on the Windows GPU
validation host; its MATLAB/Simulink and LTspice acceptance is the `matlab`
lane receipt. Do not attach its artifacts or CI requirements to the CPU lane.

The solver candidate lane never installs, uninstalls, repoints or reconnects
radia-mcp. MCP and Cubit retain their independent release lanes.

CPU `ci-verify` resolves exact-SHA Actions workflow ownership before filtering
independent distribution checks. Shared and CPU checks remain mandatory,
including `radia-native-release-build`; missing ownership metadata fails closed.

## ===
## optuna_candidate — independent radia-optuna release
## ===

`radia-optuna` is a separate distribution with its own repository, CI and
release workflow. `tools/release_quad.py` has no radia-optuna commands and the
Radia solver release never builds, installs or publishes it. Its GPU and MATLAB
qualification follows the same host rule: CPU build on the build host, Linux GPU
on the supercomputer, Windows GPU and MATLAB on the Windows GPU validation host.

## ===
## preflight_gates — impact-scoped pre-push validation
## ===

`tools/ci_preflight.py` is the single local entry point. By default it compares
the candidate and working tree with `origin/main`, then runs only the compact
gates owned by the affected paths:

```powershell
python tools/ci_preflight.py
```

Policy and version consistency always run. A `packages/radia-mcp/` change also
selects the publication-boundary lint, live catalog contracts, affected server
selftests, and impact-selected MCP tests. The default path does not collect or
run the solver validation corpus.

Use explicit escalation only when its evidence is needed:

```powershell
python tools/ci_preflight.py --full                 # lightweight tests/
python tools/ci_preflight.py --validation           # validation collect-only
python tools/ci_preflight.py --full --validation    # run validation_test/
```

Routine CI mirrors this separation: the CPU-only build host owns the fast repository lane,
`radia-mcp` uses its change selector, and native/release validation is manual or
tag-scoped. A release candidate advances only after its selected gates pass for
the exact candidate SHA; a tag must not repeat unrelated package matrices.

## ===
## mcp_quality_review — how to interpret a green radia-mcp matrix
## ===

As of the 2026-06-26 release-candidate review, the radia-mcp public
surface is considered **healthy and practical** when these gates are
green:

* radia-mcp pytest matrix passes (review evidence: 229 passed).
* policy lint passes.
* version consistency passes.
* the checked meta catalog and affected servers' live `tools/list` contracts pass.
* impact-selected package tests and affected server selftests pass.

This is strong evidence for:

* tool definition and inventory consistency,
* docs / policy baseline consistency,
* radia-mcp package test health,
* package version management,
* consolidation into `packages/radia-mcp` instead of the old
  `public-safe curated corpus` layout.

Operational completion for radia-mcp requires its independent release-dual
gates:

* package CI, published-wheel verification, and dependency resolution,
* a verified non-editable wheel and fresh import on the validation host, plus
  an editable development registration and fresh import on the build host,
* validation-host live-source and harmless affected-tool verification.

Existing build-host clients can be `next-launch-pending`; they are not a release
blocker. Numerical methods still require their owning validation evidence, but
the solver lanes, Simulink, and all-user client restart are not universal MCP
publication gates.

Fast MCP fleet regressions live under `packages/radia-mcp/tests/`, including
`test_meta_health.py`; release-specific runtime evidence is emitted by the
release gate rather than tracked as a stale package-local JSON snapshot.

## ===
## ci_failure_modes — known historical CI failures + cause + fix
## ===

| Symptom (CI log line)                                                   | Root cause                                                       | Local gate that catches it     | First seen   |
|-------------------------------------------------------------------------|------------------------------------------------------------------|--------------------------------|--------------|
| `ModuleNotFoundError: No module named 'radia_mcp.elf'`                  | stale test importing extracted subpackage                        | Gate 2 (collect-only sweep)    | v4.27.0      |
| `AssertionError: ... TOOLS.md is stale`                                 | TOOLS.md generator hardcodes a subpackage that no longer exists  | Gate 1 (regen + diff)          | v4.27.1      |
| `FAIL: examples/<dir>/ missing README.md`                               | Historical policy from the retired examples tier                 | do not add new `examples/`; classify as tests, validation, or docs | (multiple)|
| `--selftest: ModuleNotFoundError: No module named 'radia_mcp.<sub>'`    | .yml matrix lists removed subpackage                             | Gate 4 (--selftest import)     | (multiple)   |
| `gh release download ... no assets match the file pattern`              | binaries-release upload race vs tag-CI start                     | (CI has 6-attempt retry as of v4.26)| v4.25.1, v4.27.0|
| `Basic tests failed` w/ collected 0 tests                               | pytest config conflict (pyproject.toml vs pytest.ini)            | Gate 2 (collect-only sweep)    | (rare)       |
| `_radia_pybind import failure: DLL load failed`                         | Radia native extension built against the wrong Python/NGSolve ABI | clean solver rebuild in the owning environment | solver native CI |
| `cubit_mesh_curver import failure: DLL load failed`                    | cubit-mesh-export native extension built against the wrong Python ABI | independent Cubit native provenance/build gate | cubit-mesh-export CI |
| PyPI version appears although its release-tag CI did not pass           | an earlier successful main CI looked up tags by SHA after tags were added | exact per-run `ci-release-context` artifact + workflow contract tests | 4.95.30 / 1.4.24 |

## ===
## recovery — when CI on a tag fails AFTER push
## ===

CI on a tag ref running on a broken commit cannot be rescued by
pushing a fix to main.  Each CI run uploads an immutable
`ci-release-context` artifact containing its ref type, ref name, SHA, run ID,
and the release tags present on that SHA when CI started. The `Release`
workflows automatically reject a successful branch CI or a package tag added
after CI began; manual/PR CI events are ineligible. Normal publication is gated
on `workflow_run` of an exact push-triggered tag-ref CI, and the tag still
points at the broken commit.

`radia-optuna` has one narrow administrative-cancellation recovery lane. A
maintainer may manually select an immutable CI run only when that run itself
was a successful in-repository `push` execution of `build-test.yml` for the
exact commit now named by `radia-optuna-v*`, and both `build-test` and the
installed-wheel MATLAB/Simulink E2E job passed. The release workflow downloads
that run's wheel artifact and reverifies its contents, source version, and tag
before trusted publishing. This does not accept a manual CI run, a local wheel,
or a failed/cancelled tag run; it recovers when GitHub cancels redundant tag CI
after an equivalent exact-SHA push CI already completed successfully.

When the exact source CI actually fails, **the only path forward is patch bump
and re-tag**, per the skill's
"PyPI is immutable" policy.  Do NOT delete + recreate a pushed tag
unless the user explicitly authorizes it (and even then, only if
the version has not gone to PyPI. Always check PyPI before choosing the next
version: releases predating the exact ref-context gate could publish from a
racing successful branch CI even when a sibling tag CI later failed.

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
# 2. Bump only the failed distribution's metadata/runtime version pair.
# 3. Compose that package's release commit on top of the fix:

git add <package-pyproject> <package-runtime-version-file> <package-changelog>
git commit -m "Release <package> <new-version> after <failed gate>"

# 4. Tag at the new commit:
git tag -a <owning-prefix><new-version> -m "..."

# 5. Push main (already up-to-date) + new tags:
git push origin main
git push origin refs/tags/<owning-prefix><new-version>

# 6. Monitor with ci-monitor.
```

Old broken tags stay in the repository as evidence. Never move or reuse a tag;
always verify PyPI state before choosing the next package version.

## ===
## lab_lock_release — package-owned lock handling
## ===

Routine radia-mcp editable source updates do not replace locked native binaries
or console launchers and therefore require no mass process termination. Reconnect
affected clients through their supported controls; build-host clients may remain
`next-launch-pending` for a release-dual update.

Native solver and Cubit artifacts use their owning release skills. Inspect the
specific locked artifact and owning process, coordinate active work, and stop
only the confirmed owner when explicitly authorized. Never kill every MCP,
Python, Cubit, MATLAB, Codex, or Claude process as a generic release phase.

## ===
## version_pairs — independent package metadata/runtime pairs
## ===

| File | Field | Read by |
|------|-------|---------|
| `pyproject.toml` | `version = "X.Y.Z"` | setuptools, pip, PyPI metadata |
| `src/radia/__init__.py` | `__version__ = "X.Y.Z"` | runtime introspection, --version flags, tests |
| `packages/radia-mcp/pyproject.toml` | `version = "X.Y.Z"` | radia-mcp wheel metadata |
| `packages/radia-mcp/src/radia_mcp/__init__.py` | `__version__ = "X.Y.Z"` | runtime, --selftest --version |

Each row is an internal pair for one independently versioned distribution. Rows
do not have to share the same version. A mismatch within a pair
causes:
* Cryptic wheel-install bugs (pip caches the wrong version).
* `tests/test_*_version_consistency.py` failures on CI.
* Users seeing `radia.__version__` not match `pip show radia`.

The owning package CI checks its pair and exits non-zero if violated. Solver
preflight may also inspect these pairs, but it does not authorize bumping or
publishing unrelated distributions.

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

This is the standard companion to release-quad Phase 7 -- after
`git push --tags` the agent should immediately call ci-monitor with
the 3 (or 2) new run IDs and not declare the release "done" until
the monitor exits 0.
"""


_TOPICS = (
    "overview",
    "mcp_release",
    "phases",
    "simulink_candidate",
    "optuna_candidate",
    "preflight_gates",
    "mcp_quality_review",
    "ci_failure_modes",
    "recovery",
    "patch_bump_protocol",
    "lab_lock_release",
    "version_pairs",
    "ci_monitor_skill",
)


def get_release_workflow_documentation(topic: str = "") -> str:
    """Return package-scoped release workflow knowledge.

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
