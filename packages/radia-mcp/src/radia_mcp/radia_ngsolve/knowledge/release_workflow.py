"""
Package release workflows for the Radia monorepo.

The Radia monorepo ships four independent PyPI distributions and one
versioned Simulink library package from one git tree. The numerical Radia
solver uses its multi-host QUAD gate. radia-mcp and cubit-mesh-export publish
independently and deploy editable sources to LAB and 100 through DUAL lanes;
radia-optuna has its own exact-wheel lane. These release boundaries must not
be collapsed into one coupled version or deployment requirement.

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
# Package-scoped releases: Radia QUAD and MCP/Cubit DUAL

radia-mcp releases independently of Radia. Its release-dual deployment targets
are LAB and 100 only, using editable installs, like cubit-mesh-export.
Never deploy radia-mcp to hibino, mdx1 or mdx2. Isolated CI is not deployment.
For MCP-only updates use update source -> reconnect -> check an affected live
tool, following packages/radia-mcp/docs/maintenance.md, not the QUAD installer.
release-dual is the two-host scope, not a separate release_dual.py command.
The `mcp_release` topic is the authoritative MCP release contract. Other
sections describe the numerical solver, Simulink, or their named package lane;
they do not add MCP deployment targets or completion gates.

This document is the AI-readable canonical reference for the Radia
release flow.  Its canonical local orchestrator is
`tools/release_quad.py`; the former triple-machine workflow is retired.
Radia PyPI publication is explicitly held after tag CI. Accept the exact tag
run's wheel on LAB and hibino and merge its tracked acceptance.json/full6.json/
focused.xml evidence to main under validation_test/esrf_three_engine/results/
candidate_<source-sha-first-9>/. Dispatch the Release workflow on main with
ci_run_id, wheel_sha256, and full acceptance_commit. Its read-only gate checks
the GitHub run/workflow/repository, tag context and peeled SHA, main reachability,
both hosts and all mandatory tests, and wheel/native hashes. Publication uses
the same artifact without rebuilding. Changed bytes require new acceptance;
four-machine QUAD done and Simulink gates remain separate obligations.
Topics: overview, mcp_release, phases, simulink_candidate, optuna_candidate, preflight_gates, mcp_quality_review,
ci_failure_modes, recovery, patch_bump_protocol, lab_lock_release,
version_pairs, ci_monitor_skill.

## ===
## overview — package ownership and verification targets
## ===

The Radia monorepo (`ksugahar/Radia`) ships four independent distributions
to PyPI. Each distribution follows its own declared verification targets:

| Package           | Tag prefix              | What it ships                                     |
|-------------------|-------------------------|---------------------------------------------------|
| radia             | `v`                     | C++ core (.pyd) + Python/MATLAB solver support    |
| cubit-mesh-export | `cubit-mesh-export-v`   | Cubit plugin .ccm/.pyd + check-vol CLI            |
| radia-mcp         | `radia-mcp-v`           | MCP guidance and workflows, including mixed omega;|
|                   |                         | Cubit MCP is distributed by cubit-mesh-export     |
| radia-optuna      | `radia-optuna-v`        | Standalone MATLAB Optuna + lightweight MEX and    |
|                   |                         | generic Simulink optimization support             |
| Radia Simulink library | Radia GitHub Release asset | `.slx`, MATLAB support, Level-2 MATLAB S-Functions, standalone MEX handles, runtime DLLs, manifest and checksums |

Packages are independently released. radia-mcp and cubit-mesh-export use
release-dual (LAB and 100); Radia's solver QUAD covers LAB, 100, mdx1 and mdx2.
None of hibino/mdx1/mdx2 receives MCP deployment. Compatibility is still checked
where an affected tool uses another package:
* radia-mcp imports radia at runtime for several tools, and a
  schema mismatch is silent until a tool crashes.
* cubit-mesh-export owns Cubit meshing, export, and its C++ plugin binaries.
  Radia consumes the checked `.vol` boundary; it does not couple solver release
  acceptance to the exporter version.
* LAB and 100号機 use NAS editable installs for `radia`,
  `cubit-mesh-export`, and `radia-mcp`.
  On 100号機 the editable path must be the mapped drive
  configured on the target, not the self-referential UNC path, because
  Windows can crash while loading `_radia_pybind.pyd` from that UNC form.
* mdx1/mdx2 install the exact accepted `radia` wheel for solver release
  acceptance. They do not receive Cubit or radia-mcp deployments.
* hibino does not receive radia-mcp deployments; it remains a computation host,
  not an MCP release-dual target.
Solver compatibility is checked in the owning numerical lane. It does not make
radia-mcp part of solver Phase 9 or require lock-step package versions.

`radia-optuna` is independently versioned and does not enter the Radia solver,
Cubit, or radia-mcp deployment lanes. Its release gate downloads the exact wheel
artifact from a successful `main` push CI run, verifies the same wheel on all
its declared targets, and binds the CI run, source commit, package version, and
SHA-256 before its tag, PyPI upload, and GitHub Release are authorized.

## ===
## mcp_release — independent MCP publication and release-dual
## ===

radia-mcp publishes through its own reviewed main commit, `radia-mcp-v<VERSION>`
tag, package CI, verified wheel, and PyPI trusted-publishing job. Its wheel must
not contain the retired `radia_mcp.cubit` namespace or `mcp-server-cubit` entry
point. Mixed omega remains in radia-mcp.

After publication, update the maintained editable sources on LAB and 100 and
verify registration plus a fresh import on both. LAB must also confirm the
affected live source and one harmless affected tool. Existing 100 clients may
remain `next-launch-pending` until their normal restart; this does not block
release completion. Immediate all-user reconnection is a separate explicit
operation. Do not deploy radia-mcp to hibino/mdx1/mdx2, invoke the solver QUAD
installer, restore an older tree, or couple its version to another package.

## ===
## phases — numerical solver release pipeline
## ===

These phases are the internal contract of `tools/release_quad.py`. The table is
the AI-readable summary; independent package lanes are excluded.

| Phase | Action | Mandatory? | Notes |
|------:|--------|-----------|-------|
| Pre | `sync-main`: fetch -> twin-aware rebase (`--empty=drop` skips commits whose patches landed on origin as rebased twins) -> impact preflight -> push | when NAS main diverged from origin/main | added 2026-08-07 after the recurring 40-min manual rebase-archaeology sessions; refuses on a dirty tree, leaves genuine conflicts in place with instructions |
| Pre | `evidence-motor [--check|--force]`: rebuild the MEX, ship the snapshot closure to HIBINO over scp, run the MATLAB generator SYNCHRONOUSLY over ssh, fetch the artifact, verify SHA pins, align the pytest test-count | whenever `src/matlab/radia_mex.cpp` / `matlab/+radia/setup.m` / the generator changed | Windows OpenSSH reaps detached children on session exit — the run MUST stay synchronous; snapshot closure = matlab/, src/matlab/, tests/matlab/, validation_test/{radia_mcp,maglev}/, docs/maglev/demos/team28/, pyproject.toml |
| 1 | Decide the solver minor vs patch increment | always | inspect solver changes since the last Radia tag |
| 2 | Bump the Radia solver's own metadata/runtime version pair | always | independent distributions keep their own version pairs and release lanes |
| **2.5** | **Pre-flight CI validation (4 gates) — ADDED 2026-05-03** | **always** | local equivalent of CI; saves 2-3 round-trips |
| 3 | Stage exactly the release files | always | NO `git add -A`; user has WIP |
| 4 | Commit the reviewed solver release files | always | do not mix unrelated package bumps into the commit |
| 5 | Create the Radia solver tag | always | package tags are created only in their owning independent lane |
| 6 | Push approved main and the solver tag | after solver gates | Radia tag CI uploads `ci-release-context` but leaves PyPI on HOLD; other packages retain independent publication gates |
| 7 | Monitor the exact solver CI run | always | use ci-monitor skill |
| 8 | Deploy/verify the Radia solver on its current QUAD targets | always | this phase does not install radia-mcp or cubit-mesh-export |
| 8S | Verify the exact versioned Simulink ZIP on LAB / 100号機 / mdx1 / mdx2 | for every Simulink revision | `simulink-candidate --package <zip> --target all` |
| 9 | Cross-machine solver consistency probe | always | package-independent MCP state is outside this probe |
| DoD | Re-run preflight, exact-source/editable checks, and Phase 9; bind the exact ZIP hash to the same HEAD without changing the verified editable pointers | always | `done --simulink-package <zip>`; only exit 0 authorizes GitHub Release publication. `done` requires the active LAB source to be the exact tracked-clean release SHA, requires that SHA to equal the peeled `v<radia-version>` tag, verifies LAB/100号機 editable metadata and import origins, and refuses while NAS main != origin/main. Any later repoint must select and verify an explicitly intended current source; never restore an older tree merely because its path was once canonical. |

## ===
## simulink_candidate — exact MEX + SLX publication gate
## ===

The production human interface is released as a full Radia Simulink library
ZIP. Build it only from the final pushed commit. The package contains
`radia_simulink_library.slx`, application/support `.m` and `.slx` files,
standalone `radia_mex`, readable IH Level-2 MATLAB S-Functions, required runtime DLLs,
`manifest.json`, and an external `SHA256SUMS.txt`.

```powershell
python tools/package_simulink_release.py `
  --full-library --mex-dir matlab --output-dir dist/simulink
python tools/verify_simulink_release.py `
  dist/simulink/radia-simulink-library-vX.Y.Z.zip `
  --matlab "C:\\Program Files\\MATLAB\\R2026a\\bin\\matlab.exe"
python tools/release_quad.py all
python tools/release_quad.py simulink-candidate `
  --package dist/simulink/radia-simulink-library-vX.Y.Z.zip --target all
python tools/release_quad.py done `
  --simulink-package dist/simulink/radia-simulink-library-vX.Y.Z.zip
```

The candidate state is keyed by the ZIP SHA-256. Rebuilding or modifying the
archive invalidates the four-machine evidence. `done` also requires the
manifest commit to equal repository `HEAD`, so validation from another commit
cannot authorize publication. A successful `done` leaves LAB and 100号機 on
the exact editable source it verified; it never swaps a release worktree for a
possibly older path as a side effect. Advance the explicitly intended editable
source to current `main` and verify it in place; do not restore a prior tree by
historical path name. Upload the ZIP,
external `manifest.json`, and `SHA256SUMS.txt` to the matching Radia GitHub
Release only after `done` exits 0.

The standalone IH preview remains supported by the same packager without
`--full-library`; it does not replace the production full-library gate.

The solver candidate lane never installs, uninstalls, repoints, or reconnects
radia-mcp. MCP uses its own LAB/100 release-dual and records its intended source
there. A solver deployment must preserve that state without requiring an
override environment variable. Likewise, Cubit publication and plugin
deployment belong to the cubit-mesh-export release lane.

## ===
## optuna_candidate — exact standalone wheel publication gate
## ===

`radia-optuna` is a separate distribution even though its canonical sources
remain in the Radia monorepo. It must not trigger either the Radia solver or
Cubit deployment merely to publish the MATLAB Optuna component.

```powershell
python tools/release_quad.py optuna-candidate `
  --ci-run-id <successful-main-push-CI-run> --target all
python tools/release_quad.py optuna-done --wheel <retained-wheel-path>
git tag -a radia-optuna-vX.Y.Z <verified-commit> -m "radia-optuna X.Y.Z"
git push origin refs/tags/radia-optuna-vX.Y.Z
gh workflow run release-radia-optuna.yml `
  -f ci_run_id=<successful-main-push-CI-run> `
  -f candidate_sha256=<verified-wheel-sha256>
```

`optuna-candidate` accepts only a successful `main` push run of the repository
CI workflow whose `build-test` and installed-wheel MATLAB/Simulink jobs both
succeeded. It downloads that run's `radia-optuna-wheel` artifact, reruns the
wheel verifier, and executes the installed-wheel MATLAB/Simulink contract on
LAB, 100号機, mdx1, and mdx2. `optuna-done` rejects a rebuilt wheel, a source
commit other than `origin/main`, a version mismatch, or any missing machine
result. The manual release workflow re-downloads the same CI artifact, checks
the QUAD SHA-256, and publishes that exact file to both PyPI and the matching
GitHub Release.

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

Routine CI mirrors this separation: mdx owns the fast repository lane,
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
* editable registration and fresh imports on LAB and 100,
* LAB live-source and harmless affected-tool verification.

Existing 100 clients can be `next-launch-pending`; they are not a release
blocker. Numerical methods still require their owning validation evidence, but
solver QUAD, Simulink, and all-user client restart are not universal MCP
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
affected clients through their supported controls; 100 clients may remain
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
| `packages/cubit-mesh-export/pyproject.toml` | `version = "X.Y.Z"` | cubit-mesh-export wheel metadata |
| `packages/cubit-mesh-export/src/cubit_mesh_export/__init__.py` | `__version__ = "X.Y.Z"` | runtime, COMPAT constants |

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
