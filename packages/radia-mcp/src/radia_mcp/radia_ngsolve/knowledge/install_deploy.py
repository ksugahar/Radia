"""Current Radia installation, machine-role, and deployment knowledge."""

from __future__ import annotations


INSTALL_DEPLOY = """\
# Radia install and deployment contract

Radia uses isolated build and release artifacts. Development checkouts are
editable where agents work; CI and release verification never depend on native
binaries copied from another machine.

Available topics: overview, development, ci_compute, release, mcp_release, cubit, failures.

## overview

| Machine | Primary role | Installation rule |
|---|---|---|
| LAB | development and operations | editable `radia`, `cubit-mesh-export`, `radia-mcp` |
| 100号機 | development and execution | editable packages from its mapped repository path |
| mdx | CI priority and MATLAB compute | isolated per-run environment; validated workloads |
| hibino | long optimization and validation when available | release or job-specific environment |

Regular GitHub Actions CI uses the `mdx` runner label. LAB is not a CI runner.
Validation studies run on hibino first, or on mdx only when hibino is
unavailable and the mdx CI queue is idle. They are not part of every
source-change CI run.

## development

LAB and 100号機 use editable installs so Python source changes are visible on
the next import. A running MCP process still owns already-imported modules and
registered tool objects. Reconnect affected clients; compatible same-root Python
edits may use a reviewed safe reload. Root, dependency, entry-point and schema
changes require reconnecting, not just reloading Python modules.

Native extensions are built locally against the selected Python environment.
That environment owns NGSolve, Netgen, pybind11, and pip `mkl-devel`. Do not
copy `.pyd`, `.dll`, or `.mex*` files between machines to update an editable
installation.

## ci_compute

mdx owns routine CI. Each native job creates a run-local virtual environment,
installs the exact NGSolve/Netgen pins and `mkl-devel`, builds from the checked
commit, and retains the resulting artifact as CI evidence. Normal CI runs only
compact impact-selected regression tests from changed package paths. Solver studies,
benchmarks, paper data, and machine comparisons belong to `validation_test/`
and run on hibino first, with mdx reserved as an idle-CI fallback.

## release

For the Radia solver, use `tools/release_quad.py` and the `release-quad` skill. The release candidate
must be built from one immutable commit, pass its required CI and package gates,
and pass the four-machine `done` gate before GitHub Release publication. The
release workflow consumes the accepted artifact for that exact commit; it does
not rebuild from a mutable checkout or upload a developer-machine binary.

This is not the radia-mcp release lane; use topic `mcp_release` for that package.
After an interrupted release, inspect the selected source and fix forward;
do not automatically restore an old editable tree. Do not repair deployment by manually dropping
native files into `site-packages`.

## mcp_release

radia-mcp publishes independently through `radia-mcp-v<VERSION>` and its package
CI/PyPI workflow. Its release-dual updates editable installations on LAB and
100 only, not hibino/mdx1/mdx2. Do not run the solver's four-host release gate.
Mixed omega remains in radia-mcp; Cubit MCP belongs to cubit-mesh-export.

Completion requires passing package checks, verified publication, both hosts'
editable registration and fresh imports, plus LAB live source and harmless
affected-tool verification. Existing 100 clients may remain next-launch-pending
until their normal restart; they do not block release completion. Failed
installation/import still blocks deployment completion. Immediate all-user
reconnection is a separate explicit request, not a reason to force restarts
during every release. Never claim an unqueried client is live-verified.

## cubit

`cubit-mesh-export` independently owns the Cubit `.ccm` backend and the toolbar
that runs inside Cubit's bundled PySide runtime. Install or verify it with
`cubit-plugin-install`; use `check-vol` after export. Normal Radia Python does
not depend on PySide or Qt. SAT and STEP are CAD interchange inputs; checked
`.vol` is the solver boundary.

## failures

When imports resolve the wrong checkout, inspect `module.__file__`, distribution
metadata, and editable `direct_url.json` before reinstalling. When a native
module fails to load, compare the Python ABI, NGSolve/Netgen pins, MKL runtime,
and artifact hash with the build manifest. When an MCP tool list is stale after
an editable source update, reload the server code or reconnect the client; an
editable install alone cannot replace objects already held by a running process.
"""


_TOPICS = (
    "overview",
    "development",
    "ci_compute",
    "release",
    "mcp_release",
    "cubit",
    "failures",
)


def get_install_deploy_documentation(topic: str = "") -> str:
    """Return the full deployment contract or one named section."""

    if not topic:
        return INSTALL_DEPLOY
    if topic not in _TOPICS:
        return (
            f"Unknown topic: {topic!r}. Available topics: {', '.join(_TOPICS)}\n\n"
            "Pass an empty string for the full document."
        )

    marker = f"## {topic}\n"
    start = INSTALL_DEPLOY.find(marker)
    if start < 0:
        return f"Topic {topic!r} declared but not found in document."
    next_section = INSTALL_DEPLOY.find("\n## ", start + len(marker))
    end = len(INSTALL_DEPLOY) if next_section < 0 else next_section
    return INSTALL_DEPLOY[start:end].rstrip() + "\n"
