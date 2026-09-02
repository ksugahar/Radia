"""Current Radia installation, machine-role, and deployment knowledge."""

from __future__ import annotations


INSTALL_DEPLOY = """\
# Radia install and deployment contract

Radia uses isolated build and release artifacts. Development checkouts are
editable where agents work; CI and release verification never depend on native
binaries copied from another machine.

Available topics: overview, development, ci_compute, release, cubit, failures.

## overview

| Machine | Primary role | Installation rule |
|---|---|---|
| LAB | development and operations | editable `radia`, `cubit-mesh-export`, `radia-mcp` |
| 100号機 | development and execution | editable packages from its mapped repository path |
| mdx | CI priority and MATLAB compute | isolated per-run environment; validated workloads |
| hibino | long optimization and validation when available | release or job-specific environment |

Regular GitHub Actions CI uses the `mdx` runner label. LAB is not a CI runner.
Validation studies are dispatched deliberately to mdx or hibino and are not
part of every source-change CI run.

## development

LAB and 100号機 use editable installs so Python source changes are visible on
the next import. A running MCP process still owns already-imported modules and
registered tool objects; invoke the server's reload tool or reconnect once when
the reload tool itself has changed.

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
and run explicitly on mdx or hibino.

## release

Use `tools/release_quad.py` and the `release-quad` skill. The release candidate
must be built from one immutable commit, pass its required CI and package gates,
and pass the four-machine `done` gate before GitHub Release publication. The
release workflow consumes the accepted artifact for that exact commit; it does
not rebuild from a mutable checkout or upload a developer-machine binary.

After an interrupted release, use the release tool's documented recovery or
editable-restore operation. Do not repair deployment by manually dropping
native files into `site-packages`.

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
