#!/usr/bin/env python
"""The Radia source policies as ONE reusable module.

Single source of truth for BOTH:
  - .github/workflows/policy-lint.yml   (the CI gate)
  - tools/ci_preflight.py               (the local pre-push gate)

so the two can never drift (previously ci_preflight re-implemented the
policies inline -- if the workflow gained a policy, the local gate would
silently miss it).  Mirrors the historical inline bash greps; uses
`git grep` / `git ls-files` so it sees TRACKED working-tree content (what
CI checks out on a fresh runner).

    python tools/policy_lint.py        # run them all; exit 1 on any violation
    python tools/policy_lint.py --quiet

Policies (see CLAUDE.md):
  1 no FldUnits() in Radia source        5 no generated files at repo root
  2 no tracked binaries                  6 no legacy src/python import path
  3 no Helmholtz WAVE kernel in core     7 no tracked retired examples/ tier
  4 no CblasColMajor in core (allowlist)
  8 HDiv geometry/field pairs use the central capability table
  9 no scrubbed commercial-solver brand in code or stored results
"""
from __future__ import annotations

import os
import re
import subprocess
import sys

# UNC-safe repo root (see tools/ci_preflight.py for the rationale: never
# resolve() on the LAB S: drive -- it canonicalises to a UNC form Python
# file reads reject).
_THIS = os.path.abspath(__file__)
REPO = os.path.dirname(os.path.dirname(_THIS))
_norm = REPO.replace("\\", "/")
if "192.168.11.100" in _norm and "/Radia/01_GitHub" in _norm:
    REPO = "S:" + _norm[_norm.index("/Radia/01_GitHub"):]

# Policy 9: brands whose scrub is complete, so a reappearance is a regression
# rather than unfinished work.  COMSOL, FEMM, JMAG, CST and ELF are NOT here
# yet -- they still appear widely in code and join this tuple as each one's
# scrub lands.  The formulation those names once carried is Simkin and
# Trowbridge's mixed total/reduced scalar potential (``simkin1980three``).
SCRUBBED_BRANDS = ("TOSCA",)

# Policy 10: solvers whose MCP servers must never be published, whatever the
# repository they live in.  This is a stricter rule than Policy 9 and does not
# wait for a brand's prose scrub: a catalog entry carries an install path and a
# repository URL, so it publishes the server rather than merely naming a tool.
COMMERCIAL_SOLVERS = ("COMSOL", "FEMM", "JMAG", "CST", "TOSCA", "Opera")

# Policy 9 allowlist: surfaces that name a brand for a legitimate reason --
# a citation of published work, or a declaration that something is absent.
BRAND_CITATION_ALLOW = (
    # the cited paper's own title, and the verbatim abstract of the paper the
    # formulation comes from
    "packages/radia-mcp/src/radia_mcp/accelerator/bibliography_index_knowledge.py",
    "packages/radia-mcp/src/radia_mcp/radia_ngsolve/bibliography_index_knowledge.py",
    "packages/radia-mcp/src/radia_mcp/bibliography/data/references.bib",
    # conference-positioning advice, where naming the incumbent tools is the
    # substance of the advice rather than a validation provenance claim
    "packages/radia-mcp/src/radia_mcp/presentation/talk_feedback.py",
    # this file, which has to spell the brand out to forbid it
    "tools/policy_lint.py",
)

# Policy 4 allowlist: genuine LAPACK / HACApK column-major interop.
CBLAS_COLMAJOR_ALLOW = (
    "rad_mmm_matrices.cpp", "rad_relaxation_methods.cpp",
    "rad_hacapk.cpp", "rad_stream_function.cpp",
)


def _sh(cmd):
    if cmd[0] == "git":
        # CI service and interactive checkout owners may differ. Trust only
        # this invocation's repository, as ci_preflight's diff helper does.
        cmd = [cmd[0], "-c", f"safe.directory={REPO}", *cmd[1:]]
    p = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    return p.returncode, (p.stdout or "") + (p.stderr or "" if p.returncode not in (0, 1) else "")


class LintExecutionError(RuntimeError):
    """The lint could not run a check, which is not the same as passing it."""


def _git_grep(pattern, pathspecs, extra=()):
    # git grep exits 0 with matches and 1 without; anything else means the
    # search itself failed (not a repository, bad pathspec, git missing).
    # That used to return None, which every caller read as "no violation" --
    # a broken git made all policies pass.  A check that could not run has
    # to fail the run, so this raises and main() turns it into exit 1.
    rc, out = _sh(["git", "grep", "-n", *extra, pattern, "--", *pathspecs])
    if rc not in (0, 1):
        raise LintExecutionError(
            f"git grep exited {rc} for pattern {pattern!r}; the policy "
            f"could not be checked and is therefore not passed: {out.strip()}")
    return [ln for ln in out.splitlines() if ln.strip()]


def _git_ls_files(*patterns):
    # Same contract as _git_grep: a listing that could not be produced is
    # not an empty listing.  Three policies read "no tracked files" as
    # "no violation", so a failing git used to pass all three.
    rc, out = _sh(["git", "ls-files", *patterns])
    if rc != 0:
        raise LintExecutionError(
            f"git ls-files exited {rc} for {patterns!r}; the policy could "
            f"not be checked and is therefore not passed")
    return [ln for ln in out.splitlines() if ln.strip()]


def check_all():
    """Return list of (policy_name, ok, detail) for all 8 policies."""
    results = []

    # 1: the removed unit-switching API must not return to executable Radia.
    hits = _git_grep("FldUnits", ["src/radia/*.py", "src/radia/**/*.py"])
    results.append(("Policy 1: no FldUnits() in Radia source", not hits,
                    hits[0] if hits else "ok"))

    # 2: no tracked binaries
    bins = _git_ls_files("*.pyd", "*.dll", "*.so", "*.lib", "*.exe", "*.obj")
    results.append(("Policy 2: no tracked binaries", not bins,
                    str(bins[:3]) if bins else "ok"))

    # 3: no Helmholtz wave kernel in src/core (hodge / decomposition excepted)
    hits = _git_grep("helmholtz", ["src/core/*.cpp", "src/core/*.h",
                                   "src/core/**/*.cpp", "src/core/**/*.h"],
                     extra=["-i"]) or []
    bad = [h for h in hits if not any(
        x in h.lower() for x in ("helmholtz-hodge", "helmholtz hodge",
                                 "helmholtz-decomposition", "helmholtz decomposition"))]
    results.append(("Policy 3: no Helmholtz wave kernel in core", not bad,
                    bad[0] if bad else "ok"))

    # 4: no CblasColMajor in src/core (LAPACK/HACApK allowlisted)
    hits = _git_grep("CblasColMajor", ["src/core/*.cpp", "src/core/*.h",
                                       "src/core/**/*.cpp", "src/core/**/*.h"]) or []
    bad = [h for h in hits if not any(a in h for a in CBLAS_COLMAJOR_ALLOW)]
    results.append(("Policy 4: no CblasColMajor in core", not bad,
                    bad[0] if bad else "ok"))

    # 5: no TRACKED generated files at repo root
    root_gen = [f for f in _git_ls_files("*.msh", "*.vtu", "*.vtk", "*.vol",
                                         "*.vts") if "/" not in f]
    results.append(("Policy 5: no tracked generated files at root", not root_gen,
                    str(root_gen) if root_gen else "ok"))

    # 6: forbid the legacy import path (use src/radia).  The needle is built
    # from fragments below so this module never self-matches -- Policy 6 greps
    # for that token on a line that also references a sys.path insertion.
    _legacy = "src/" + "python"
    hits = _git_grep(_legacy, ["*.py", "**/*.py"]) or []
    bad = [h for h in hits if "sys.path" in h]
    results.append((f"Policy 6: no legacy {_legacy} import path", not bad,
                    bad[0] if bad else "ok"))

    # 7: the examples tier is retired. Public demonstrations belong in docs,
    # numerical evidence in validation_test, and regressions in tests.
    examples = _git_ls_files("examples")
    results.append(("Policy 7: no tracked retired examples tier", not examples,
                    str(examples[:3]) if examples else "ok"))

    # 8: Piola-mapped HDiv field and geometry orders are dimension-dependent.
    # Keep one explicit table instead of reintroducing a tempting but incorrect
    # global relation such as ``geometry_order <= hdiv_order + 1``.
    cap = os.path.join(REPO, "src", "radia", "vim", "_capabilities.py")
    required = {
        "_vim.py": "validate_hdiv_configuration",
        "_vim2d.py": "validate_hdiv_configuration",
        "_solve.py": "validate_hdiv_configuration",
    }
    bad = []
    try:
        with open(cap, encoding="utf-8") as f:
            cap_text = f.read()
        for token in ('HDivCapability(2, "quad", 2, (1, 2, 3), 3)',
                      'HDivCapability(3, "hex", 2, (1, 2), 2)',
                      'HDivCapability(3, "wedge", 2, (1, 2), 2)'):
            if token not in cap_text:
                bad.append(f"_capabilities.py missing {token}")
    except OSError as exc:
        bad.append(f"cannot read _capabilities.py: {exc}")

    vim_dir = os.path.join(REPO, "src", "radia", "vim")
    relation = re.compile(
        r"(?:geometry_order|curve_order)\s*(?:<=|>=|<|>)\s*"
        r"(?:self\.)?(?:order|p)(?:\s*[+-]\s*\d+)?|"
        r"(?:self\.)?(?:order|p)(?:\s*[+-]\s*\d+)?\s*"
        r"(?:<=|>=|<|>)\s*(?:geometry_order|curve_order)")
    for name in sorted(os.listdir(vim_dir)):
        if not name.endswith(".py") or name == "_capabilities.py":
            continue
        path = os.path.join(vim_dir, name)
        try:
            with open(path, encoding="utf-8") as f:
                text = f.read()
        except OSError as exc:
            bad.append(f"cannot read {name}: {exc}")
            continue
        if name in required and required[name] not in text:
            bad.append(f"{name} does not use the central HDiv capability validator")
        match = relation.search(text)
        if match:
            line = text.count("\n", 0, match.start()) + 1
            bad.append(f"{name}:{line}: ad-hoc HDiv geometry/order relation")
    results.append(("Policy 8: central HDiv geometry/order capabilities", not bad,
                    bad[0] if bad else "ok"))

    # 9: a scrubbed commercial-solver brand must not return to Radia's own
    # identifiers, labels, source, tests or stored result JSONs.  Prose that
    # CITES published work keeps its brand words, which is why this checks
    # code and result data rather than Markdown.
    hits = _git_grep("|".join(SCRUBBED_BRANDS),
                     ["src/*", "tests/*", "matlab/*", "tools/*",
                      "validation_test/*.py", "validation_test/**/*.py",
                      "validation_test/**/*.json", "packages/**/*.py"],
                     extra=["-E", "-w", "-i"]) or []
    bad = [h for h in hits
           if not any(h.startswith(a + ":") for a in BRAND_CITATION_ALLOW)]
    results.append((f"Policy 9: no scrubbed commercial brand "
                    f"({'/'.join(SCRUBBED_BRANDS)}) in Radia code or results",
                    not bad, bad[0] if bad else "ok"))

    # 10: the published radia-mcp catalog must not advertise an MCP server for
    # a commercial solver.  Naming the reference tools is one thing; shipping
    # a public directory entry with their install path and repository URL
    # publishes the servers themselves, which the boundary forbids outright.
    catalog = os.path.join(REPO, "packages", "radia-mcp", "src", "radia_mcp",
                           "meta", "catalog.py")
    bad = []
    try:
        with open(catalog, encoding="utf-8") as f:
            entries = re.findall(r'^    "([A-Za-z0-9_-]+)": \{\n(.*?)^    \},$',
                                 f.read(), re.S | re.M)
        for name, body in entries:
            named = [b for b in COMMERCIAL_SOLVERS
                     if re.search(rf"\b{b}\b", name + body, re.I)]
            if named:
                bad.append(f"catalog entry {name!r} advertises {'/'.join(named)}")
    except OSError as exc:
        bad.append(f"cannot read the radia-mcp catalog: {exc}")
    results.append(("Policy 10: the public MCP catalog advertises no "
                    "commercial-solver server", not bad, bad[0] if bad else "ok"))

    return results


def main(argv=None):
    quiet = "--quiet" in (argv if argv is not None else sys.argv[1:])
    try:
        results = check_all()
    except (LintExecutionError, OSError) as exc:
        # A lint that could not run has not passed.  Say so and fail the
        # run; silently reporting PASS here is how a broken git once made
        # every policy green.
        print(f"FAIL  policy lint could not run: {exc}", file=sys.stderr)
        return 1
    nfail = 0
    for name, ok, detail in results:
        if ok:
            if not quiet:
                print(f"PASS  {name}")
        else:
            nfail += 1
            print(f"FAIL  {name}: {detail}")
    if nfail:
        print(f"\n{nfail}/{len(results)} policies FAILED", file=sys.stderr)
        return 1
    if not quiet:
        print(f"\nall {len(results)} policies pass")
    return 0


if __name__ == "__main__":
    sys.exit(main())
