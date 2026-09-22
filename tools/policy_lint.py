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
  9 no product name on a Radia formulation (citations allowlisted)
 10 commercial-solver interchange stays in its non-public mcp-server
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

# Policy 4 allowlist: genuine LAPACK / HACApK column-major interop.
CBLAS_COLMAJOR_ALLOW = (
    "rad_mmm_matrices.cpp", "rad_relaxation_methods.cpp",
    "rad_hacapk.cpp", "rad_stream_function.cpp",
)


def _sh(cmd):
    p = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    return p.returncode, (p.stdout or "")


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
            f"could not be checked and is therefore not passed")
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
    """Return list of (policy_name, ok, detail) for every policy."""
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

    # 9: a formulation is not named after a product. Reference solvers are
    # instruments, and a public artifact must not carry one as the name of
    # Radia's own method -- it both advertises a brand and says what the work
    # was scored against. Citing a published paper is a different act and is
    # allowed, so the allowlist below holds back the places where the word is
    # somebody else's name and not ours: a paper's title or quoted abstract,
    # the catalog entry citing it, an external input-deck file name, a path a
    # historical run recorded, a list of products an audience uses, and the
    # policy text that forbids commercial artifacts in the first place.
    allowed = {
        "docs/api/EARLY_TIMES_CPP_API_DESIGN.md",
        "packages/radia-mcp/src/radia_mcp/__init__.py",
        "packages/radia-mcp/src/radia_mcp/accelerator/bibliography_index_knowledge.py",
        "packages/radia-mcp/src/radia_mcp/accelerator/knowledge.py",
        "packages/radia-mcp/src/radia_mcp/accelerator/server.py",
        "packages/radia-mcp/src/radia_mcp/presentation/skill.md",
        "packages/radia-mcp/src/radia_mcp/presentation/talk_feedback.py",
        "packages/radia-mcp/src/radia_mcp/radia_ngsolve/bibliography_index_knowledge.py",
        "validation_test/omega_quadrature/candidate_de7feea0/README.md",
        "validation_test/omega_quadrature/omega_quadrature_embedding_p12_20260911.json",
        "validation_test/omega_quadrature/omega_quadrature_small_20260911.json",
        # The audit record of the rename itself: it has to name the label it
        # replaced, or restoring it to check the pre-rename hash is not
        # possible and the record proves nothing.
        "validation_test/c_type_three_engine/results/"
        "c_type_20260903_nonlinear_bdm2_mesh_convergence_certificate.json",
        # The same audit, carried out in a test rather than recorded in a
        # certificate: it has to hold the old label to substitute it back.
        "tests/test_esrf6_repaired_acceptance_record.py",
        # Detectors have to name what they detect. Renaming a detector's
        # needle disarms it silently, which is the failure mode this policy
        # exists to prevent in the first place.
        "tests/test_label_rename_audit.py",
        "tests/test_formulation_attribution.py",
        "tools/policy_lint.py",
    }
    branded = []
    for token in ("tosca", "opera-3d"):
        for line in _git_grep(rf"\b{token}\b", (), extra=("-i",)) or ():
            path = line.split(":", 1)[0]
            if path not in allowed:
                branded.append(line[:160])
    results.append(("Policy 9: no product name on a Radia formulation",
                    not branded, branded[0] if branded else "ok"))

    # 10: a converter between a commercial solver and Radia/NGSolve, and any
    # MCP server exposing one, lives in that tool's non-public mcp-server
    # home.  What makes it non-public is the commercial format it reads and
    # writes, not who wrote it -- "we wrote it ourselves" is not a licence to
    # publish a reader for somebody's product.  Two halves are checked: the
    # converter must not be HERE, and the public catalog must not ADVERTISE a
    # commercial-solver server, which is the same boundary crossed by
    # description rather than by code.
    tools = ("comsol", "femm", "jmag", "cst", "elf")
    converter = re.compile(
        r"\b(?:" + "|".join(tools) + r")_converter\b"
        r"|\b(?:mph|fem|jmag|cst|elf)-to-radia\b"
        r"|\bradia-to-(?:mph|fem|jmag|cst|elf)\b", re.IGNORECASE)
    leaked = []
    hits = _git_grep(r"_converter\|-to-radia\|radia-to-", (), extra=("-i",)) or []
    # The packaging guard in packages/radia-mcp already rejects a commercial
    # converter wired into public scripts; its test builds a synthetic repo
    # containing exactly that to prove it. A detector's fixture has to name
    # what it detects, so it is held back here rather than disarmed.
    converter_allowed = {
        "tools/policy_lint.py",
        "CLAUDE.md",
        "packages/radia-mcp/tests/test_policy_lint.py",
        "packages/radia-mcp/tools/policy_lint.py",
    }
    for line in hits:
        path, _, rest = line.partition(":")
        if path in converter_allowed:
            continue
        if converter.search(rest):
            leaked.append(line[:160])

    # The catalog is data, so it is read rather than grepped: an entry for a
    # commercial tool is a violation however its description is worded.  The
    # documentation-only ELF server is the one carve-out -- no solver, no
    # converter, no vendor-derived numbers.
    catalog_path = os.path.join(
        REPO, "packages", "radia-mcp", "src", "radia_mcp", "meta",
        "catalog.py")
    advertised = []
    try:
        with open(catalog_path, encoding="utf-8") as f:
            catalog_text = f.read()
        for name in re.findall(r'^    "([A-Za-z0-9_\-]+)": \{', catalog_text,
                               re.M):
            stem = name.replace("mcp-server-", "")
            if stem in tools and stem != "elf":
                advertised.append(f"catalog.py advertises {name!r}")
    except OSError as exc:
        advertised.append(f"cannot read the catalog: {exc}")

    bad = leaked + advertised
    results.append(("Policy 10: commercial interchange stays non-public",
                    not bad, bad[0] if bad else "ok"))

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
