#!/usr/bin/env python
"""The 8 Radia source policies as ONE reusable module.

Single source of truth for BOTH:
  - .github/workflows/policy-lint.yml   (the CI gate)
  - tools/ci_preflight.py               (the local pre-push gate)

so the two can never drift (previously ci_preflight re-implemented the
policies inline -- if the workflow gained a policy, the local gate would
silently miss it).  Mirrors the historical inline bash greps; uses
`git grep` / `git ls-files` so it sees TRACKED working-tree content (what
CI checks out on a fresh runner).

    python tools/policy_lint.py        # run all 8; exit 1 on any violation
    python tools/policy_lint.py --quiet

Policies (see CLAUDE.md):
  1 no FldUnits() in examples            5 no generated files at repo root
  2 no tracked binaries                  6 no legacy src/python import path
  3 no Helmholtz WAVE kernel in core     7 README.md in every examples/ dir
  4 no CblasColMajor in core (allowlist)
  8 HDiv geometry/field pairs use the central capability table
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


def _git_grep(pattern, pathspecs, extra=()):
    rc, out = _sh(["git", "grep", "-n", *extra, pattern, "--", *pathspecs])
    if rc not in (0, 1):
        return None
    return [ln for ln in out.splitlines() if ln.strip()]


def check_all():
    """Return list of (policy_name, ok, detail) for all 8 policies."""
    results = []

    # 1: no FldUnits() in examples/*.py
    hits = _git_grep("FldUnits", ["examples/*.py", "examples/**/*.py"])
    results.append(("Policy 1: no FldUnits() in examples", not hits,
                    hits[0] if hits else "ok"))

    # 2: no tracked binaries
    rc, out = _sh(["git", "ls-files", "*.pyd", "*.dll", "*.so",
                   "*.lib", "*.exe", "*.obj"])
    bins = [b for b in out.splitlines() if b.strip()]
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
    rc, out = _sh(["git", "ls-files", "*.msh", "*.vtu", "*.vtk", "*.vol", "*.vts"])
    root_gen = [f for f in out.splitlines() if f.strip() and "/" not in f.strip()]
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

    # 7: every examples/*/ has README.md
    ex = os.path.join(REPO, "examples")
    missing = []
    if os.path.isdir(ex):
        for d in sorted(os.listdir(ex)):
            dp = os.path.join(ex, d)
            if os.path.isdir(dp) and not os.path.exists(os.path.join(dp, "README.md")):
                missing.append(f"examples/{d}/")
    results.append(("Policy 7: README.md in every example dir", not missing,
                    str(missing[:3]) if missing else "ok"))

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

    return results


def main(argv=None):
    quiet = "--quiet" in (argv if argv is not None else sys.argv[1:])
    results = check_all()
    nfail = 0
    for name, ok, detail in results:
        if ok:
            if not quiet:
                print(f"PASS  {name}")
        else:
            nfail += 1
            print(f"FAIL  {name}: {detail}")
    if nfail:
        print(f"\n{nfail}/8 policies FAILED", file=sys.stderr)
        return 1
    if not quiet:
        print("\nall 8 policies pass")
    return 0


if __name__ == "__main__":
    sys.exit(main())
