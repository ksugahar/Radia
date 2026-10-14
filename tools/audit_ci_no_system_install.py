#!/usr/bin/env python3
"""Audit: self-hosted CI jobs must NOT pip-install Radia packages into system Python.

The self-hosted runner executes on mdx and shares its service-account Python.
In a job on that runner, any of:

  * ``pip install -e .``            (radia editable)
  * ``pip install -e packages/...`` (radia-mcp / cubit-mesh-export editable)
  * ``pip install --force-reinstall <our-wheel>.whl``

rewrites the host's installed package state to the CI checkout (or replaces it
with a wheel snapshot). CI must use a run-local virtual environment or import
the checkout through ``PYTHONPATH`` instead.

Jobs on GitHub-hosted runners (``ubuntu-latest``, ``windows-latest``, ...) are
isolated, so the same install there is harmless and is NOT flagged.

This audit fails (exit 1) if a self-hosted job reintroduces such an install,
so the regression cannot sneak back.  Third-party tools (``pip install mcp``,
``pip install build``) and isolated installs (``--target``) are allowed.

Run:  python tools/audit_ci_no_system_install.py     # exit 0 = clean
"""
from __future__ import annotations

import glob
import os
import re
import sys

import yaml

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_WF_DIR = os.path.join(_REPO, ".github", "workflows")

# `pip install ... -e .`  or  `pip install ... -e packages/<anything>`
_EDITABLE = re.compile(r"pip\s+install\b[^\n#]*\s-e\s+(\.|packages/)", re.I)
# `pip install ... --force-reinstall ... <something>.whl`
_WHEEL_FORCE = re.compile(
    r"pip\s+install\b[^\n#]*--force-reinstall\b[^\n#]*\.whl", re.I)

# GitHub-hosted runner labels are isolated (a fresh VM per run).  Anything
# else (custom label / 'self-hosted') may share a persistent host interpreter.
_GH_HOSTED = re.compile(r"^(ubuntu|windows|macos)-(latest|\d[\w.-]*)$", re.I)


def _is_self_hosted(runs_on) -> bool:
    labels = runs_on if isinstance(runs_on, list) else [runs_on]
    # Self-hosted if ANY label is not a github-hosted image name.
    return any(not _GH_HOSTED.match(str(lbl).strip()) for lbl in labels if lbl)


def audit():
    violations = []
    files = (sorted(glob.glob(os.path.join(_WF_DIR, "*.yml"))) +
             sorted(glob.glob(os.path.join(_WF_DIR, "*.yaml"))))
    for wf in files:
        try:
            with open(wf, encoding="utf-8") as fh:
                doc = yaml.safe_load(fh)
        except Exception:
            continue  # YAML validity is checked elsewhere
        if not isinstance(doc, dict):
            continue
        for jobname, job in (doc.get("jobs", {}) or {}).items():
            if not isinstance(job, dict):
                continue
            if not _is_self_hosted(job.get("runs-on", "")):
                continue  # GitHub-hosted -> isolated -> safe
            for step in (job.get("steps", []) or []):
                run = (step.get("run", "") or "") if isinstance(step, dict) else ""
                sname = step.get("name", "?") if isinstance(step, dict) else "?"
                for ln in run.splitlines():
                    if ln.lstrip().startswith("#"):
                        continue
                    if _EDITABLE.search(ln):
                        violations.append((wf, jobname, sname,
                                           "editable install into persistent "
                                           "system Python", ln.strip()))
                    elif _WHEEL_FORCE.search(ln):
                        violations.append((wf, jobname, sname,
                                           "force-reinstall wheel into persistent "
                                           "system Python", ln.strip()))
    return violations


def main():
    if not os.path.isdir(_WF_DIR):
        print(f"[SKIP] no workflows dir at {_WF_DIR}")
        return 0
    violations = audit()
    if not violations:
        print("[OK] no self-hosted CI job installs a Radia package into a "
              "persistent system Python (use a venv or PYTHONPATH).")
        return 0
    print("::error::self-hosted CI job(s) would pip-install a Radia package into "
          "persistent system Python, which drifts the execution host:")
    for wf, jobname, sname, why, line in violations:
        print(f"  {os.path.relpath(wf, _REPO)}  job={jobname}  step='{sname}'")
        print(f"      [{why}]  {line}")
    print("\nFix: use a run-local virtual environment or import the checkout via "
          "PYTHONPATH instead of pip-installing it.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
