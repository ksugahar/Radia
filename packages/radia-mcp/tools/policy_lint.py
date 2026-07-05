#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""policy_lint -- enforce the Sugahara-lab PUBLISH BOUNDARY on radia-mcp.

Two complementary guards live here:

A. STRUCTURE guard (``check_packaging``): MCP servers that TARGET / WRAP a
   commercial tool (COMSOL / FEMM / JMAG) must NOT be published -- regardless of
   who wrote the code. "It's my own code" is not an exemption. Catches a wrapper
   wired into ``[project.scripts]`` or shipped in the wheel. It also keeps
   Optuna study/trial operation external to radia-mcp: use the official public
   ``optuna/optuna-mcp`` server, while radia-mcp keeps CAE objectives, design
   spaces, and analysis helpers.

B. PROVENANCE guard (``scan_text_tree``): public artifacts must NOT attribute
   verification to commercial / licensed tools, and must not leak internal
   filesystem paths. Per the lab boundary policy (2026-06-06): cite the open /
   published basis, but never write "vs FEMM <number>", "REAL COMSOL",
   "ONELAB vs JMAG", a tool-attributed benchmark table, any lab-local
   ``S:`` / ``W:`` absolute path, a user-profile ``C:\\Users\\Administrator``
   path, an ``_crossval`` path, etc. A commercial-derived verification
   *number* may REMAIN, but reworded as a neutral "stored regression reference"
   WITHOUT the tool name. This scanner reports ``file:line: <reason>`` so each
   such leak can be scrubbed (number kept, attribution removed).

   ALLOWED (never flagged): open-source tool citations (NGSolve, ONELAB/GetDP,
   Gmsh), published method / author citations (e.g. "Abdel-Razek air-gap
   element", "Dowell", "Hollaus MSFEM"), generic "FEM" / "finite element", and
   formulation names ("A-Omega", "T-Omega"). Only tool-ATTRIBUTED verification
   and internal paths are flagged -- bare capability mentions are fine.  Open
   tools may be cited, but their lab-local clone / NAS paths must be replaced
   by repo-relative or generic corpus descriptions before publication.

This is the executable form of "失敗からも lint に学ぶ": it turns the near-miss of
attributing public results to a commercial benchmark into a permanent,
CI-enforceable guard. Run it before any public commit / push / publish.

Usage:
  python tools/policy_lint.py [--root <radia-mcp dir>] [--strict] [--selftest]
      (default: run BOTH guards on the real package; exit 1 on any finding.)
  python tools/policy_lint.py --root <repo-root> --tracked-only
      (provenance-scan tracked src/examples/tests files from a public repo root;
       untracked private scratch files are ignored.)
  python tools/policy_lint.py --selftest
      (self-test the provenance patterns on crafted strings; no repo needed.)
"""
from __future__ import annotations

import argparse
import ast
import re
import subprocess
import sys
from pathlib import Path

try:
    import tomllib                       # Python 3.11+ (stdlib)
except ModuleNotFoundError:              # Python 3.10: tomllib is not in stdlib
    try:
        import tomli as tomllib          # drop-in backport (CI installs it on 3.10)
    except ModuleNotFoundError:          # provenance scan does not need toml
        tomllib = None

# ===========================================================================
# A. STRUCTURE guard -- packaging / wheel boundary (commercial wrappers)
# ===========================================================================
# Subpackages that WRAP a commercial tool -> must never be published.
HARD_DENY = {
    "comsol_converter",   # COMSOL .mph -> NGSolve/Java converter (named in policy)
    "cst_converter",      # CST .cst -> ngsolve.bem converter (non-public; lives in S:\CST)
}
# Subpackages that CONTAIN commercial-tool content and need a human call on
# whether they belong in the public package (warn, don't hard-fail).
REVIEW = {
    "interop",            # radia<->COMSOL/MATLAB LiveLink interop tips
}
# Public modules whose FILENAME is dedicated to a commercial tool warrant a
# human review for content / bench-number leakage. Filename-based on purpose:
# a free-text token scan flags the legitimate benchmark *mentions* (".mph",
# "JMAG-Designer", "COMSOL multilingual RAG") that pervade an open, FEMM-parity
# tool -- those are allowed; only their content / models / bench numbers are not.
COMMERCIAL_NAME_RE = re.compile(r"(comsol|jmag|cst_)", re.IGNORECASE)
OPTUNA_EXTERNAL_DEPS = {"optuna", "optuna-mcp", "optuna-dashboard"}


def _requirement_name(req: str) -> str:
    """Return normalized distribution name from a PEP 508-ish requirement."""
    return re.split(r"[\[<>=!~;\s]", req, 1)[0].strip().lower().replace("_", "-")


def _subpackage_of(target: str) -> str | None:
    """'radia_mcp.comsol_converter.server:main' -> 'comsol_converter'."""
    mod = target.split(":", 1)[0]
    parts = mod.split(".")
    if len(parts) >= 2 and parts[0] == "radia_mcp":
        return parts[1]
    return None


def _optuna_imports(src: Path) -> list[str]:
    """Find direct imports of Optuna runtime packages inside public source."""
    findings = []
    if not src.exists():
        return findings
    for path in src.rglob("*.py"):
        try:
            tree = ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
        except SyntaxError as exc:
            findings.append(f"{path}: cannot parse Python source for policy lint: {exc}")
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    root = alias.name.split(".", 1)[0].replace("_", "-").lower()
                    if root in OPTUNA_EXTERNAL_DEPS:
                        findings.append(
                            f"{path}:{node.lineno}: imports '{alias.name}'. "
                            "Optuna runtime/MCP operation must stay external to radia-mcp.")
            elif isinstance(node, ast.ImportFrom) and node.module:
                root = node.module.split(".", 1)[0].replace("_", "-").lower()
                if root in OPTUNA_EXTERNAL_DEPS:
                    findings.append(
                        f"{path}:{node.lineno}: imports from '{node.module}'. "
                        "Optuna runtime/MCP operation must stay external to radia-mcp.")
    return findings


def check_packaging(root: Path) -> tuple[list[str], list[str]]:
    """Structure guard: returns (errors, warnings) for the packaging boundary."""
    pyproject = root / "pyproject.toml"
    src = root / "src" / "radia_mcp"
    errors: list[str] = []
    warns: list[str] = []

    if tomllib is None or not pyproject.is_file():
        warns.append("[packaging] pyproject.toml unreadable (tomllib missing or "
                     "file absent) -- structure guard skipped.")
        return errors, warns

    with open(pyproject, "rb") as f:
        cfg = tomllib.load(f)

    scripts = cfg.get("project", {}).get("scripts", {})
    deps = cfg.get("project", {}).get("dependencies", [])
    optional_deps = cfg.get("project", {}).get("optional-dependencies", {})
    find = cfg.get("tool", {}).get("setuptools", {}).get("packages", {}).get("find", {})
    excludes = find.get("exclude", [])

    def _is_excluded(sub: str) -> bool:
        return any(sub in pat for pat in excludes)

    # 1) public console-script entry points must not point at a wrapper
    for name, target in scripts.items():
        sub = _subpackage_of(target)
        if sub in HARD_DENY:
            errors.append(
                f"[scripts] '{name} = {target}' publishes commercial-wrapper "
                f"subpackage '{sub}'. Remove from [project.scripts].")
        elif sub in REVIEW:
            warns.append(
                f"[scripts] '{name}' exposes REVIEW subpackage '{sub}' "
                f"(contains commercial content -- confirm it may be public).")
        if "optuna" in name.lower() or sub == "optuna":
            errors.append(
                f"[scripts] '{name} = {target}' publishes Optuna operation from "
                "radia-mcp. Use the official external optuna/optuna-mcp server "
                "instead; keep CAE objectives/design spaces in radia-mcp.")

    # 1b) radia-mcp must not depend on Optuna runtime packages; users install
    #     optuna-mcp separately when they need Study/Trial/Dashboard operation.
    reqs: list[tuple[str, str]] = [("project.dependencies", r) for r in deps]
    for group, values in optional_deps.items():
        reqs.extend((f"project.optional-dependencies.{group}", r) for r in values)
    for group, req in reqs:
        if _requirement_name(req) in OPTUNA_EXTERNAL_DEPS:
            errors.append(
                f"[dependencies] {group} contains '{req}'. Optuna runtime/MCP "
                "packages are external to radia-mcp; install optuna-mcp "
                "separately.")

    # 1c) the old in-package radia_mcp.optuna server must not return.
    if (src / "optuna").is_dir():
        errors.append(
            "[source] radia_mcp.optuna exists in the public source tree. "
            "Optuna Study/Trial operation belongs to the external "
            "optuna/optuna-mcp package.")
    errors.extend(f"[source] {finding}" for finding in _optuna_imports(src))

    # 2) wrapper subpackages must be excluded from the wheel (packages.find)
    for sub in sorted(HARD_DENY):
        if (src / sub).is_dir():
            if not _is_excluded(sub):
                errors.append(
                    f"[packaging] subpackage '{sub}' is shipped in the wheel "
                    f"(not in [tool.setuptools.packages.find] exclude). Add "
                    f"'radia_mcp.{sub}*' to exclude, or relocate it out of src/.")
            # 3) even if excluded from the wheel, source in the public repo leaks on GitHub
            warns.append(
                f"[source] commercial-wrapper '{sub}' source still lives under "
                f"the public src/ tree -> relocate to a lab-private package for "
                f"full public-GitHub compliance.")

    # 4) dedicated commercial-tool modules in PUBLIC subpackages -> review for
    #    content / bench-number leakage (benchmark *mentions* are fine, so this
    #    is filename-based, not a free-text token scan).
    public_subs = set()
    for target in scripts.values():
        s = _subpackage_of(target)
        if s and s not in HARD_DENY and s not in REVIEW:
            public_subs.add(s)
    for sub in sorted(public_subs):
        d = src / sub
        if not d.is_dir():
            continue
        for py in sorted(d.rglob("*.py")):
            if COMMERCIAL_NAME_RE.search(py.stem):
                warns.append(
                    f"[content] {py.relative_to(root)} is a commercial-tool "
                    f"module in public subpackage '{sub}' -- review for content "
                    f"/ bench-number leakage (benchmark mentions are OK).")

    return errors, warns


# ===========================================================================
# B. PROVENANCE guard -- commercial-tool ATTRIBUTION + internal-path scan
# ===========================================================================
# Each rule is (compiled regex, human reason). A line that matches ANY rule is a
# finding. Rules are deliberately narrow: they fire on tool-ATTRIBUTED
# verification (a commercial tool adjacent to a number / another solver / a
# ground-truth claim) or on an internal filesystem path -- NOT on a bare
# capability mention, an open-source tool, or a published author/method.
_TOOL = r"(?:FEMM|JMAG|COMSOL|CST)"          # commercial / licensed tools
# Tools whose INFLUENCE on a public capability must be hidden ("何から学んだかも隠す").
# FEMM is deliberately EXCLUDED: radia-ngsolve's "FEMM-parity" is the project's STATED open goal
# (allowlisted), so a FEMM design-reference is allowed; COMSOL / JMAG / CST design-references are not.
_REF_TOOL = r"(?:COMSOL|JMAG|CST)"
# A "value" token that signals a benchmark NUMBER (a MEASUREMENT), NOT a
# capability label ("prob1-4big") nor a VERSION ("FEMM 4.2", "COMSOL 6.2").
# A measurement is: signed; or %-suffixed; or a physical unit (T/H/us/...); or a
# decimal with >=2 fractional digits (0.2498, 269.62, 0.78) -- versions have a
# single fractional digit (4.2, 6.2) and are excluded. Must follow a value
# context char (space/=/() so a digit glued to letters ("prob1") does not count.
_NUM = (r"(?<=[\s=(])"
        r"(?:[-+]\d|"                              # signed value
        r"\d+(?:\.\d+)?\s*%|"                      # percent
        r"\d+(?:\.\d+)?\s*(?:[TH]|[µu][Hs])\b|"    # tesla / henry / us / uH
        r"\d+\.\d{2,}|"                            # decimal, >=2 frac (not 4.2)
        r"rad\s*=\s*\d)")

PROVENANCE_RULES: list[tuple[re.Pattern, str]] = [
    # --- tool-attributed VERIFICATION ------------------------------------
    (re.compile(r"\b(?:REAL|real)\s+" + _TOOL + r"\b"),
     "commercial tool named as ground-truth ('REAL <tool>')"),
    (re.compile(r"\bvs\.?\s+" + _TOOL + r"\b"),
     "verification attributed to a commercial tool ('vs <tool>')"),
    (re.compile(r"\b" + _TOOL + r"\b[^\n|]{0,24}?" + _NUM),
     "commercial tool adjacent to a benchmark number"),
    # markdown VALUE-ROW: a table row led by a commercial tool name whose later
    # cell carries a number (e.g. '| COMSOL (A-form) | 0.664661 T | ...').
    (re.compile(r"^\s*\|\s*" + _TOOL + r"\b[^\n]*\|\s*[-+]?\d"),
     "commercial tool leads a benchmark table value-row"),
    (re.compile(r"\b" + _TOOL + r"\b\s*(?:<?-+>|==|<->|<=>)\s*"
                r"(?:NGSolve|analytic|radia)", re.IGNORECASE),
     "tool-to-result agreement attributed to a commercial tool"),
    (re.compile(r"(?:NGSolve|analytic|radia)\s*(?:<?-+>|==|<->|<=>)\s*\b"
                + _TOOL + r"\b", re.IGNORECASE),
     "tool-to-result agreement attributed to a commercial tool"),
    (re.compile(r"\bONELAB\b[^\n]{0,12}?\bvs\.?\b[^\n]{0,12}?\bJMAG\b",
                re.IGNORECASE),
     "open-vs-commercial comparison names JMAG ('ONELAB vs JMAG')"),
    (re.compile(r"\bJMAG\b[^\n]{0,12}?\bvs\.?\b[^\n]{0,12}?\bONELAB\b",
                re.IGNORECASE),
     "open-vs-commercial comparison names JMAG ('JMAG vs ONELAB')"),
    (re.compile(r"\breference\s+commercial\s+(?:code|tool)\b", re.IGNORECASE),
     "result attributed to 'the reference commercial code'"),
    (re.compile(r"\bcommercial\s+(?:code|tool|reference)\b.{0,40}" + _TOOL,
                re.IGNORECASE),
     "commercial tool named as the reference code"),
    (re.compile(_TOOL + r".{0,40}\bcommercial\s+(?:code|reference)\b",
                re.IGNORECASE),
     "commercial tool named as the reference code"),
    (re.compile(_TOOL + r"[^\n]{0,40}\bground[\s-]*truth\b", re.IGNORECASE),
     "commercial tool named as ground truth"),
    (re.compile(r"\bground[\s-]*truth\b[^\n]{0,40}" + _TOOL, re.IGNORECASE),
     "commercial tool named as ground truth"),
    # --- commercial tool named as the DESIGN / CAPABILITY reference -------
    # "公開物では何から学んだかも隠す": a public capability described as "<tool>-modeller",
    # "<tool>-style", "<tool>-equivalent", "like <tool>", "what <tool> calls X", or "<tool>'s X verb"
    # leaks the commercial tool as the SOURCE even with no number. (FEMM-parity is the project's
    # stated open goal, so _REF_TOOL excludes FEMM -- these fire on COMSOL / JMAG / CST only.)
    (re.compile(r"\b" + _REF_TOOL + r"[\s-]modell?(?:er|ing|ed)\b", re.IGNORECASE),
     "commercial tool's modeller named as a reference ('<tool>-modeller / -modelling')"),
    (re.compile(r"\b" + _REF_TOOL + r"[\s-](?:equivalent|style|like)\b", re.IGNORECASE),
     "public capability described as commercial-tool-equivalent ('<tool>-equivalent / -style')"),
    (re.compile(r"\b(?:like|as in|mimic\w*|emulat\w*|modell?ed after|equivalent to|inspired by|"
                r"akin to)\s+" + _REF_TOOL + r"\b", re.IGNORECASE),
     "public capability attributed to a commercial tool as its model ('like <tool>')"),
    (re.compile(r"\bwhat\s+" + _REF_TOOL + r"\s+calls\b", re.IGNORECASE),
     "public capability described by what a commercial tool calls it"),
    (re.compile(r"\b" + _REF_TOOL + r"(?:'s)?\s+\w{2,14}\s+"
                r"(?:verb|feature|operation|primitive)\b", re.IGNORECASE),
     "commercial tool named as the source of a modelling verb / feature"),
    # --- internal / commercial-harness paths & identifiers ----------------
    # Public radia-mcp must not reveal lab-local drive topology.  Open-source
    # tools and literature corpora may be cited by name, but not by S:/ or W:/
    # clone/NAS paths.  Windows installation examples under C:/Program Files
    # and temporary examples under C:/temp are intentionally not matched here.
    (re.compile(r"(?<![A-Za-z])(?:file:///)?[SW]:[\\/]", re.IGNORECASE),
     "lab-local absolute drive path (S:/ or W:/)"),
    (re.compile(r"\\\\(?:\d{1,3}\.){3}\d{1,3}[\\/]", re.IGNORECASE),
     "lab-local UNC path with private network address"),
    (re.compile(r"\b192\.168\.\d{1,3}\.\d{1,3}\b"),
     "private LAN address leaked into public text"),
    (re.compile(r"\b00_CAE\b", re.IGNORECASE),
     "lab-local CAE workspace name leaked into public text"),
    (re.compile(r"\b[SW]:\s*(?:drive|scripts?|share)\b", re.IGNORECASE),
     "bare lab drive topology leaked into public text"),
    (re.compile(r"C:[\\/]Users[\\/]Administrator\b", re.IGNORECASE),
     "user-local absolute profile path"),
    (re.compile(r"<lab-local path>|private local path", re.IGNORECASE),
     "temporary scrub placeholder leaked into public text"),
    (re.compile(r"private\s+correspondence", re.IGNORECASE),
     "private correspondence/provenance note"),
    (re.compile(r"(?<![A-Za-z])[SW]:[\\/](?:FEMM|JMAG|COMSOL|CST)\b", re.IGNORECASE),
     "internal commercial-tool drive path (S:/W: \\<tool>)"),
    (re.compile(r"(?<![A-Za-z])[SW]:[\\/](?:ELF_MAGIC|ELF)\b", re.IGNORECASE),
     "internal ELF/MAGIC drive path"),
    (re.compile(r"(?<![A-Za-z])W:[\\/]00_CAE[\\/](?:FEMM|JMAG|COMSOL|CST)\b", re.IGNORECASE),
     "internal commercial-tool path (W:\\00_CAE\\<tool>)"),
    (re.compile(r"(?<![A-Za-z])W:[\\/]00_CAE[\\/](?:ELF_MAGIC|ELF)\b", re.IGNORECASE),
     "internal ELF/MAGIC path"),
    (re.compile(r"(?<![A-Za-z])[SW]:[\\/][^\n\"']*[\\/](?:FEMM|JMAG|COMSOL|CST)[\\/]",
                re.IGNORECASE),
     "internal path with a commercial-tool directory component"),
    (re.compile(r"(?<![A-Za-z])[SW]:[\\/][^\n\"']*[\\/](?:ELF_MAGIC|ELF)[\\/]",
                re.IGNORECASE),
     "internal path with an ELF/MAGIC directory component"),
    (re.compile(r"\b_crossval\b", re.IGNORECASE),
     "internal cross-validation path token (_crossval)"),
    (re.compile(r"\bcomsol_xval\b", re.IGNORECASE),
     "commercial cross-validation identifier (comsol_xval)"),
    (re.compile(r"[Cc]:[\\/]te?mp[\\/](?:cc_lab|comsol_xval|cc_)", re.IGNORECASE),
     "commercial-tool verification harness path"),
    (re.compile(r"\bcc_lab\b"),
     "commercial LiveLink harness identifier (cc_lab)"),
    (re.compile(r"reference_comsol_rdp_license"),
     "internal memory reference naming a commercial tool"),
]

# Lines containing any of these substrings are EXEMPT (allowlist): published
# repo / package citations and policy commentary that legitimately name a tool
# without attributing verification or leaking a path.
_ALLOW_SUBSTR = (
    "wjc9011/COMSOL_Multiphysics_MCP",      # cited open pattern source (GitHub)
    "COMSOL_Multiphysics_MCP",
    "EED/CST London",                       # quoted paper affiliation (CST = a London
                                            # institution, NOT the CST software) in a
                                            # bibliography abstract -- not an attribution
)

# File globs to scan under each tree.
_SCAN_GLOBS = ("*.py", "*.md", "*.m", "*.wls")
_SCAN_SUFFIXES = tuple(
    glob[1:] for glob in _SCAN_GLOBS if glob.startswith("*.") and len(glob) > 1
)
# Subtrees of the package to scan.
_SCAN_DIRS = ("src", "examples", "tests")
# Files exempt from the provenance scan: this guard's OWN test, which holds
# synthetic boundary-pattern fixtures (e.g. "vs COMSOL .mph") on purpose to
# exercise the STRUCTURE guard -- they are not real attributions.
_EXEMPT_FILES = ("test_policy_lint.py",)


def _line_is_allowed(line: str) -> bool:
    return any(tok in line for tok in _ALLOW_SUBSTR)


def scan_text(text: str, label: str = "<string>") -> list[tuple[str, int, str]]:
    """Scan text; return list of (label, lineno, reason) provenance findings."""
    out: list[tuple[str, int, str]] = []
    for i, line in enumerate(text.splitlines(), 1):
        if _line_is_allowed(line):
            continue
        for rx, reason in PROVENANCE_RULES:
            if rx.search(line):
                out.append((label, i, reason))
                break          # one finding per line is enough to flag it
    return out


def scan_file(path: Path, root: Path | None = None) -> list[tuple[str, int, str]]:
    if path.name in _EXEMPT_FILES:
        return []
    try:
        text = path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        return []
    label = path.relative_to(root).as_posix() if root else path.as_posix()
    return scan_text(text, label)


def _git_tracked_files(root: Path) -> list[Path] | None:
    """Return tracked files under the current git worktree, or None outside git."""
    try:
        top = subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "--show-toplevel"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        out = subprocess.check_output(
            ["git", "-C", top, "ls-files"],
            text=True,
            stderr=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError):
        return None

    top_path = Path(top)
    tracked: list[Path] = []
    try:
        root_resolved = root.resolve()
    except OSError:
        root_resolved = root
    for rel in out.splitlines():
        candidate = (top_path / rel).resolve()
        try:
            candidate.relative_to(root_resolved)
        except ValueError:
            continue
        tracked.append(candidate)
    return tracked


def scan_text_tree(root: Path, *, tracked_only: bool = False) -> list[tuple[str, int, str]]:
    """Scan src/ + examples/ + tests/ files matching _SCAN_GLOBS for findings."""
    findings: list[tuple[str, int, str]] = []
    seen: set[Path] = set()
    try:
        scan_root = root.resolve()
    except OSError:
        scan_root = root
    tracked = _git_tracked_files(root) if tracked_only else None

    if tracked_only and tracked is not None:
        scan_bases = [
            (scan_root / sub).resolve()
            for sub in _SCAN_DIRS
            if (scan_root / sub).is_dir()
        ]
        for f in sorted(tracked):
            if "__pycache__" in f.parts or f.suffix not in _SCAN_SUFFIXES:
                continue
            if not any(_is_relative_to(f, base) for base in scan_bases):
                continue
            if f in seen:
                continue
            seen.add(f)
            findings.extend(scan_file(f, scan_root))
        return findings

    for sub in _SCAN_DIRS:
        base = root / sub
        if not base.is_dir():
            continue
        for glob in _SCAN_GLOBS:
            for f in sorted(base.rglob(glob)):
                if f in seen or "__pycache__" in f.parts:
                    continue
                seen.add(f)
                findings.extend(scan_file(f, root))
    return findings


def _is_relative_to(path: Path, base: Path) -> bool:
    try:
        path.relative_to(base)
    except ValueError:
        return False
    return True


# ===========================================================================
# Self-test (no repo required): assert the provenance scanner's precision.
# ===========================================================================
def selftest() -> int:
    must_flag = [
        "vs FEMM rad=100.mat 0.78% of peak",
        "**+ real COMSOL 3-way (b2 0.0005 %)**",
        "ALL of #5-#12 directly confirmed against REAL COMSOL",
        "| COMSOL  (A-form)                  | 0.664661 T | 0.30 %          |",
        "NGSolve F_z = -4.842 N vs COMSOL -4.930 N (~1.8 %)",
        "=> COMSOL <-> NGSolve agree to 0.11 %.",
        "every one a three-way agreement COMSOL == analytic == NGSolve",
        "comparing ONELAB vs JMAG torque waveforms",
        "Liu Xinyao 2025 JMAG vs ONELAB comparison notes",
        "the lab's JMAG corpus (S:\\JMAG).",
        r"FEMM_DIR = r'W:\00_CAE\FEMM\2020_05_05'",
        "Recipes/harnesses in `S:\\COMSOL\\_crossval\\`",
        "against the COMSOL values in comsol_xval, so the",
        "harness lives in C:\\temp\\cc_lab (lab machine)",
        "trip-wires in [[reference_comsol_rdp_license]]",
        "a convergence study vs CST plus mesh-dependence",       # CST attribution (HF tool)
        "the proposed method matches CST to 0.5 %",
        "validated against REAL CST Studio Suite",
        r"the lab's CST corpus (S:\CST).",                        # CST internal path
        r"private validation mesh at S:\ELF_MAGIC\_crossval\case.meg",
        r"licensed reference path W:\00_CAE\ELF_MAGIC\_crossval\run.json",
        r"YANO_PATH = r'W:\00_CAE\NGSolve\矢野'",
        "Distilled from S:/ONELAB/ElectricMachines/",
        "distilled from s:/local/private/notes",
        r"absorbed from W:\03_文献・論文\00_電磁界解析 corpus",
        r"see S:\cubit\knowledge\license.py",
        r"C:\Users\Administrator\.claude\skills\pdf2ppt-pdfgear",
        r"pip install -e \\192.168.11.100\work\00_CAE\Radia\01_GitHub",
        "ssh 192.168.11.100",
        "self-hosted CI runner cannot access S: drive",
        "the scrub left <lab-local path> in public docs",
        "Dave Meeker, private correspondence",
        # design / capability reference (no number) -- the leak class this guard now catches:
        'The CST-modeller "rotate with copies" verb',            # <tool>-modeller
        "a CST-style near-field equivalent source",              # <tool>-style
        "NGSolve recipes for a COMSOL-equivalent workflow",      # <tool>-equivalent
        "the sweep is modelled after COMSOL's loft tool",        # 'modelled after <tool>'
        "this reproduces what CST calls a Blend",                # 'what <tool> calls'
        "COMSOL's sweep feature is reproduced here",             # "<tool>'s X feature"
        r"FEMM itself as the ground truth",
        "**radia** (0.2512 vs FEMM 0.2498)",
        "validated against ... the reference commercial code",
        r"distilled from S:\FEMM\等価定理の基礎原理\, S:\FEMM\2015_05_21",
        r"S:/COMSOL/2022_05_12_XFEM/An eXtended Finite Element Method",
    ]
    must_not_flag = [
        "# FEMM-parity on NGSolve -- executable, tested",         # bare framing
        "# Magnetics API (FEMM prob1-4big analogs)",              # capability label
        "Cross-check with the virtual work method (Coulomb 1983)",
        "the de-facto standard since Davat-Ren-Lajoie-Mangin (1985)",
        "This is the Arkkio formula (Arkkio 1987 PhD).",
        "an independent open-boundary axisymmetric reference, 0.78% of peak",
        "ONELAB GetDP is the canonical open-source reference",
        "validated to <2 % against a closed-form / exact reference",
        "Hollaus MSFEM multiscale formulation",
        "the A-Omega / T-Omega formulation names",
        "NGSolve reproduces the analytic 2/3 Br to 0.41 %",       # open tool + number
        "Pattern adopted from wjc9011/COMSOL_Multiphysics_MCP (2026-05-24).",
        "see CHANGELOG.md for per-release notes",
        "a temp file under C:\\temp\\foo.py is fine",             # generic C:\temp ok
        "| Capability | FEMM | JMAG | radia-ngsolve | Strongest |",  # capability header
        "Distilled from the ONELAB/GetDP ElectricMachines templates",
        "absorbed from the lab literature corpus",
        "see the Cubit license knowledge module",
        # VERSION strings are not benchmark numbers:
        "in FEMM 4.2 (build dates 2018-2023).",
        "Created in COMSOL Multiphysics 6.2",
        "licensed under the COMSOL Software License Agreement 6.2.",
        "matches FEMM 4.2 production",
        "Created in CST Studio Suite 2026.0",                    # CST version, not a benchmark
        "ELF/MAGIC product documentation server is public-safe when solver data is absent",
        "EMF is at the EED/CST London, UK, e-mail 10016.3130",   # quoted paper affiliation (allowlisted)
        "high-frequency tools (HFSS, FEKO) exist",               # competitor landscape (not in _TOOL)
        "candidate_doi='HTTPS://DOI.ORG/10.1109/tap.2006.880747'",
        "public-safe curated corpus",
        # FEMM design-references are ALLOWED (FEMM-parity is the stated open goal):
        "radia-ngsolve is FEMM-like -- the FEMM-parity goal",    # FEMM excluded from _REF_TOOL
        "the FEMM modeller workflow, reimplemented openly",      # FEMM, not COMSOL/JMAG/CST
        # generic capability prose (no commercial tool named) stays clean:
        "a generic loft / sweep / shell modelling operation",
        "Created in COMSOL Multiphysics 6.2 with the AC/DC module",  # version + module, no feature-ref
    ]
    fails = []
    for s in must_flag:
        if not scan_text(s):
            fails.append(f"MISSED (should flag): {s!r}")
    for s in must_not_flag:
        hits = scan_text(s)
        if hits:
            fails.append(f"FALSE POSITIVE (should NOT flag): {s!r} -> {hits[0][2]}")
    print("=" * 70)
    print("policy_lint provenance self-test")
    print("=" * 70)
    if fails:
        for f in fails:
            print("  [X]", f)
        print(f"\nSELFTEST FAIL ({len(fails)} case(s))")
        return 1
    print(f"  OK -- {len(must_flag)} flagged, {len(must_not_flag)} clean, all correct.")
    return 0


# ===========================================================================
# Entry point
# ===========================================================================
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="radia-mcp publish-boundary lint")
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[1]),
                    help="radia-mcp package dir (contains pyproject.toml)")
    ap.add_argument("--strict", action="store_true",
                    help="promote structure warnings to errors")
    ap.add_argument("--tracked-only", action="store_true",
                    help="provenance-scan only files tracked by git under root")
    ap.add_argument("--selftest", action="store_true",
                    help="run the provenance pattern self-test and exit")
    args = ap.parse_args(argv if argv is not None else None)

    if args.selftest:
        return selftest()

    root = Path(args.root)

    # --- A. structure guard ---
    errors, warns = check_packaging(root)

    # --- B. provenance guard ---
    findings = scan_text_tree(root, tracked_only=args.tracked_only)

    print("=" * 70)
    print("radia-mcp publish-boundary lint (CLAUDE.md 公開境界)")
    print("=" * 70)

    # report provenance findings (file:line: reason)
    suffix = " (git-tracked only)" if args.tracked_only else ""
    print(f"\nPROVENANCE (commercial-tool attribution / internal paths) in "
          f"src+examples+tests{suffix}:")
    if findings:
        for label, lineno, reason in findings:
            print(f"  {label}:{lineno}: {reason}")
        print(f"  -> {len(findings)} finding(s).")
    else:
        print("  OK -- no commercial-tool attribution or internal-path leaks.")

    # report structure findings
    if errors:
        print(f"\nSTRUCTURE ERRORS ({len(errors)}) -- block public publish:")
        for e in errors:
            print(f"  [X] {e}")
    if warns:
        print(f"\nSTRUCTURE WARNINGS ({len(warns)}):")
        for w in warns:
            print(f"  [!] {w}")
    if not errors and not warns:
        print("\nSTRUCTURE: OK -- no packaging-boundary issues.")

    fail = bool(errors) or bool(findings) or (args.strict and bool(warns))
    print("\nRESULT:", "FAIL" if fail else "PASS")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
