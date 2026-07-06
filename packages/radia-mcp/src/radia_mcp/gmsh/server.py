"""
GMSH MCP Server for Radia Project

Provides tools for:
- GMSH visualization and post-processing documentation
- .msh file format reference (v2.2, v4.1)
- Linting Python scripts for GMSH policy violations
- High-order element display guidance

In the Radia project, GMSH is used for VISUALIZATION ONLY,
not mesh generation. Mesh generation uses Netgen or Cubit.

Usage:
    mcp-server-gmsh              # Start MCP server (stdio transport)
    mcp-server-gmsh --selftest   # Run self-test
    mcp-server-gmsh --selftest --audit-repo
                                  # Run self-test plus durable repo-lane audit
"""

from collections import Counter
import os
import errno
import sys
from pathlib import Path, PurePosixPath

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

from .rules import ALL_RULES
from .gmsh_knowledge import get_gmsh_documentation
from .gmsh_reference import get_gmsh_reference
from .gmsh_examples import get_gmsh_examples
from .post_display import (
    build_gmsh_post_display_contract,
    gmsh_post_display_manifest_gate,
)

_RULE_REMEDIATIONS = {
    "numsubedges-missing": (
        "For high-order curved display, add a companion .geo with "
        "Mesh.NumSubEdges = 4 or launch gmsh with -numsubedges 4."
    ),
    "gmsh-mesh-generation": (
        "Keep GMSH as visualization/post-processing only. Generate meshes with "
        "Netgen or Cubit, then open/export files for display."
    ),
    "pip-gmsh-import": (
        "Keep the gmsh Python API out of computation/mesh-generation scripts. "
        "Use the `gmsh` launcher on PATH for viewer/inspection workflows, and "
        "prefer .msh v4.1 data plus a .geo launch file for Radia post."
    ),
}

# Create MCP server
mcp = FastMCP("gmsh-lint")

# Resolve relative paths against current working directory
PROJECT_ROOT = Path.cwd()


def _lint_file(filepath: str) -> list[dict]:
    """Run all lint rules on a single file."""
    try:
        with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
    except (OSError, IOError) as e:
        return [{'line': 0, 'severity': 'ERROR', 'rule': 'read-error',
                 'message': f'Cannot read file: {e}'}]

    findings = []
    for rule_fn in ALL_RULES:
        findings.extend(rule_fn(filepath, lines))

    severity_order = {'CRITICAL': 0, 'HIGH': 1, 'MODERATE': 2,
                      'LOW': 3, 'INFO': 4, 'ERROR': -1}
    findings.sort(key=lambda f: (severity_order.get(f['severity'], 9),
                                 f['line']))
    return findings


def _format_findings(filepath: str, findings: list[dict]) -> str:
    """Format findings for display."""
    if not findings:
        return f"[OK] {filepath}: No issues found."

    lines = [f"[{len(findings)} issue(s)] {filepath}:"]
    for f in findings:
        lines.append(
            f"  L{f['line']:>4d} [{f['severity']}] {f['rule']}: {f['message']}"
        )
    return '\n'.join(lines)


def _lint_directory_summary(directory: str = "examples", top_n: int = 10) -> dict:
    """Machine-readable directory lint summary for loop/audit bookkeeping."""
    d = Path(directory)
    if not d.is_absolute():
        d = PROJECT_ROOT / d
    if not d.exists():
        return {
            "ok": False,
            "error": f"Directory not found: {d}",
            "directory": str(d),
        }

    py_files = sorted(d.rglob("*.py"))
    limit = max(0, min(int(top_n), 50))
    by_severity = Counter({"CRITICAL": 0, "HIGH": 0, "MODERATE": 0, "LOW": 0})
    by_rule: Counter[str] = Counter()
    top_files = []
    total_findings = 0

    for py_file in py_files:
        findings = _lint_file(str(py_file))
        if not findings:
            continue
        total_findings += len(findings)
        try:
            rel_path = py_file.relative_to(PROJECT_ROOT).as_posix()
        except ValueError:
            rel_path = str(py_file)
        top_files.append({"path": rel_path, "findings": len(findings)})
        for finding in findings:
            by_severity[str(finding.get("severity", "UNKNOWN"))] += 1
            by_rule[str(finding.get("rule", "unknown"))] += 1

    top_rules = [
        {
            "rule": rule,
            "count": count,
            "action": _RULE_REMEDIATIONS.get(
                rule,
                "Inspect representative findings and add a specific remediation note.",
            ),
        }
        for rule, count in sorted(by_rule.items(), key=lambda item: (-item[1], item[0]))[:limit]
    ]

    return {
        "ok": True,
        "directory": str(d),
        "files_scanned": len(py_files),
        "files_with_findings": len(top_files),
        "total_findings": total_findings,
        "clean": total_findings == 0,
        "by_severity": dict(by_severity),
        "top_rules": top_rules,
        "dominant_rule": top_rules[0] if top_rules else None,
        "top_files": sorted(top_files, key=lambda item: (-item["findings"], item["path"]))[:limit],
    }


def _relative_to_project(path: Path) -> str:
    try:
        return path.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def _numsubedges_triggers(lines: list[str]) -> list[str]:
    triggers = []
    if any("mesh.Curve(" in line or "Curve(" in line.split("#")[0] for line in lines):
        triggers.append("high_order_curve")
    if any("GmshPostExport" in line for line in lines):
        triggers.append("gmsh_post_export")
    return triggers


def _directory_numsubedges_companion(directory: str) -> dict:
    companion = str(PurePosixPath(directory) / "_gmsh_display.geo")
    return {
        "geo_companion": companion,
        "geo_template": (
            f"// Shared GMSH display companion for {directory}\n"
            "// Use with any high-order .msh output from this directory.\n"
            "Mesh.NumSubEdges = 4;\n"
            "// Merge \"<result>.msh\";\n"
        ),
    }


def _line_excerpt(lines: list[str], line_number: int) -> str:
    if line_number <= 0 or line_number > len(lines):
        return ""
    return lines[line_number - 1].strip()


# ============================================================
# Tools
# ============================================================

@mcp.tool()
def lint_gmsh_script(filepath: str) -> str:
    """
    Lint a single Python script for GMSH policy violations.

    Checks:
    - GMSH Python API used for mesh generation (CRITICAL)
    - GmshBuilder import (removed) (CRITICAL)
    - pip gmsh package import (HIGH)
    - meshio import (removed) (HIGH)
    - MSH version mismatch with NGSolve (HIGH)
    - Missing NumSubEdges for high-order elements (MODERATE)

    Args:
        filepath: Absolute or relative path to the Python file.
    """
    p = Path(filepath)
    if not p.is_absolute():
        p = PROJECT_ROOT / p
    if not p.exists():
        return f"Error: File not found: {p}"
    if not p.suffix == '.py':
        return f"Error: Not a Python file: {p}"

    findings = _lint_file(str(p))
    return _format_findings(str(p), findings)


@mcp.tool()
def lint_gmsh_directory(directory: str = "examples") -> str:
    """
    Lint all Python scripts in a directory for GMSH policy violations.

    Args:
        directory: Directory path relative to project root (default: "examples").
    """
    d = Path(directory)
    if not d.is_absolute():
        d = PROJECT_ROOT / d
    if not d.exists():
        return f"Error: Directory not found: {d}"

    py_files = sorted(d.rglob("*.py"))
    if not py_files:
        return f"No Python files found in {directory}."

    total_findings = 0
    file_results = []
    summary = {'CRITICAL': 0, 'HIGH': 0, 'MODERATE': 0, 'LOW': 0}

    for py_file in py_files:
        findings = _lint_file(str(py_file))
        if findings:
            total_findings += len(findings)
            file_results.append(_format_findings(str(py_file), findings))
            for f in findings:
                sev = f['severity']
                if sev in summary:
                    summary[sev] += 1

    if total_findings == 0:
        return f"[OK] All {len(py_files)} files clean."

    header = (
        f"GMSH Lint: {total_findings} issue(s) in "
        f"{len(file_results)}/{len(py_files)} files\n"
        f"  CRITICAL: {summary['CRITICAL']}  HIGH: {summary['HIGH']}  "
        f"MODERATE: {summary['MODERATE']}  LOW: {summary['LOW']}\n"
    )
    return header + "\n".join(file_results)


@mcp.tool()
def gmsh_audit_summary(directory: str = "examples", top_n: int = 10) -> dict:
    """
    Return a machine-readable GMSH lint audit summary.

    Use this for learning-loop bookkeeping and dashboards. It reports files
    scanned, total findings, severity counts, top rule counts, and top files
    without printing the long per-file audit body.
    """
    return _lint_directory_summary(directory, top_n)


@mcp.tool()
def gmsh_numsubedges_remediation_plan(directory: str = "examples",
                                      limit: int = 20) -> dict:
    """
    List scripts that need high-order GMSH display settings.

    This is the actionable companion to the `numsubedges-missing` audit rule.
    It does not edit files; it returns a bounded list of affected scripts,
    trigger reasons, CLI hints, and a companion `.geo` template.
    """
    d = Path(directory)
    if not d.is_absolute():
        d = PROJECT_ROOT / d
    if not d.exists():
        return {
            "ok": False,
            "error": f"Directory not found: {d}",
            "directory": str(d),
        }

    max_items = max(0, min(int(limit), 200))
    affected = []
    by_directory: Counter[str] = Counter()
    total = 0
    for py_file in sorted(d.rglob("*.py")):
        findings = _lint_file(str(py_file))
        if not any(f.get("rule") == "numsubedges-missing" for f in findings):
            continue
        total += 1
        rel_for_group = _relative_to_project(py_file)
        by_directory[str(PurePosixPath(rel_for_group).parent)] += 1
        if len(affected) >= max_items:
            continue
        try:
            lines = py_file.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            lines = []
        rel = rel_for_group
        geo_name = f"{py_file.stem}_display.geo"
        affected.append({
            "script": rel,
            "triggers": _numsubedges_triggers(lines),
            "cli_hint": "gmsh <result>.msh -numsubedges 4",
            "geo_companion": str(PurePosixPath(rel).with_name(geo_name)),
            "geo_template": (
                f"// Display companion for outputs from {rel}\n"
                "Mesh.NumSubEdges = 4;\n"
                "// Merge \"<result>.msh\";\n"
            ),
        })

    return {
        "ok": True,
        "directory": str(d),
        "rule": "numsubedges-missing",
        "total_affected": total,
        "returned": len(affected),
        "truncated": total > len(affected),
        "action": _RULE_REMEDIATIONS["numsubedges-missing"],
        "directory_groups": [
            {
                "directory": directory,
                "count": count,
                "directory_companion": _directory_numsubedges_companion(directory),
            }
            for directory, count in sorted(
                by_directory.items(),
                key=lambda item: (-item[1], item[0]),
            )
        ],
        "affected": affected,
    }


@mcp.tool()
def gmsh_mesh_generation_remediation_plan(directory: str = "examples",
                                          limit: int = 20) -> dict:
    """
    List scripts that still use GMSH as a mesh generator.

    This is the actionable companion to the `gmsh-mesh-generation` audit rule.
    It reports affected files, line snippets, directory grouping, and a
    public-safe migration hint toward Netgen/Cubit `.vol` mesh generation plus
    GMSH-only visualization.
    """
    d = Path(directory)
    if not d.is_absolute():
        d = PROJECT_ROOT / d
    if not d.exists():
        return {
            "ok": False,
            "error": f"Directory not found: {d}",
            "directory": str(d),
        }

    max_items = max(0, min(int(limit), 200))
    affected = []
    by_directory: Counter[str] = Counter()
    total = 0
    total_findings = 0

    for py_file in sorted(d.rglob("*.py")):
        findings = [
            finding
            for finding in _lint_file(str(py_file))
            if finding.get("rule") == "gmsh-mesh-generation"
        ]
        if not findings:
            continue
        total += 1
        total_findings += len(findings)
        rel = _relative_to_project(py_file)
        by_directory[str(PurePosixPath(rel).parent)] += len(findings)
        if len(affected) >= max_items:
            continue
        try:
            lines = py_file.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            lines = []
        affected.append({
            "script": rel,
            "findings": [
                {
                    "line": finding.get("line", 0),
                    "message": finding.get("message", ""),
                    "snippet": _line_excerpt(lines, int(finding.get("line", 0))),
                }
                for finding in findings
            ],
            "migration_hint": (
                "Replace GMSH geometry/mesh generation with Netgen OCC or "
                "Cubit/Coreform export netgen .vol. Keep GMSH only for opening "
                "or post-processing existing .msh/.geo visualization files."
            ),
            "mesh_output_hint": "Prefer Mesh('model.vol') for NGSolve inputs.",
        })

    return {
        "ok": True,
        "directory": str(d),
        "rule": "gmsh-mesh-generation",
        "total_affected": total,
        "total_findings": total_findings,
        "returned": len(affected),
        "truncated": total > len(affected),
        "action": _RULE_REMEDIATIONS["gmsh-mesh-generation"],
        "directory_groups": [
            {"directory": directory, "findings": count}
            for directory, count in sorted(
                by_directory.items(),
                key=lambda item: (-item[1], item[0]),
            )
        ],
        "affected": affected,
    }


@mcp.tool()
def gmsh_usage(topic: str = "all") -> str:
    """
    Get GMSH documentation for visualization and post-processing.

    Topics: policy, overview, cli, shortcuts, options, msh_format,
            geo, high_order, workflow, onelab, pitfalls

    Args:
        topic: Documentation topic (default: "all").
    """
    return get_gmsh_documentation(topic)


@mcp.tool()
def gmsh_reference(topic: str = "all") -> str:
    """
    Get GMSH technical reference (options, algorithms, fields, formats).

    Topics: algorithms, formats, fields, python_api_field, python_api_postproc,
            mesh_options, transfinite, boolean, extrusion, periodic,
            mesh_commands, plugins, view_options, geometry_options

    Lab policy (2026-04-19): gmsh は Python API 経由の軽量**ポスト処理**専用。
    mesh 生成は Netgen/Cubit が担当。Netgen は Tcl 実装でポストが弱いため、
    gmsh Python API 駆動でポストを埋める構成。
    - `python_api_postproc` — gmsh.view.*, gmsh.plugin.*, I/O（**研究室の本題**）
    - `python_api_field` — gmsh.model.mesh.field.* (reference only, 研究室では不使用)
    - `fields` — .geo syntax (legacy 参照用)

    Args:
        topic: Reference topic (default: "all").
    """
    return get_gmsh_reference(topic)


@mcp.tool()
def gmsh_examples(topic: str = "all") -> str:
    """
    Get GMSH tutorial and example documentation.

    Includes ONELAB tutorials (t1-t21), electric machine models
    (PMSM, IM, SRM), magnet/inductor models, and common patterns.

    Topics: tutorials, machines, patterns, magnets

    Args:
        topic: Example topic (default: "all").
    """
    return get_gmsh_examples(topic)


@mcp.tool()
def gmsh_post_display_contract(msh_path: str,
                               output_base: str | None = None,
                               camera_preset: str = "z_up_xz_from_positive_y") -> dict:
    """
    Return the shared .geo/.geo.opt/.msh.opt contract for Gmsh post artifacts.

    Use this before writing radia-acoustic or Gypsilab-derived display files:
    the raw data stay in one MSH v4.1 file, users open the sibling .geo,
    and the exact sidecars preserve camera, cut-plane, and view state.
    """
    return build_gmsh_post_display_contract(
        msh_path,
        output_base=output_base,
        camera_preset=camera_preset,
        views=[{"index": 0, "name": "primary_post_view", "kind": "scalar"}],
    )


@mcp.tool()
def gmsh_post_display_gate(manifest: dict) -> dict:
    """
    Validate a Radia/Gypsilab Gmsh post-display manifest.

    The gate checks MSH v4.1, .geo launch target, exact .geo.opt/.msh.opt
    autoload naming, Z-up camera metadata, cut-plane metadata, and named views.
    """
    return gmsh_post_display_manifest_gate(manifest)


@mcp.tool()
def get_gmsh_lint_rules() -> str:
    """List all available GMSH lint rules with descriptions."""
    lines = [f"GMSH Lint Rules ({len(ALL_RULES)} rules):", ""]
    for rule_fn in ALL_RULES:
        doc = rule_fn.__doc__ or "No description"
        lines.append(f"  {rule_fn.__name__}: {doc.strip()}")
    return "\n".join(lines)


# ============================================================
# Self-Test
# ============================================================

def _selftest(audit_repo: bool = False):
    """Run fixture self-test; optionally audit durable repository lanes."""
    print("=" * 70)
    print("GMSH Lint Self-Test")
    print("=" * 70)

    # --- Fixtures validation ---
    fixtures_dir = (
        Path(__file__).parent.parent.parent.parent.parent / "tests"
        / "mcp_server" / "fixtures"
    )
    if not fixtures_dir.exists():
        fixtures_dir = Path(__file__).parent / "fixtures"

    if fixtures_dir.exists():
        bad_file = fixtures_dir / "bad_gmsh_script.py"
        clean_file = fixtures_dir / "clean_gmsh_script.py"
        if bad_file.exists():
            findings = _lint_file(str(bad_file))
            gmsh_findings = [f for f in findings
                             if f['rule'].startswith(('gmsh-', 'pip-gmsh',
                                                      'meshio-', 'msh-',
                                                      'numsubedges',
                                                      'readgmsh'))]
            print(f"  bad_gmsh_script.py: {len(gmsh_findings)} finding(s)")
            if not gmsh_findings:
                print("  WARNING: bad_gmsh_script.py has no GMSH findings")
        if clean_file.exists():
            findings = _lint_file(str(clean_file))
            gmsh_findings = [f for f in findings
                             if f['rule'].startswith(('gmsh-', 'pip-gmsh',
                                                      'meshio-', 'msh-',
                                                      'numsubedges',
                                                      'readgmsh'))]
            print(f"  clean_gmsh_script.py: {len(gmsh_findings)} finding(s)")
            if gmsh_findings:
                for f in gmsh_findings:
                    print(f"    L{f['line']} [{f['severity']}] {f['rule']}: {f['message']}")
                print("  FAIL: clean script should have zero GMSH findings")
                sys.exit(1)
        print("  fixture validation: PASSED")
        print()

    audit_dirs = [
        "docs",
        "validation_test",
        "src/radia/panels/samples",
    ]
    existing_audit_dirs = [PROJECT_ROOT / d for d in audit_dirs if (PROJECT_ROOT / d).exists()]
    if not existing_audit_dirs:
        if not fixtures_dir.exists():
            print("No durable audit lanes found")
        print("PASSED")
        return

    if not audit_repo:
        print("  repo audit: SKIPPED (run --selftest --audit-repo)")
        print("PASSED")
        return

    py_files = []
    for audit_dir in existing_audit_dirs:
        py_files.extend(sorted(audit_dir.rglob("*.py")))
    total = 0
    issues = 0

    for py_file in py_files:
        findings = _lint_file(str(py_file))
        total += 1
        gmsh_findings = [f for f in findings
                         if f['rule'].startswith(('gmsh-', 'pip-gmsh',
                                                  'meshio-', 'msh-',
                                                  'numsubedges', 'gmsh-builder'))]
        if gmsh_findings:
            issues += len(gmsh_findings)
            print(_format_findings(str(py_file), gmsh_findings))

    print()
    print(f"Scanned: {total} files")
    print(f"GMSH issues: {issues}")
    print("PASSED")




register_status_tool(
    mcp,
    server_name='mcp-server-gmsh',
    description='GMSH MSH v4.1 inspect/validate/convert/write_node_data',
    subpackage='radia_mcp.gmsh',
    related_servers=["cubit"],
    optional_deps=["gmsh"],
    audit_command="mcp-server-gmsh --selftest --audit-repo",
)


def _is_closed_stdout_error(exc: BaseException) -> bool:
    """Windows raises EINVAL when a downstream PowerShell pipe closes early."""
    return isinstance(exc, BrokenPipeError) or (
        isinstance(exc, OSError)
        and getattr(exc, "errno", None) in {errno.EPIPE, errno.EINVAL}
    )


def main():
    """Entry point for mcp-server-gmsh console script."""
    if '--selftest' in sys.argv[1:]:
        from radia_mcp.common.utf8_stdout import use_utf8_stdout
        use_utf8_stdout()
        try:
            _selftest(audit_repo='--audit-repo' in sys.argv[1:])
        except (BrokenPipeError, OSError) as exc:
            if _is_closed_stdout_error(exc):
                return
            raise
    else:
        mcp.run(transport="stdio")


if __name__ == '__main__':
    main()
