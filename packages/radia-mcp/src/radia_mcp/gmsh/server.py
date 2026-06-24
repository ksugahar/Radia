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
    mcp-server-gmsh --selftest --audit-examples
                                  # Run self-test plus repo-wide examples audit
"""

import os
import errno
import sys
from pathlib import Path

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

from .rules import ALL_RULES
from .gmsh_knowledge import get_gmsh_documentation
from .gmsh_reference import get_gmsh_reference
from .gmsh_examples import get_gmsh_examples

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

def _selftest(audit_examples: bool = False):
    """Run fixture self-test; optionally run the repo-wide examples audit."""
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

    # --- Examples audit ---
    examples_dir = PROJECT_ROOT / "examples"
    if not examples_dir.exists():
        if not fixtures_dir.exists():
            print(f"Examples directory not found: {examples_dir}")
        print("PASSED")
        return

    if not audit_examples:
        print("  examples audit: SKIPPED (run --selftest --audit-examples)")
        print("PASSED")
        return

    py_files = sorted(examples_dir.rglob("*.py"))
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
            _selftest(audit_examples='--audit-examples' in sys.argv[1:])
        except (BrokenPipeError, OSError) as exc:
            if _is_closed_stdout_error(exc):
                return
            raise
    else:
        mcp.run(transport="stdio")


if __name__ == '__main__':
    main()
