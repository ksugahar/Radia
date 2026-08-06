"""
GMSH MCP Server for Radia Project

Provides tools for:
- Inspecting and validating MSH v4.1 files (structure, NodeData, Jacobians)
- Validating .geo launch files (Merge targets, invalid options)
- Writing the shared .geo/.geo.opt/.msh.opt post-display launch artifact
- GMSH visualization and post-processing documentation
- .msh file format reference (v4.1)
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

import errno
import sys
from collections import Counter
from pathlib import Path, PurePosixPath

from mcp.server.fastmcp import FastMCP

from ..common import register_status_tool
from .detect import detect_capabilities
from .gmsh_examples import get_gmsh_examples
from .gmsh_knowledge import get_gmsh_documentation
from .gmsh_reference import get_gmsh_reference
from .msh_inspect import (
    audit_msh_directory,
    diff_msh,
    field_stats,
    inspect_msh,
    mesh_quality,
    probe_options,
    validate_geo,
    validate_msh,
)
from .post_display import (
    build_gmsh_post_display_contract,
    gmsh_post_display_manifest_gate,
    write_gmsh_post_launch_artifact,
)
from .post_process import (
    curve_profile,
    cut_plane_extract,
    derived_field,
    export_view_csv,
    extract_skin,
    field_histogram,
    flux_lines,
    harmonic_to_time,
    integrate_view,
    isosurface,
    line_profile,
    math_eval,
    mirror_expand,
    modulus_phase,
    point_history,
    probe_field,
    resample_grid,
    smooth_to_nodes,
    streamlines,
    streamlines_2d,
    threshold,
    transform_view,
    view_min_max,
    warp_view,
)
from .render import (
    export_animation,
    render_montage,
    render_png,
)
from .rules import ALL_RULES
from .session import (
    session_exec,
    session_run_file,
    session_shutdown,
    session_status,
)
from .verify import verify_artifact

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
    "invalid-gmsh-option": (
        "Replace option names that do not exist in GMSH 4.x: Mesh.Volumes -> "
        "Mesh.VolumeEdges/VolumeFaces, Mesh.Surfaces -> Mesh.SurfaceEdges/"
        "SurfaceFaces, General.GraphicsSizeX/Y -> General.GraphicsWidth/Height."
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
    except OSError as e:
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
    if p.suffix != '.py':
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
            geo, high_order, workflow, onelab, pitfalls, animation,
            paraview (ParaView-filter -> gmsh-tool correspondence
            matrix with measured semantics and honest gaps)

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
def gmsh_inspect_msh(msh_path: str) -> dict:
    """
    Inspect an MSH v4.1 file and summarize its structure.

    Pure-Python parser (no gmsh dependency): reports version, physical
    names, entity counts, node/element counts, element types with order,
    bounding box, NodeData/ElementData views (name, components, time
    steps), plus display hints (e.g. NumSubEdges=4 for order>=2 meshes).

    Args:
        msh_path: Path to the .msh file (absolute or relative to cwd).
    """
    p = Path(msh_path)
    if not p.is_absolute():
        p = PROJECT_ROOT / p
    return inspect_msh(p)


@mcp.tool()
def gmsh_validate_msh(msh_path: str, check_jacobians: bool = False,
                      quadrature: str = "Gauss2") -> dict:
    """
    Validate MSH v4.1 structural consistency; optional Jacobian check.

    Structural checks (pure Python, always available): v4.1 ASCII format,
    balanced sections, header counts vs actual, unique node/element tags,
    element node references exist, known element types, NodeData declared
    counts vs data rows, data tags exist, numComponents in {1,3,9},
    per-view component consistency across time steps.

    With check_jacobians=True the gmsh Python API runs in a SUBPROCESS
    and evaluates getJacobians on all 3D element types (repo policy for
    high-order export verification): negative determinants = inverted
    elements from wrong node ordering. Also returns the integrated
    volume per element type for cross-checks against CAD volume.

    Args:
        msh_path: Path to the .msh file.
        check_jacobians: Run the gmsh getJacobians subprocess check.
        quadrature: Integration rule for the Jacobian check (e.g. Gauss2,
                    Gauss4).
    """
    p = Path(msh_path)
    if not p.is_absolute():
        p = PROJECT_ROOT / p
    return validate_msh(p, check_jacobians=check_jacobians,
                        quadrature=quadrature)


@mcp.tool()
def gmsh_field_stats(msh_path: str, view_name: str | None = None) -> dict:
    """
    Per-view, per-time-step field statistics for an MSH v4.1 file.

    Answers "what is |B| max?" / "did the solve produce NaN?" without
    opening a GUI: scalars report signed min/max/mean/rms per step,
    vectors/tensors report Euclidean-magnitude stats plus pooled
    per-component min/max. NaN/Inf are counted (they render silently
    wrong in GMSH) and excluded from the statistics. Pure Python.

    Args:
        msh_path: Path to the .msh file.
        view_name: Restrict to one view (error lists available names).
    """
    p = Path(msh_path)
    if not p.is_absolute():
        p = PROJECT_ROOT / p
    return field_stats(p, view_name=view_name)


@mcp.tool()
def gmsh_validate_geo(geo_path: str, deep: bool = True) -> dict:
    """
    Validate a .geo launch/companion file BEFORE opening it in GMSH.

    Catches the classic "Open GMSH doesn't work / window is black" bugs:
    missing Merge targets, invalid GMSH 4.x option names (Mesh.Volumes,
    Mesh.Surfaces, General.GraphicsSizeX/Y). With deep=True (default)
    the merged .msh files are scanned for their view count, so
    out-of-range View[N] references (silently ignored by GMSH, field
    stays invisible) are caught, and the exact-autoload sidecars
    (.geo.opt / .msh.opt) are reported. Also reports Mesh.NumSubEdges.

    Args:
        geo_path: Path to the .geo file.
        deep: Scan Merge targets for view counts and sidecars.
    """
    p = Path(geo_path)
    if not p.is_absolute():
        p = PROJECT_ROOT / p
    return validate_geo(p, deep=deep)


@mcp.tool()
def gmsh_audit_msh_directory(directory: str,
                             check_jacobians: bool = False,
                             pattern: str = "**/*.msh") -> dict:
    """
    Validate every .msh under a directory and summarize the health.

    The .msh companion of gmsh_audit_summary (which lints Python
    scripts): one call answers "are the repository's mesh artifacts
    structurally sound?" -- per-file status, failed checks for
    problem files, node/element/view counts, high-order flags.
    check_jacobians=True adds the getJacobians inverted-element gate
    per file (slower; runs gmsh in subprocesses).

    Args:
        directory: Directory to scan (recursively).
        check_jacobians: Also run the per-file Jacobian gate.
        pattern: Glob pattern relative to the directory.
    """
    d = Path(directory)
    if not d.is_absolute():
        d = PROJECT_ROOT / d
    return audit_msh_directory(d, check_jacobians=check_jacobians,
                               pattern=pattern)


@mcp.tool()
def gmsh_diff_msh(msh_a: str, msh_b: str, rel_tol: float = 1e-9) -> dict:
    """
    Compare two MSH v4.1 files: structure + field statistics.

    Built for before/after verification (re-export, node-order fix,
    solver change): reports node/element-type/physical-name/view
    structure differences, bbox drift, and per-common-view relative
    drift of the min/max statistics. identical_structure and
    fields_match give the one-glance verdict; `differences` lists every
    deviation in plain language. Field values are compared through
    statistics, not tag-by-tag. Pure Python.

    Args:
        msh_a: First .msh file (reference).
        msh_b: Second .msh file (candidate).
        rel_tol: Relative tolerance for field min/max drift.
    """
    pa = Path(msh_a)
    if not pa.is_absolute():
        pa = PROJECT_ROOT / pa
    pb = Path(msh_b)
    if not pb.is_absolute():
        pb = PROJECT_ROOT / pb
    return diff_msh(pa, pb, rel_tol=rel_tol)


@mcp.tool()
def gmsh_write_post_launch_artifact(msh_path: str,
                                    output_base: str | None = None,
                                    title: str = "Gmsh post launch artifact",
                                    camera_preset: str = "z_up_xz_from_positive_y",
                                    cut_plane: dict | None = None,
                                    views: list | None = None,
                                    mesh: dict | None = None,
                                    animation: dict | None = None) -> dict:
    """
    Write the shared .geo/.geo.opt/.msh.opt post-display launch artifact.

    File-writing companion to gmsh_post_display_contract: emits case.geo
    (launch target with Merge + display options), the exact autoload
    sidecars case.geo.opt and case.msh.opt, plus case.display.json, and
    returns the gated manifest. The .msh itself is NOT modified.

    Args:
        msh_path: Existing MSH v4.1 data file the artifact points at.
        output_base: Base path for .geo/.opt outputs (default: msh stem).
        title: Human-readable artifact title written into the .geo.
        camera_preset: One of z_up_xz_from_positive_y, positive_y_oblique,
                       front_xz, custom.
        cut_plane: Optional {enabled, normal:[3], offset, whole_elements,
                   only_volume} dict.
        views: Optional list of view dicts ({index, name, kind:
               scalar|vector|displacement, time_step, range, ...}).
        mesh: Optional mesh display overrides ({surface_faces,
              surface_edges, volume_faces, volume_edges, num_sub_edges}).
        animation: Optional {delay, cycle, step, link_time_steps} dict.
    """
    p = Path(msh_path)
    if not p.is_absolute():
        p = PROJECT_ROOT / p
    if not p.is_file():
        return {"ok": False, "error": f"msh file not found: {p}"}
    base = None
    if output_base is not None:
        base = Path(output_base)
        if not base.is_absolute():
            base = PROJECT_ROOT / base
    return write_gmsh_post_launch_artifact(
        p,
        output_base=base,
        title=title,
        camera_preset=camera_preset,
        cut_plane=cut_plane,
        views=views,
        mesh=mesh,
        animation=animation,
    )


@mcp.tool()
def gmsh_render(path: str,
                png_out: str | None = None,
                width: int = 1000, height: int = 800,
                numsubedges: int = 4,
                camera_preset: str | None = None,
                time_step: int | None = None,
                cut_plane: dict | None = None,
                options: dict | None = None,
                string_options: dict | None = None,
                adapt_views: bool = True,
                smooth_normals: bool = True) -> dict:
    """
    Render a .msh or .geo file to PNG headlessly (gmsh subprocess).

    High-order aware by default: Mesh.NumSubEdges=4 renders curved
    mesh edges curved, and every post view gets AdaptVisualizationGrid=1
    so >8-node elements (TET10, HEX20, ...) are not silently skipped.
    Opening a .geo auto-loads its .geo.opt sidecar, so the PNG shows
    exactly what a user double-click would show. Uses -noconfig plus
    explicit window geometry (immune to the stale off-screen-monitor
    window position pitfall). Meshes with no views auto-enable
    SurfaceFaces + physical-group coloring.

    Args:
        path: .msh or .geo file to render.
        png_out: Output PNG path (default: alongside input).
        width: Requested window width (exported PNG can be narrower by
               the FLTK sidebar; the result reports the actual size).
        height: Window height in pixels.
        numsubedges: Subdivisions for curved high-order mesh display.
        camera_preset: Optional preset (z_up_xz_from_positive_y,
                       positive_y_oblique, front_xz).
        time_step: Optional time step applied to all views before render.
        cut_plane: Optional structured cut plane ({enabled, normal:[3],
                   offset, whole_elements, only_volume}); views get
                   Clip=1 automatically.
        options: Numeric Gmsh options applied after opening the input.
        string_options: String-valued Gmsh options applied after open.
        adapt_views: Enable adaptive visualization for high-order views.
        smooth_normals: Average surface normals for smooth shading.
    """
    p = Path(path)
    if not p.is_absolute():
        p = PROJECT_ROOT / p
    out = None
    if png_out is not None:
        out = Path(png_out)
        if not out.is_absolute():
            out = PROJECT_ROOT / out
    return render_png(p, out, width=width, height=height,
                      numsubedges=numsubedges, camera_preset=camera_preset,
                      time_step=time_step, cut_plane=cut_plane,
                      options=options, string_options=string_options,
                      adapt_views=adapt_views,
                      smooth_normals=smooth_normals)


@mcp.tool()
def gmsh_export_animation(path: str,
                          gif_out: str | None = None,
                          keep_frames: bool = False,
                          num_steps: int | None = None,
                          delay_ms: int = 40,
                          width: int = 1000, height: int = 800,
                          camera_preset: str | None = None,
                          orbit_axis: str | None = None,
                          orbit_degrees: float = 360.0,
                          orbit_frames: int = 36,
                          time_step: int | None = None,
                          cut_plane: dict | None = None,
                          options: dict | None = None,
                          string_options: dict | None = None,
                          adapt_views: bool = True,
                          smooth_normals: bool = True,
                          link_views: bool = True) -> dict:
    """
    Export a time-stepped post-view animation as GIF (gmsh subprocess).

    Codifies the lab animation recipe: links all views
    (PostProcessing.Link=1, AnimationCycle=0), steps every view's
    TimeStep explicitly, writes one PNG frame per step, and assembles
    the GIF with Pillow. High-order views get AdaptVisualizationGrid=1
    automatically. Works on a .geo launch artifact (auto-loads
    .geo.opt) or directly on a time-stepped .msh.

    orbit_axis switches to a CAMERA-ORBIT fly-around instead (ParaView
    camera path): the data stays at time_step while the camera sweeps
    orbit_degrees in orbit_frames frames -- works on mesh-only files.

    Args:
        path: .geo or .msh with time-stepped NodeData/ElementData views.
        gif_out: Output GIF path (default: alongside input).
        keep_frames: Keep per-step PNGs in a <gif stem>_frames dir.
        num_steps: Number of steps (default: max NbTimeStep over views).
        delay_ms: Per-frame delay in the GIF.
        width: Requested window width in pixels.
        height: Window height in pixels.
        camera_preset: Optional camera preset (see gmsh_render).
        orbit_axis: x | y | z enables the camera-orbit mode.
        orbit_degrees: Total camera sweep angle.
        orbit_frames: Number of orbit frames.
        time_step: Data step shown during an orbit.
        cut_plane: Optional structured cut plane; views are clipped.
        options: Numeric Gmsh options applied after opening the input.
        string_options: String-valued Gmsh options applied after open.
        adapt_views: Enable adaptive visualization for high-order views.
        smooth_normals: Average surface normals for smooth shading.
        link_views: Step all compatible views together.
    """
    p = Path(path)
    if not p.is_absolute():
        p = PROJECT_ROOT / p
    out = None
    if gif_out is not None:
        out = Path(gif_out)
        if not out.is_absolute():
            out = PROJECT_ROOT / out
    return export_animation(p, out, keep_frames=keep_frames,
                            num_steps=num_steps, delay_ms=delay_ms,
                            width=width, height=height,
                            camera_preset=camera_preset,
                            orbit_axis=orbit_axis,
                            orbit_degrees=orbit_degrees,
                            orbit_frames=orbit_frames,
                            time_step=time_step, cut_plane=cut_plane,
                            options=options,
                            string_options=string_options,
                            adapt_views=adapt_views,
                            smooth_normals=smooth_normals,
                            link_views=link_views)


@mcp.tool()
def gmsh_exec(code: str, timeout_s: float = 120.0) -> dict:
    """
    Execute Python code in a PERSISTENT gmsh session (stateful evaluate).

    matlab-mcp-core-server-style engine session: the first call lazily
    starts a headless gmsh worker subprocess; later calls reuse it, so
    models, options, and views persist (open a big .msh once, then
    interrogate it across many calls). Assign to a variable named
    `result` to return a JSON value; stdout is captured and returned.
    The worker is import-isolated: a crash or hang kills only the
    worker (the call raises loudly) and the next call starts fresh.

    Do NOT use gmsh.fltk here (no GUI in a persistent server-owned
    process) -- use gmsh_render / gmsh_export_animation for screenshots.
    One-shot gating stays with gmsh_inspect_msh / gmsh_validate_msh.

    Example:
        gmsh_exec("gmsh.open(r'C:/models/case.msh')")
        gmsh_exec("result = gmsh.model.getBoundingBox(-1, -1)")
        gmsh_exec("result = [gmsh.view.getTags(),
                             gmsh.option.getNumber('Mesh.NumSubEdges')]")

    Args:
        code: Python source executed in the session globals (gmsh is
              pre-imported and initialized with -noconfig).
        timeout_s: Hard per-call timeout; on expiry the worker is killed
                   and the call fails loudly.
    """
    return session_exec(code, timeout_s=timeout_s)


@mcp.tool()
def gmsh_detect() -> dict:
    """
    Detect gmsh capabilities on this machine (detect_matlab_toolboxes twin).

    One call reports: gmsh/Pillow package presence, gmsh version, build
    features (OpenCASCADE, FLTK, MED, ...), whether an FLTK graphics
    context can ACTUALLY be created (probed in a subprocess -- decides
    upfront if gmsh_render / gmsh_export_animation will work), the
    session state, and the tool-lane map (gating / rendering / session).
    """
    return detect_capabilities()


@mcp.tool()
def gmsh_run_file(path: str, timeout_s: float = 120.0) -> dict:
    """
    Open/run a file in the persistent gmsh session (run_matlab_file twin).

    gmsh.open semantics: a .geo executes as a script (its .geo.opt
    autoloads), .msh/.pos load models and views, .step loads geometry.
    Starts the session lazily like gmsh_exec; returns the post-open
    status (models, current model, view count, bounding box) so the
    follow-up gmsh_exec calls can interrogate the loaded state.

    Args:
        path: File to open in the session (.geo/.msh/.pos/.step).
        timeout_s: Hard per-call timeout.
    """
    p = Path(path)
    if not p.is_absolute():
        p = PROJECT_ROOT / p
    return session_run_file(p, timeout_s=timeout_s)


@mcp.tool()
def gmsh_verify(path: str, check_jacobians: bool = True,
                check_options: bool = False) -> dict:
    """
    Run ALL applicable gates on an artifact (run_matlab_test_file twin).

    The one-call "test runner" for GMSH artifacts: a .msh gets the
    structural + NaN/Inf + (default) Jacobian gates plus its sibling
    .geo deep check; a .geo gets the deep launch check plus every
    merged .msh. Returns a structured pass/fail report (passed /
    failed gate lists, failed checks and first errors per gate).

    Args:
        path: .msh or .geo artifact.
        check_jacobians: Run the gmsh getJacobians gate per .msh.
        check_options: Also probe every .geo option assignment against
                       the gmsh option database.
    """
    p = Path(path)
    if not p.is_absolute():
        p = PROJECT_ROOT / p
    return verify_artifact(p, check_jacobians=check_jacobians,
                           check_options=check_options)


@mcp.tool()
def gmsh_probe(msh_path: str, points: list,
               view_name: str | None = None,
               step: int = -1) -> dict:
    """
    Probe post-processing views at arbitrary points (interpolated).

    "What is B at the gap center?" without a GUI: returns per-view,
    per-point interpolated values (all time steps with step=-1, split
    per step). Points outside the mesh report found=false plus the
    distance to the nearest element.

    Args:
        msh_path: .msh with NodeData/ElementData views.
        points: List of [x, y, z] probe points.
        view_name: Restrict to one view by name (default: all views).
        step: Time step (-1 = all steps).
    """
    p = Path(msh_path)
    if not p.is_absolute():
        p = PROJECT_ROOT / p
    return probe_field(p, points, view=view_name, step=step)


@mcp.tool()
def gmsh_line_profile(msh_path: str, start: list, end: list, n: int = 100,
                      view_name: str | None = None,
                      plot_png: str | None = None) -> dict:
    """
    Sample a view along a straight line; optionally plot a PNG graph.

    The classic "Bz along the axis" post plot in one call: n
    interpolated samples between start and end, distances included,
    and (with plot_png) a matplotlib graph -- scalars plot the value,
    vectors plot the magnitude, one curve per time step.

    Args:
        msh_path: .msh with views.
        start: Line start [x, y, z].
        end: Line end [x, y, z].
        n: Number of samples (>= 2).
        view_name: Restrict to one view (default: first/all).
        plot_png: Optional output PNG path for the graph.
    """
    p = Path(msh_path)
    if not p.is_absolute():
        p = PROJECT_ROOT / p
    out = None
    if plot_png is not None:
        out = Path(plot_png)
        if not out.is_absolute():
            out = PROJECT_ROOT / out
    return line_profile(p, start, end, n, view=view_name, plot_png=out)


@mcp.tool()
def gmsh_integrate(msh_path: str, view_name: str | None = None,
                   dimension: int = 3) -> dict:
    """
    Integrate a view over its elements (per time step).

    Total loss, total flux, stored energy: Plugin(Integrate) with the
    dimension pinned (default 3 = volume elements only). CAUTION:
    dimension=-1 SUMS integrals of every element dimension present
    (volume + surface + line) -- verified behavior, rarely what a
    physical quantity means. Accuracy note (measured): the plugin
    integrates the view at piecewise-LINEAR accuracy even on
    high-order elements, so nonlinear integrands carry O(h^2) error;
    exact FE integrals belong to NGSolve Integrate on the solver side.

    Args:
        msh_path: .msh with views.
        view_name: View to integrate (default: first).
        dimension: Element dimension to integrate over (3, 2, 1; -1=all).
    """
    p = Path(msh_path)
    if not p.is_absolute():
        p = PROJECT_ROOT / p
    return integrate_view(p, view=view_name, dimension=dimension)


@mcp.tool()
def gmsh_math_eval(msh_path: str, expressions: list,
                   view_name: str | None = None,
                   other_view_name: str | None = None,
                   result_name: str = "math_eval",
                   out_file: str | None = None) -> dict:
    """
    Create a derived view with Plugin(MathEval) and save it.

    |B| from components, differences between two solutions, scaled
    fields: expressions use v0..v8 (view components) and w0..w8
    (other view). NOTE: expressions apply to NODAL values and the
    result interpolates f(node values) -- standard FEM post semantics
    (interp(T^2) is the mean of squared vertex values, not
    (interp T)^2). Output defaults to <stem>_math.pos next to the
    input; feed it back to gmsh_probe / gmsh_field_stats / gmsh_render.

    Args:
        msh_path: Input .msh with views.
        expressions: 1-9 expressions (one per output component), e.g.
                     ["Sqrt(v0^2+v1^2+v2^2)"] or ["v0-w0"].
        view_name: Source view (default: first).
        other_view_name: Optional second view bound to w0..w8.
        result_name: Name of the generated view.
        out_file: Output .pos/.msh path (default: <stem>_math.pos).
    """
    p = Path(msh_path)
    if not p.is_absolute():
        p = PROJECT_ROOT / p
    out = None
    if out_file is not None:
        out = Path(out_file)
        if not out.is_absolute():
            out = PROJECT_ROOT / out
    return math_eval(p, expressions, view=view_name,
                     other_view=other_view_name, result_name=result_name,
                     out_file=out)


@mcp.tool()
def gmsh_isosurface(msh_path: str, value: float,
                    view_name: str | None = None,
                    recur_level: int = 0, target_error: float = 1e-4,
                    out_file: str | None = None) -> dict:
    """
    Extract the isosurface of a scalar view (e.g. the saturation front).

    Plugin(Isosurface): returns the piece counts (triangles on volume
    elements, lines on surface elements) and saves the extracted
    surface. recur_level > 0 enables ADAPTIVE extraction on high-order
    data (order-2 GmshPostExport output): elements subdivide along the
    actual high-order interpolant, so the surface follows the curved
    field instead of the P1 chord (measured: radial error 0.21 ->
    0.008 at level 4 on a quadratic field). Render the result with
    gmsh_render (smooth shading is on by default) and use
    View[i].ColormapAlpha < 1 for nested transparent surfaces.

    Args:
        msh_path: Input .msh with a scalar view.
        value: Iso value.
        view_name: Source view (default: first).
        recur_level: Adaptive subdivision depth, 0..6 (0 = plain P1 cut).
        target_error: Adaptive refinement error target.
        out_file: Output path (default: <stem>_iso.pos).
    """
    p = Path(msh_path)
    if not p.is_absolute():
        p = PROJECT_ROOT / p
    out = None
    if out_file is not None:
        out = Path(out_file)
        if not out.is_absolute():
            out = PROJECT_ROOT / out
    return isosurface(p, value, view=view_name, recur_level=recur_level,
                      target_error=target_error, out_file=out)


@mcp.tool()
def gmsh_cut_plane_extract(msh_path: str, normal: list, offset: float,
                           view_name: str | None = None,
                           out_file: str | None = None) -> dict:
    """
    Cut a view with the plane A*x+B*y+C*z+D=0 and save the section DATA.

    Unlike the render-time clip (visual only), Plugin(CutPlane)
    extracts the section as data: probe it, integrate it, render it,
    or stats it downstream.

    Args:
        msh_path: Input .msh with views.
        normal: Plane normal [A, B, C].
        offset: Plane offset D.
        view_name: Source view (default: first).
        out_file: Output path (default: <stem>_cut.pos).
    """
    p = Path(msh_path)
    if not p.is_absolute():
        p = PROJECT_ROOT / p
    out = None
    if out_file is not None:
        out = Path(out_file)
        if not out.is_absolute():
            out = PROJECT_ROOT / out
    return cut_plane_extract(p, normal, offset, view=view_name,
                             out_file=out)


@mcp.tool()
def gmsh_harmonic_to_time(msh_path: str, n_steps: int = 20,
                          view_name: str | None = None,
                          real_step: int = 0, imag_step: int = 1,
                          out_file: str | None = None) -> dict:
    """
    Expand a complex (re/im) view into a time-domain animation view.

    AC field post: v(t_k) = re*cos(2*pi*k/n) - im*sin(2*pi*k/n)
    (Plugin HarmonicToTime). The output feeds gmsh_export_animation
    directly -- eddy-current phasor solutions become rotating-field
    GIFs in two calls.

    Args:
        msh_path: .msh whose view carries re/im as two time steps.
        n_steps: Number of time samples over one period.
        view_name: Source view (default: first).
        real_step: Step index holding the real part.
        imag_step: Step index holding the imaginary part.
        out_file: Output path (default: <stem>_time.pos).
    """
    p = Path(msh_path)
    if not p.is_absolute():
        p = PROJECT_ROOT / p
    out = None
    if out_file is not None:
        out = Path(out_file)
        if not out.is_absolute():
            out = PROJECT_ROOT / out
    return harmonic_to_time(p, view=view_name, real_step=real_step,
                            imag_step=imag_step, n_steps=n_steps,
                            out_file=out)


@mcp.tool()
def gmsh_streamlines(msh_path: str, seed_start: list, seed_end: list,
                     n_seeds: int = 10, view_name: str | None = None,
                     step_size: float | None = None, max_steps: int = 400,
                     both_directions: bool = True,
                     adaptive: bool = True, closure: bool = True,
                     max_turn_deg: float = 10.0, arrows_every: int = 0,
                     return_points: bool = False,
                     out_file: str | None = None) -> dict:
    """
    Trace field lines of a vector view from seeds on a line segment.

    Probe-driven arc-length RK4 with curvature ADAPTIVITY (step_size is
    the MAXIMUM step; it halves wherever the turn per step exceeds
    max_turn_deg) and CLOSED-LOOP detection: a magnetic field line that
    returns to its seed closes exactly instead of overdrawing or
    stopping mid-loop. One merged polyline per seed (backward+forward),
    |v| as line color, per-line termination reasons (left_data | closed
    | stagnation | max_steps) for solution debugging, and optional
    direction arrows as a companion VP view. (This gmsh build's
    Plugin(StreamLines) only re-emits seed points.)

    Args:
        msh_path: .msh with a vector view.
        seed_start: Seed segment start [x, y, z].
        seed_end: Seed segment end [x, y, z].
        n_seeds: Number of seeds along the segment.
        view_name: Vector view (default: first).
        step_size: Maximum arc-length step (default: bbox diag / 200).
        max_steps: Max integration steps per direction.
        both_directions: Trace backward as well as forward.
        adaptive: Halve the step where the line turns too fast.
        closure: Detect and exactly close returning loops.
        max_turn_deg: Adaptive turn-angle bound per step.
        arrows_every: Emit a direction arrow every k-th point (0=off).
        return_points: Include polyline coordinates in the result.
        out_file: Output path (default: <stem>_stream.pos).
    """
    p = Path(msh_path)
    if not p.is_absolute():
        p = PROJECT_ROOT / p
    out = None
    if out_file is not None:
        out = Path(out_file)
        if not out.is_absolute():
            out = PROJECT_ROOT / out
    return streamlines(p, seed_start, seed_end, n_seeds=n_seeds,
                       view=view_name, step_size=step_size,
                       max_steps=max_steps, both_directions=both_directions,
                       adaptive=adaptive, closure=closure,
                       max_turn_deg=max_turn_deg, arrows_every=arrows_every,
                       return_points=return_points, out_file=out)


def _abs_path(path_str: str | None) -> Path | None:
    """Resolve a possibly-relative tool path against PROJECT_ROOT."""
    if path_str is None:
        return None
    p = Path(path_str)
    return p if p.is_absolute() else PROJECT_ROOT / p


@mcp.tool()
def gmsh_flux_lines(msh_path: str, n_levels: int = 20,
                    levels: list | None = None,
                    view_name: str | None = None,
                    recur_level: int = 0, target_error: float = 1e-4,
                    out_file: str | None = None) -> dict:
    """
    Equally spaced isolines of a scalar view, merged into ONE view.

    THE 2D-magnetics field-line tool: on a planar A_z view (or the
    axisymmetric flux function psi = r*A_theta) the equally spaced
    levels ARE the exact field lines with EQUAL FLUX between adjacent
    lines -- the FEMM-style motor flux plot with physically correct
    density, no integration, no seeds, no tuning. On a 3D scalar the
    same call stacks isosurfaces at the given levels (recur_level > 0
    for smooth adaptive extraction on order-2 data; render with
    View[i].ColormapAlpha < 1 for nested transparency).

    Args:
        msh_path: .msh/.pos with the scalar view (A_z, psi, |B|, T...).
        n_levels: Number of interior levels between view min and max.
        levels: Explicit level values (overrides n_levels).
        view_name: Source view (default: first).
        recur_level: Adaptive subdivision depth, 0..6 (0 = plain P1 cut).
        target_error: Adaptive refinement error target.
        out_file: Output path (default: <stem>_flux.pos).
    """
    return flux_lines(_abs_path(msh_path), n_levels=n_levels,
                      levels=levels, view=view_name,
                      recur_level=recur_level, target_error=target_error,
                      out_file=_abs_path(out_file))


@mcp.tool()
def gmsh_streamlines_2d(msh_path: str, origin: list, u_point: list,
                        v_point: list, d_sep: float | None = None,
                        view_name: str | None = None,
                        step_size: float | None = None,
                        first_seed: list | None = None,
                        max_lines: int = 200,
                        arrows_every: int = 0,
                        return_points: bool = False,
                        out_file: str | None = None) -> dict:
    """
    Evenly spaced streamlines on a plane slice (Jobard-Lefer placement).

    The uniform-density streamline picture ParaView only offers for
    native-2D datasets, here on ANY plane cut of a 3D field: the
    in-plane projection of the vector view is traced, new seeds spawn
    automatically d_sep away from accepted lines, and lines stop at
    d_sep/2 from their neighbors -- no manual seeding, no bunching,
    closed loops detected. Exact field lines on symmetry planes
    (B.n = 0); elsewhere the standard projected-field portrait. The
    patch follows the resample_grid convention (origin + u/v edge
    endpoints).

    Args:
        msh_path: .msh/.pos with a vector view.
        origin: Patch corner [x, y, z].
        u_point: End of the patch U edge.
        v_point: End of the patch V edge.
        d_sep: Line separation (default: patch diagonal / 30).
        view_name: Vector view (default: first).
        step_size: Integration step (default: d_sep / 4).
        first_seed: Optional [u, v] in-plane start seed.
        max_lines: Cap on the number of lines.
        arrows_every: Emit a direction arrow every k-th point (0=off).
        return_points: Include line coordinates in the result.
        out_file: Output path (default: <stem>_stream2d.pos).
    """
    return streamlines_2d(_abs_path(msh_path), origin, u_point, v_point,
                          d_sep=d_sep, view=view_name,
                          step_size=step_size, first_seed=first_seed,
                          max_lines=max_lines, arrows_every=arrows_every,
                          return_points=return_points,
                          out_file=_abs_path(out_file))


@mcp.tool()
def gmsh_derived_field(msh_path: str, operation: str,
                       view_name: str | None = None,
                       check_point: list | None = None,
                       out_file: str | None = None) -> dict:
    """
    Derive gradient, curl, divergence, or tensor eigenvalues of a view.

    The ParaView Gradient-filter analog: grad(phi) for equipotential
    checks, curl(A) = B, div(B) = 0 sanity maps, and min/mid/max
    eigenvalues of a Maxwell stress tensor view (written as three
    views in one file). Values are exact derivatives of the P1
    interpolant, element-wise constant, discontinuous across elements.

    Args:
        msh_path: .msh/.pos with the source view.
        operation: gradient | curl | divergence | eigenvalues.
        view_name: Source view (default: first).
        check_point: Optional [x,y,z]; probes the result there.
        out_file: Output path (default: <stem>_<operation>.pos).
    """
    return derived_field(_abs_path(msh_path), operation, view=view_name,
                         check_point=check_point,
                         out_file=_abs_path(out_file))


@mcp.tool()
def gmsh_threshold(msh_path: str, min_val: float, max_val: float,
                   view_name: str | None = None, time_step: int = 0,
                   dimension: int = -1,
                   out_file: str | None = None) -> dict:
    """
    Keep only elements whose MEAN value lies in [min_val, max_val].

    The ParaView Threshold analog (Plugin ExtractElements): isolate
    saturated iron (|B| > 1.8 T), loss hot spots, or any value band as
    a new view. Selection uses the ELEMENT MEAN at time_step (measured
    semantics) -- threshold |v| of vector fields via gmsh_math_eval
    first.

    Args:
        msh_path: .msh/.pos with a scalar view.
        min_val: Lower bound of the kept band.
        max_val: Upper bound of the kept band.
        view_name: Source view (default: first).
        time_step: Step whose values select elements.
        dimension: Restrict to one element dimension (-1 = all).
        out_file: Output path (default: <stem>_thresh.pos).
    """
    return threshold(_abs_path(msh_path), min_val, max_val, view=view_name,
                     time_step=time_step, dimension=dimension,
                     out_file=_abs_path(out_file))


@mcp.tool()
def gmsh_extract_skin(msh_path: str, view_name: str | None = None,
                      from_mesh: bool = False,
                      out_file: str | None = None) -> dict:
    """
    Extract the boundary skin of a volume view (surface + field).

    The ParaView ExtractSurface analog: the boundary triangles/quads
    of the view's volume elements with the field interpolated on them
    -- surface |B| maps from volume solutions without re-export.

    Args:
        msh_path: .msh/.pos with a volume view.
        view_name: Source view (default: first).
        from_mesh: Skin the model mesh instead of the view data.
        out_file: Output path (default: <stem>_skin.pos).
    """
    return extract_skin(_abs_path(msh_path), view=view_name,
                        from_mesh=from_mesh, out_file=_abs_path(out_file))


@mcp.tool()
def gmsh_mirror_expand(msh_path: str, planes: list, parity: str = "scalar",
                       view_name: str | None = None,
                       origin: list | None = None,
                       result_name: str = "mirrored",
                       out_file: str | None = None) -> dict:
    """
    Expand a half/quarter/eighth symmetric model view by mirroring.

    The ParaView Reflect analog with the field physics done right:
    every subset of planes (["x"], ["x","y"], ...) adds a mirrored
    copy and ALL copies merge into one view. parity controls the data
    transform under a mirror M: "scalar" copies values; "vector"
    (polar: A, J, force) maps v' = M v; "pseudovector" (axial: B and
    H!) maps v' = det(M) M v. Element orientation is repaired, so a
    quarter magnet model becomes the full-field picture in one call.

    Args:
        msh_path: .msh/.pos with the symmetric-sector view.
        planes: Mirror planes, subset of ["x","y","z"].
        parity: scalar | vector | pseudovector (B/H are pseudovectors).
        view_name: Source view (default: first).
        origin: Mirror-plane crossing point (default: [0,0,0]).
        result_name: Name of the merged output view.
        out_file: Output path (default: <stem>_full.pos).
    """
    return mirror_expand(_abs_path(msh_path), planes, parity=parity,
                         view=view_name, origin=origin,
                         result_name=result_name,
                         out_file=_abs_path(out_file))


@mcp.tool()
def gmsh_transform_view(msh_path: str, matrix: list,
                        translation: list | None = None,
                        view_name: str | None = None,
                        value_expressions: list | None = None,
                        swap_orientation: bool = False,
                        out_file: str | None = None) -> dict:
    """
    Apply an affine transform x' = A x + t to a COPY of a view.

    The ParaView Transform analog: move a rotor view by its rotation
    angle, offset a coil view for assembly pictures. matrix is
    row-major 3x3. Plugin(Transform) does NOT rotate the data
    components -- pass value_expressions (v0..v8) to rewrite vectors
    consistently, and swap_orientation=True when det(A) < 0.

    Args:
        msh_path: .msh/.pos with the source view.
        matrix: Row-major 3x3 matrix (9 numbers).
        translation: [tx, ty, tz] (default zero).
        view_name: Source view (default: first).
        value_expressions: Optional data rewrite during the copy.
        swap_orientation: Repair element orientation (det < 0).
        out_file: Output path (default: <stem>_xform.pos).
    """
    return transform_view(_abs_path(msh_path), matrix,
                          translation=translation, view=view_name,
                          value_expressions=value_expressions,
                          swap_orientation=swap_orientation,
                          out_file=_abs_path(out_file))


@mcp.tool()
def gmsh_warp(msh_path: str, factor: float = 1.0,
              view_name: str | None = None, time_step: int = 0,
              out_file: str | None = None) -> dict:
    """
    Displace a vector view's geometry by factor * its own vectors.

    The ParaView WarpByVector analog: exaggerated deformation display
    for displacement or force-density fields (magnetostriction,
    magnet-pull visualization). The output is a self-contained .pos
    at the displaced coordinates.

    Args:
        msh_path: .msh/.pos with a vector view.
        factor: Displacement scale factor.
        view_name: Source view (default: first).
        time_step: Step supplying the displacement vectors.
        out_file: Output path (default: <stem>_warp.pos).
    """
    return warp_view(_abs_path(msh_path), factor, view=view_name,
                     time_step=time_step, out_file=_abs_path(out_file))


@mcp.tool()
def gmsh_smooth_to_nodes(msh_path: str, view_name: str | None = None,
                         out_file: str | None = None) -> dict:
    """
    Average element-wise data to nodes (CellDataToPointData analog).

    Each node receives the mean of its adjacent elements' values
    (measured: elements 10/20 -> shared node 15). Run this on
    per-element views (loss density, |J| per cell) before probing,
    isosurfacing, or contour plots, which need nodal continuity.

    Args:
        msh_path: .msh with an ElementData view.
        view_name: Source view (default: first).
        out_file: Output path (default: <stem>_nodal.pos).
    """
    return smooth_to_nodes(_abs_path(msh_path), view=view_name,
                           out_file=_abs_path(out_file))


@mcp.tool()
def gmsh_view_min_max(msh_path: str, view_name: str | None = None,
                      over_time: bool = False) -> dict:
    """
    Locate the min and max of a scalar view WITH their coordinates.

    "Where is the hottest point / peak |B|?" in one call: returns
    {"min": {point, values}, "max": {point, values}} per time step
    (Plugin MinMax with Argument=1). For vector views build |v| with
    gmsh_math_eval first.

    Args:
        msh_path: .msh/.pos with a scalar view.
        view_name: Source view (default: first).
        over_time: Reduce over all time steps as well.
    """
    return view_min_max(_abs_path(msh_path), view=view_name,
                        over_time=over_time)


@mcp.tool()
def gmsh_modulus_phase(msh_path: str, view_name: str | None = None,
                       real_step: int = 0, imag_step: int = 1,
                       out_file: str | None = None) -> dict:
    """
    Convert a complex re/im two-step view to modulus and phase steps.

    AC post companion to gmsh_harmonic_to_time: step 0 becomes
    sqrt(re^2 + im^2), step 1 becomes atan2(im, re) -- amplitude and
    phase maps of eddy-current phasor solutions in one call.

    Args:
        msh_path: .msh whose view carries re/im as two time steps.
        view_name: Source view (default: first).
        real_step: Step index holding the real part.
        imag_step: Step index holding the imaginary part.
        out_file: Output path (default: <stem>_modphase.pos).
    """
    return modulus_phase(_abs_path(msh_path), view=view_name,
                         real_step=real_step, imag_step=imag_step,
                         out_file=_abs_path(out_file))


@mcp.tool()
def gmsh_curve_profile(msh_path: str, x_expr: str, y_expr: str,
                       z_expr: str, u_min: float, u_max: float,
                       n: int = 100, view_name: str | None = None,
                       plot_png: str | None = None,
                       csv_out: str | None = None,
                       out_file: str | None = None) -> dict:
    """
    Sample a view along a parametric curve x(u), y(u), z(u).

    gmsh_line_profile generalized to curves (Plugin CutParametric,
    MathEx syntax). THE air-gap tool: B(theta) on a circle is
    x_expr="0.05*Cos(u)", y_expr="0.05*Sin(u)", z_expr="0", u in
    [0, 2*Pi]. Returns u, points, per-step values; optionally writes
    a PNG graph, a CSV, and an SL line view for rendering.

    Args:
        msh_path: .msh/.pos with the source view.
        x_expr: MathEx expression for x(u).
        y_expr: MathEx expression for y(u).
        z_expr: MathEx expression for z(u).
        u_min: Parameter start.
        u_max: Parameter end.
        n: Number of samples.
        view_name: Source view (default: first).
        plot_png: Optional PNG graph path (value / |v| vs u).
        csv_out: Optional CSV path.
        out_file: Optional SL line view output (.pos).
    """
    return curve_profile(_abs_path(msh_path), x_expr, y_expr, z_expr,
                         u_min, u_max, n=n, view=view_name,
                         plot_png=_abs_path(plot_png),
                         csv_out=_abs_path(csv_out),
                         out_file=_abs_path(out_file))


@mcp.tool()
def gmsh_resample_grid(msh_path: str, origin: list, u_point: list,
                       v_point: list, w_point: list, nu: int, nv: int,
                       nw: int, view_name: str | None = None,
                       csv_out: str | None = None,
                       out_file: str | None = None) -> dict:
    """
    Resample a view on a regular grid spanned by three box edges.

    The ParaView ResampleToImage analog (Plugin CutBox): origin plus
    the endpoints of the U/V/W edges, nu x nv x nw samples (W varies
    fastest). Turns any unstructured result into a uniform grid for
    numpy/MATLAB via the CSV export.

    Args:
        msh_path: .msh/.pos with the source view.
        origin: Box origin [x,y,z].
        u_point: End of the U edge.
        v_point: End of the V edge.
        w_point: End of the W edge.
        nu: Samples along U.
        nv: Samples along V.
        nw: Samples along W.
        view_name: Source view (default: first).
        csv_out: Optional CSV path.
        out_file: Optional point-view output (.pos).
    """
    return resample_grid(_abs_path(msh_path), origin, u_point, v_point,
                         w_point, nu, nv, nw, view=view_name,
                         csv_out=_abs_path(csv_out),
                         out_file=_abs_path(out_file))


@mcp.tool()
def gmsh_export_csv(msh_path: str, csv_out: str,
                    view_name: str | None = None,
                    kind: str = "auto") -> dict:
    """
    Dump view data to CSV (the ParaView spreadsheet / SaveData analog).

    Pure Python on the radia-mcp MSH v4.1 parser -- no gmsh needed.
    kind="nodes" writes tag,x,y,z + one column per view step/component
    from $NodeData; kind="elements" writes element CENTROIDS +
    $ElementData columns; "auto" prefers nodes. List-based .pos files
    carry no node table -- use gmsh_resample_grid/gmsh_curve_profile
    with csv_out for those.

    Args:
        msh_path: .msh with NodeData/ElementData sections.
        csv_out: Output CSV path.
        view_name: Restrict to one view (default: all).
        kind: auto | nodes | elements.
    """
    return export_view_csv(_abs_path(msh_path), _abs_path(csv_out),
                           view=view_name, kind=kind)


@mcp.tool()
def gmsh_field_histogram(msh_path: str, view_name: str | None = None,
                         step: int | None = None,
                         component: int | None = None, bins: int = 32,
                         value_range: list | None = None,
                         plot_png: str | None = None) -> dict:
    """
    Histogram of a view's values (the ParaView Histogram analog).

    Value-distribution checks without a GUI: how much of the iron
    sits above 1.8 T, is the loss density long-tailed. Scalars bin
    the value, vectors bin |v| unless component selects one; step=None
    pools all time steps. Pure Python (no gmsh); optional PNG chart.

    Args:
        msh_path: .msh with NodeData/ElementData sections.
        view_name: Source view (default: first found).
        step: Time step (None = pool all steps).
        component: Component index (None = magnitude for vectors).
        bins: Number of bins.
        value_range: [lo, hi] bin range (default: data min/max).
        plot_png: Optional PNG bar chart path.
    """
    return field_histogram(_abs_path(msh_path), view=view_name, step=step,
                           component=component, bins=bins,
                           value_range=value_range,
                           plot_png=_abs_path(plot_png))


@mcp.tool()
def gmsh_point_history(msh_path: str, point: list,
                       view_name: str | None = None,
                       plot_png: str | None = None) -> dict:
    """
    Value of a view at one point across ALL time steps.

    The ParaView PlotDataOverTime analog for a probe point: transient
    solutions and gmsh_harmonic_to_time outputs become per-step value
    lists plus the recorded step TIMES ($NodeData headers), with an
    optional value-vs-time PNG.

    Args:
        msh_path: .msh with a multi-step view.
        point: Probe point [x, y, z].
        view_name: Source view (default: first).
        plot_png: Optional PNG graph path.
    """
    return point_history(_abs_path(msh_path), point, view=view_name,
                         plot_png=_abs_path(plot_png))


@mcp.tool()
def gmsh_render_montage(images: list, out_png: str,
                        cols: int | None = None,
                        labels: list | None = None) -> dict:
    """
    Compose rendered PNGs into one comparison grid (side-by-side).

    The ParaView comparative-views analog for static output:
    before/after, per-frequency, or per-design renders lined up in a
    single image with optional per-cell labels.

    Args:
        images: PNG paths (at least 2).
        out_png: Output montage path.
        cols: Grid columns (default: 3, or 2 for 4 images).
        labels: Optional per-image labels.
    """
    return render_montage([_abs_path(p) for p in images],
                          _abs_path(out_png), cols=cols, labels=labels)


@mcp.tool()
def gmsh_mesh_quality(msh_path: str, threshold: float = 0.1,
                      quadrature: str = "Gauss4") -> dict:
    """
    Gmsh minSICN shape-quality distribution for all 3D elements.

    Complements the sign-only Jacobian gate of gmsh_validate_msh: a
    non-inverted affine or high-order element can still be nearly
    degenerate. Gmsh's signed inverse condition number detects both
    aspect-ratio and shape degradation. The sampled min(detJ)/max(detJ)
    ratio is also reported as a separate curvature diagnostic.
    ok=True only when no element is inverted or below the threshold.

    Args:
        msh_path: Path to the .msh file.
        threshold: Minimum acceptable minSICN quality (default 0.1).
        quadrature: Integration rule for sampling detJ (default Gauss4).
    """
    p = Path(msh_path)
    if not p.is_absolute():
        p = PROJECT_ROOT / p
    return mesh_quality(p, threshold=threshold, quadrature=quadrature)


@mcp.tool()
def gmsh_probe_options(names: list) -> dict:
    """
    Ask gmsh ITSELF whether option names exist (subprocess probe).

    Complements the static invalid-option lint: any typo or removed
    option is caught, existing options report kind (number/string/
    color) and default value. View[N].X is normalized to the View.X
    template. ok=True only when every requested name exists.

    Args:
        names: Option names to verify (e.g. ["Mesh.NumSubEdges",
               "View[0].Visible"]).
    """
    return probe_options([str(n) for n in names])


@mcp.tool()
def gmsh_session_status() -> dict:
    """
    Report the persistent gmsh session state WITHOUT starting one.

    Returns running flag, worker pid, gmsh/python versions, uptime,
    call count, open model names, current model, and view count.
    """
    return session_status()


@mcp.tool()
def gmsh_session_shutdown() -> dict:
    """
    Shut down the persistent gmsh session (idempotent).

    Graceful shutdown request first; the worker is killed only if it
    does not exit within the timeout. The next gmsh_exec starts a
    fresh session.
    """
    return session_shutdown()


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
    # server.py -> gmsh -> radia_mcp -> src -> radia-mcp (package root)
    fixtures_dir = (
        Path(__file__).parents[3] / "tests" / "mcp_server" / "fixtures"
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

    # --- MSH v4.1 inspect/validate smoke (inline, no fixture needed) ---
    import tempfile
    tiny_msh = (
        "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n"
        "$PhysicalNames\n1\n3 1 \"block\"\n$EndPhysicalNames\n"
        "$Entities\n0 0 0 1\n1 0 0 0 1 1 1 1 1 0\n$EndEntities\n"
        "$Nodes\n1 4 1 4\n3 1 0 4\n1\n2\n3\n4\n"
        "0 0 0\n1 0 0\n0 1 0\n0 0 1\n$EndNodes\n"
        "$Elements\n1 1 1 1\n3 1 4 1\n1 1 2 3 4\n$EndElements\n"
        "$NodeData\n1\n\"T\"\n1\n0.0\n3\n0\n1\n4\n"
        "1 1.0\n2 2.0\n3 3.0\n4 4.0\n$EndNodeData\n"
    )
    with tempfile.TemporaryDirectory() as td:
        msh = Path(td) / "selftest.msh"
        msh.write_text(tiny_msh, encoding="utf-8")
        info = inspect_msh(msh)
        vres = validate_msh(msh)
        geo = Path(td) / "selftest.geo"
        geo.write_text('Merge "selftest.msh";\nMesh.NumSubEdges = 4;\n',
                       encoding="utf-8")
        gres = validate_geo(geo)
        print(f"  msh smoke: inspect ok={info['ok']}, "
              f"validate={vres['status']}, geo={gres['status']}")
        if not (info["ok"] and vres["ok"] and gres["ok"]):
            print("  FAIL: inline MSH inspect/validate smoke failed")
            for err in vres.get("errors", []) + gres.get("errors", []):
                print(f"    {err}")
            sys.exit(1)
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

    # Structural .msh audit over the same durable lanes (no Jacobians:
    # keep the repo audit light and gmsh-free).
    msh_total = 0
    msh_issues = 0
    for audit_dir in existing_audit_dirs:
        audit = audit_msh_directory(audit_dir)
        if not audit.get("ok"):
            continue
        msh_total += audit["files_scanned"]
        for issue in audit["issues"]:
            msh_issues += 1
            print(f"  [MSH] {audit_dir.name}/{issue['path']}: "
                  f"{issue['failed_checks']}")

    print()
    print(f"Scanned: {total} py files, {msh_total} msh files")
    print(f"GMSH issues: {issues} lint, {msh_issues} msh")
    print("PASSED")




register_status_tool(
    mcp,
    server_name='mcp-server-gmsh',
    description='GMSH MSH v4.1 inspect/validate + ParaView-parity post '
                'verbs (derive/threshold/mirror/resample/CSV) + '
                'post-display launch artifacts + visualization policy lint',
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
