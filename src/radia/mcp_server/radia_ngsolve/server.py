"""
Radia + NGSolve Unified MCP Server

Provides tools for both Radia (MMM/MSC/PEEC C++ core) and NGSolve (FEM/BEM):
- Unified linting (33 rules: Radia API + NGSolve FEM + BEM + PEEC)
- Radia C++ library usage (MMM, MSC, field computation, materials, solver)
- NGSolve FEM usage (22 topics: EM formulations, axisymmetric, materials)
- ngsolve.bem (BEM operators, inductance extraction)
- ngsolve.la (Compact AMS/COCR/ICCG preconditioners)
- Kelvin transformation for open boundary FEM
- md2html converter documentation

App-specific knowledge is in separate servers:
- mcp-server-ih: Induction heating (SIBC, ESIM, Karl iteration)
- mcp-server-cubit: Cubit scripting, mesh export
- mcp-server-gmsh: GMSH post-processing

Usage:
    mcp-server-radia-ngsolve              # Start MCP server (stdio transport)
    mcp-server-radia-ngsolve --selftest   # Run self-test
"""

import os
import sys
from pathlib import Path

from mcp.server.fastmcp import FastMCP

from .rules import ALL_RULES
from .radia_knowledge import get_radia_documentation
from .md2html_knowledge import get_md2html_documentation
from .ngsolve_knowledge import get_ngsolve_documentation
from .sparsesolv_knowledge import get_sparsesolv_documentation
from .kelvin_knowledge import get_kelvin_documentation
from .ngsbem_inductance_knowledge import get_ngsbem_inductance_documentation
from .panel_gui_pitfalls_knowledge import get_panel_gui_pitfalls
from .gmsh_post_spec import get_gmsh_post_spec
from .panel_describer import (
    find_panel_file as _find_panel_file,
    parse_panel_file as _parse_panel_file,
    describe_panel_jp as _describe_panel_jp,
    widget_locations as _widget_locations,
)

# NOTE: induction_heating_knowledge is in mcp-server-ih (not here)

mcp = FastMCP("mcp-server-radia-ngsolve")

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

    severity_order = {'CRITICAL': 0, 'HIGH': 1, 'MODERATE': 2, 'LOW': 3, 'INFO': 4, 'ERROR': -1}
    findings.sort(key=lambda f: (severity_order.get(f['severity'], 9), f['line']))
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


@mcp.tool()
def lint_radia_script(filepath: str) -> str:
    """
    Lint a single Python script for Radia + NGSolve convention violations.

    Radia checks:
    - ObjBckg called with list instead of callable (CRITICAL)
    - Missing UtiDelAll cleanup (HIGH)
    - Removed APIs: FldUnits, FldBatch, old solver params (HIGH)
    - MSC sign convention, eval point (HIGH/MODERATE)

    NGSolve checks:
    - BEM on HDivSurface without .Trace() (CRITICAL)
    - Circular SIBC using jv instead of iv (CRITICAL)
    - HCurl magnetostatics without nograds=True (HIGH)
    - Eddy current FE space missing complex=True (HIGH)
    - EFIE V_LL term with wrong (minus) sign (HIGH)
    - PEEC P/(jw) low-frequency breakdown (HIGH)
    - BDDC preconditioner registered after assembly (MODERATE)
    - Overwriting x/y/z coordinate variables (MODERATE)
    - Direct .vec assignment without .data (MODERATE)
    - 2D OCC geometry without dim=2 (MODERATE)
    - CG on A-Omega saddle-point system (MODERATE)
    - Kelvin domain without bonus_intorder (MODERATE)
    - VectorH1 for electromagnetic fields (MODERATE)
    - PINVIT/LOBPCG without gradient projection (MODERATE)
    - Joule heat missing Conj() for complex fields (MODERATE)
    - PEEC n_seg too low for coupling accuracy (MODERATE)
    - Classical EFIE 1/kappa^2 low-frequency breakdown (MODERATE)
    - BEM GenerateMesh without curvaturesafety (MODERATE)
    - TaskManager with BEM non-determinism (LOW)

    Args:
        filepath: Absolute or relative path to the Python file to check.
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
def lint_radia_directory(directory: str = ".") -> str:
    """
    Lint all Python scripts in a directory for NGSolve convention violations.

    Recursively scans .py files and reports findings grouped by file.

    Args:
        directory: Directory path (default: current directory).
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
    summary_by_severity = {'CRITICAL': 0, 'HIGH': 0, 'MODERATE': 0, 'LOW': 0, 'INFO': 0}

    for py_file in py_files:
        findings = _lint_file(str(py_file))
        if findings:
            total_findings += len(findings)
            rel_path = py_file.relative_to(PROJECT_ROOT) if py_file.is_relative_to(PROJECT_ROOT) else py_file
            file_results.append(_format_findings(str(rel_path), findings))
            for f in findings:
                sev = f['severity']
                if sev in summary_by_severity:
                    summary_by_severity[sev] += 1

    output_parts = [
        f"NGSolve Lint Report: {len(py_files)} files scanned, {total_findings} issues found.",
        "",
        f"Summary: {summary_by_severity['CRITICAL']} CRITICAL, "
        f"{summary_by_severity['HIGH']} HIGH, "
        f"{summary_by_severity['MODERATE']} MODERATE, "
        f"{summary_by_severity['LOW']} LOW",
        "",
    ]

    if file_results:
        output_parts.append("=" * 70)
        output_parts.extend(file_results)
    else:
        output_parts.append("All files passed!")

    return '\n'.join(output_parts)


@mcp.tool()
def get_radia_lint_rules() -> str:
    """
    List all available NGSolve lint rules with descriptions.

    Returns a summary of each rule, its severity, and what it checks for.
    """
    rules_info = [
        {
            'rule': 'ngsolve-missing-trace-bem',
            'severity': 'CRITICAL',
            'description': (
                'BEM operators (LaplaceSL/HelmholtzSL) on HDivSurface require '
                '.Trace() on trial/test functions. Without it, boundary-edge '
                'DOFs get corrupted, causing wildly wrong results.'
            ),
            'fix': 'Use j_trial.Trace()*ds(...) instead of j_trial*ds(...).',
        },
        {
            'rule': 'bessel-jv-not-iv',
            'severity': 'CRITICAL',
            'description': (
                'Circular wire SIBC must use modified Bessel functions iv (I0, I1), '
                'NOT regular jv (J0, J1). jv gives correct R_ac/R_dc but wrong sign '
                'on internal inductance Im(Z).'
            ),
            'fix': 'from scipy.special import jv -> from scipy.special import iv',
        },
        {
            'rule': 'hcurl-missing-nograds',
            'severity': 'HIGH',
            'description': (
                'HCurl space for magnetostatics should use nograds=True to '
                'remove gradient null space. Without it, the curl-curl system '
                'is singular.'
            ),
            'fix': 'Add nograds=True: HCurl(mesh, order=2, nograds=True)',
        },
        {
            'rule': 'eddy-current-missing-complex',
            'severity': 'HIGH',
            'description': (
                'HCurl/H1 space in eddy current context without complex=True. '
                'Frequency-domain analysis requires complex-valued FE spaces.'
            ),
            'fix': 'Add complex=True: HCurl(mesh, order=2, complex=True)',
        },
        {
            'rule': 'efie-v-minus-sign',
            'severity': 'HIGH',
            'description': (
                'EFIE system (Zs*M_LL + jw*mu_0*V_LL)*I = -jw*b must use '
                'POSITIVE sign on V_LL term. A minus sign violates Lenz\'s law.'
            ),
            'fix': 'Change minus to plus: Zs*M_LL + jw*mu_0*V_LL (not minus).',
        },
        {
            'rule': 'peec-p-over-jw',
            'severity': 'HIGH',
            'description': (
                'PEEC Loop-Star P/(jw) causes low-frequency breakdown (40-340% error). '
                'Use reformulated Schur complement or stabilized EFIE.'
            ),
            'fix': 'Precompute P^{-1}@M_LS, multiply by jw. Or use mode="stabilized".',
        },
        {
            'rule': 'ngsolve-precond-after-assemble',
            'severity': 'MODERATE',
            'description': (
                'BDDC Preconditioner must be registered BEFORE .Assemble() '
                'to access element matrices.'
            ),
            'fix': 'Move Preconditioner(a, "bddc") BEFORE a.Assemble().',
        },
        {
            'rule': 'ngsolve-overwrite-xyz',
            'severity': 'MODERATE',
            'description': (
                'Loop variable x/y/z overwrites NGSolve coordinate '
                'CoefficientFunction. After the loop, the variable is a '
                'scalar, not a coordinate.'
            ),
            'fix': 'Use different loop variable: "for xi in ..." instead of "for x in ...".',
        },
        {
            'rule': 'ngsolve-vec-assign',
            'severity': 'MODERATE',
            'description': (
                'Direct .vec = assignment creates symbolic expression, '
                'not evaluated result. Must use .vec.data = to evaluate.'
            ),
            'fix': 'Use gfu.vec.data = ... instead of gfu.vec = ...',
        },
        {
            'rule': 'ngsolve-dim2-occ',
            'severity': 'MODERATE',
            'description': (
                '2D OCC geometry (Rectangle, Face) requires dim=2 parameter '
                'in OCCGeometry(). Without it, a 3D surface mesh is generated.'
            ),
            'fix': 'Add dim=2: OCCGeometry(shape, dim=2)',
        },
        {
            'rule': 'ngsolve-cg-on-saddle-point',
            'severity': 'MODERATE',
            'description': (
                'CG solver on A-Omega mixed formulation (saddle-point system). '
                'The system is indefinite, CG may diverge.'
            ),
            'fix': 'Replace solvers.CG() with solvers.GMRes() or MinRes().',
        },
        {
            'rule': 'ngsolve-kelvin-missing-bonus-intorder',
            'severity': 'MODERATE',
            'description': (
                'Kelvin domain integration without bonus_intorder. The varying '
                'Jacobian requires higher quadrature for accurate results.'
            ),
            'fix': 'Add bonus_intorder=4: dx("Kelvin", bonus_intorder=4)',
        },
        {
            'rule': 'ngsolve-vectorh1-for-em',
            'severity': 'MODERATE',
            'description': (
                'VectorH1 used in electromagnetic context. VectorH1 enforces full '
                'C^0 continuity on ALL components, which is wrong for EM fields.'
            ),
            'fix': 'Replace VectorH1 with HCurl (for E, A) or HDiv (for B, J).',
        },
        {
            'rule': 'ngsolve-pinvit-no-projection',
            'severity': 'MODERATE',
            'description': (
                'PINVIT/LOBPCG eigenvalue solver on HCurl without gradient '
                'projection. Curl-curl null space produces spurious zero eigenvalues.'
            ),
            'fix': 'Build gradient projection via fes.CreateGradient().',
        },
        {
            'rule': 'joule-heat-missing-conj',
            'severity': 'MODERATE',
            'description': (
                'Joule heat computed as InnerProduct(E, E) instead of '
                'InnerProduct(E, Conj(E)). Complex E*E != |E|^2.'
            ),
            'fix': 'Use: 0.5 * sigma * InnerProduct(E, Conj(E)).real',
        },
        {
            'rule': 'peec-low-nseg',
            'severity': 'MODERATE',
            'description': (
                'Circular coil PEEC with n_seg < 32 may give poor coupling accuracy.'
            ),
            'fix': 'Increase n_seg to 64 or higher.',
        },
        {
            'rule': 'classical-efie-breakdown',
            'severity': 'MODERATE',
            'description': (
                'Classical EFIE using 1/kappa^2 has O(kappa^{-2}) '
                'condition number blow-up at low frequency.'
            ),
            'fix': 'Use stabilized EFIE: [A_k, Q_k; Q_k^T, kappa^2*V_k].',
        },
    ]

    lines = ["NGSolve Lint Rules", "=" * 50, ""]
    for r in rules_info:
        lines.append(f"[{r['severity']}] {r['rule']}")
        lines.append(f"  {r['description']}")
        lines.append(f"  Fix: {r['fix']}")
        lines.append("")

    return '\n'.join(lines)


@mcp.tool()
def ngsolve_usage(topic: str = "all") -> str:
    """
    Get NGSolve finite element library usage documentation.

    NGSolve is a high-performance FEM library for electromagnetic simulation.
    This tool provides API patterns, best practices, and common pitfalls
    gathered from official tutorials, documentation, and community forums.

    Sources:
      - https://docu.ngsolve.org/latest/i-tutorials/
      - https://forum.ngsolve.org/

    Args:
        topic: Documentation topic. Options:
            "all"              - Complete documentation
            "overview"         - Installation, workflow, direct solvers
            "spaces"           - FE spaces (H1, HCurl, HDiv, HDivSurface, SurfaceL2)
            "maxwell"          - Maxwell/magnetostatics (A-formulation, BDDC, materials)
            "solvers"          - Direct & iterative solver selection guide
            "preconditioners"  - BDDC, multigrid, Jacobi, AMG configuration
            "bem"              - Boundary element method (ngsolve.bem, LaplaceSL, FEM-BEM coupling)
            "mesh"             - Mesh generation (OCC geometry, STEP import, surface mesh)
            "nonlinear"        - Newton's method for nonlinear problems
            "pitfalls"         - Common mistakes and how to avoid them (40 items)
            "linalg"           - Vector/matrix operations, NumPy interop
            "formulations"     - EM formulations: A, Omega, A-Phi, T-Omega, Kelvin (EMPY)
            "adaptive"         - Adaptive mesh refinement with ZZ error estimator (EMPY)
            "darwin"           - Darwin approximation, Surface Impedance BC, Extended Darwin
            "esim"             - ESIM: nonlinear Zs(H,w) Robin BC for any FEM formulation
            "treecotree"       - Tree-cotree splitting, low-freq stability, field-circuit coupling
            "pml"              - Perfectly Matched Layers for open boundary (full-wave)
            "decomposition"    - Domain decomposition: FETI-DP, BDDC, DFDD, AWE/SSP
            "material"         - Material modeling: anisotropy, BH curves, Fixed-Point method
            "ironloss"         - Iron loss estimation: decomposition, FEM computation, steel grades
            "practical"        - Practical techniques: voltage source, force/torque, rotation, coupling
            "team7"            - TEAM Problem 7: eddy current benchmark (A-formulation, OCC geometry, BDDC/AMS solver)
    """
    return get_ngsolve_documentation(topic)


@mcp.tool()
def sparsesolv(topic: str = "all") -> str:
    """
    Get sparsesolv documentation and code examples (now in ngsolve.la).

    Since v3.1.0, sparsesolv types are unified into ngsolve.la module.
    Import: from ngsolve.la import CompactAMSPreconditioner, COCRSolver, etc.

    Repository: https://github.com/ksugahar/ngsolve-sparsesolv

    Args:
        topic: Documentation topic. Options:
            "all"              - Complete documentation
            "overview"         - Library overview, add-on positioning, features
            "api"              - Python API reference (solvers, preconditioners)
            "examples"         - Usage examples (Poisson, curl-curl, complex, etc.)
            "abmc"             - ABMC ordering: parallel triangular solve optimization
            "compact_ams"      - Compact AMS: theory, benchmarks, COCR solver
            "best_practices"   - Preconditioner selection, complex systems, tips
            "build"            - Build and installation instructions
            "example_poisson"  - Ready-to-run: 2D Poisson with ICCG
            "example_curlcurl" - Ready-to-run: 3D curl-curl with auto-shift IC
            "example_eddy"     - Ready-to-run: Complex eddy current problem
            "example_precond"  - Ready-to-run: IC/SGS with NGSolve CGSolver
            "example_divergence" - Ready-to-run: Divergence detection
            "example_compact_ams" - Ready-to-run: Compact AMS + COCR eddy current
    """
    return get_sparsesolv_documentation(topic)


@mcp.tool()
def kelvin_transformation(topic: str = "all") -> str:
    """
    Get Kelvin transformation documentation for open boundary FEM problems.

    The Kelvin transformation maps an unbounded exterior domain to a bounded
    computational domain, enabling FEM solutions without artificial truncation.

    Args:
        topic: Documentation topic. Options:
            "all"            - Complete documentation
            "overview"       - Mathematical foundation and key principles
            "h_formulation"  - H-field perturbation potential formulation
            "a_formulation"  - Vector potential formulation (coils)
            "3d"             - 3D sphere/solid examples (H1)
            "hcurl_3d"       - 3D HCurl A-formulation (calc_fem_kelvin.py)
            "adaptive"       - Adaptive mesh refinement with Kelvin
            "identify"       - Periodic boundary Identify() best practices
            "tips"           - Common mistakes and performance tips
            "robustness"     - Robustness checklist (mesh copy, material scaling,
                               FreeDofs verification, symmetry models, GND)
            "verification"   - Numerical verification (single-domain approach)
            "periodic_wedge" - 1/n sector (symmetry model) with Periodic BC
    """
    return get_kelvin_documentation(topic)


@mcp.tool()
def radia_usage(topic: str = "all") -> str:
    """
    Get Radia C++ library usage documentation.

    Covers: MMM/MSC solvers, field computation, materials (MatLin, MatSatIsoTab,
    hysteresis), background fields, solver configuration (LU/BiCGSTAB/HACApK),
    NGSolve integration (RadiaField CF), memory management, IMA.

    Args:
        topic: Documentation topic. Options:
            "all"            - Complete documentation
            "overview"       - Architecture, MMM vs MSC vs BEM
            "elements"       - ObjRecMag, ObjHexahedron, ObjTetrahedron, ObjWedge
            "materials"      - MatLin, MatSatIsoTab, hysteresis, permanent magnets
            "solver"         - rad.Solve, SolverConfig, LU/BiCGSTAB/HACApK
            "field"          - rad.Fld, batch evaluation, A field
            "removed_apis"   - Removed APIs reference (FldEnr/ObjDivMag/UtiDmp/...)
            "ngsolve"        - RadiaField CF, netgen_mesh_to_radia
            "background"     - ObjBckg, Biot-Savart source
            "ima"            - Image Method of Analysis
    """
    return get_radia_documentation(topic)


@mcp.tool()
def md2html_usage() -> str:
    """Get md2html converter documentation (MathJax, reference links, styled HTML)."""
    return get_md2html_documentation()


@mcp.tool()
def ngsbem_inductance(topic: str = "all") -> str:
    """
    Get ngsolve.bem boundary element method documentation for inductance extraction.

    ngsolve.bem is NGSolve's native boundary element module. Combined with Cubit
    mesh export and SetGeomInfo, it enables accurate inductance extraction
    on high-order curved surface elements.

    Key workflow:
      Cubit mesh -> export_NGSolveCurvedMesh(cubit, order=N) -> LaplaceSL BEM -> L extraction

    Sources:
      - https://docu.ngsolve.org/latest/how_to/ngsbem.html
      - https://github.com/Weggler/docu-ngsbem/ (stabilized BEM)
      - https://github.com/ksugahar/Radia (cubit_mesh_export / export plugin)

    Args:
        topic: Documentation topic. Options:
            "all"            - Complete documentation
            "overview"       - What is ngsolve.bem, comparison with Radia PEEC
            "api"            - LaplaceSL operator, HDivSurface, matrix extraction
            "cubit_workflow" - Cubit -> SetGeomInfo -> Curve -> BEM pipeline
            "curve_order"    - Curve order convergence study for BEM accuracy
            "stabilized"     - Weggler's stabilized BEM for low-frequency
            "examples"       - Runnable examples (circular loop, Cubit torus)
            "best_practices" - Common pitfalls, validation, performance tips
            "known_limitations" - curvaturesafety, TaskManager, QUAD hang, grad(G) gap
    """
    return get_ngsbem_inductance_documentation(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def new_ngsolve_simulation(description: str, formulation: str = "magnetostatic") -> str:
    """Create a new NGSolve electromagnetic simulation script."""
    return (
        f"Create an NGSolve simulation script: {description}\n"
        f"Formulation: {formulation}\n\n"
        "Follow these conventions:\n"
        "1. Use appropriate FE spaces from the de Rham complex:\n"
        "   - HCurl for vector potential A, electric field E\n"
        "   - HDiv for magnetic flux B, current J\n"
        "   - H1 for scalar potential Phi, temperature T\n"
        "   - Do NOT use VectorH1 for EM fields\n"
        "2. Magnetostatics: use HCurl(mesh, order=2, nograds=True)\n"
        "3. Eddy current: use complex=True on HCurl/H1 spaces\n"
        "4. BDDC: register Preconditioner BEFORE .Assemble()\n"
        "5. Do not overwrite x/y/z variables in loops\n"
        "6. Use .vec.data = for vector assignment\n"
        "7. 2D OCC: OCCGeometry(shape, dim=2)\n"
        "8. Saddle-point systems: use GMRes/MinRes, not CG\n"
        "9. BEM with HDivSurface: use .Trace() on trial/test functions\n"
    )


@mcp.prompt()
def ngsolve_eddy_current(geometry: str) -> str:
    """Set up an NGSolve eddy current / induction heating simulation."""
    return (
        f"Set up an NGSolve eddy current simulation for: {geometry}\n\n"
        "Use the mcp-server-ih ih_sibc tool for SIBC method selection.\n"
        "Key points:\n"
        "1. A-Phi formulation: HCurl(complex=True) * H1(complex=True)\n"
        "2. P_total: use BEM (ScalarBIESIBCSolver), not FEM BND integral\n"
        "3. Thermal: transient theta-scheme with H1 space (real)\n"
        "4. Kelvin transform for open boundary (bonus_intorder=4)\n"
    )


# ============================================================
# MCP Resources
# ============================================================

@mcp.resource("ngsolve://spaces")
def ngsolve_spaces_reference() -> str:
    """NGSolve FE space selection quick reference."""
    return (
        "# NGSolve FE Space Selection\n\n"
        "## de Rham Complex\n"
        "```\n"
        "H1 --grad--> HCurl --curl--> HDiv --div--> L2\n"
        "```\n\n"
        "| Space | Continuity | Use For |\n"
        "|-------|-----------|----------|\n"
        "| H1 | Full C^0 | Scalar potential Phi, temperature T |\n"
        "| HCurl | Tangential | Vector potential A, electric field E |\n"
        "| HDiv | Normal | Magnetic flux B, current density J |\n"
        "| HDivSurface | Normal (surface) | BEM surface currents |\n"
        "| SurfaceL2 | None (surface) | BEM charges |\n"
        "| VectorH1 | Full C^0 (all) | Elasticity (NOT for EM!) |\n\n"
        "## Common Parameters\n"
        "- `order=2`: polynomial order (default 1)\n"
        "- `nograds=True`: remove gradient null space (magnetostatics)\n"
        "- `complex=True`: complex-valued (eddy current, time-harmonic)\n"
        "- `dirichlet='bnd'`: essential BC on named boundary\n"
    )


@mcp.resource("ngsolve://solvers")
def ngsolve_solvers_reference() -> str:
    """NGSolve solver selection quick reference."""
    return (
        "# NGSolve Solver Selection\n\n"
        "## Direct Solvers\n"
        "| Solver | Strengths | When to Use |\n"
        "|--------|-----------|------------|\n"
        "| UMFPACK | Default, robust | N < 50K DOFs |\n"
        "| PARDISO | Fast, parallel | N > 50K, Intel MKL available |\n"
        "| MUMPS | Distributed, out-of-core | Very large, MPI |\n\n"
        "## Iterative Solvers\n"
        "| Solver | System Type | Preconditioner |\n"
        "|--------|------------|----------------|\n"
        "| CG | SPD only | BDDC, Jacobi, AMG |\n"
        "| MinRes | Symmetric indefinite | BDDC |\n"
        "| GMRes | General | Any |\n\n"
        "## Preconditioners\n"
        "- BDDC: domain decomposition, must register BEFORE .Assemble()\n"
        "- local (Jacobi/GS): simple, good for smoothing\n"
        "- multigrid: geometric or algebraic, best for H1\n"
    )


# ============================================================
# Panel Registry Tools
# ============================================================

def _load_panel_registry():
    """Load panel_registry.json from panels directory."""
    import json
    registry_path = (Path(__file__).parent.parent.parent / "panels"
                     / "panel_registry.json")
    if not registry_path.exists():
        return None
    with open(registry_path, "r", encoding="utf-8") as f:
        return json.load(f)


@mcp.tool()
def panel_schema(panel_name: str = "") -> str:
    """
    Show Radia-NGSolve panel definitions with Japanese labels and physics.

    When called without arguments, lists all available panels.
    When called with a panel name, shows detailed parameter definitions
    including Japanese names, physical meaning, CLI flags, and defaults.

    This enables natural language <-> CLI parameter mapping:
      "周波数を50kHzに" -> --frequency 50000
      "銅のワークピース" -> --material copper --sigma 5.8e7

    Args:
        panel_name: Panel ID (e.g. "inductance", "fem_kelvin").
                    Empty string returns overview of all panels.
    """
    reg = _load_panel_registry()
    if reg is None:
        return ("Error: panel_registry.json not found. "
                "Run: python src/radia/panels/sync_registry.py")

    panels = reg.get("panels", {})

    if not panel_name:
        lines = ["# Radia-NGSolve Panels\n"]
        for pid, p in panels.items():
            n = len(p.get("params", []))
            lines.append(f"## {pid}: {p['ja_name']}")
            lines.append(f"  {p['ja_description']}")
            lines.append(f"  Script: {p['script']} | Method: {p['method']}")
            lines.append(f"  Parameters: {n}")
            lines.append("")
        lines.append("Use panel_schema(panel_name) for parameter details.")
        return "\n".join(lines)

    if panel_name not in panels:
        return (f"Unknown panel: {panel_name}. "
                f"Available: {', '.join(panels.keys())}")

    p = panels[panel_name]
    lines = [
        f"# {p['ja_name']} ({panel_name})",
        f"Script: `{p['script']}` | Function: `{p['function']}`",
        f"Method: {p['method']}",
        f"Description: {p['ja_description']}",
        "",
        "## Parameters",
        "",
        "| CLI | 日本語 | Type | Default | Physics |",
        "|-----|--------|------|---------|---------|",
    ]
    for param in p.get("params", []):
        cli = param.get("cli", "")
        ja = param.get("ja", "")
        typ = param.get("type", "str")
        default = param.get("default", "")
        if param.get("required"):
            default = "**required**"
        physics = param.get("physics", "")
        choices = param.get("choices", [])
        if choices:
            physics = f"{physics} [{'/'.join(str(c) for c in choices)}]"
        lines.append(f"| `{cli}` | {ja} | {typ} | {default} | {physics} |")

    if p.get("command_builder"):
        lines.append(f"\nCommand builder: `{p['command_builder']}`")

    return "\n".join(lines)


@mcp.tool()
def panel_add_param(panel_name: str, param_name: str, param_type: str = "float",
                    cli_flag: str = "", default: str = "",
                    ja: str = "", physics: str = "",
                    help_text: str = "") -> str:
    """
    Plan where to add a new parameter to a Radia-NGSolve panel.

    Does NOT modify code. Returns a checklist of files and locations
    that need to be updated, so the LLM can make precise edits.

    Args:
        panel_name: Panel ID (e.g. "fem_kelvin", "inductance")
        param_name: Python parameter name (e.g. "coil_sigma")
        param_type: "float", "int", "str", "bool"
        cli_flag: CLI flag (e.g. "--coil-sigma"). Auto-generated if empty.
        default: Default value as string
        ja: Japanese label (e.g. "コイル導電率")
        physics: Physics description (e.g. "R = L/(sigma*A)")
        help_text: English help text for argparse
    """
    reg = _load_panel_registry()
    if reg is None:
        return "Error: panel_registry.json not found."

    panels = reg.get("panels", {})
    if panel_name not in panels:
        return f"Unknown panel: {panel_name}. Available: {', '.join(panels.keys())}"

    p = panels[panel_name]
    if not cli_flag:
        cli_flag = "--" + param_name.replace("_", "-")

    # Check if param already exists
    existing = [x["cli"] for x in p.get("params", [])]
    if cli_flag in existing:
        return f"Parameter {cli_flag} already exists in {panel_name}."

    script = p["script"]
    function = p["function"]
    builder = p.get("command_builder", "")

    lines = [
        f"# Add `{param_name}` to {panel_name} ({p['ja_name']})",
        f"  日本語: {ja}",
        f"  Physics: {physics}",
        "",
        "## Checklist (4 locations):",
        "",
        f"### 1. `panels/{script}` — argparse",
        f"  Add: `parser.add_argument(\"{cli_flag}\", type={param_type}, "
        f"default={default}, help=\"{help_text}\")`",
        "",
        f"### 2. `panels/{script}` — function `{function}()`",
        f"  Add parameter: `{param_name}: {param_type} = {default}`",
        f"  Wire: `args.{param_name.replace('-', '_')}` -> function call",
        "",
    ]

    if builder:
        mod, method = builder.split(":")
        lines.extend([
            f"### 3. `{mod}` — `{method}()`",
            f"  Add: `cmd += [\"{cli_flag}\", self.val(\"{param_name}\")]`",
            "",
            f"### 4. `{mod}` — widget definition",
            f"  Add QLineEdit/QSpinBox for `{param_name}`",
            f"  Label: \"{ja}\" (displayed in Qt panel)",
            "",
        ])
    else:
        lines.extend([
            "### 3-4. No command builder (standalone script)",
            "",
        ])

    lines.extend([
        "### 5. Update registry",
        f"  Run: `python panels/sync_registry.py`",
        "",
        "### 6. Update MCP knowledge (if physics-relevant)",
        f"  File: mcp_server knowledge related to {panel_name}",
        "",
        "### 7. Avoid the GUI pitfalls",
        "  Before committing, call `panel_gui_pitfalls()` and check",
        "  that the new param does not regress any of the listed",
        "  bugs (combo state save/restore, hidden-widget read in",
        "  build_command, mode-switch widget visibility, GMSH viz,",
        "  subprocess argparse choices, Cubit .jou id capture, ...).",
    ])

    return "\n".join(lines)


@mcp.tool()
def panel_describe_jp(panel_name: str) -> str:
    """
    現在のパネルソースを AST 解析して日本語で詳細に説明する。

    Reads the actual ``radia_<panel_name>.py`` source file (NOT the
    cached panel_registry.json which can be stale) and returns a
    Japanese hierarchical description of:

      - all widgets (key, label, type, default, combo items)
      - mode-switch visibility logic per handler
      - subprocess command builders with their CLI flag mapping

    Use this to:
      1. Confirm what the panel actually looks like before editing
      2. Generate a "spec" for the user to confirm in plain Japanese
      3. Diff against panel_registry.json to find drift

    Args:
        panel_name: Panel id (e.g. "ih", "em", "pcb"). Resolved to
                    ``src/radia/radia_<panel_name>.py``.

    Returns markdown text. Combine with panel_gui_pitfalls() output
    when planning a panel modification — first describe the current
    state, then check the relevant pitfalls.
    """
    path = _find_panel_file(panel_name)
    if path is None:
        return (f"Panel file not found for {panel_name!r}. "
                f"Expected at src/radia/radia_{panel_name}.py.")
    try:
        info = _parse_panel_file(path)
    except SyntaxError as e:
        return f"SyntaxError in {path}:{e.lineno}: {e.msg}"
    return _describe_panel_jp(info)


@mcp.tool()
def panel_widget_locations(panel_name: str, widget_key: str) -> str:
    """
    Return file:line locations for everything that touches a widget.

    For a given widget key (e.g. ``"half_thickness"``), returns:

      - **Definition** location: which add_line/add_combo/add_spin
        call created the widget, with the line number, default
        value, and combo items.
      - **Visibility rules**: every ``self._set_row_visible(key, ...)``
        call across all _on_*_changed handlers, with the conditional
        branch and the visibility expression.
      - **Command builder uses**: every ``cmd += ["--flag",
        self.val("key")]`` line in _build_*_command methods.

    Use this BEFORE editing a widget so you can update every
    location that references it in one consistent commit. The MCP
    output is JSON-pretty so the LLM can structure follow-up
    edits programmatically.

    Args:
        panel_name:  Panel id (e.g. "ih")
        widget_key:  Internal widget key (e.g. "half_thickness",
                     "wp_sigma", "method")
    """
    import json
    path = _find_panel_file(panel_name)
    if path is None:
        return f"Panel file not found for {panel_name!r}."
    try:
        info = _parse_panel_file(path)
    except SyntaxError as e:
        return f"SyntaxError in {path}:{e.lineno}: {e.msg}"
    locs = _widget_locations(info, widget_key)
    return json.dumps(locs, indent=2, ensure_ascii=False)


@mcp.tool()
def gmsh_post_spec() -> str:
    """
    GMSH post-processing specification for Radia panels.

    Returns the SINGLE SOURCE OF TRUTH for what the GMSH output must
    look like: file format, physical groups, mesh curving, display
    options (.msh.opt), NodeData requirements, and the "kirei"
    reference from v3.6.1.

    Read this BEFORE writing any GMSH export code in calc_*.py.
    Every item is mandatory — no exceptions.
    """
    return get_gmsh_post_spec()


@mcp.tool()
def panel_gui_pitfalls(topic: str = "") -> str:
    """
    Pitfalls and lessons learned from Radia GUI / Cubit panel development.

    Read this BEFORE adding a new parameter, mode, or method to a
    `radia_*.py` panel, BEFORE renaming a combo item, and BEFORE
    writing a new sample .jou. Each pitfall is paired with a "rule"
    that prevents it from coming back.

    Topics:
      combo_state             -- save/restore by text, not index
      mode_switch             -- hidden widgets must not feed build_command
      layout_unification      -- shared widget set across solver methods
      gmsh_viz                -- companion .geo, hide volume mesh, vector only
      gmsh_arrow_size         -- ArrowSizeMin/Max=20 — without this the
                                 field arrows are functionally invisible
      subprocess_args         -- calc_*.py choices must match GUI combos
      cubit_jou               -- subtract id semantics, surface id renumbering
      sample_jou              -- one .jou per (panel, method) pair
      silent_action           -- menu actions must produce visible feedback
      silent_except           -- never bare-except; always log type+traceback
                                 tail; always provide a fallback path
      result_keys             -- subprocess result dict is an API contract
      regression_blast_radius -- run BOTH panels after touching shared
                                 helpers; opaque casts (PointId) bite
      panel_qt_testing        -- use tests/panels/test_*_qt.py headless
                                 PySide6 tests as regression guards;
                                 string-grep tests do not catch behaviour
      learn_edition_cap       -- ignore the 50k warning, radia_export bypasses it

    Args:
        topic: Empty for the full document, or one of the topic
               keywords above for a single section.
    """
    return get_panel_gui_pitfalls(topic)


def _selftest():
    """Run lint on built-in fixtures to verify rules work correctly."""
    print("=" * 70)
    print("NGSolve Lint Self-Test")
    print("=" * 70)
    print()

    fixtures_dir = Path(__file__).parent.parent.parent.parent.parent / "tests" / "mcp_server" / "fixtures"
    if not fixtures_dir.exists():
        fixtures_dir = Path(__file__).parent / "fixtures"

    if fixtures_dir.exists():
        print(f"Using fixtures: {fixtures_dir}")
        print()
        py_files = sorted(fixtures_dir.glob("*.py"))
        total_findings = 0
        total_files = 0
        for py_file in py_files:
            findings = _lint_file(str(py_file))
            if findings:
                print(_format_findings(str(py_file), findings))
                print()
            total_findings += len(findings)
            total_files += 1

        print("=" * 70)
        print(f"Summary: {total_findings} finding(s) in {total_files} fixture file(s)")

        # Verify: bad_ngsolve should have findings, clean_ngsolve should not
        ngsolve_bad = fixtures_dir / "bad_ngsolve_script.py"
        if ngsolve_bad.exists():
            findings = _lint_file(str(ngsolve_bad))
            if len(findings) == 0:
                print(f"  WARNING: {ngsolve_bad.name} expected findings but got none")

        ngsolve_clean = fixtures_dir / "clean_ngsolve_script.py"
        if ngsolve_clean.exists():
            findings = _lint_file(str(ngsolve_clean))
            if len(findings) > 0:
                print(f"  FAIL: {ngsolve_clean.name} should be clean but got {len(findings)} finding(s)")
                sys.exit(1)

        print("Self-test PASSED")
    else:
        print("SKIP: No fixtures found.")


def main():
    """Entry point for mcp-server-ngsolve console script."""
    if len(sys.argv) > 1 and sys.argv[1] == '--selftest':
        _selftest()
    else:
        mcp.run(transport="stdio")


if __name__ == '__main__':
    main()
