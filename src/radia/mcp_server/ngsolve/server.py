"""
NGSolve FEM MCP Server

Provides tools for:
- Linting Python scripts against NGSolve conventions (18 rules)
- NGSolve FEM usage documentation (21 topics incl. EM formulations, axisymmetric, material modeling)
- ngsolve.bem (BEM operators: LaplaceSL/DL, HelmholtzSL/DL, Maxwell, BiotSavart)
- ngsolve.la (Compact AMS/COCR/ICCG, formerly ngsolve-sparsesolv) solver documentation
- Kelvin transformation reference for open boundary FEM
- Induction heating workflow (EM -> Joule heat -> transient thermal, 7 topics)
- ngsolve.bem BEM inductance extraction (Cubit -> SetGeomInfo -> Curve -> LaplaceSL, 7 topics)

Usage:
    mcp-server-ngsolve              # Start MCP server (stdio transport)
    mcp-server-ngsolve --selftest   # Run self-test on fixtures
"""

import os
import sys
from pathlib import Path

from mcp.server.fastmcp import FastMCP

from .rules import ALL_RULES
from .ngsolve_knowledge import get_ngsolve_documentation
from .sparsesolv_knowledge import get_sparsesolv_documentation
from .kelvin_knowledge import get_kelvin_documentation
from .induction_heating_knowledge import get_induction_heating_documentation
from .ngsbem_inductance_knowledge import get_ngsbem_inductance_documentation

mcp = FastMCP("ngsolve-lint")

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
def lint_ngsolve_script(filepath: str) -> str:
    """
    Lint a single Python script for NGSolve convention violations.

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
def lint_ngsolve_directory(directory: str = ".") -> str:
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
def get_ngsolve_lint_rules() -> str:
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
            "3d"             - 3D sphere/solid examples
            "adaptive"       - Adaptive mesh refinement with Kelvin
            "identify"       - Periodic boundary Identify() best practices
            "tips"           - Common mistakes and performance tips
    """
    return get_kelvin_documentation(topic)


@mcp.tool()
def induction_heating(topic: str = "all") -> str:
    """
    Get induction heating simulation documentation for NGSolve.

    Complete workflow for electromagnetic induction heating analysis:
    EM eddy current solve (A-Phi) -> Joule heat computation -> transient
    thermal analysis, including Gmsh mesh loading, rotating workpiece,
    VTK output, and post-processing patterns.

    Args:
        topic: Documentation topic. Options:
            "all"           - Complete documentation
            "overview"      - Physics overview, parameters, skin depth
            "gmsh_mesh"     - ReadGmsh loading, physical groups, boundary layers
            "eddy_current"  - A-Phi formulation (production code with Gmsh mesh)
            "thermal"       - Transient heat equation (theta-scheme, material props)
            "rotating"      - Rotating workpiece (mesh.SetDeformation)
            "postprocess"   - VTK output, point/line evaluation, CSV/MAT export
            "pitfalls"      - Common mistakes in induction heating simulation
    """
    return get_induction_heating_documentation(topic)


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
      - https://github.com/ksugahar/Radia (cubit_mesh_export integrated in src/radia/)

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
        "Use the induction_heating tool for workflow documentation.\n"
        "Key points:\n"
        "1. A-Phi formulation: HCurl(complex=True) * H1(complex=True)\n"
        "2. Joule heat: Q = 0.5 * sigma * InnerProduct(E, Conj(E)).real\n"
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


def _selftest():
    """Run lint on built-in fixtures to verify rules work correctly."""
    print("=" * 70)
    print("NGSolve Lint Self-Test")
    print("=" * 70)
    print()

    fixtures_dir = Path(__file__).parent.parent.parent.parent / "tests" / "fixtures"
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
