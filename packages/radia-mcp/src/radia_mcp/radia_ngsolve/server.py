"""
Radia + NGSolve Unified MCP Server

Provides tools for both Radia (MMM/MSC/PEEC C++ core) and NGSolve (FEM/BEM):
- Unified linting (33 rules: Radia API + NGSolve FEM + BEM + PEEC)
- Radia C++ library usage (MMM, MSC, field computation, materials, solver)
- NGSolve FEM usage (22 topics: EM formulations, axisymmetric, materials)
- ngsolve.bem (BEM operators, inductance extraction)
- radia.sparsesolv_ngsolve (Compact AMS/COCR/ICCG preconditioners)
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
from ..common import register_status_tool

from .rules import ALL_RULES
from .knowledge.radia import get_radia_documentation
from .knowledge.md2html import get_md2html_documentation
from .knowledge.ngsolve import get_ngsolve_documentation
from .knowledge.sparsesolv import get_sparsesolv_documentation
from .knowledge.kelvin import get_kelvin_documentation
from .knowledge.kelvin_identify_post_hoc import (
    get_post_hoc_documentation as _get_kelvin_identify_post_hoc_doc)
from .knowledge.axifem import get_axifem_documentation
from .knowledge.ngsbem_inductance import get_ngsbem_inductance_documentation
from .knowledge.peec_inductance import get_peec_inductance_documentation
from .knowledge.esim import get_esim_documentation
from .knowledge.panel_gui_pitfalls import get_panel_gui_pitfalls
from .knowledge.analytical_formulas import get_analytical_formulas_documentation
from .knowledge.force_validation import get_force_validation_documentation
from .knowledge.install_deploy import get_install_deploy_documentation
from .knowledge.release_workflow import get_release_workflow_documentation
from .knowledge.standalone_panels import get_standalone_panels_documentation
from .knowledge.loop_learning import get_loop_learning_documentation
from .knowledge.basis_functions import get_basis_functions_documentation
from .knowledge.taskmanager import get_taskmanager_knowledge
from .knowledge.cln_sibc_orthogonal import (
    get_cln_sibc_orthogonal_documentation,
    get_cln_sibc_orthogonal_section,
)
from .knowledge.cln_3d import (
    get_cln_3d_documentation,
    get_cln_3d_notebook,
)
from .knowledge.bem_cln import get_bem_cln_documentation
from .knowledge.cln_sphere_dd import get_cln_sphere_dd_documentation
from .knowledge.mmm_core import get_mmm_core_documentation
from .knowledge.hdiv_vim import get_hdiv_vim_documentation
from .knowledge.femm_parity import get_femm_parity_documentation
from .knowledge.fem_bem_schur import get_fem_bem_schur_documentation
from .knowledge.airgap_motor_workflow import get_airgap_motor_workflow_documentation
from .knowledge.dtn_coarse_mesh import get_dtn_coarse_mesh_documentation
from .knowledge.urn import get_urn_documentation, urn_fit_from_csv
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
            'rule': 'scattered-eddy-missing-a0',
            'severity': 'HIGH',
            'description': (
                'joule_loss_density() called WITHOUT A0= in a file that sets a '
                'background/applied field. The FE unknown gfA is the SCATTERED '
                'potential; total E = -jw (A0 + gfA + grad(Phi)), so dropping A0 '
                'overestimates the loss ~10x. (Coil/total-field sources where gfA '
                'is already total correctly use A0=None and are not flagged.)'
            ),
            'fix': 'Pass A0=<background>: joule_loss_density(gfA, gfPhi, sigma, omega, A0=A0).',
        },
        {
            'rule': 'ngsolve-set-definedon-in-loop',
            'severity': 'HIGH',
            'description': (
                'GridFunction.Set(..., definedon=...) called inside a loop. Set() '
                'zeros the whole vector first, so looping over boundaries/materials '
                'leaves only the LAST region (prior ones wiped -> e.g. conductor '
                'potentials collapse to 0, energy 0).'
            ),
            'fix': 'Set all regions in ONE call: gf.Set(mesh.BoundaryCF({...}, default=0), definedon=mesh.Boundaries("a|b")).',
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
def taskmanager(topic: str = "overview") -> str:
    """
    NGSolve TaskManager parallelism — usage, MKL interaction, audit, C++.

    ★ AGENTS.md policy: follow the NGSolve-native execution model. Use
    NGSolve TaskManager for Radia thread-level parallelisation (NOT raw
    OpenMP / std::thread / private pools).  External threaded kernels
    such as MKL dense LU are guarded with SuspendTaskManager.  This tool documents the
    `with TaskManager():` pattern, `SetNumThreads()` + `--nthreads`
    CLI convention, MKL nesting trap, and the C++ `ngcore::ParallelFor`
    equivalent.  Includes a 2026-05-27 audit of the IH panel solver
    scripts (`calc_inductance.py` / `calc_fem_kelvin.py` /
    `calc_fem_coilmesh.py` / `calc_heat*.py`) — all PASS.

    Args:
        topic: One of:
            "overview"          - Policy + what TaskManager does (DEFAULT)
            "usage"             - `with TaskManager():` pattern, when to wrap
            "set_num_threads"   - `SetNumThreads()` + `--nthreads` CLI convention
            "mkl_interaction"   - TaskManager vs MKL thread pool nesting
                                  (alias: "mkl", "pardiso")
            "cpp_kernels"       - `ngcore::ParallelFor` in custom C++ TUs
                                  (alias: "cpp", "parallel_for")
            "audit_radia_ih"    - 2026-05-27 IH panel solver audit
                                  (alias: "audit", "radia_ih")
            "common_mistakes"   - Wrap-outside-Assemble, missing
                                  --nthreads, OMP_NUM_THREADS inherited
                                  from Cubit, etc.
                                  (alias: "mistakes", "pitfalls")
            "all"               - Everything concatenated
    """
    return get_taskmanager_knowledge(topic)


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
            "ngsolve_bem_50"    - .vol visualization + 50-case NGSolve.BEM/MATLAB comparison lane
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
            "multiphysics"     - COMSOL-class couplings: induction heating EM->thermal (joule_loss_density + solve_heat_steady), the scattered-field A0 gotcha
            "cross_validation_registry"
                               - Reusable validation scripts/summary JSONs and the
                                 public-safe MCP knowledge hooks that learned from them
    """
    return get_ngsolve_documentation(topic)


@mcp.tool()
def sparsesolv(topic: str = "all") -> str:
    """
    Get sparsesolv documentation and code examples.

    Since 2026-05-08, sparsesolv ships inside the radia wheel as the
    submodule `radia.sparsesolv_ngsolve` (built from src/ext/sparsesolv/).
    The legacy standalone `ngsolve-sparsesolv` PyPI package is retired.
    Import: from radia.sparsesolv_ngsolve import CompactAMSPreconditioner, COCRSolver, ...

    Source: src/ext/sparsesolv/ in the Radia monorepo.

    See also (theory + decision tree, not code-usage):
        - radia_mcp.matrix_solvers — solver theory + genealogy + which-to-pick
            * matrix_solvers_overview('decision_tree')
            * matrix_solvers_preconditioners('ams_hiptmair_xu') — CompactAMS theory
            * matrix_solvers_krylov('cocg_cocr') — COCR Sogabe-Zhang 2007 origin
            * matrix_solvers_em_specific('eddy_current_stab') — shifted-prec rationale

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
def esim(topic: str = "all") -> str:
    """
    Get ESIM (Effective Surface Impedance Method) general documentation.

    ESIM extends linear SIBC to nonlinear magnetic materials by solving a 1D
    cell problem through the conductor depth.  Returns a field-dependent
    surface impedance Z_s(H_t) for use in BEM/FEM with surface impedance
    boundary conditions.

    This tool documents the GENERAL technique (cell problem mathematics,
    Karl iteration, module API) without coupling to any specific
    application.  For application-specific use of ESIM (induction heating
    workpieces with steel BH curves), see `mcp-server-ih.ih_sibc(topic="esim")`.

    Args:
        topic: Documentation topic. Options:
            "all"             - Complete documentation
            "overview"        - When ESIM vs linear SIBC; nonlinear conductors
            "cell_problem"    - 1D BVP, BCs, geometries (slab/cylinder/finite_slab)
            "karl_iteration"  - Picard relaxation, convergence pitfalls,
                                per-element vs per-node Z_s
            "module_api"      - radia.esim_cell_problem.ESIMFiniteSlabSolver
                                + BEM-SIBC / FEM-SIBC coupling examples
    """
    return get_esim_documentation(topic)


@mcp.tool()
def urn(topic: str = "all") -> str:
    """
    Universal Relaxation Network (URN): causal/passive rational fitting of a
    complex frequency response, with direct time-domain (relaxation-network /
    SPICE / auxiliary-ODE) synthesis.

    URN decomposes Z(omega) (impedance, dispersive eps/mu, or an open-boundary
    DtN symbol G_n(omega)) into a SPARSE sum of physical relaxation mechanisms
    (Debye, Cole-Cole, Cole-Davidson, Havriliak-Negami, CPE, Warburg, Gerischer,
    RLC, skin-effect; series + admittance branch) with KAN-style adaptive tau and
    attention.  Every basis is passive, so the fit is causal/passive BY
    CONSTRUCTION and maps to one first-order auxiliary ODE per pole (fractional
    terms -> short RC/RL ladder) -- the local-in-time operator an FETD /
    Newmark-beta solver needs.  Beats Vector Fitting on fractional/Cole-Cole data
    (avg ~22.8% lower NRMSE on NASA battery + TDK ferrite datasets).

    Use to turn a frequency-domain absorbing-BC / dispersive-layer response into
    a stable broadband time-domain model.  Run the fit with the urn_fit tool.
    Ref: Sugahara & Sato, IEEE Access 2026; impl docs/universal_relaxation_network.

    Args:
        topic: all | overview | method | api | timedomain | application
    """
    return get_urn_documentation(topic)


@mcp.tool()
def urn_fit(data_csv: str, freq_col: int = 0, real_col: int = 1,
            imag_col: int = 2, delimiter: str = ",", skip_rows: int = 0,
            n_debye: int = 3, n_cole_cole: int = 2, n_warburg: int = 1,
            n_cole_davidson: int = 0, sparsity_weight: float = 0.01,
            n_epochs: int = 2000, n_restarts: int = 3, spice_out: str = "") -> str:
    """
    Fit a complex frequency response with a Universal Relaxation Network and
    return the discovered relaxation mechanisms, the fit NRMSE, and a SPICE
    netlist (== the auxiliary-ODE ladder for a time-domain / Newmark-beta solver).

    Input is a CSV with columns (frequency_Hz, Re(Z), Im(Z)).  Requires torch;
    training is iterative -- lower n_epochs / n_restarts for a faster, rougher
    fit (defaults ~2000/3 are a responsive compromise; the paper uses 6000/10).

    Args:
        data_csv: path to a CSV of the frequency response.
        freq_col/real_col/imag_col: 0-based column indices.
        delimiter, skip_rows: CSV parsing.
        n_debye/n_cole_cole/n_warburg/n_cole_davidson: basis counts to try.
        sparsity_weight: L1-ish penalty that prunes unused bases.
        n_epochs/n_restarts: Adam iterations / random restarts.
        spice_out: if set, also write the SPICE netlist to this path.
    """
    return urn_fit_from_csv(
        data_csv, freq_col, real_col, imag_col, delimiter, skip_rows,
        n_debye, n_cole_cole, n_warburg, n_cole_davidson, sparsity_weight,
        n_epochs, n_restarts, spice_out)


@mcp.tool()
def cln_3d(topic: str = "all") -> str:
    """
    Get 3D Cauer Ladder Network (CLN) / Kameari-Tanimoto iteration
    documentation for eddy current analysis with NGSolve.

    Captures Tanimoto-Kameari iterative methods from the master's thesis
    + production code (W:/00_CAE/NGSolve/谷本/, ~25 notebooks):
      - A-T formulation (primary)
      - T-Ω formulation (H1 confined to conductor)
      - A-Φ formulation (HCurl + H1 mixed)
      - Constraint variants: penalty stabilization, explicit Coulomb gauge
      - Production solvers: SparseSolvPy ICCG, accICCG, NGSolve CG, direct

    Each formulation produces a Cauer-II ladder {R_n, L_n} via Schmidt
    orthogonalization on impressed J source. Validated against
    cylindrical TM-mode analytical R/L for n=0..9.

    Open research: Kameari + Kelvin transformation combination (3D
    HCurl A-formulation hits ~25× discrepancy with mpmath BEM Foster
    target due to A_ext gauge being unbounded at infinity; future
    work includes T-Ω with reduced-Ω = -H_0·z + Ω_r).

    Args:
        topic: Documentation topic. Options:
            "all"           - Complete documentation
            "overview"      - Mathematical foundation, three formulations
            "notebooks"     - Index of 修論 / 定式_誤差検証 / 静止器回転機用
            "formulas"      - Cauer-II synthesis, drift diagnostic,
                              bonus_intorder critical setting
    """
    full_doc = get_cln_3d_documentation()
    if topic == "all":
        return full_doc
    sections = {
        "overview": "OVERVIEW",
        "notebooks": "NOTEBOOK_INDEX",
        "formulas": "KEY_FORMULAS",
    }
    if topic in sections:
        # Split on H2 markdown headers as section breaks
        from .knowledge import cln_3d
        if topic == "overview":
            return cln_3d.CLN_3D_OVERVIEW
        if topic == "notebooks":
            return cln_3d.CLN_3D_NOTEBOOK_INDEX
        if topic == "formulas":
            return cln_3d.CLN_3D_KEY_FORMULAS
    return full_doc


@mcp.tool()
def bem_cln(topic: str = "all") -> str:
    """
    Get BEM-CLN (per-element multipole CLN with Schur-F termination)
    documentation: multi-conductor extension of single-conductor
    Schur-F CLN, using polarizability alpha(s) and integral-equation
    Green's function coupling.

    Backs Sugahara, Nagamine, Hane (2026) IEEE Trans Mag submission,
    sections V.G (DOF accounting) and V.H (verification).

    Key features:
      - polarizability alpha(s) = V - Y_cln(s) / sigma (DC = 0, PEC = V built in)
      - 2D coupling: 1/D^2, 3D coupling: mu_0 / (4 pi D^3)
      - bounded alpha -> no phenomenological saturation factor needed
      - per-element DOF = N_Cauer + 1; total = N (N_Cauer + 1)

    Args:
        topic: Documentation section. Options:
            "all"           - Complete documentation
            "overview"      - Framework summary, DOF accounting
            "2d_rigorous"   - Phase 2.5 canonical 2D rigorous
            "3d"            - Phase 3 B rigorous 3D cuboid extension
            "scripts"       - Index of Mathematica verification scripts
    """
    if topic == "all":
        return get_bem_cln_documentation()
    from .knowledge import bem_cln as bcln
    if topic == "overview":
        return bcln.BEM_CLN_OVERVIEW
    if topic == "2d_rigorous":
        return bcln.BEM_CLN_2D_RIGOROUS
    if topic == "3d":
        return bcln.BEM_CLN_3D
    if topic == "scripts":
        return bcln.BEM_CLN_NOTEBOOK_INDEX
    return get_bem_cln_documentation()


@mcp.tool()
def cln_sibc_orthogonal(section: str = "all") -> str:
    """
    Get CLN expansion-point + SIBC orthogonal-residual theory documentation.

    New theory established 2026-05-16 in IGTE 2026 Sugahara work:
    a small number N of Cauer ladder stages plus the SIBC analytical
    asymptote span the full eddy-current frequency response via
    Foster eigenmode L^2-orthogonality:

        L^2(Ω) = span{φ_0, ..., φ_{N-1}} ⊕ span{φ_N, φ_{N+1}, ...}
        Y_exact(jω) = Y_CLN^N(jω) + Y_SIBC^⊥(jω)
        Y_SIBC^⊥(jω) -> K_SIBC (jω)^{-1/2}  as  ω -> ∞

    K_SIBC = √(σ/μ) × geometric factor (slab=1, 2D prism=perimeter,
    3D body=surface integral with edge/corner corrections).  This
    avoids the precision frontier of high-stage Cauer extraction
    (Nagamine et al. 2026 verified interval arithmetic shows
    ~60×/stage interval growth even with 192-bit MPFR + affine).

    Verified analytically with Mathematica on:
      - 1D slab (c = 1 mm Cu)
      - 2D rectangular prism (17.72 × 2 mm, Dirichlet A_z = 0)

    Open: 2D square (Nagamine geometry, degeneracy bundling) and
    3D cuboid (HDivDivFreeHex + Green's theorem on surface integrals).

    Companion to Nagamine, Yamaguchi, Sugahara, Hiruma, Mifune,
    Matsuo (JJIAM 2026 submitted) verified 2D square prism Cauer
    extraction — Nagamine = high-N verified rungs, this theory =
    small-N + SIBC asymptote.

    Args:
        section: Documentation section to return:
            "list"           - list available sections
            "overview"       - one-page summary of the construction
            "matsuo"         - relation to Matsuo SA-26-014 expansion point
            "kuriyama"       - Kuriyama 2019 multi-expansion canonical method
            "math"           - derivation, orthogonality, asymptote
            "verification"   - Mathematica results (1D slab + 2D rectangle)
            "nagamine"       - link to Nagamine 2026 verified extraction
            "outlook"        - 2D Nagamine square + 3D cuboid completion path
            "xfem_vs_sibc"   - XFEM (Hiruma 2023) / classical SIBC /
                               augmented CLN decision framework, with
                               port-driven scope and stacking strategy
                               for volume-source problems
            "all"            - full documentation (default)

    Returns:
        Markdown text of the requested section.
    """
    if section == "all":
        return get_cln_sibc_orthogonal_documentation()
    return get_cln_sibc_orthogonal_section(section)


@mcp.tool()
def cln_3d_notebook(name: str = "list") -> str:
    """
    Retrieve Tanimoto's raw 3D CLN notebook Python code.

    Provides direct access to the canonical Python code from Tanimoto's
    master's thesis + production notebooks at W:/00_CAE/NGSolve/谷本/.
    Use this when you need to see the actual implementation details
    (HCurl space construction, Kameari iteration loop, ICCG solver
    invocation, output extraction, etc.).

    Args:
        name: Notebook identifier. Options:
            "list"       - List available notebooks with file sizes
            "AT"         - A-T formulation (primary 修論 reference,
                           cylinder, 10-stage Kameari, SparseSolvPy ICCG)
            "T_Omega"    - T-Ω formulation (HCurl × H1 with Ω confined
                           to conductor)
            "APhi"       - A-Φ formulation (HCurl + H1, body current
                           via σ∇Φ)
            "2D"         - 2D scalar reference (pedagogical, Kameari
                           formula validation)
            "production" - 2024-09-17 production: A + ICCG with inline
                           gauge correction, accICCG params, type1 HCurl

    Returns:
        Full Python script content (~3-9 KB each), or list of available.
    """
    return get_cln_3d_notebook(name)


@mcp.tool()
def cln_sphere_dd_pipeline() -> str:
    """
    Get the Sphere DD (double-double, ~32 digit) VIM Cauer Ladder Network
    extraction pipeline reference.

    Verified-arithmetic CLN extractor demonstrated on the canonical Cu
    sphere benchmark (R=10mm, sigma=5.8e7 S/m, uniform B_z=1T). Pure
    Python/CuPy implementation of the entire VIM-CLN chain in DD
    precision: kernel evaluation (mpmath elliptic K(m), E(m) at 40 digit
    via dd_axisym_kernel.py), DD K/M/b assembly with multiprocessing
    (dd_sphere_axisym_mp.py), mpmath Cholesky-based generalized eigh
    at dps=35-50, and verified-interval Hankel-Pade Cauer extraction
    (mpmath.iv at 80 digit).

    Reaches DD precision floor K_lo/K_hi = 1.1e-16 across all matrix
    entries (the entire 14 trailing decimal digits below FP64 are
    correctly captured), with verified-interval relative width < 1e-30
    on all extracted Cauer rungs.

    Returns:
        Markdown documentation covering algorithm steps, sigma-rescaling
        for canonical normalization, multiprocessing speedup, file
        inventory, production scaling estimates, key implementation
        lessons, and cross-references to Nagamine, Stoll, Hiruma,
        Sugahara TEAM 28.
    """
    return get_cln_sphere_dd_documentation()


@mcp.tool()
def mmm_core(topic: str = "chubar_1998") -> str:
    """
    MMM/MSC (Magnetic Moment Method / Magnetic Surface Charge) -- how
    to BUILD models, the interaction-matrix structure, the near-null
    eigenvalue behavior, plus Radia heritage.

    Practical "how to make an MMM/MSC model" recipes: topic
    "build_msc_mmm". Matrix/system probes (MMM dense N vs current
    multipole-moment MMM MSC system): "matrix_structure". Near-null "loop"
    modes / conditioning / beautiful->ugly BiCGSTAB (CEFC 2026 study):
    "eigenvalue_nullspace". Heritage (Chubar 1998, Wakao 2007 ACA,
    Janet MMM+RNM, Le-Duc PEEC+MMM, Weddemann FEM-BEM): remaining
    topics.

    Args:
        topic: One of:
            "build_msc_mmm"      - HOW TO BUILD an MMM/MSC model
                                   (elements, materials, solve, results)
            "matrix_structure"   - Matrix probes: GetInteractMatrix is
                                   MMM/dense legacy; multipole-moment MMM MSC uses
                                   BuildMomentSystem / MomentSystemDenseRaw /
                                   dense-only diagnostics
            "eigenvalue_nullspace" - Near-null loop modes, cond ~ mu_r,
                                   beautiful->ugly BiCGSTAB (CEFC 2026)
            "multipole_modes"    - What field the 6-DoF MSC creates:
                                   mono+dipole+2 quadrupole (SVD/eig + .wls
                                   proof); cond = multipole field-strength
                                   ratio; aspect ratio (not size) sets the
                                   iteration count; distortion only rotates
                                   the quadrupole (Sugahara-lab 2026-06-22)
            "chubar_1998"        - Original Radia paper (ESRF 1998)
            "takahashi_2007_aca" - Wakao group MMM + ACA H-matrix
                                   (Sugahara lab heritage)
            "pradhan_2007"       - Kolkata cyclotron Radia/TOSCA validation
            "janet_2005_mmm_rnm" - MMM + Reluctance Network mixed (Schneider)
            "le_duc_peec_mmm"    - PEEC + MMM coupling (G2Elab)
            "weddemann_fembem"   - FEM-BEM hybrid alternative
            "mmm_vs_other"       - Decision matrix: when to use MMM vs FEM
                                   vs PEEC vs BEM
            "all"                - Everything
    """
    return get_mmm_core_documentation(topic)


@mcp.tool()
def hdiv_vim(topic: str = "overview") -> str:
    """
    HDiv-type VIM (Volume Integral Method) demag operator -- the lab's FEEC H(div) RT
    alternative/complement to the canonical multipole-moment MMM MSC kernel.  Canonical reference:
    docs/hdiv_vim/README.md.

    Key idea: SYMMETRIC demag operator N = B^T G B with the loop modes FIELD-NULL BY
    CONSTRUCTION (loops = ker B) -> mu_r-INDEPENDENT convergence + NO hand-crafted loop-star.
    VALIDATED: linear demag (sphere/spheroid/triaxial EXACT vs analytic), NONLINEAR (damped
    Newton; cube & C-yoke <1-3% vs shipped Radia in ~6 iters), distorted-mesh mu_r-independence,
    CURVED + high-order (demag exact; field accuracy-per-DOF ~10-30x vs flat Radia), and SYMMETRY
    models 1/2,1/4,1/8 (loops automatic + image-method demag).  Production entry
    radia.vim.hdiv_demag_solve; the C++ _ChargeGramHMatrix kernel is the SOLE demag operator (the
    dense Python Gram path + the analytic_gram/wilton_surface kwargs were REMOVED 2026-06-23).
    TaskManager is assumed: shared HACApK build paths and long C++ solve loops stand up or reuse
    NGSolve RegionTaskManager; direct diagnostic `.matvec()` calls plus Python/NGSolve assembly
    follow caller-wraps `with ng.TaskManager():`.  The C++ HDiv CG/MINRES/Picard kernels use
    ParallelFor/ParallelForRange for charge gather, dot products, preconditioner/vector updates,
    and AtomicAdd for sparse face-vector scatters.

    Args:
        topic: One of:
            "overview"       - what it is + why (symmetric, loops field-null, mu_r-independent) [DEFAULT]
            "implementation" - C++/pybind/Python files + APIs (rad_hdiv_vim, _ChargeGramHMatrix, ...)
            "scaling"        - the C++ charge-Gram H-matrix (exact analytic near AND far)
            "verification"   - golden tests (tests/feec/) + the verify-first bug catches
            "nonlinear"      - damped Newton on the exact C++ charge Gram; C-yoke vs Radia;
                               fail-loud on non-convergence; the honest reference distinctions
            "curved"         - curved + high-order demag (ngsolve.bem single-layer): sphere/spheroid/
                               triaxial exact; accuracy-per-DOF ~10-30x vs flat Radia; curved x nonlinear
            "symmetry"       - 1/2, 1/4, 1/8 models: loops automatic (ker B) + image-method demag value
            "cross_method"   - the demag tensor cross-validated by 3 independent discretizations
                               (FEEC surface-charge VIM == volume-FE A-formulation == BEM == Osborn)
            "reference_audit"- HOW TO DEBUG a disagreement with a FEM cross-validation reference:
                               the audit ladder (solution scalar -> evaluator invariance ->
                               drive-equivalence probe -> split tests -> closed-form arbiter) +
                               the reference-side trap catalog (coil-polygon current deficit,
                               frozen-edge truncation, conjugate-potential sign/branch cuts, ...)
            "status"         - done / open summary
            "all"            - everything
    """
    return get_hdiv_vim_documentation(topic)


@mcp.tool()
def kelvin_identify_post_hoc(topic: str = "all", vol_path: str = "") -> str:
    """
    Add Kelvin Periodic Identifications to an existing NGSolve mesh
    AFTER load (no Cubit launcher / OCC Identify needed).

    Use this when you have a `.vol` file with `kelvin_int` /
    `kelvin_ext` boundary labels but no `Identifications` section yet
    -- typical when the Cubit C++ exporter's all-or-nothing vertex
    matching skipped due to a single tolerance miss (1/8-octant
    geometry), or when the mesh came from outside the Cubit panel
    pipeline.

    Public Python API (call this from user code, NOT from the AI):
        from radia import add_kelvin_identification
        info = add_kelvin_identification(mesh)
        # -> dict(n_pairs, n_unmatched, kelvin_offset, max_dist, ...)

    Args:
        topic: Documentation topic. Options:
            "all"           - Full documentation (overview + API + workflow + caveats)
            "overview"      - Why this exists, when to use, the 3-line solution
            "api"           - add_kelvin_identification() signature + return dict
            "workflow"      - 3-step usage: detect_offset -> add -> verify
            "snippet"       - Copy-paste minimal example code
            "caveats"       - 1:1 mesh requirement, tolerance, idnr handling
            "verify"        - If vol_path given: load .vol and run helper
                              (returns the helper's diagnostic dict as text)
        vol_path: optional path to a .vol file to actually run the helper
                  against (only used when topic="verify").

    Companion tool: `kelvin_transformation` provides the underlying
    theory (Kelvin inversion, material scaling, formulations).  This
    tool is specifically about the post-hoc identification helper.
    """
    if topic == "verify":
        if not vol_path:
            return ("topic='verify' requires vol_path argument.  "
                    "Pass the path to a .vol file with kelvin_int / "
                    "kelvin_ext labels.")
        import json
        import os
        if not os.path.exists(vol_path):
            return f"vol_path not found: {vol_path}"
        try:
            from ngsolve import Mesh
            from radia import (add_kelvin_identification,
                                detect_kelvin_offset,
                                has_kelvin_identification)
        except ImportError as e:
            return (f"Cannot import ngsolve / radia in this process: {e}.  "
                    f"This tool needs radia + ngsolve installed in the MCP "
                    f"server's Python environment.")
        mesh = Mesh(vol_path)
        report = {
            "vol_path": vol_path,
            "materials": list(mesh.GetMaterials()),
            "boundaries": list(mesh.GetBoundaries()),
            "has_existing_identifications": has_kelvin_identification(mesh),
            "detected_offset": detect_kelvin_offset(mesh),
        }
        try:
            info = add_kelvin_identification(mesh)
            report["add_kelvin_identification"] = info
            report["status"] = "ok"
        except ValueError as e:
            report["status"] = "error"
            report["error"] = str(e)
        return json.dumps(report, indent=2, default=str)
    return _get_kelvin_identify_post_hoc_doc(topic)


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
            "verified_recipe" - Verified A-formulation + Periodic Kelvin
                                + BDDC cookbook (5 mandatory elements,
                                failure-mode -> root-cause table).
                                **Read this FIRST when debugging an
                                HCurl eddy-current script that is off
                                by O(10x) from analytical reference.**
            "adaptive"       - Adaptive mesh refinement with Kelvin
            "mesh_control"   - WHERE to spend elements (consolidated): the
                               Gamma-conforming constraint + six MEASURED
                               pillars (exterior volume is free, floor=Curve
                               order, p>=n & p-vs-h, optimal R/a~3, corner hp,
                               DoF-cost 1/45; holds in 2D on the -n/R ladder)
                               + the CLOSED-FORM mesh-adequacy criterion
                               (source order p*=ceil(ln eps/ln(d_max/R)),
                               geometry floor (h/R)^2k, eccentric/multi-body,
                               the A-form centre, the apparatus design calc)
                               + the non-separable build + DtN->CLN arc
                               (square C4v / cube O_h, 2D conformal disk)
            "identify"       - Periodic boundary Identify() best practices
            "tips"           - Common mistakes and performance tips
            "robustness"     - Robustness checklist (mesh copy, material scaling,
                               FreeDofs verification, symmetry models, GND)
            "verification"   - Numerical verification (single-domain approach)
            "periodic_wedge" - 1/n sector (symmetry model) with Periodic BC
            "material_exterior" - Provenance of material/conducting exterior
                               Kelvin: Freeman-Lowther (air exterior) vs
                               transformation-optics sigma (Ward-Pendry) vs
                               Sugahara 2022's validated sigma-conformal ECT
                               fusion (conductor crosses the truncation); the
                               formulation basis for radia.open_boundary.
                               kelvin_dtn's (a/r)^4 sigma / (a/r)^2 mu weights.
                               Read before claiming "material Kelvin" novel.
    """
    return get_kelvin_documentation(topic)


@mcp.tool()
def axifem_documentation(topic: str = "all") -> str:
    """
    Get radia-core axifem documentation: Henrotte axisymmetric FE
    add-on for NGSolve (registered FESpace name: "axihenrotte").

    Use this when designing or reviewing axisymmetric eddy-current /
    magnetostatic FEM problems with axis-touching elements, or when
    comparing axihenrotte to standard NGSolve H1 elements.

    Canonical API: ``FESpace("axihenrotte", mesh, order=k)`` or
    ``H1Henrotte(mesh, order=k)``.  order=1 dispatches to P1 triangles
    or Q1 quads; order=2 dispatches to P2 triangles or Q2 quads.  P2
    triangles are curved-mesh aware after ``mesh.Curve(2)``.  Curved Q2
    quads are production opt-in via ``H1Henrotte(mesh, order=2,
    curvedquad=True)``.

    Args:
        topic: Documentation section. Options:
            "all"             - all sections concatenated
            "overview"        - what it is and why NGSolve doesn't already have it
            "support_matrix"  - exact P1/P2/P2 curved/Q1/Q2/Q2 curved status
            "api"             - FESpace("axihenrotte", mesh, order=k) usage
            "hodge_geometry"  - differential-geometry/Hodge view of Henrotte
                                axisymmetric reduction
            "taskmanager"     - NGSolve-native parallel execution contract
            "basis_p1"        - order=1 P1 triangle + Q1 quad basis details
            "basis_p2"        - order=2 P2 triangle + Q2 quad basis details
            "curved_geometry" - P2 curved triangle and opt-in Q2 curved quad support
            "vs_standard_h1"  - 6-property comparison table vs H1 order=2
            "validation"      - cross-validation references; Hessian-of-W convention
            "kelvin"          - Phase B3 z-offset Kelvin recipe (Periodic + H1Henrotte,
                                sphere -0.001 % vs Stoll, mu-vs-nu factor & Curve(2) gotchas)
            "magnet"          - permanent-magnet source term (FEMM prob3big.cpp port):
                                weak-form RHS + magnetized-sphere validation (-0.05%)
            "file_layout"     - where each piece lives (C++, Mathematica, tests)
            "why_dropped_p3"  - why p=3 was attempted and reverted (Vandermonde cond ~ 1e30)
    """
    return get_axifem_documentation(topic)


@mcp.tool()
def femm_parity_documentation(topic: str = "all") -> str:
    """
    Get FEMM-parity documentation: which FEMM (Finite Element Method Magnetics,
    D. Meeker) analyses are reproduced as EXECUTABLE + TESTED NGSolve capability
    in radia-ngsolve, with the function to call and the analytical benchmark each
    was validated against (all <2%, most <0.2%).

    Use this when designing an FEM analysis that a FEMM user would run in FEMM 4.2,
    to find the equivalent radia-ngsolve function and its API.

    Args:
        topic: Documentation section. Options:
            "all"        - all sections concatenated
            "overview"   - design rule: build capability not a number
            "matrix"     - planar / axisymmetric capability matrix (21 analyses)
            "magnetics"  - magnetostatic / nonlinear / eddy / circuit + axi API
            "scalar"     - electrostatic / heat / current-flow (csolv/hsolv) API
            "lamination" - laminated steel (anisotropic + complex-mu), multi-
                           conductor proximity circuit, FEMM open-bdry cross-check
            "validation" - regression test list and per-test error bounds
    """
    return get_femm_parity_documentation(topic)


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
            "elements"       - ObjRecMag, ObjHexahedron, ObjTetrahedron, ObjWedge, ObjPyramid
            "materials"      - MatLin, MatSatIsoTab, hysteresis, permanent magnets
            "solver"         - rad.Solve, SolverConfig, LU/BiCGSTAB/HACApK
            "field"          - rad.Fld, batch evaluation, A field
            "removed_apis"   - Removed APIs reference (FldEnr/ObjDivMag/UtiDmp/...)
            "ngsolve"        - RadiaField CF, netgen_mesh_to_radia
            "background"     - ObjBckg, Biot-Savart source
            "ima"            - Image Method of Analysis
            "best_practices" - lab best practices, incl. #11 storing
                               numerical-experiment results as JSON
                               (reproducible) + delete-by-folder cleanup
    """
    return get_radia_documentation(topic)


@mcp.tool()
def md2html_usage() -> str:
    """Get md2html converter documentation (MathJax, reference links, styled HTML)."""
    return get_md2html_documentation()


@mcp.tool()
def analytical_formulas(topic: str = "all") -> str:
    """
    Get documentation for radia.analytical_formulas (closed-form reference layer).

    The package collects nine modules of closed-form expressions taken from
    the Wakao-Igarashi-Fujiwara-Kameari review series (IEE Japan, 2002-2004).
    These are the trusted-baseline analytical results that the rest of Radia
    (MMM, MSC, PEEC, FEM panels, ngsolve.bem) is sanity-checked against.

    For any new analysis the FIRST QUESTION to ask is "is there a closed form
    here that I can validate against?". Use the validation_use_cases topic
    for the practical mapping from analysis type -> applicable formula.

    Args:
        topic: Documentation topic. Options:
            "all"                  - Complete documentation
            "overview"             - Module index, conventions, quick API surface
            "ellipsoid"            - Rotational ellipsoid demag factor + torque
                                     (Part 5 §5, eq 38-44; Osborn references included)
            "ac_locus"             - B_max / B_min of an AC vector phasor's locus
                                     (Part 5 §4, eq 29-37)
            "shielding"            - Magnetic-shell shielding factor S
                                     (Part 1 §5, eq 23-24)
            "rect_magnet_2d"       - 2D rectangular bar A_z, B_x, B_y
                                     (Part 2 §2, eq 2-3)
            "plate_eddy"           - Thin rectangular-plate eddy current
                                     (Part 1 §6.1, eq 26-27)
            "solenoid_central"     - Fabri form factor + closed-form axial field
                                     (Part 4 §4, eq 26-27)
            "three_phase_line"     - Triangle / planar / helical line field
                                     (Part 4 §5, Part 5 §3)
            "elliptic_integrals"   - K(k), E(k) Hastings polynomial approximation
                                     (Part 3 §3, Tables 1-2)
            "gauss_legendre"       - Gauss-Legendre nodes / weights to n=24
                                     (Part 3 §4, Table 3)
            "validation_use_cases" - "I have X, sanity-check it with Y" mapping

    Sources:
        - src/radia/analytical_formulas/   (Python modules, ~40 functions across Part 1-9)
        - tests/analytical_formulas/       (170+ pytest tests, < 1 s)
        - docs/analytical_formulas/        (analytical_formulas.ipynb: 12 demos + PNGs)
        - docs/analytical_formulas.md      (PDF -> code cross-reference)
        - PDFs themselves: lab-internal, not redistributed with the repo.
    """
    return get_analytical_formulas_documentation(topic)


@mcp.tool()
def force_validation(topic: str = "all") -> str:
    """
    EM force extraction in NGSolve + independent <-> NGSolve cross-validation.

    Records how the radia-ngsolve FEM path chooses and computes electromagnetic
    force/torque: weighted Maxwell stress, Maxwell surface stress, Lorentz
    conductor force, air-gap pressure, energy/coenergy checks, dq torque, and
    the cross-validation results that make the NGSolve magnetostatic path
    trustworthy. Reference values are kept as stored regression references.

    Validated (linear magnetostatics, A-form, HCurl order 2):
      * uniformly magnetized sphere: reference == NGSolve to 0.11 %, both <0.5 %
        of the analytic (2/3)Br interior field.
      * coil + linear-iron sphere force: reference == NGSolve to ~3 %, both near
        the dipole-in-gradient analytic.

    Read this when: implementing/validating an EM force computation, deciding
    whether to trust an NGSolve magnetostatic result, or reviewing the
    reference-solve note (the reference_note topic summarises the shared
    geometry / material / source / force-method contract of the cross-check).

    Args:
        topic: all (default) | method_map | eggshell | cross_validation | reference_note
    """
    return get_force_validation_documentation(topic)


@mcp.tool()
def ngsbem_inductance(topic: str = "all") -> str:
    """
    Get ngsolve.bem boundary element method documentation for inductance extraction.

    ngsolve.bem is NGSolve's native boundary element module. Combined with Cubit
    mesh export and SetGeomInfo, it enables accurate inductance extraction
    on high-order curved surface elements.

    Key workflow:
      Cubit mesh -> export netgen "mesh.vol" order N -> LaplaceSL BEM -> L extraction

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


@mcp.tool()
def peec_inductance(topic: str = "all") -> str:
    """
    Get documentation for the Radia PEEC-inductance (coil only, STEP) panel mode.

    Lightest mode in the Radia panel family: STEP / Cubit .jou → perimeter
    filaments → Biot-Savart + Loop-bundle PEEC solve → L_coil, R_coil.
    No workpiece, no BEM, no FEM mesh.  Use when the question is "what
    is this coil's L / R at this frequency?" and you do not need
    heating or workpiece physics.

    Centerline extraction is **STEP-only** (v4.48.1, 2026-05-15):
    no caller-supplied JSON, no ``--path-points-m`` flag.
    Classification-based single dispatch picks ONE of 5 predicates;
    if none cover the conductor bbox, raises with a HINT pointing at
    CAD regeneration or BEM-A switch (NEVER a silent fallback).

    Five centerline predicates (positive-match, no try/except cascade):
      1. Loft of profiles (>= 5 planar end-caps)  -> NN-chain centroids
      2. United multi-turn pancake (CIRCLE-edge stations) -> NN-chain
         circle centres
      3. Single-loop revolution sweep (TORUS / CYLINDER / CONE /
         REVOLUTION + PLANE caps) -> analytical arc
      4. OPEN coil with caps -> longest open lateral rim edge
         (handles "arc + leads" e.g. keiko outsideline.step)
      5. CLOSED full revolution (no caps) -> coil_topology spine

    Sibling-.jou auto-preference: if the user picks ``foo.step`` and
    ``foo.jou`` (case-insensitive stem match) coexists in the same
    directory AND contains the PEEC explicit-centerline pattern,
    the calc switches to the .jou parser automatically.

    Unicode / Japanese path: verified end-to-end (Python + Cubit
    plugin via utf8_path.hpp + PySide6 QProcess).

    Args:
        topic: Options:
            "all"             - Complete documentation
            "overview"        - What it solves, when to use, perimeter placement
            "centerline"      - STEP centerline extraction (5 classification predicates)
            "filament_dispatch" - n_peri filament placement paths (single-dispatch)
            "step_authoring"  - Cubit / build123d recipes for auto-detect-friendly STEPs
            "jou"             - .jou explicit centerline parser
            "sibling_jou"     - Auto-prefer sibling .jou when co-located
            "japanese_path"   - Unicode / Japanese path support
    """
    return get_peec_inductance_documentation(topic)


@mcp.tool()
def fem_bem_schur(topic: str = "all") -> str:
    """
    Get FEM-BEM Schur coupling documentation -- exact open boundary for interior FEM.

    The exterior Laplace DtN assembled from ngsolve.bem replaces the unbounded
    exterior domain with a boundary operator on the coupling surface, giving an
    EXACT transparent BC without mesh truncation, PML, or Kelvin transformation.

    Coupled system:  (K_FEM + P^T Λ_ext P) u = f
    where Λ_ext = V⁻¹ (−½M + K)  [exterior BIE, NGSolve sign convention]
    and P is the H1 volume → SurfaceL2 boundary projection.

    Verification: spherical shell with inner Dirichlet u=cosθ, BEM DtN at outer
    sphere; exact solution u = R_inner² cosθ / r²; L2 rel_err < 2% at order=2.

    Companion code: radia_ngsolve.fem_bem_coupling
      - laplace_fem_bem_schur():  coupled solve (dense, for verification)
      - sphere_shell_fem_bem():   verification on spherical shell

    Args:
        topic: Documentation topic. Options:
            "all"          - Complete documentation
            "overview"     - Concept, DtN operator, coupling matrix
            "api"          - laplace_fem_bem_schur and sphere_shell_fem_bem usage
            "applications" - Open BC magnetostatics, motor models, AGE + BEM
    """
    return get_fem_bem_schur_documentation(topic)


@mcp.tool()
def airgap_motor_workflow(topic: str = "all") -> str:
    """
    Get AGE rotating machine workflow documentation -- nonlinear iron + AGE coupling.

    Wraps the linear AGE core (airgap_machine.py) in a Picard fixed-point iteration
    for machines with saturating iron B-H characteristics.  The air-gap element stays
    analytic (no gap mesh); only the iron FE stiffness is updated each iteration.

    Key validation:
      - nu_cf_fn=None (or constant nu) → 1 iteration, identical to direct airgap_solve
      - Froelich saturation model → converges in 5-15 iterations
      - age_motor_rotation_sweep → torque-angle curve, all angles converged

    Companion code: radia_ngsolve.airgap_motor_workflow
      - age_motor_nonlinear_solve(): Picard loop for one rotor position
      - age_motor_rotation_sweep():  torque–angle curve

    Args:
        topic: Documentation topic. Options:
            "all"               - Complete documentation
            "overview"          - AGE core, Picard iteration, linear limit
            "api"               - age_motor_nonlinear_solve and rotation_sweep usage
            "saturation_models" - Froelich, piecewise linear, NGSolve BSpline B-H
            "validated"         - Test patterns: linear limit, saturation, sweep
    """
    return get_airgap_motor_workflow_documentation(topic)


@mcp.tool()
def dtn_coarse_mesh(topic: str = "all") -> str:
    """
    Why open-boundary methods stay accurate on COARSE meshes -- a spectral
    (DtN-matrix) explanation of the Kelvin-transformation coarse-mesh accuracy.

    Kameari demonstrated the Kelvin transform's coarse-mesh accuracy empirically,
    by mesh refinement.  This reframes it as a PROPERTY OF THE DtN MATRIX across
    mesh sizes: every open-BC closure (Kelvin / BEM / PML / Robin) approximates
    the one exterior DtN operator Λ_ext, whose sphere eigenvalues are the
    mesh-independent ladder λ_n = −(n+1)/R.  The discrete matrix Λ_h reproduces
    the LOW-degree eigenvalues accurately and almost independently of h even on
    the coarsest mesh; the per-mode error grows with degree n.  Since a compact
    source's field is dominated by low multipoles, the coarse mesh already
    resolves everything that matters -- Kameari's result, read off the spectrum.

    Companion code (both sides MEASURED -- two discretisations of the one Λ_ext):
      - bem_integral.exterior_dtn_spectrum():  eigenvalues of the BEM Λ_h matched
                                  to −(n+1)/R (boundary operator spectrum)
      - bem_integral.dtn_spectrum_vs_mesh():   per-degree error vs mesh size,
                                  with coarse-low-mode and accurate-band summary
      - fem_bem_coupling.kelvin_dtn_eigenvalue(): the Kelvin closure's effective
                                  DtN (volume FEM); order≥n kills the polynomial
                                  error, then a curved-geometry floor (~5-6 digits
                                  in 3D = Kameari's result; the dominant dipole
                                  inverts to a linear field, accurate at order 1)

    Args:
        topic: Documentation topic. Options:
            "all"          - Complete documentation
            "overview"     - One operator behind every open BC; the spectral argument
            "numerics"     - Measured DtN spectrum vs mesh table + how to read it
            "api"          - exterior_dtn_spectrum / dtn_spectrum_vs_mesh usage
            "applications" - Air-box sizing, method choice, trusting coarse Kelvin, debugging
            "p_method"     - Kelvin is a p-method not an h-method (measured p-vs-h,
                             ~20-80x DOF gap); polyhedron faceting error scales with
                             multipole degree (dipole robust); (R,p) design rule
            "formulation"  - Differential-geometry view: Omega vs A as Hodge-dual,
                             complementary dual bracketing (certified bounds), DtN
                             gradient block is formulation-independent, conformal
                             pullback material (scalar 0-form / tensor 1-form),
                             infinity = one-point conformal compactification, FEEC
            "datasheet"    - Problem-INDEPENDENT performance: open-BC error factors
                             into source multipoles x method eigenvalue-defect; the
                             Kelvin closure has an analytic, universal "datasheet"
                             (certify once, predict any problem by multipole content).
                             Includes COST measured by DoF INCREMENT (not time): the
                             closure adds a Gamma-scale coarse ball (measured dDoF=58,
                             closure error ~1/45 of the interior FE error) -> Kelvin is
                             cheap; keep it SPARSE (condensing to a dense DtN BC is
                             N^4/3, 10-20x nnz).
                             Includes the C/L dual: capacitance C <- n=0 (monopole,
                             exact) and external inductance L_ext <- n=1 (dipole) are
                             two Steklov modes of the SAME exterior DtN (W_ext =
                             1/2 mu0 (n+1)/R oint phi^2); certified via the magnetic
                             POTENTIAL exterior (A no-cut / Omega+cut), NOT the ngsbem
                             vector single-layer L=mu0 J^T(LaplaceSL)J (a different
                             operator with no -(n+1)/R ladder)
            "method_map"   - UNIFIED open-boundary map: Kelvin / BEM / PML / IABC /
                             CLN on the three selection axes (frequency / geometry /
                             space-vs-time) + the no-free-lunch modal axis + the
                             selection table. Audit-verified anchors.
                             doc: docs/open_boundary/OPEN_BOUNDARY_MAP.md
    """
    return get_dtn_coarse_mesh_documentation(topic)


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
      panel_qt_testing        -- retired PySide note; current gate is
                                 validation_test/panels/test_notebook_workbench.py
      learn_edition_cap       -- ignore the 50k warning, export bypasses it

    Args:
        topic: Empty for the full document, or one of the topic
               keywords above for a single section.
    """
    return get_panel_gui_pitfalls(topic)


@mcp.tool()
def install_deploy(topic: str = "") -> str:
    """
    Radia install / deploy policy and recipes — 2-tier configuration
    (LAB + 100号機 editable / mdx + hibino PyPI), reversible migration steps, and the
    non-obvious gotchas that cause silent breakage.

    Read this when:
      * Setting up a new lab machine.
      * Migrating a machine between editable / PyPI install.
      * Diagnosing "import works but pip says wrong version" or
        "DLL load failed" on a freshly-deployed machine.

    Topics:
      two_tier                    -- current LAB-editable / PyPI-consumer policy
      lab_editable                -- LAB editable install
      hyaku_editable              -- 100号機 NAS editable install
      mdx_pypi                    -- mdx PyPI install
      hibino_pypi                 -- hibino PyPI install via `ssh hibino`
      editable_to_pypi_migration  -- e.g. 100号機 NAS-editable -> PyPI
      pypi_to_editable_migration  -- e.g. mdx PyPI -> editable
      metadata_sync               -- pip metadata vs radia.__version__
      pyd_dll_bootstrap           -- cubit_mesh_curver requires `import radia` first
      cubit_plugin_layers         -- Cubit plugin lives in TWO places
      common_failure_modes        -- symptoms and fixes table

    Args:
        topic: Empty for the full document, or one of the topics above.
    """
    return get_install_deploy_documentation(topic)


@mcp.tool()
def release_workflow(topic: str = "") -> str:
    """
    Release-QUD workflow for the Radia monorepo
    (3 packages / 4 machines: LAB, 100号機, mdx, hibino). Documents the 9-phase
    pipeline, the 4 pre-flight gates added 2026-05-03, the historical
    CI failure modes + their root causes, and the patch-bump recovery
    protocol when a tag CI fails.

    Read this when:
      * The user asks for a release / version bump / PyPI publish.
      * CI on a tag ref fails after `git push --tags`.
      * "Release" workflow shows skipped in `gh run list` (CI never
        went green).
      * A user reports "is X.Y.Z on PyPI yet?" and propagation is
        stuck.

    Topics:
      overview               -- what gets released and why atomically
      phases                 -- the 9-phase pipeline (table)
      preflight_gates        -- Phase 2.5 4-gate pre-push validation
      ci_failure_modes       -- known CI failures + cause + fix table
      recovery               -- when CI on a tag fails AFTER push
      patch_bump_protocol    -- exact steps for retry after CI failure
      lab_lock_release       -- pre-deploy: stop processes that hold .pyd
      monorepo_lockstep      -- 4-6 version files that must stay in sync
      ci_monitor_skill       -- companion skill for Phase 7

    Args:
        topic: Empty for the full document, or one of the topics above.
    """
    return get_release_workflow_documentation(topic)


@mcp.tool()
def standalone_panels(topic: str = "") -> str:
    """
    Retired standalone PySide panel topic.  The canonical Radia panel surface
    is now the Jupyter notebook workbench (`radia_<app>.ipynb` +
    `radia.<app>_notebook`).  This tool remains as a compatibility redirect.

    Read this when:
      * A user asks about the old standalone panel entry-points.
      * An MCP client still calls the historical `standalone_panels` topic.
      * You need the post-migration notebook route and no-PySide boundary.

    Topics:
      quick_start      -- current notebook route
      four_panels      -- active notebook workbenches
      build_notebook_gui -- construction recipe for new notebook GUIs
      cubit_panels_migration -- examples/cubit_panels promotion route
      vol_sources      -- Cubit / Netgen-OCC / build123d / etc.
      vs_cubit         -- notebook route vs Cubit export/plugin boundary
      ih_methods       -- IH through `radia_ih.ipynb`
      troubleshooting  -- common post-migration issues

    Args:
        topic: Empty for the full document, or one of the topics above.
    """
    return get_standalone_panels_documentation(topic)


@mcp.tool()
def loop_learning(topic: str = "overview") -> str:
    """
    Public-safe CAE loop learning rules distilled from repeated validation
    rotations. Use this after a multi-tool loop has produced artifacts and the
    user asks whether the MCP server has actually learned from them.

    Topics:
      overview                -- artifact -> MCP learning workflow
      dual_lane               -- split one artifact into public/source-tool lanes
      mesh_geometry_vol       -- .vol, tri/tet, block registration, geometry checks
      force_moment            -- Lorentz, Maxwell traction, coenergy, moments
      motor_airgap_torque     -- Br/Bt harmonic torque phase and sign gates
      rf_acoustic_passivity   -- acoustic impedance and two-port passivity gates
      artifact_feedback       -- JSON/notebook/result artifact -> MCP knowledge
      mcp_closure             -- collected/distilled/encoded/verified/learned labels
      all                     -- complete document

    Args:
        topic: Topic name, or "all".
    """
    return get_loop_learning_documentation(topic)


@mcp.tool()
def basis_functions(topic: str = "") -> str:
    """
    Finite-element basis function library — Mathematica-canonical
    reference for H1, HCurl, HDiv, L2 on triangle / tet / quad / hex
    / prism / pyramid, arbitrary order.

    Phase 1 (radia-mcp 0.38.0) covers triangle + tetrahedron with
    H1 Lagrange P1-P5 (triangle) / P1-P3 (tet), HDiv RT₀ (= RWG used
    by BEM-A), and L2 P0-P3.  Mathematica notation is the canonical
    form; SymPy / NumPy translations are CI-verified in
    ``tests/basis/test_basis_functions.py``.

    Read this when:
      * Implementing a new BEM / FEM panel and need the canonical
        basis formula + a verified NumPy translation.
      * Cross-checking NGSolve's hierarchical output against a
        textbook Lagrange definition.
      * Debugging RWG sign / divergence issues in BEM-A.
      * Generating a code-gen template (NumPy / MATLAB / SymPy /
        Maple) from the Mathematica canonical form.

    Topics:
      theory_overview               -- what each space (H1/HCurl/HDiv/L2) means
      code_gen_pattern              -- Mathematica → NumPy / SymPy translation rules
      triangle_h1_lagrange          -- Triangle P1, P2, P3, P4, P5 (Lagrange nodal)
      triangle_h1_hierarchical      -- Triangle hierarchical (NGSolve compat)
      triangle_rwg                  -- Triangle HDiv RT₀ (= RWG, BEM-A workhorse)
      triangle_l2                   -- Triangle L2 P0-P3, Dubiner orthogonal
      tet_h1_lagrange               -- Tetrahedron P1, P2, P3
      tet_rwg                       -- Tetrahedron HDiv RT₀ (3D RWG, future)
      verification_recipes          -- partition of unity, div=σ/A, etc.

    Companion files:
      packages/radia-mcp/src/radia_mcp/mathematica/basis_functions/RadiaBasis.m
        -- canonical Mathematica package
      tests/basis/test_basis_functions.py
        -- CI-verified NumPy ports

    Args:
        topic: Empty for the full document, or one of the topics above.
    """
    return get_basis_functions_documentation(topic)


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




register_status_tool(
    mcp,
    server_name='mcp-server-radia-ngsolve',
    description='Radia + NGSolve: Kelvin / sparsesolv / CLN / PEEC / analytical formulas / lint',
    subpackage='radia_mcp.radia_ngsolve',
    related_servers=["fem", "bem", "matrix-solvers"],
    optional_deps=["radia", "ngsolve"],
)


def main():
    """Entry point for mcp-server-ngsolve console script."""
    if len(sys.argv) > 1 and sys.argv[1] == '--selftest':
        from radia_mcp.common.utf8_stdout import use_utf8_stdout
        use_utf8_stdout()
        _selftest()
    else:
        mcp.run(transport="stdio")


if __name__ == '__main__':
    main()
