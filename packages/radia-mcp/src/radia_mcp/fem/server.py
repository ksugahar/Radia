"""MCP Server: radia_mcp.fem

FEM formulations theory + decision layer for EM analysis.

Covers all major potential formulations (A-Omega, T-Omega, H, Reduced,
Darwin), element technology (edge, high-order, XFEM, isogeometric, DG),
gauging (tree-cotree, ungauged + AMS), open boundary (Kelvin, ABC, PML),
domain-specific (axisym Henrotte, TD-FEM, harmonic balance, HF),
circuit coupling, and large-scale + multi-scale (Hollaus MSFEM).

This is the **theory/genealogy** layer.  For code usage:
- NGSolve API: see `radia_mcp.radia_ngsolve`
- Axisymmetric: see `radia_mcp.radia_ngsolve.axifem_documentation`
- Hollaus lamination: see `radia_mcp.motor.hollaus_eddy`
- Solver+preconditioner: see `radia_mcp.matrix_solvers`
- BEM/MoM: see `radia_mcp.bem`

Usage:
    mcp-server-fem              # stdio
    mcp-server-fem --selftest   # self-test
"""

import asyncio
import json
import sys

from mcp.server.fastmcp import FastMCP
from mcp.types import ToolAnnotations
from ..common import register_status_tool
from ..common.mcp_contract import apply_tool_contract

from .overview_knowledge import get_overview_knowledge
from .potential_formulations_knowledge import get_potential_formulations_knowledge
from .elements_knowledge import get_elements_knowledge
from .gauge_open_boundary_knowledge import get_gauge_open_boundary_knowledge
from .time_domain_axisym_knowledge import get_time_domain_axisym_knowledge
from .large_scale_special_knowledge import get_large_scale_special_knowledge
from .ngsolve_hierarchy_knowledge import get_ngsolve_hierarchy_knowledge
from .xfem_em_hiruma_knowledge import get_em_xfem_knowledge
# Non-conforming mesh coupling (mortar / Nitsche / FETI-DP / DG /
# HFSS Mesh Fusion analogs).  2026-05-26: distilled from
# public-safe curated corpus (33 PDFs across 6 subfolders).
from .nonconforming_mesh_coupling_knowledge import (
    get_nonconforming_mesh_coupling_documentation,
)
# Equivalence-theorem near-field source (near-field-source equivalent).
# 2026-05-26: distilled from lab background on the equivalence theorem
# (axisymmetric basics, the 2015 open-domain generalization, a WPT
# reconstruction study).  IABC / SDI content intentionally excluded --
# Kelvin transformation is the lab BC.
from .equivalence_source_knowledge import (
    get_equivalence_source_knowledge,
)
from ..radia_ngsolve.profile2d_handoff import profile2d_handoff_gate
from ..radia_ngsolve.vol2d_transient_runtime import execute_transient_runtime
from ..radia_ngsolve.validation_evidence import validate_evidence_bundle
from .uninstall_safety import validate_solver_uninstall_safety_evidence
from .axifem_retirement import validate_axifem_element_evidence
from .legacy_corpus_absorption import validate_legacy_corpus_evidence
from .legacy_signature_migration import (
    validate_legacy_signature_migration_evidence,
)
from .axifem_signature_execution import validate_axifem_signature_execution


mcp = FastMCP("mcp-server-fem")


def _decode_worker_json(stdout: bytes) -> dict:
    """Decode the final JSON object while tolerating native diagnostics."""

    text = stdout.decode("utf-8", errors="replace").strip()
    if not text:
        raise RuntimeError("scalar worker returned no JSON")
    for line in reversed(text.splitlines()):
        candidate = line.strip()
        if not candidate.startswith("{"):
            continue
        try:
            result = json.loads(candidate)
        except json.JSONDecodeError:
            continue
        if isinstance(result, dict):
            return result
    raise RuntimeError("scalar worker output did not end with a JSON object")


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True,
        destructiveHint=False,
        idempotentHint=True,
        openWorldHint=False,
    )
)
def fem_axifem_element_evidence_gate(evidence_json: str) -> str:
    """Validate content-addressed evidence for P1/Q1/P2/Q2 curved axifem paths."""

    try:
        evidence = json.loads(evidence_json)
        result = validate_axifem_element_evidence(evidence)
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "schema": "radia.axifem-element-evidence-gate.v1",
            "status": "invalid_input",
            "pass": False,
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True,
        destructiveHint=False,
        idempotentHint=True,
        openWorldHint=False,
    )
)
def fem_solver_uninstall_safety_gate(evidence_json: str) -> str:
    """Validate reversible solver-uninstall evidence without local path access."""

    try:
        evidence = json.loads(evidence_json)
        result = validate_solver_uninstall_safety_evidence(evidence)
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "schema": "radia.solver-uninstall-safety-gate.v1",
            "status": "invalid_input",
            "pass": False,
            "ready_for_explicit_uninstall_approval": False,
            "solver_uninstall_performed": False,
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True,
        destructiveHint=False,
        idempotentHint=True,
        openWorldHint=False,
    )
)
def fem_legacy_corpus_absorption_gate(evidence_json: str) -> str:
    """Gate solver-neutral model, automation, document, and topic coverage."""

    try:
        evidence = json.loads(evidence_json)
        result = validate_legacy_corpus_evidence(evidence)
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "schema": "radia.legacy-solver-corpus-gate.v1",
            "status": "invalid_input",
            "pass": False,
            "live_dependency_required": False,
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True,
        destructiveHint=False,
        idempotentHint=True,
        openWorldHint=False,
    )
)
def fem_legacy_signature_migration_gate(evidence_json: str) -> str:
    """Gate solver-neutral signature routing and false-ready rejection evidence."""

    try:
        evidence = json.loads(evidence_json)
        result = validate_legacy_signature_migration_evidence(evidence)
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "schema": "radia.legacy-signature-migration-gate.v1",
            "status": "invalid_input",
            "pass": False,
            "live_dependency_required": False,
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True,
        destructiveHint=False,
        idempotentHint=True,
        openWorldHint=False,
    )
)
def fem_axifem_signature_execution_gate(evidence_json: str) -> str:
    """Gate axisymmetric Henrotte executions separately from retirement readiness."""

    try:
        evidence = json.loads(evidence_json)
        result = validate_axifem_signature_execution(evidence)
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "schema": "radia.axifem-signature-execution-gate.v1",
            "status": "invalid_input",
            "pass": False,
            "retirement_ready": False,
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True,
        destructiveHint=False,
        idempotentHint=True,
        openWorldHint=False,
    )
)
async def fem_vol2d_scalar_analysis(analysis_json: str) -> str:
    """Solve or replay a portable 2-D scalar-PDE ``.vol`` artifact.

    The closed-world physics choices are electrostatics, current flow (DC or
    harmonic lossy dielectric), and steady heat conduction.  Planar models
    carry an explicit depth; axisymmetric models use the full ``2*pi*r``
    measure.  One owned worker returns terminal/balance observables and exact
    JSON/CSV/Gmsh sidecars without touching a shared CAE session or a
    caller-selected output path.
    """

    process = None
    try:
        json.loads(analysis_json)
        process = await asyncio.create_subprocess_exec(
            sys.executable,
            "-m",
            "radia_mcp.radia_ngsolve.vol2d_scalar_worker",
            stdin=asyncio.subprocess.PIPE,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        stdout, stderr = await asyncio.wait_for(
            process.communicate(analysis_json.encode("utf-8")),
            timeout=180.0,
        )
        if process.returncode != 0:
            message = stderr.decode("utf-8", errors="replace").strip()
            raise RuntimeError(message[-1000:] or "scalar worker failed")
        result = _decode_worker_json(stdout)
    except asyncio.TimeoutError:
        if process is not None and process.returncode is None:
            process.kill()
            await process.wait()
        result = {
            "schema": "radia.vol2d-scalar-analysis.v1",
            "status": "timeout",
            "error": "scalar worker exceeded 180 seconds",
        }
    except asyncio.CancelledError:
        if process is not None and process.returncode is None:
            process.kill()
            await process.wait()
        raise
    except (json.JSONDecodeError, OSError, TypeError, ValueError, RuntimeError) as exc:
        result = {
            "schema": "radia.vol2d-scalar-analysis.v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True,
        destructiveHint=False,
        idempotentHint=True,
        openWorldHint=False,
    )
)
def fem_profile2d_handoff_gate(packet_json: str) -> str:
    """Gate a source-neutral 2-D CAD/mesh/solver handoff packet.

    The packet uses metre coordinates, explicit line and signed-arc topology,
    and separate material/boundary/conductor semantics.  STEP is accepted as
    geometry only and therefore requires a matching semantic JSON sidecar.
    Optional MATLAB MEX or Simulink S-function I/O must be fixed width, use
    explicit SI units, and contain no dynamic paths.
    """

    try:
        packet = json.loads(packet_json)
        if not isinstance(packet, dict):
            raise ValueError("packet_json must decode to an object")
        result = profile2d_handoff_gate(packet)
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "schema": "radia.profile2d-handoff-gate.v1",
            "status": "invalid_input",
            "pass": False,
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool(
    title="Validation evidence bundle gate",
    annotations=ToolAnnotations(
        readOnlyHint=True,
        destructiveHint=False,
        idempotentHint=True,
        openWorldHint=False,
    ),
)
def fem_validation_evidence_bundle(packet_json: str) -> str:
    """Validate balanced, content-addressed solver-replacement evidence.

    Artifact content is supplied inline, so this public tool never reads a
    caller-selected local path.  It requires verified public and source lanes,
    full commit identities, explicit capability claims, and exact required
    coverage.  An incomplete bundle remains useful but cannot authorize a
    dependency retirement.
    """

    try:
        packet = json.loads(packet_json)
        if not isinstance(packet, dict):
            raise ValueError("packet_json must decode to an object")
        result = validate_evidence_bundle(packet)
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "schema": "radia.validation-evidence-bundle.v1",
            "status": "invalid_input",
            "contract_valid": False,
            "retirement_ready": False,
            "solver_uninstall_performed": False,
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool(
    title="Kelvin open-boundary validation",
    annotations=ToolAnnotations(
        readOnlyHint=True,
        destructiveHint=False,
        idempotentHint=True,
        openWorldHint=False,
    ),
)
async def fem_kelvin_open_boundary_validation(request_json: str = "{}") -> str:
    """Run static 2-D/3-D Kelvin validation in one owned native worker.

    This is a magnetostatic/Poisson gate.  It deliberately rejects PML and any
    attempt to infer a wave-radiation boundary from static Kelvin evidence.
    The genuine 3-D two-sphere periodic solve returns mesh, operator, residual,
    and result identities together with analytic multipole errors.
    """

    process = None
    try:
        request = json.loads(request_json)
        if not isinstance(request, dict):
            raise ValueError("request_json must decode to an object")
        process = await asyncio.create_subprocess_exec(
            sys.executable,
            "-m",
            "radia_mcp.radia_ngsolve.kelvin_open_boundary_worker",
            stdin=asyncio.subprocess.PIPE,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        stdout, stderr = await asyncio.wait_for(
            process.communicate(json.dumps(request).encode("utf-8")), timeout=180.0
        )
        if process.returncode != 0:
            message = stderr.decode("utf-8", errors="replace").strip()
            raise RuntimeError(message[-1000:] or "Kelvin validation worker failed")
        result = _decode_worker_json(stdout)
    except asyncio.TimeoutError:
        if process is not None and process.returncode is None:
            process.kill()
            await process.wait()
        result = {
            "schema": "radia.kelvin-open-boundary-validation.v1",
            "status": "timeout",
            "pass": False,
            "error": "Kelvin validation worker exceeded 180 seconds",
        }
    except asyncio.CancelledError:
        if process is not None and process.returncode is None:
            process.kill()
            await process.wait()
        raise
    except (json.JSONDecodeError, OSError, TypeError, ValueError, RuntimeError) as exc:
        result = {
            "schema": "radia.kelvin-open-boundary-validation.v1",
            "status": "invalid_input",
            "pass": False,
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool(
    title="Harmonic magnetic validation",
    annotations=ToolAnnotations(
        readOnlyHint=True,
        destructiveHint=False,
        idempotentHint=True,
        openWorldHint=False,
    ),
)
async def fem_harmonic_magnetic_validation(request_json: str = "{}") -> str:
    """Run planar and weighted-axisymmetric complex magnetic solves."""

    process = None
    try:
        request = json.loads(request_json)
        if not isinstance(request, dict):
            raise ValueError("request_json must decode to an object")
        process = await asyncio.create_subprocess_exec(
            sys.executable,
            "-m",
            "radia_mcp.radia_ngsolve.harmonic_magnetic_worker",
            stdin=asyncio.subprocess.PIPE,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        stdout, stderr = await asyncio.wait_for(
            process.communicate(json.dumps(request).encode("utf-8")), timeout=120.0
        )
        if process.returncode != 0:
            message = stderr.decode("utf-8", errors="replace").strip()
            raise RuntimeError(message[-1000:] or "harmonic magnetic worker failed")
        result = _decode_worker_json(stdout)
    except asyncio.TimeoutError:
        if process is not None and process.returncode is None:
            process.kill()
            await process.wait()
        result = {
            "schema": "radia.harmonic-magnetic-validation.v1",
            "status": "timeout",
            "pass": False,
            "error": "harmonic magnetic worker exceeded 120 seconds",
        }
    except asyncio.CancelledError:
        if process is not None and process.returncode is None:
            process.kill()
            await process.wait()
        raise
    except (json.JSONDecodeError, OSError, TypeError, ValueError, RuntimeError) as exc:
        result = {
            "schema": "radia.harmonic-magnetic-validation.v1",
            "status": "invalid_input",
            "pass": False,
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool(
    title="2-D transient runtime gate",
    annotations=ToolAnnotations(
        readOnlyHint=True,
        destructiveHint=False,
        idempotentHint=True,
        openWorldHint=False,
    ),
)
def fem_vol2d_transient_runtime(packet_json: str) -> str:
    """Run one digest-bound fixed-width field/thermal lifecycle operation.

    Operations are ``initialize``, ``step``, ``reset``, and ``terminate``.
    The returned state token is portable across Python, MATLAB MEX, and a
    Simulink S-function: it contains no file path or hidden server handle.
    Optional thermal states receive conductivity loss from the field step.
    """

    try:
        packet = json.loads(packet_json)
        if not isinstance(packet, dict):
            raise ValueError("packet_json must decode to an object")
        result = execute_transient_runtime(packet)
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "schema": "radia.vol2d-transient-runtime.v1",
            "status": "invalid_input",
            "pass": False,
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def fem_overview(topic: str = "decision_tree") -> str:
    """
    FEM landscape: lab stack, decision tree, genealogy.

    Args:
        topic: One of:
            "lab_stack"      - Production NGSolve + sparsesolv + Radia stack
            "decision_tree"  - Which formulation for which problem (DEFAULT)
            "history"        - Genealogy 1980 Nedelec -> 2025
            "all"            - Everything
    """
    return get_overview_knowledge(topic)


@mcp.tool()
def fem_potential_formulations(topic: str = "catalog") -> str:
    """
    Potential formulations: A-Omega, T-Omega, H, Reduced, Darwin.

    Args:
        topic: One of:
            "catalog"        - All formulations + decision tree (DEFAULT)
            "a_omega"        - A-Omega (mixed vector-scalar)
            "t_omega"        - T-Omega (electric T + magnetic Omega)
            "h_formulation"  - H-formulation (HCurl direct, SC / nonlinear)
            "reduced"        - ★ Reduced potential (accelerator magnet)
            "darwin"         - Darwin model (DC to MHz transition)
            "all"            - Everything
    """
    return get_potential_formulations_knowledge(topic)


@mcp.tool()
def fem_elements(topic: str = "catalog") -> str:
    """
    Element technology: edge (Nedelec), high-order, XFEM, isogeometric, DG.

    Args:
        topic: One of:
            "catalog"        - H1/HCurl/HDiv/L2 + hierarchical (DEFAULT)
            "edge"           - Edge elements (Nedelec) - foundation of HCurl
            "high_order"     - p-version (hierarchical basis)
            "xfem"           - eXtended FEM (cracks, interfaces)
            "isogeometric"   - Isogeometric Analysis (NURBS)
            "dg"             - Discontinuous Galerkin
            "all"            - Everything
    """
    return get_elements_knowledge(topic)


@mcp.tool()
def fem_gauge_open_boundary(topic: str = "gauge") -> str:
    """
    Gauging + open boundary techniques.

    Args:
        topic: One of:
            "gauge"             - Gauge choices for HCurl A (DEFAULT)
            "kelvin_transform"  - ★ Kelvin transformation (lab default)
            "abc_pml"           - Asymptotic BC and PML alternatives
            "all"               - Everything
    """
    return get_gauge_open_boundary_knowledge(topic)


@mcp.tool()
def fem_time_domain_axisym(topic: str = "henrotte_axisym") -> str:
    """
    Time-domain, axisymmetric (Henrotte), harmonic balance, HF, circuit coupling.

    Args:
        topic: One of:
            "time_domain_fem"    - TD-FEM time-stepping
            "henrotte_axisym"    - ★ Henrotte basis for axisym (lab CORE)
            "harmonic_balance"   - HB method for periodic steady-state
            "high_frequency"     - HF Helmholtz (reference, out of lab scope)
            "circuit_coupling"   - FE-circuit coupling (Hanser 2025 Schur)
            "all"                - Everything
    """
    return get_time_domain_axisym_knowledge(topic)


@mcp.tool()
def fem_ngsolve_hierarchy(topic: str = "decision_tree") -> str:
    """
    NGSolve hierarchical H(curl) bases - Zaglmayr / nograds / tree-cotree.

    Migrated 2026-05-23 from mcp_server_document/mathematica/notes_fem_hcurl.md
    per "all FEM knowledge lives in radia-mcp" policy.

    Args:
        topic: One of:
            "decision_tree"            - Which strategy for which question (DEFAULT)
            "zaglmayr_split"           - Hierarchical gradient + rotational split
            "nograds_option"           - NGSolve `nograds=True` behavior
            "tree_cotree"              - Tree-cotree gauge + generalizations
            "mathematica_verification" - Symbolic verification recipes (Wolfram)
            "all"                      - Everything
    """
    return get_ngsolve_hierarchy_knowledge(topic)


@mcp.tool()
def fem_xfem_em_hiruma(topic: str = "overview") -> str:
    """
    EM-XFEM (Hiruma 2023): electromagnetic XFEM for eddy-current
    skin-depth resolution.  Lab-verified implementation in
    docs/hiruma_xfem_comparison/ (beta-1 Phase 1-4, 2026-05-25).

    Reference:
        Hiruma et al., "Extended Finite Element Method for Eddy-Current
        Problems with Surface Skin Effect", IEEE Trans. Magn. 59(5), 2023.

    Independent NGSolve reproduction:
        88-DOF compound H1 x H1 + manual product-rule gradient gives
        0.14% error at r/delta = 15 on the Hiruma cylinder (Phase 2).

    Args:
        topic: One of:
            "overview"               - High-level summary, vs solid-mech XFEM (DEFAULT)
            "enrichment_function"    - psi = exp(-gamma*xi); frequency dependence;
                                       freezing strategy for ROM use
            "ngsolve_implementation" - NGSolve compound H1 x H1 pattern,
                                       three pitfalls (bonus_intorder,
                                       product-rule, enrichment scope)
            "cylinder_validation"    - Phase 2: 0.14% at r/delta=15, 88 DOFs
            "volume_source_scope"    - Phase 3: cures FE residual that breaks
                                       augmented CLN on volume-source
            "cln_stacking_negative"  - Phase 4: XFEM does NOT extend canonical
                                       CLN Hankel-QD stability (rejected)
            "decision_table"         - When to use EM-XFEM: layer + use case
            "ngsxfem_relation"       - vs the ngsxfem library (cut-FEM XFEM,
                                       Lehrenfeld et al.); when ngsxfem WOULD
                                       be appropriate (it is NOT for Hiruma
                                       PUFEM)
            "reproduction"           - Run order for Phase 1-4 benchmarks
            "lab_files"              - Paths to scripts and figures
            "all"                    - Everything (concatenated)
    """
    return get_em_xfem_knowledge(topic)


@mcp.tool()
def fem_large_scale_special(topic: str = "large_scale") -> str:
    """
    Large-scale, error theory, multi-scale (Hollaus MSFEM), misc techniques.

    Args:
        topic: One of:
            "large_scale"       - Parallel scaling, AMS/AMG (DEFAULT)
            "error_theory"      - A-posteriori estimators, adaptivity
            "hollaus_msfem"     - ★ Multi-scale FEM for laminated cores
            "misc_techniques"   - Fixed-point, polarization, Helmholtz-Weyl
            "all"               - Everything
    """
    return get_large_scale_special_knowledge(topic)


@mcp.tool()
def fem_nonconforming_mesh_coupling(topic: str = "overview") -> str:
    """
    Non-conforming mesh coupling: mortar / Nitsche / FETI-DP / BDDC / DG /
    HFSS Mesh Fusion analogs in NGSolve.

    Distilled from 33 PDFs across 6 subfolders at
    public-safe curated corpus (Sugahara Lab curated).

    Independent meshing of multiple regions coupled weakly at the
    interface -- the academic foundation for Ansys HFSS Mesh Fusion.
    Distinct from XFEM (single mesh + trial enrichment, see
    fem_xfem_em_hiruma).

    Args:
        topic: One of:
            "overview"          - When to use; decision matrix vs
                                  XFEM/PML/Kelvin (DEFAULT)
            "mortar"            - Lagrange-multiplier coupling
                                  (Bernardi-Maday-Patera 1994 +
                                  Buffa-Maday-Rapetti 2001 EM lineage)
            "nitsche"           - Penalty-based, no LM space
                                  (Nitsche 1971, Becker-Hansbo-
                                  Stenberg 2003)
            "feti_bddc"         - Iterative substructuring for HPC
                                  (Farhat 1991, Dohrmann 2003,
                                  Hofer-Langer IETI-DP for IGA)
            "dg"                - Discontinuous Galerkin as mesh-
                                  coupling (Arnold-Brezzi-Cockburn-
                                  Marini 2002)
            "hfss_mesh_fusion"  - Ansys product + academic analogs
                                  (Zhang-Liang 2015/2020 transfinite
                                  mortar)
            "ngsolve_recipe"    - Concrete NGSolve code: H1 mortar /
                                  Nitsche / HCurl Maxwell mortar
            "em_applications"   - EM-specific: rotor-stator, harmonic
                                  mortar (Egger 2020/2021), magnetic
                                  brake
            "lab_scenarios"     - Sugahara Lab adoption priority
                                  (rotor-stator > accelerator shim >
                                  IH coil > BEM-CLN)
            "bibliography"      - 33-PDF catalog organized by subfolder
            "all"               - Everything (~30 KB)

    Aliases: feti_dp/bddc -> feti_bddc, hfss/mesh_fusion ->
             hfss_mesh_fusion, ngsolve/recipe -> ngsolve_recipe,
             em/applications -> em_applications, lab/scenarios ->
             lab_scenarios, bib/papers -> bibliography.
    """
    return get_nonconforming_mesh_coupling_documentation(topic)


@mcp.tool()
def fem_equivalence_source(topic: str = "overview") -> str:
    """
    Equivalence-theorem near-field source (Schelkunoff/Love -- Stratton-Chu).

    Solve an FEM problem with whatever BC is convenient (typically
    Kelvin), then "record" E, H on a closed surface around the source
    and "replay" it via the Stratton-Chu surface integral to
        (a) probe the exterior field at any external point, or
        (b) re-radiate as a source in a downstream simulation
            (ngsolve.bem, Radia rad.Fld, external MoM via Nastran
            Near_Field_Area_*.dat).

    Production: ``radia.equivalence_source.NearFieldSource`` with C++
    accelerator at ``src/core/rad_equivalence_source.cpp`` (Phase A
    static, Phase B harmonic with full dyadic Green's function).

    See also:
        * docs/equivalence_source/USER_GUIDE.md
        * docs/equivalence_source/CPP_DESIGN.md
        * docs/equivalence_source/FMM_DESIGN.md (Phase D acceleration)
        * bem_low_freq("equivalence_source_lf") -- low-frequency rule:
          Weggler EFIE stabilization does NOT apply here (no matrix
          to invert); Laplace ML routing is the right tool.
        * bem_low_freq("weggler") -- the canonical ngsolve.bem LF
          stabilization (Lucy Weggler's product-space EFIE), for the
          SEPARATE case when the user computes (E, H) on the surface
          from a BEM solve instead of an FEM solve.

    Distilled from lab background on the equivalence theorem:
        - 2007-2008 axisymmetric basics
        - 2015 generalization (open-domain via the equivalence theorem)
        - a related WPT reconstruction study

    IABC / SDI content intentionally excluded -- Radia's Kelvin
    transformation is the canonical inner BC, and the 2015 lab work
    proved the equivalence reconstruction OUTSIDE is insensitive to
    the inner BC anyway.

    Args:
        topic: One of:
            "overview"               - scope + production module pointers (DEFAULT)
            "schelkunoff_love"       - Stratton-Chu math + static reduction
                                       (the rho_m = mu_0 n.H trap explained)
            "ngsolve_recipe"         - NGSolve extraction pattern
                                       (mesh.Boundaries + per-face sampling)
            "radia_recipe"           - Radia extract from rad.Fld + how to
                                       re-radiate via rad.ObjCnt
            "validation_use_cases"   - Golden tests + reference datasets
            "all"                    - Everything concatenated

    Aliases: schelkunoff/love/stratton_chu/math -> schelkunoff_love;
             ngsolve/extract -> ngsolve_recipe;
             radia/rad_fld/reradiate -> radia_recipe;
             validation/use_cases/tests -> validation_use_cases.
    """
    return get_equivalence_source_knowledge(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def pick_a_fem_formulation(problem_class: str) -> str:
    """Recommend FEM formulation for a class of EM problem."""
    guidance = {
        "accelerator_magnet": (
            "Accelerator electromagnet (★ lab core use case):\n"
            "1. Coil source: Biot-Savart from CoilBuilder (no coil mesh!)\n"
            "   → fem_potential_formulations('reduced')\n"
            "2. Iron yoke: Reduced Omega + nonlinear iron\n"
            "   → fem_potential_formulations('reduced')\n"
            "3. Open boundary: Kelvin transform\n"
            "   → fem_gauge_open_boundary('kelvin_transform')\n"
            "4. Solver: CG + CompactAMS (HCurl) or direct LU (small)\n"
            "   → matrix_solvers MCP\n"
            "5. Code recipe via radia\n"
            "   → radia_ngsolve.radia_usage, electromagnet MCP\n"
        ),
        "induction_heating_workpiece": (
            "IH workpiece (★ lab IH path):\n"
            "1. Workpiece: bulk eddy current with SIBC\n"
            "   → fem_potential_formulations('a_omega')\n"
            "   → radia_mcp.ih.sibc\n"
            "2. Coil: PEEC filaments (not FEM)\n"
            "   → radia_mcp.peec\n"
            "3. Open boundary: Kelvin\n"
            "   → fem_gauge_open_boundary('kelvin_transform')\n"
            "4. Frequency-domain: complex symmetric → COCR + ComplexCompactAMS\n"
            "   → matrix_solvers_krylov('cocg_cocr')\n"
            "   → matrix_solvers_preconditioners('ams_hiptmair_xu')\n"
            "5. Code: calc_fem_kelvin.py\n"
        ),
        "motor_transient": (
            "Motor transient (PMSM/IM/SRM):\n"
            "1. Vector A with HCurl + Picard nonlinear (iron BH)\n"
            "   → fem_potential_formulations('a_omega')\n"
            "2. Time-stepping: implicit Euler / Crank-Nicolson\n"
            "   → fem_time_domain_axisym('time_domain_fem')\n"
            "3. Rotor rotation: rotor PM via Br(theta) CoefficientFunction\n"
            "   → radia_mcp.motor\n"
            "4. Circuit coupling: Hanser 2025 Schur voltage drive\n"
            "   → fem_time_domain_axisym('circuit_coupling')\n"
            "5. Laminated stator: Hollaus MSFEM effective material\n"
            "   → fem_large_scale_special('hollaus_msfem')\n"
            "6. Code: calc_motor_transient.py, calc_motor_lamination.py\n"
        ),
        "axisymmetric_magnetic": (
            "Axisymmetric magnetic A_phi (★ Henrotte CORE):\n"
            "1. Use Henrotte basis {1, r^2, z}\n"
            "   → fem_time_domain_axisym('henrotte_axisym')\n"
            "2. radia.axifem.H1Henrotte\n"
            "3. Standard H1 for scalar T, V, etc. (NOT Henrotte)\n"
            "4. Code: radia_ngsolve.axifem_documentation\n"
        ),
        "open_boundary_static": (
            "Open boundary magnetostatic (★ lab default):\n"
            "1. Kelvin transform (radia + NGSolve, pure FEM)\n"
            "   → fem_gauge_open_boundary('kelvin_transform')\n"
            "2. NOT FEM-BEM hybrid (extra complexity)\n"
            "3. NOT PML (fails at DC)\n"
            "4. Code: calc_fem_kelvin.py, kelvin_transformation MCP\n"
        ),
        "eddy_current_mqs": (
            "Eddy current MQS (frequency domain):\n"
            "1. Bulk conductor, multiply-connected → A-V Biro-Preis\n"
            "   → fem_potential_formulations('catalog')\n"
            "2. Surface skin only → A + SIBC Robin BC (lab IH)\n"
            "   → radia_mcp.ih.sibc\n"
            "3. Solver: COCR + ComplexCompactAMS\n"
            "   → matrix_solvers_krylov('cocg_cocr')\n"
            "4. Low-freq stabilization: shifted preconditioner\n"
            "   → matrix_solvers_preconditioners('shifted')\n"
        ),
    }
    return guidance.get(problem_class, (
        f"Unknown problem_class '{problem_class}'. Available:\n"
        "  accelerator_magnet, induction_heating_workpiece, motor_transient,\n"
        "  axisymmetric_magnetic, open_boundary_static, eddy_current_mqs.\n"
        "For overview, see fem_overview('decision_tree').\n"
    ))




register_status_tool(
    mcp,
    server_name='mcp-server-fem',
    description='FEM formulations theory layer (A-Omega / T-Omega / H / Reduced / Darwin, edge / HO / XFEM / IGA / DG, gauging + Kelvin, MSFEM, Schur...',
    subpackage='radia_mcp.fem',
    related_servers=["radia-ngsolve", "bem", "differential-forms"],
)

apply_tool_contract(
    mcp,
    server_name="mcp-server-fem",
    version="1.4.19",
)


def main():
    """Entry point for mcp-server-fem."""
    if "--selftest" in sys.argv:
        print("fem MCP server self-test:")
        print(f"  overview: {len(get_overview_knowledge('all'))} chars")
        print(f"  potential_formulations: {len(get_potential_formulations_knowledge('all'))} chars")
        print(f"  elements: {len(get_elements_knowledge('all'))} chars")
        print(f"  gauge_open_boundary: {len(get_gauge_open_boundary_knowledge('all'))} chars")
        print(f"  time_domain_axisym: {len(get_time_domain_axisym_knowledge('all'))} chars")
        print(f"  large_scale_special: {len(get_large_scale_special_knowledge('all'))} chars")
        print(f"  ngsolve_hierarchy: {len(get_ngsolve_hierarchy_knowledge('all'))} chars")
        print(f"  xfem_em_hiruma: {len(get_em_xfem_knowledge('all'))} chars")
        print(f"  nonconforming_mesh_coupling: {len(get_nonconforming_mesh_coupling_documentation('all'))} chars")
        print(f"  equivalence_source: {len(get_equivalence_source_knowledge('all'))} chars")
        print("OK")
        return
    mcp.run()


if __name__ == "__main__":
    main()
