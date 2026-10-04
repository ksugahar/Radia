"""Keep executable MCP guidance aligned with the current solver boundary."""

import re

from radia_mcp.bem.h_matrix_knowledge import get_h_matrix_knowledge
from radia_mcp.bem.overview_knowledge import get_overview_knowledge
from radia_mcp.electromagnet.em_knowledge import get_electromagnet_documentation
from radia_mcp.matrix_solvers.direct_solvers_knowledge import (
    get_direct_solvers_knowledge,
)
from radia_mcp.matrix_solvers.overview_knowledge import (
    get_overview_knowledge as get_matrix_overview_knowledge,
)
from radia_mcp.radia_ngsolve.knowledge.radia import get_radia_documentation


_RETIRED_SOLVE_CALL = re.compile(
    r"rad\.Solve\([^\n]*,\s*[12](?:\s*[,)]|\s*#)"
)


def _current_solver_guidance() -> str:
    topics = [
        get_radia_documentation("solving"),
        get_radia_documentation("parallelization"),
        get_radia_documentation("best_practices"),
        get_radia_documentation("play_models"),
        get_radia_documentation("hysteresis"),
        get_direct_solvers_knowledge("lu_radia"),
        get_overview_knowledge("decision_tree"),
        get_h_matrix_knowledge("hacapk"),
        get_electromagnet_documentation("ima"),
    ]
    return "\n".join(topics)


def test_urn_manual_keeps_private_inputs_and_noise_aware_basis_selection():
    from radia_mcp.radia_ngsolve.knowledge.urn import get_urn_documentation
    guidance = get_urn_documentation('overview')
    assert 'RADIA_NASA_EIS_CSV' in guidance
    assert '34-basis' in guidance
    assert 'uncertainty_ohm' in guidance
    assert 'down to ONE basis' in guidance
    assert 'retained basis count' in guidance
    assert 'Do not force an arbitrary 1e-3' in guidance
    assert 'Fewer bases alone do not prove physical truth' in guidance
    assert 'model-assigned' in guidance


def test_mcp_solver_guidance_does_not_call_retired_relaxation_methods():
    guidance = _current_solver_guidance()

    assert _RETIRED_SOLVE_CALL.search(guidance) is None
    assert "HDivSolver" in guidance
    assert "dense LU" in guidance


def test_mcp_hysteresis_guidance_uses_mesh_backed_history_solver():
    guidance = get_radia_documentation("hysteresis")

    assert "vim.SolveHysteresis" in guidance
    assert "vim.HDivSolver" in guidance
    assert "Mesh-less soft-iron `rad.Solve` is retired" in guidance
    assert "rad.MatApl(" not in guidance


def test_mcp_result_ownership_distinguishes_docs_from_validation():
    guidance = get_radia_documentation("best_practices")
    flat_guidance = " ".join(guidance.split())

    assert "stores its result and WebGUI output in the notebook itself" in flat_guidance
    assert "does not require a result JSON sidecar" in flat_guidance
    assert "`validation_test/` owns" in guidance


def test_legacy_direct_topic_explains_migration_without_enabling_pardiso():
    guidance = get_direct_solvers_knowledge("pardiso")
    flat_guidance = " ".join(guidance.split())

    assert "ngsolve-openblas" in guidance
    assert "mkl>=2026,<2027" in guidance
    assert "does not select PARDISO" in flat_guidance
    assert "Use SparseCholesky explicitly" in flat_guidance
    assert 'inverse="pardiso"' not in guidance


def test_matrix_solver_guidance_preserves_the_direct_solve_contract():
    direct = get_direct_solvers_knowledge("sparsecholesky")
    all_direct = get_direct_solvers_knowledge("all")
    decision = get_matrix_overview_knowledge("decision_tree")

    assert direct in all_direct
    assert 'inverse="sparsecholesky"' in direct
    assert "Projector(fes.FreeDofs(), True)" in direct
    assert "f.vec - a.mat * gfu.vec" in direct
    assert "isfinite(relative_residual)" in direct
    assert "same native" in direct
    assert 'inverse="pardiso"' not in all_direct
    assert "BDDC + AMS" in decision
    assert "Do not silently fall back" in decision
    assert "PARDISO complex direct" not in decision
    assert "radia.Solve method=0" not in decision


def test_radia_guidance_uses_current_demo_and_validation_ownership():
    geometry = get_radia_documentation("geometry")
    esim = get_radia_documentation("esim")

    assert "docs/background_fields/background_fields_results.json" not in geometry
    assert "results are embedded in the notebook" in geometry
    assert "docs/esim/WORKFLOWS.md" in esim
    assert "validation_test/induction_heating/cubit_panels_legacy" in esim
    assert "examples/.../" not in esim


def test_bem_inductance_guidance_requires_physical_current_constraints():
    from radia_mcp.radia_ngsolve.knowledge.ngsbem_inductance import NGBEM_OVERVIEW
    from radia_mcp.radia_ngsolve.knowledge.ngsolve import NGSOLVE_BEM

    assert "rank 269/269" in NGBEM_OVERVIEW
    assert "transforming both A and e" in NGBEM_OVERVIEW
    assert "Do NOT fix this" not in NGBEM_OVERVIEW
    assert "GetDofNrs" in NGBEM_OVERVIEW
    assert "c(J)=1 A" in NGBEM_OVERVIEW
    assert "flat-mesh convergence" in NGBEM_OVERVIEW
    assert "including at order=0 (RT0)" in NGSOLVE_BEM
    assert "only valid for order=0" not in NGSOLVE_BEM


def test_capacity_is_retrievable_through_existing_public_tools():
    from radia_mcp.matrix_solvers.server import matrix_solvers_direct, pick_a_solver
    from radia_mcp.radia_ngsolve.server import ngsolve_usage, sparsesolv, ngsolve_solvers_reference
    from radia_mcp.ih.server import induction_heating

    for tool in (matrix_solvers_direct, ngsolve_usage, sparsesolv, induction_heating):
        for topic in ("memory", "capacity", "diagnostics", "bddc"):
            text = tool(topic)
            assert "Unknown topic" not in text
            assert "factor-entry storage" in text
            assert "commit limit minus committed bytes" in text
            assert "outer solver -> preconditioner -> coarse solver" in text
            assert "not yet implemented" in text
            assert "memory-confounded" in text
            assert "no silent" in text
        assert "Hermitian positive-definite" in tool("memory")
        assert "Maxwell AMS is not an H1 prescription" in tool("memory")
    for problem in ("magnetostatic_h1", "magnetostatic_hcurl", "eddy_current_mqs",
                    "frequency_sweep", "non_symmetric", "indefinite_saddle_point"):
        assert "matrix_solvers_direct('memory')" in pick_a_solver(problem)
    assert "ngsolve_usage('memory')" in ngsolve_solvers_reference()
    assert "Residual checks do not protect allocation" in ngsolve_solvers_reference()
    h1 = pick_a_solver("magnetostatic_h1")
    assert "N<" not in h1
    assert "method=0" not in h1


def test_default_solver_guidance_reports_capacity_without_size_heuristics():
    from radia_mcp.radia_ngsolve.knowledge.ngsolve import get_ngsolve_documentation
    from radia_mcp.radia_ngsolve.knowledge.sparsesolv import get_sparsesolv_documentation
    for text in (get_direct_solvers_knowledge(),
                 get_matrix_overview_knowledge("decision_tree"),
                 get_ngsolve_documentation("solvers"),
                 get_sparsesolv_documentation("overview")):
        assert "Direct solvers can be fast when the factor fits" in text
        assert "An iterative outer solve with BDDC can still use a direct" in text
        assert "true linear residual gate remains 1e-6" in text
        assert "<10K DOF" not in text and ">100K DOF" not in text
    practices = get_sparsesolv_documentation("best_practices")
    assert "Krylov workspace (not total memory)" in practices
    assert "ICCG (fallback)" not in practices
