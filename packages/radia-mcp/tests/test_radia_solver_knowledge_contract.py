"""Keep executable MCP guidance aligned with the current solver boundary."""

import re

from radia_mcp.bem.h_matrix_knowledge import get_h_matrix_knowledge
from radia_mcp.bem.overview_knowledge import get_overview_knowledge
from radia_mcp.electromagnet.em_knowledge import get_electromagnet_documentation
from radia_mcp.matrix_solvers.direct_solvers_knowledge import (
    get_direct_solvers_knowledge,
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


def test_pardiso_guidance_distinguishes_ngsolve_openblas_from_radia_mkl():
    guidance = get_direct_solvers_knowledge("pardiso")
    flat_guidance = " ".join(guidance.split())

    assert "ngsolve-openblas" in guidance
    assert "mkl>=2026,<2027" in guidance
    assert "does not bundle or select MKL" in flat_guidance
    assert "a compatible MKL runtime is a separate environment dependency" in flat_guidance
    assert "included with NGSolve PyPI wheel" not in guidance
    assert "NGSolve binary distributions automatically" not in guidance
