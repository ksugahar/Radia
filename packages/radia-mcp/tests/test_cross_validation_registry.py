# -*- coding: utf-8 -*-
"""Cross-validation registry knowledge should stay queryable by MCP tools."""

from radia_mcp.radia_ngsolve.knowledge.ngsolve import get_ngsolve_documentation
from radia_mcp.radia_ngsolve.server import ngsolve_usage


def test_ngsolve_cross_validation_registry_topic_records_reusable_artifacts():
    doc = get_ngsolve_documentation("cross_validation_registry")

    assert "validation_*.py" in doc
    assert "validation_*_summary.json" in doc
    assert "examples/acoustic_bem/validation_*_summary.json" in doc
    assert "examples/cubit_mesh_export/validation_vol_*_summary.json" in doc
    assert "examples/build123d_netgen_gmsh_flow/validation_*_summary.json" in doc
    assert "examples/electric_machine/validation_*_summary.json" in doc
    assert "examples/rf_waveguide/validation_*_summary.json" in doc
    assert 'force_validation("cross_validation")' in doc
    assert 'fem_bem_schur("api")' in doc


def test_ngsolve_cross_validation_registry_is_exposed_through_mcp_tool():
    doc = ngsolve_usage("cross_validation_registry")

    assert "Cross-validation registry" in doc
    assert "official external `optuna/optuna-mcp`" in doc
    assert "`radia_mcp.optuna` package" in doc
    assert "policy_lint.py" in doc
    assert "Unknown topic" not in doc[:120]

