import copy
import json
import math

import pytest

from radia_mcp.differential_forms.bibliography_knowledge import (
    get_bibliography_documentation,
)
from radia_mcp.differential_forms.server import (
    differential_forms_geometry_gate,
    differential_forms_starter,
    differential_forms_visual_geometry,
)
from radia_mcp.differential_forms.visual_geometry_gate import (
    NORMALIZATION,
    SCHEMA,
    evaluate_visual_geometry,
)
from radia_mcp.differential_forms.visual_geometry_knowledge import (
    get_visual_geometry_documentation,
)


def mapped_em_summary() -> dict:
    return {
        "schema": SCHEMA,
        "profile": "mapped_em",
        "normalization_policy": NORMALIZATION,
        "de_rham": {
            "curl_grad_relative_residual": 2.0e-14,
            "div_curl_relative_residual": 3.0e-14,
        },
        "pullback": {
            "exterior_derivative_commutator_relative_residual": 4.0e-13,
        },
        "hodge": {
            "symmetry_relative_error": 5.0e-14,
            "minimum_eigenvalue": 0.25,
        },
        "maxwell": {"dF_relative_residual": 6.0e-13},
        "tolerances": {
            "de_rham_relative": 1.0e-12,
            "pullback_relative": 1.0e-11,
            "hodge_symmetry_relative": 1.0e-12,
            "hodge_positive_min": 1.0e-9,
            "maxwell_relative": 1.0e-11,
        },
    }


def surface_summary() -> dict:
    return {
        "schema": SCHEMA,
        "profile": "surface",
        "normalization_policy": NORMALIZATION,
        "orientation_convention": "outward normal; positive boundary orientation",
        "cartan_sign_convention": "T=dtheta+omega^theta; Omega=domega+omega^omega",
        "hodge": {
            "symmetry_relative_error": 2.0e-14,
            "minimum_eigenvalue": 0.5,
        },
        "connection": {
            "metric_compatibility_relative_residual": 3.0e-13,
            "torsion_relative_residual": 3.5e-13,
            "curvature_form_relative_residual": 4.0e-13,
            "bianchi_relative_residual": 4.5e-13,
        },
        "surface": {
            "curvature_integral_rad": 0.5 * math.pi,
            "boundary_geodesic_curvature_rad": 1.5 * math.pi,
            "corner_turning_rad": 0.0,
            "euler_characteristic": 1,
        },
        "holonomy": {
            "angle_rad": 0.5 * math.pi,
            "curvature_integral_rad": 0.5 * math.pi,
            "orientation_sign": 1,
        },
        "tolerances": {
            "hodge_symmetry_relative": 1.0e-12,
            "hodge_positive_min": 1.0e-9,
            "connection_relative": 1.0e-11,
            "gauss_bonnet_rad": 1.0e-10,
            "holonomy_rad": 1.0e-10,
        },
    }


def full_summary() -> dict:
    mapped = mapped_em_summary()
    surface = surface_summary()
    mapped["profile"] = "full"
    for key in (
        "orientation_convention",
        "cartan_sign_convention",
        "connection",
        "surface",
        "holonomy",
    ):
        mapped[key] = surface[key]
    mapped["tolerances"].update(surface["tolerances"])
    return mapped


def test_visual_geometry_topics_connect_geometry_to_em_spaces() -> None:
    doc = get_visual_geometry_documentation("all")
    for phrase in (
        "The five-act map",
        "Intrinsic versus extrinsic",
        "Gauss-Bonnet",
        "Cartan's structure equations",
        "dF = 0",
        "H1/HCurl/HDiv/L2",
        "Levi-Civita holonomy",
        "de Rham period",
        "gauge-theoretic",
    ):
        assert phrase in doc
    assert "W:\\" not in doc


def test_visual_geometry_aliases_and_server_wrapper() -> None:
    assert "Hodge" in get_visual_geometry_documentation("metric")
    assert "Gauss-Bonnet" in differential_forms_visual_geometry("holonomy")
    unknown = get_visual_geometry_documentation("not-a-topic")
    assert "Unknown topic" in unknown
    assert "source_scope" in unknown


def test_visual_geometry_source_is_registered_without_private_path() -> None:
    bibliography = get_bibliography_documentation()
    assert "T. Needham" in bibliography
    assert "978-4-621-31240-7" in bibliography
    assert "five\nacts: space, metric, curvature" in bibliography
    assert "W:\\" not in bibliography


def test_geometry_starter_exposes_geometry_first_path() -> None:
    starter = differential_forms_starter("geometry")
    assert "five_acts" in starter
    assert "curvature_holonomy" in starter
    assert "radia_workflow" in starter
    assert "topology, metric" in starter


def test_mapped_em_gate_passes_and_server_returns_json() -> None:
    summary = mapped_em_summary()
    result = evaluate_visual_geometry(summary)
    assert result["status"] == "ok"
    assert all(result["checks"].values())
    wrapped = json.loads(differential_forms_geometry_gate(json.dumps(summary)))
    assert wrapped["status"] == "ok"
    assert wrapped["schema"] == SCHEMA


def test_surface_gate_checks_cartan_gauss_bonnet_and_wrapped_holonomy() -> None:
    result = evaluate_visual_geometry(surface_summary())
    assert result["status"] == "ok"
    assert result["checks"]["connection_is_metric_compatible"] is True
    assert result["checks"]["cartan_first_structure_equation_is_torsion_free"] is True
    assert result["checks"]["cartan_second_structure_equation_holds"] is True
    assert result["checks"]["bianchi_identity_holds"] is True
    assert result["checks"]["gauss_bonnet_closes"] is True
    assert result["checks"]["holonomy_matches_integrated_curvature"] is True


def test_full_gate_runs_every_geometry_family_together() -> None:
    result = evaluate_visual_geometry(full_summary())
    assert result["status"] == "ok"
    assert len(result["checks"]) == 12
    assert all(result["checks"].values())


def test_gate_reports_independent_structural_failures() -> None:
    bad = copy.deepcopy(mapped_em_summary())
    bad["de_rham"]["div_curl_relative_residual"] = 1.0e-3
    bad["pullback"]["exterior_derivative_commutator_relative_residual"] = 2.0e-3
    bad["hodge"]["minimum_eigenvalue"] = -0.1
    bad["maxwell"]["dF_relative_residual"] = 3.0e-3
    result = evaluate_visual_geometry(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["div_curl_is_zero"] is False
    assert result["checks"]["pullback_commutes_with_exterior_derivative"] is False
    assert result["checks"]["hodge_is_positive"] is False
    assert result["checks"]["homogeneous_maxwell_form_is_closed"] is False


@pytest.mark.parametrize(
    ("key", "message"),
    [
        ("orientation_convention", "orientation_convention"),
        ("cartan_sign_convention", "cartan_sign_convention"),
    ],
)
def test_surface_gate_requires_declared_conventions(key: str, message: str) -> None:
    bad = surface_summary()
    del bad[key]
    with pytest.raises(ValueError, match=message):
        evaluate_visual_geometry(bad)


def test_gate_rejects_wrong_schema_nonfinite_values_and_zero_positive_floor() -> None:
    bad_schema = mapped_em_summary()
    bad_schema["schema"] = "wrong"
    with pytest.raises(ValueError, match="schema"):
        evaluate_visual_geometry(bad_schema)

    nonfinite = mapped_em_summary()
    nonfinite["maxwell"]["dF_relative_residual"] = math.nan
    with pytest.raises(ValueError, match="finite"):
        evaluate_visual_geometry(nonfinite)

    zero_floor = mapped_em_summary()
    zero_floor["tolerances"]["hodge_positive_min"] = 0.0
    with pytest.raises(ValueError, match="must be positive"):
        evaluate_visual_geometry(zero_floor)

    with pytest.raises(ValueError, match="valid JSON"):
        differential_forms_geometry_gate("{broken")
