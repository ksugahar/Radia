import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

from validation_test.maglev.team28_hcurl_eddy_bubble import run


def test_team28_p6_hcurl_eddy_bubble_structural_acceptance():
    result = run()

    assert result["structural_and_reference_acceptance_passed"] is True
    assert all(result["checks"].values())
    assert result["surface_model"]["selected_model"] == "volumetric"
    assert result["surface_model"]["selected_sibc_face_count"] == 0
    assert result["p6_spatial_reduction"]["parent_order"] == 6
    assert result["reduced_mode_summary"]["parent_ndof"] > 20_000
    assert result["reduced_mode_summary"]["estimated_total_modes"] < 200
    assert result["reduced_mode_summary"]["estimated_reduction_ratio"] < 0.01
    assert result["cln_reference_acceptance"]["cln_vs_full_fem"][
        "max_abs_error_N"
    ] < 5.0e-6
    assert result["hcurl_vim_force_acceptance_complete"] is True
    assert result["hcurl_vim_force_acceptance"]["validation_host"] == "mdx"
    assert result["hcurl_vim_force_acceptance"][
        "maximum_force_relative_error"
    ] < 0.01
