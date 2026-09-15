"""Cubit geometry-query terminology contract."""
from cubit_mesh_export.mcp.api_reference import get_api_reference


def test_cubit_geometry_reference_distinguishes_center_from_mass_centroid():
    docs = get_api_reference("geometry_queries")
    assert "representative geometric center" in docs
    assert "volume_id).centroid()" in docs
    assert "CenterOf.MASS" in docs
