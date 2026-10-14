import json

from radia_mcp.radia_ngsolve.hmatrix_scaling_gate import hmatrix_compression_scaling_gate
from radia_mcp.radia_ngsolve.server import hmatrix_compression_scaling_gate as mcp_gate


ROWS = [
    {"pointCount": 60, "maxRank": 5, "lowRankBlocks": 1, "storedEntries": 600,
     "denseEntries": 3600, "compressionRatio": 1/6, "matvecRelativeError": 3.95e-13,
     "buildSeconds": 0.117, "matvecSeconds": 0.0022, "denseReferenceSeconds": 0.0018},
    {"pointCount": 120, "maxRank": 5, "lowRankBlocks": 1, "storedEntries": 1200,
     "denseEntries": 14400, "compressionRatio": 1/12, "matvecRelativeError": 6.25e-13,
     "buildSeconds": 0.0088, "matvecSeconds": 0.0012, "denseReferenceSeconds": 0.0007},
    {"pointCount": 240, "maxRank": 5, "lowRankBlocks": 1, "storedEntries": 2400,
     "denseEntries": 57600, "compressionRatio": 1/24, "matvecRelativeError": 2.44e-13,
     "buildSeconds": 0.0033, "matvecSeconds": 0.0007, "denseReferenceSeconds": 0.0010},
]


def test_hmatrix_scaling_accepts_linear_storage_and_dense_accuracy():
    result = hmatrix_compression_scaling_gate(ROWS)
    assert result["status"] == "ok"
    assert result["metrics"]["max_ranks"] == [5, 5, 5]
    assert result["metrics"]["storage_growth_exponents"] == [1.0, 1.0]


def test_hmatrix_scaling_rejects_rank_and_quadratic_storage_regression():
    bad = [dict(row) for row in ROWS]
    bad[2].update(maxRank=40, storedEntries=40000, compressionRatio=40000/57600)
    result = hmatrix_compression_scaling_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["rank_is_bounded"] is False
    assert result["checks"]["compression_improves_with_size"] is False
    assert result["checks"]["stored_entry_growth_is_subquadratic"] is False


def test_hmatrix_scaling_rejects_stale_ratio_and_inaccurate_matvec():
    bad = [dict(row) for row in ROWS]
    bad[1].update(compressionRatio=0.5, matvecRelativeError=1e-3)
    result = hmatrix_compression_scaling_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["compression_ratio_matches_storage"] is False
    assert result["checks"]["matvec_matches_dense_reference"] is False


def test_hmatrix_scaling_mcp_dispatches_json_rows():
    result = json.loads(mcp_gate(json.dumps(ROWS)))
    assert result["status"] == "ok"
    assert result["policy"] == "hmatrix_compression_scaling_gate_v1"
