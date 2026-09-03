"""Keep HDiv-MMM environment controls finite, classified, and observable."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]

DIAGNOSTIC = {
    "RADIA_HDIV_BLOCK_CACHE_STATS",
    "RADIA_HDIV_HMATVEC_STATS",
}
FAILURE_INJECTION = {"RADIA_HDIV_TEST_FAIL_FILL_AFTER"}
PERFORMANCE = {
    "RADIA_HDIV_DISABLE_TRANS_CACHE",
    "RADIA_HDIV_HEX_BLOCK_CACHE_LIMIT",
    "RADIA_HDIV_WEDGE_TRANS_CACHE",
}
NUMERICAL_PATH = {
    "RADIA_HDIV_AUTO_JACOBI_TET_NFACE",
    "RADIA_HDIV_CURVED_DIRECT",
    "RADIA_HDIV_DISABLE_HO_ANALYTIC_BLOCK",
    "RADIA_HDIV_DISABLE_HO_IMAGE_BLOCK",
    "RADIA_HDIV_HEX_DISTORTED_FAR_FACTOR",
    "RADIA_HDIV_HEX_FAR_ONESIDED",
    "RADIA_HDIV_HO_FAR_ONESIDED",
    "RADIA_HDIV_WEDGE_FAR_ONESIDED",
}


def _environment_names():
    sources = (
        ROOT / "src" / "core" / "rad_hacapk_hdiv.cpp",
        ROOT / "src" / "ext" / "HACApK" / "cHACApK_cpp_impl.c",
        ROOT / "src" / "radia" / "vim" / "_solve.py",
    )
    return set().union(*(
        set(re.findall(r"RADIA_HDIV_[A-Z0-9_]+", path.read_text(encoding="utf-8")))
        for path in sources
    ))


def test_every_hdiv_environment_control_is_classified():
    classified = DIAGNOSTIC | FAILURE_INJECTION | PERFORMANCE | NUMERICAL_PATH
    assert _environment_names() == classified
    assert sum(map(len, (DIAGNOSTIC, FAILURE_INJECTION, PERFORMANCE, NUMERICAL_PATH))) == len(classified)


def test_retired_hex_cache_stats_alias_stays_removed():
    assert "RADIA_HDIV_HEX_CACHE_STATS" not in _environment_names()


def test_native_stats_expose_release_relevant_effective_controls():
    source = (ROOT / "src" / "core" / "rad_hacapk_hdiv.cpp").read_text(encoding="utf-8")
    expected_fields = {
        "hex_block_cache_limit",
        "hex_trans_cache_enabled",
        "hex_far_one_sided_threshold",
        "hex_distorted_far_factor",
        "curved_direct_enabled",
        "wedge_far_one_sided_threshold",
        "wedge_trans_cache_scope",
        "ho_far_one_sided_enabled",
        "ho_analytic_block_enabled",
        "ho_image_block_enabled",
        "hmatvec_stats_enabled",
        "hmatvec_mkl_threads",
        "nonproduction_numerical_override_active",
        "nondefault_performance_override_active",
        "release_claim_eligible",
    }
    for field in expected_fields:
        assert '"%s"' % field in source


def test_review_records_every_hdiv_environment_control():
    review = (ROOT / "docs" / "hdiv_vim" / "HDiv-MMM_review.md").read_text(encoding="utf-8")
    for name in DIAGNOSTIC | FAILURE_INJECTION | PERFORMANCE | NUMERICAL_PATH:
        assert name in review
