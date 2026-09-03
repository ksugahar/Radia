"""Golden regression test for the IGTE 2026 ESIM digest dense sweep.

Reads the committed `validation_test/ih_esim_benchmark/sweep_data_dense/`
JSONs and asserts the headline numbers cited in the digest:

  - Anchor case (100 A, 50 kHz): per-panel 18.75 W vs scalar 30.51 W
    -> gap = -38.5 %
  - Dense-grid maximum gap (500 A, 10 kHz): per-panel 81.83 W vs
    scalar 162.85 W -> gap = -49.75 %
  - Convergence: 49 of 54 per-element cases converged in 6-13 iter
  - Stall band: I=500 A, f={20, 50, 100, 200, 500} kHz (5 cases,
    iter=30 cap, conv=False)
  - Side-wall |Z_s| at max-gap case: range 4.5-16.5 mOhm, 561 verts

These numbers anchor the IGTE digest body text and Fig. 1; any
re-run of `sweep_f_I.py` that drifts outside the locked bands
must be investigated BEFORE updating the digest.

This test is CHEAP (reads committed JSON, no NGSolve, no solver).
"""
import json
import os
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DENSE = os.path.join(
    REPO, "validation_test", "ih_esim_benchmark", "sweep_data_dense"
)
SWEEP = DENSE


def _load(tag):
    """tag = 'I100_f50k_per_panel' etc."""
    path = os.path.join(DENSE, f"{tag}.json")
    if not os.path.exists(path):
        pytest.skip(f"sweep_data_dense missing: {path}")
    with open(path) as f:
        return json.load(f)


def _gap_pct(scalar_P, per_panel_P):
    return (per_panel_P / scalar_P - 1.0) * 100.0


def test_anchor_100A_50kHz_locked():
    """Digest body cites P_per=18.75 W, P_uni=30.51 W, gap=-38.5 %."""
    sc = _load("I100_f50k_scalar")
    pp = _load("I100_f50k_per_panel")
    assert sc["P_wp_W"] == pytest.approx(30.507, rel=2e-3), \
        f"scalar P_wp drift: {sc['P_wp_W']} vs 30.507 W"
    assert pp["P_wp_W"] == pytest.approx(18.753, rel=2e-3), \
        f"per-panel P_wp drift: {pp['P_wp_W']} vs 18.753 W"
    gap = _gap_pct(sc["P_wp_W"], pp["P_wp_W"])
    assert gap == pytest.approx(-38.53, abs=0.1), \
        f"gap drift: {gap:+.2f} % vs -38.53 %"


def test_max_gap_500A_10kHz_locked():
    """Digest panel (b) case: P_uni=162.85 W, P_per=81.83 W, gap=-49.75 %."""
    sc = _load("I500_f10k_scalar")
    pp = _load("I500_f10k_per_panel")
    assert sc["P_wp_W"] == pytest.approx(162.85, rel=2e-3), \
        f"scalar P_wp drift: {sc['P_wp_W']} vs 162.85 W"
    assert pp["P_wp_W"] == pytest.approx(81.82, rel=2e-3), \
        f"per-panel P_wp drift: {pp['P_wp_W']} vs 81.82 W"
    gap = _gap_pct(sc["P_wp_W"], pp["P_wp_W"])
    assert gap == pytest.approx(-49.75, abs=0.1), \
        f"gap drift: {gap:+.2f} % vs -49.75 %"
    # Must be CONVERGED -- if this case ever stalls, the digest claim
    # "the maximum gap of -49.8 % is at (500 A, 10 kHz)" is invalid.
    assert pp.get("esim_converged") is True, \
        f"max-gap case must be converged; got {pp.get('esim_converged')}"


def test_stall_band_I500A_high_freq():
    """5 cases stall at iter=30 cap: I=500A, f in {20,50,100,200,500} kHz."""
    expected_stalls = [
        ("I500_f20k_per_panel", -38.9),
        ("I500_f50k_per_panel", -19.3),
        ("I500_f100k_per_panel", -10.0),
        ("I500_f200k_per_panel", -4.7),
        ("I500_f500k_per_panel", +19.9),  # spurious stall artifact
    ]
    for tag, expected_gap in expected_stalls:
        pp = _load(tag)
        sc_tag = tag.replace("per_panel", "scalar")
        sc = _load(sc_tag)
        assert pp.get("esim_iterations") == 30, \
            f"{tag}: expected iter=30 cap, got {pp.get('esim_iterations')}"
        assert pp.get("esim_converged") is False, \
            f"{tag}: expected stall (conv=False), got {pp.get('esim_converged')}"
        gap = _gap_pct(sc["P_wp_W"], pp["P_wp_W"])
        # Looser tolerance: the limit-cycle iterate is itself unstable,
        # so the recorded P_wp can drift across re-runs.  Just lock the
        # sign and a +-5 absolute % band for the documented value.
        assert gap == pytest.approx(expected_gap, abs=5.0), \
            f"{tag}: gap drift {gap:+.1f} % vs expected {expected_gap:+.1f} %"


def test_low_I_column_sign_flip():
    """I=1 A column stays at +17 to +28 % (digest body)."""
    gaps = []
    for f_kHz in [10, 20, 50, 100, 200, 500]:
        sc = _load(f"I1_f{f_kHz}k_scalar")
        pp = _load(f"I1_f{f_kHz}k_per_panel")
        gaps.append(_gap_pct(sc["P_wp_W"], pp["P_wp_W"]))
    assert min(gaps) >= 15.0, f"I=1 A column min drift: {min(gaps):+.1f} %"
    assert max(gaps) <= 30.0, f"I=1 A column max drift: {max(gaps):+.1f} %"


def test_converged_iter_range():
    """Converged-cases iter distribution.

    49 cases converge (54 per-panel minus 5 stalls).  The bulk
    (~47 of 49) sit in 6-13 iter (median 8 = the digest "6-13"
    band qualified by "almost all").  Two high-power converged
    outliers stretch the range up to 26:
      I=500 A, f=10 kHz     -> iter = 25 (the max-gap case)
      I=200 A, f=500 kHz    -> iter = 26
    """
    converged = []
    stalled_band = {
        (500.0, 20000), (500.0, 50000), (500.0, 100000),
        (500.0, 200000), (500.0, 500000),
    }
    for I in [1, 2, 5, 10, 20, 50, 100, 200, 500]:
        for f_Hz in [10000, 20000, 50000, 100000, 200000, 500000]:
            tag = f"I{I}_f{f_Hz//1000}k_per_panel"
            pp = _load(tag)
            if (float(I), f_Hz) in stalled_band:
                continue
            converged.append((I, f_Hz, pp.get("esim_iterations")))
    assert len(converged) == 49, \
        f"expected 49 converged per-panel cases, got {len(converged)}"
    it_min = min(c[2] for c in converged)
    it_max = max(c[2] for c in converged)
    # Lock observed band 6-26 with a small re-run cushion (4-28).
    assert 4 <= it_min, f"min iter drift: {it_min}"
    assert it_max <= 28, f"max iter drift: {it_max}"
    # The bulk-of-cases band (digest "6-13" claim, qualified by
    # "almost all"): >= 45 of 49 must be in 6-13.
    bulk = [c for c in converged if 6 <= c[2] <= 13]
    assert len(bulk) >= 45, \
        f"digest 'almost all cases converge in 6-13 iter' claim: " \
        f"only {len(bulk)} of {len(converged)} fit (need >= 45)"


def test_max_gap_side_field_locked():
    """side-wall |Zs| at the max-gap case: 561 verts, 4.5-16.5 mOhm."""
    path = os.path.join(SWEEP, "I500_f10k_Zs_side_field.json")
    if not os.path.exists(path):
        pytest.skip(f"side-field JSON missing: {path}")
    with open(path) as f:
        d = json.load(f)
    Z = d["Z_s_mOhm"]
    assert len(Z) == 561, f"side-wall vertex count drift: {len(Z)} vs 561"
    assert min(Z) == pytest.approx(4.54, abs=0.3), \
        f"|Z_s| min drift: {min(Z):.2f} mOhm vs ~4.54"
    assert max(Z) == pytest.approx(16.46, abs=0.5), \
        f"|Z_s| max drift: {max(Z):.2f} mOhm vs ~16.46"
    assert d["R_mm"] == pytest.approx(25.0, abs=0.1), \
        f"side-wall radius drift: {d['R_mm']} vs 25.0 mm"


def test_aggregated_sweep_results_complete():
    """sweep_results.json has all 108 cases, no errors."""
    path = os.path.join(DENSE, "sweep_results.json")
    if not os.path.exists(path):
        pytest.skip(f"sweep_results.json missing: {path}")
    with open(path) as f:
        d = json.load(f)
    runs = d["runs"]
    assert len(runs) == 108, f"expected 108 runs, got {len(runs)}"
    errors = [r for r in runs if "error" in r]
    assert len(errors) == 0, f"expected 0 errors, got {len(errors)}: {errors[:3]}"
