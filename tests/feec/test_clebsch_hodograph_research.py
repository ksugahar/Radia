"""Golden lock for the Clebsch hodograph RESEARCH examples
(examples/clebsch_hodograph/): the A-method (vector-potential primary) dual
and the Kelvin (exact open boundary) variant.  Imports each script's
solve() and asserts the headline accuracy + bidirectional consistency bands.

These are research demonstrations (not panel modes); the verified FORWARD
panel mode is locked separately by tests/panels/test_clebsch_golden.py.

Verified 2026-06-12:
  A-method 2-D   field_error ~3e-5  consistency ~4e-4
  Kelvin axisym  field_error ~1e-7  consistency ~2e-5  (exact open boundary)
"""
import sys
from pathlib import Path

import pytest

EXDIR = Path(__file__).resolve().parents[2] / "examples" / "clebsch_hodograph"
sys.path.insert(0, str(EXDIR))


def test_a_method_clebsch_2d():
    import a_method_clebsch_2d as am
    r = am.solve(mu_r=1000.0, order=3, maxh=0.10)   # coarser for CI speed
    # Vector-potential A-method recovers the 2-D interior field 2 mu_r/(mu_r+1).
    assert r["field_error"] < 2e-3, r
    # A_z / V accuracy vs the exact conjugate pair.
    assert r["Az_error"] < 2e-3 and r["V_error"] < 2e-3, r
    # Hodograph self-consistency B(from A_z) vs B(from V).
    assert 0.0 < r["consistency"] < 2e-3, r


def test_hodograph_kelvin_axisym():
    import hodograph_kelvin_axisym as hk
    r = hk.solve(mu_r=100.0, order=3, maxh=0.05)    # coarser for CI speed
    # Kelvin open boundary is EXACT -> field_error must be tiny (<< the
    # ~3e-3 of the far-truncated panel-mode sphere).
    assert r["field_error"] < 1e-3, r
    assert 0.0 < r["consistency"] < 5e-3, r


def test_cohomology_currentlink():
    """The hodograph's scalar coordinate is a 1st-cohomology class iff a
    current threads a hole (period != 0).  Locks the radia.cohomology
    generator + the grad-fit obstruction contrast."""
    import cohomology_hodograph_currentlink as ch
    r = ch.solve(maxh=0.025)
    assert r["b1_solid"] == 0 and r["b1_washer"] == 1, r
    assert r["curl_rel"] < 1e-6, r                       # generator is curl-free
    assert abs(abs(r["oint_hole"]) - 1.0) < 0.05, r      # unit circulation
    assert abs(r["oint_contractible"]) < 1e-2, r
    # current-linking field is NOT a gradient -> cohomology required.
    assert r["residual_cohomology_field"] > 0.5, r
    # zero-period field IS a gradient -> single-valued scalar coordinate.
    assert r["residual_gradient_field"] < 1e-3, r


def test_accel_pole_design_multipole_analyzer():
    """The design+measurement foundation for the hodograph + HDiv-MMM
    combination: the multipole analyzer returns ONLY the quad for a pure
    quad field and detects an injected octupole to machine precision."""
    import accel_pole_design as apd
    r = apd.solve()
    assert r["quad_spurious_rel"] < 1e-10, r          # pure quad -> only b_2
    assert r["qo_b4_rel_err"] < 1e-9, r               # octupole b_4 exact
    assert abs(r["pole_hyperbola_xy_const"] - 0.5) < 1e-12, r   # xy = r0^2/2


def test_accel_pole_harmonics_design_lever():
    """(A) 'iron face off the equipotential = harmonics': ideal hyperbola
    pole -> pure quad; a shim -> sextupole a_6 grows ~linearly with it."""
    import accel_pole_harmonics as aph
    r = aph.solve()
    assert r["ideal_spurious_rel"] < 1e-10, r         # ideal pole -> pure quad
    assert r["ideal_residual"] < 1e-10, r             # exact equipotential
    assert r["a6_at_shim"]["0.04"] > 1e-2, r          # shim -> nonzero a_6
    assert r["a6_slope_spread"] < 0.1, r              # a_6 linear in the shim


def test_one_turn_streamfunction_limit():
    """(B) 1-turn coil via the stream function: the single best wire (one
    contour) is the coarsest realization -- worse than the multi-turn
    stream-function current, the honest 1-turn limit."""
    import one_turn_coil_streamfunction as otc
    r = otc.solve()
    assert r["err_multiturn"] < 0.05, r               # full SF current is good
    assert r["err_one_turn"] > r["err_multiturn"], r  # 1 turn is coarser
    assert 0.0 < r["one_turn_radius"], r              # a real wire came out


def test_accel_pole_ends_3d_integrated():
    """(A) the magnet ENDS in 3-D: the INTEGRATED multipole analyzer + the
    equipotential-following end rule (DESIGN_METHODOLOGY sec 3.2).  A
    Maxwellian (equipotential-following) end -> the integrated field is a PURE
    quad (fringe pseudo-multipoles are total z-derivatives -> integrate to 0);
    a non-equipotential end defect -> spurious integrated b_6 growing linearly
    with the deviation."""
    import accel_pole_ends_3d as ae
    r = ae.solve()
    # Maxwellian end: integrated quad strength is exact, NO spurious harmonic.
    assert r["good_b2_rel_err"] < 1e-4, r              # bbar_2 = (INT G) r_ref
    assert r["good_spurious_b6_rel"] < 1e-10, r        # fringe integrates away
    # Non-equipotential end defect: spurious integrated b_6 ~ linear in deviation.
    assert r["bad_b6_rel_at_c6"]["8000"] > 1e-6, r     # defect -> nonzero bbar_6
    assert r["bad_b6_slope_spread"] < 0.05, r          # linear in the deviation


def test_accel_pole_dipole_body_2d():
    """(A) the dipole BODY lever: the transverse b_3,5 are a pole-SHAPE lever
    (width + curvature) -- the honest other half of the two-lever split, the
    lever the END chamfer cannot move.  A finite flat pole droops at its edges
    (b_3 < 0); a wider pole flattens it; a concave shim z_face = g/2-delta(x/w)^2
    drives b_3 through zero, leaving |b_5| as the residual."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import accel_pole_dipole_body_2d as bd
    r = bd.solve()
    # WIDTH lever: a finite flat pole droops (b_3 < 0), a wider pole flattens it.
    assert r["b3_flat"] < 0, r                              # edge droop
    w_b3 = [abs(row[1]) for row in r["width_sweep"]]
    assert w_b3[-1] < w_b3[0] / 10.0, r                     # widening cuts |b3| >10x
    assert r["wide_spurious_rel"] < r["narrow_spurious_rel"], r
    # CURVATURE lever: b_3 crosses zero at a finite in-range concavity, and
    # zeroing it genuinely improves field quality (residual ~ |b5| < |b3_flat|).
    assert 0.0 < r["delta_opt_m"] < 0.0024, r
    cs = r["curvature_sweep"]
    assert cs[0][1] < 0 < cs[-1][1], r                      # flat<0, max-delta>0
    assert r["spurious_at_opt_rel"] < abs(r["b3_flat"]), r  # genuine improvement


@pytest.mark.slow
@pytest.mark.filterwarnings("ignore::UserWarning")   # benign CoilBuilder gimbal-lock
def test_accel_pole_ends_fem_forward():
    """(A) the FEM rung: reduced-Omega + CoilBuilder finite-length dipole
    (netgen.occ, no Cubit) fed to the SAME integrated analyzer.  Loose
    qualitative bands (a forward FEM solve, not a precision benchmark): a clean
    flat-top dipole with a real end fringe and a sane integrated dipole."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    pytest.importorskip("radia")
    import accel_pole_ends_fem as af
    # faster-than-default knobs for CI; bands are loose so tuning won't break it.
    r = af.solve(maxh_air=0.06, maxh_iron=0.03, n_beam=61, n_theta=24)
    assert 0.08 < abs(r["bz_body_T"]) < 0.30, r            # a real flat-top dipole
    assert r["bx_over_bz_centre"] < 0.15, r                # x-symmetric -> small skew
    assert r["L_eff_m"] > 0.001 + 0.120, r                 # L_eff > iron L (fringe)
    assert r["end_overshoot"] > 0.02, r                    # pole-end flux concentration
    assert r["integrated_dipole_bbar1_Tm"] > 0.01, r       # sane integrated dipole
    assert r["integrated_spurious_rel"] < 0.25, r          # ends + finite pole width
    # (B) the equipotential end contour: body recovers the flat pole face g/2
    # (self-consistency), and the equipotential lifts past the iron end (chamfer).
    assert abs(r["z_pole_body_m"] - af.GAP / 2) < 0.001, r  # body z_p == g/2
    assert r["end_contour_lift_m"] > 0.002, r               # equipotential lift at the end


@pytest.mark.slow
@pytest.mark.filterwarnings("ignore::UserWarning")   # benign CoilBuilder gimbal-lock
def test_accel_pole_ends_fem_design_loop():
    """(C) close the §3.2 design loop: re-shape the pole END (chamfer it) ->
    re-solve -> the longitudinal pole-end enhancement is driven DOWN (through
    zero).  Honest: the integrated TRANSVERSE spurious b_3,5 is body-dominated
    and is NOT the end-shaping lever."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    pytest.importorskip("radia")
    import accel_pole_ends_fem as af
    straight = af.solve(maxh_air=0.06, maxh_iron=0.03, n_beam=61, n_theta=24,
                        chamfer_depth=0.0)
    chamfered = af.solve(maxh_air=0.06, maxh_iron=0.03, n_beam=61, n_theta=24,
                         chamfer_depth=0.012)
    # the chamfer drives the longitudinal end enhancement clearly DOWN (the
    # magnitude is mesh-dependent; the CI mesh is coarse, so assert a clear
    # >2 pt reduction -- finer meshes drive it through zero to negative).
    assert chamfered["end_overshoot"] < straight["end_overshoot"] - 0.02, \
        (straight["end_overshoot"], chamfered["end_overshoot"])
    # ... while the integrated transverse spurious stays body-dominated (small change).
    assert abs(chamfered["integrated_spurious_rel"]
               - straight["integrated_spurious_rel"]) < 0.05, \
        (straight["integrated_spurious_rel"], chamfered["integrated_spurious_rel"])
