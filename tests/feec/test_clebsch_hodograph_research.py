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


def test_hodograph_kelvin_2d():
    """rung 1 of 'Kelvin in the hodograph': the 2-D CARTESIAN Kelvin transform
    (conformal -> WEIGHT-FREE, mu'=mu0) gives an EXACT open boundary -- interior
    B = 2mu_r/(mu_r+1)B0 to a tiny field error, vastly better than a truncated
    air box (r/a=6) -- and the flux/scalar conjugate net W = A_z + i mu0 V is
    self-consistent."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    pytest.importorskip("radia")
    import hodograph_kelvin_2d as hk
    r = hk.solve(mu_r=100.0, order=2, maxh=0.07)
    assert r["field_error"] < 1e-3, r                        # Kelvin = exact open bdry
    assert r["airbox_error"] > 1e-2, r                       # truncated box is bad
    assert r["airbox_error"] > 50.0 * r["field_error"], r    # the Kelvin win
    assert abs(r["Bx_in"]) < 1e-6, r                         # no transverse field (symmetry)
    assert 0.0 < r["consistency"] < 1e-2, r                  # conjugate net consistent


def test_clebsch_kelvin_3d():
    """rung 2 of 'Kelvin in the hodograph': the 3-D CARTESIAN Kelvin two-sphere
    (mu'=(R/rho')^2 mu0 -- NOT weight-free, unlike 2-D) gives an EXACT open
    boundary -- interior H = 3/(mu_r+2)H0 to a small field error, far better
    than a truncated air box -- and the 3-D CLEBSCH potentials
    B = grad(psi) x grad(chi) (psi,chi 0-forms) reproduce the field."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    pytest.importorskip("radia")
    import clebsch_kelvin_3d as ck
    r = ck.solve(mu_r=100.0, order=3, maxh=0.08, with_airbox=True)
    assert r["field_error"] < 2e-3, r                        # interior (boundary-insensitive)
    # STRONG test: the EXTERIOR field matches the exact uniform+induced-dipole --
    # this is what actually stresses the Kelvin open boundary (mesh-limited here;
    # the example at maxh 0.045 reaches ~1.4e-3).  Confirms the Kelvin transform
    # is interpreted correctly, not just the boundary-insensitive interior.
    assert r["exterior_error"] < 1.2e-2, r
    assert r["airbox_error"] > 5e-3, r                       # truncated box is worse
    assert r["airbox_error"] > 10.0 * r["field_error"], r    # the Kelvin win
    assert abs(r["Hx_in"]) < 1e-3, r                         # no transverse field (symmetry)
    assert 0.0 < r["consistency"] < 2e-2, r                  # Clebsch net reproduces B


def test_saturation_loop_2d():
    """The nonlinear saturation LOOP (the reference for the Chaplygin rung):
    the A-formulation nu(|B|) Picard.  B_in/B0 falls from the unsaturated demag
    value 2mu_r0/(mu_r0+1) toward 1 (saturated), monotone, respecting the demag
    limit -- and the iterate is its OWN frozen re-solve to machine precision (the
    diagnostic that the B-input loop found the TRUE solution, not the spurious
    fixed point the ill-conditioned H-input mu(|H|) Picard converges to)."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import saturation_loop_2d as sl
    r = sl.solve_saturation(B0_list=(0.05, 0.5, 5.0), mur0=10.0, Bk=1.0,
                            order=2, maxh=0.12, R=3.0)
    # unsaturated low field -> the demag limit; clearly saturating at high field.
    assert 0.96 * r["demag_linear"] < r["ratio_lowfield"] <= r["demag_linear"] + 1e-2, r
    assert r["ratio_highfield"] < 1.25, r
    assert r["monotone"], r                                # B_in/B0 monotone decreasing
    assert r["respects_demag"], r                          # all in [1, 2mu_r0/(mu_r0+1)]
    # the clean-Picard diagnostic: the iterate IS the solution (true fixed point).
    assert r["max_inconsistency"] < 1e-5, r


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


@pytest.mark.slow
@pytest.mark.filterwarnings("ignore::UserWarning")   # benign CoilBuilder gimbal-lock
def test_accel_pole_ends_fem_curved_chamfer():
    """(C+) follow z_p(y) EXACTLY: a CURVED end chamfer matching the measured
    equipotential bow-out (a parameter-free shape), vs the linear taper.  The
    robust, mesh-noise-tolerant facts: the bow-out shape is CONVEX (rises faster
    than linear); following it drives the end bump from + through zero to
    negative, so the curved profile zeros the bump; the naive SINGLE-PASS depth
    OVER-corrects (shape right, depth needs one knob); and the integrated
    transverse spurious stays body-dominated (the end shape is not its lever)."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    pytest.importorskip("radia")
    pytest.importorskip("scipy")
    import numpy as np
    import accel_pole_ends_fem as af
    r = af.curved_chamfer_study(maxh_air=0.06, maxh_iron=0.03, n_beam=61, n_theta=22)
    # the equipotential bow-out is a real, convex curve (not the linear taper).
    assert np.interp(0.5, r["ghat_x"], r["ghat_y"]) > 0.6, r   # Ghat(0.5) >> 0.5
    # the straight pole has an end bump, and the curved chamfer drives it through
    # zero: at the over-deep (~natural) depth it is clearly reduced/over-corrected.
    assert r["end_overshoot_straight"] > 0.02, r
    assert r["end_overshoot_curved_overdeep"] < r["end_overshoot_straight"] - 0.04, r
    # the naive single-pass equipotential depth OVER-corrects (zero-bump depth is
    # a fraction of it): the shape is right, the depth needs one knob/iteration.
    assert 0.0 < r["depth_zerobump_m"] < r["depth_natural_m"] / 1.5, r
    # the END shape does not move the body-dominated transverse spurious.
    assert abs(r["spurious_curved_zerobump"] - r["spurious_straight"]) < 0.05, r


@pytest.mark.slow
@pytest.mark.filterwarnings("ignore::UserWarning")   # benign CoilBuilder gimbal-lock
def test_accel_pole_ends_fem_open_boundary():
    """Open boundary: for this iron-FLUX-RETURN dipole the Dirichlet air-box
    truncation is below the mesh-noise floor -- the integrated dipole bbar_1 is
    stable across box sizes, so the open boundary is NOT the limiting error
    (an exact Kelvin would not change the answer; Kelvin matters for a
    flux-return-FREE magnet, demonstrated exactly in hodograph_kelvin_axisym.py)."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    pytest.importorskip("radia")
    import accel_pole_ends_fem as af
    ob = af.open_boundary_convergence(air_halves=(0.20, 0.34), maxh_air=0.04,
                                      maxh_iron=0.02, n_beam=61, n_theta=22)
    # bbar_1 is stable across box sizes -> the air-box truncation is small
    # (below the mesh-noise floor; the example at a finer mesh shows < 1%): the
    # open boundary is adequate here.  Loose band to tolerate coarse-mesh noise.
    assert ob["bbar1_box_spread"] < 0.06, ob
    # both boxes give a sane integrated dipole (same flat-top magnet).
    for a, ne, bz, bb in ob["rows"]:
        assert 0.015 < bb < 0.035, (a, bb)


@pytest.mark.slow
def test_accel_quad_ends_fem():
    """(A) the QUADRUPOLE FEM rung: the integrated analyzer handles any
    multipole.  A real finite-length 4-pole hyperbola quad (scalar-potential
    high-mu model, netgen.occ) fed to the SAME integrated analyzer gives a CLEAN
    integrated quadrupole: main b_2; the symmetry-forbidden normals n=1,3,5
    suppressed well below the first ALLOWED spurious, the 12-pole b_6."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import accel_quad_ends_fem as aq
    r = aq.solve()
    assert r["rel_harmonics"][2] == 1.0, r                      # main = b_2 (the quad)
    # the hodograph pole face is the hyperbola xy = r0^2/2.
    assert abs(r["hyperbola_xy_const"] - aq.R0 ** 2 / 2) < 1e-12, r
    # the symmetry-forbidden normals (n=1,3,5) are suppressed near the numerical floor.
    assert r["forbidden_max_rel"] < 1.5e-3, r
    # the 12-pole b_6 is the dominant ALLOWED spurious, clearly above that floor.
    assert 2e-3 < r["b6_rel"] < 8e-3, r
    assert r["b6_rel"] > 2.0 * r["forbidden_max_rel"], r        # allowed >> forbidden
    assert r["b6_rel"] > r["b10_rel"], r                        # b_6 the leading allowed spurious


@pytest.mark.slow
def test_chaplygin_hodograph_2d():
    """Chaplygin rung 1.5b: the 2-D hodograph LINEARISES saturation.  On a
    simple-hodograph geometry (a saturable tapered flux guide) the nonlinearity
    becomes a COEFFICIENT mu(q), so a single 1-shot quadrature reproduces the
    full 2-D nonlinear FEM loop -- and the agreement tightens toward the segment
    (lubrication) limit as the throat is made more slender."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import chaplygin_hodograph_2d as ch
    r = ch.solve_chaplygin(Psi_list=(0.008, 0.016, 0.032), order=2, maxh=0.010)
    # 1-shot quadrature matches the nonlinear loop across the saturation range.
    assert r["max_relerr"] < 3e-2, r
    # the magnetomotive drive rises monotonically with flux (physical solve).
    fem = [row[3] for row in r["rows"]]
    assert all(fem[i] < fem[i + 1] for i in range(len(fem) - 1)), r
    # the SATURATION BEND of drive/flux agrees between 1-shot and loop.
    assert abs(r["bend_fem"] - r["bend_1shot"]) / r["bend_fem"] < 0.05, r
    assert r["bend_fem"] > 1.5, r                               # genuinely saturating
    # the hodograph image is a thin BAND about theta=0 (not a tautological 0,
    # not a 2-D blob): off-axis throat tilt is bounded and nonzero.
    assert 1.0 < r["hodo_theta_max_deg"] < 45.0, r
    # the segment 1-shot is the slender limit: rel.err shrinks as the throat
    # gets more gradual (larger notch radius).
    trend = ch.slenderness_trend(Rcs=(0.11, 0.30), order=2, maxh=0.010)
    assert trend[1][1] < trend[0][1], trend                    # gentler -> tighter
