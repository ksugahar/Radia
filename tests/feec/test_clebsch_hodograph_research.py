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


@pytest.mark.slow
def test_clebsch_kelvin_nonlinear_3d():
    """Rung 3: the 3-D MERGE -- one Picard handles the exact Kelvin open boundary
    (geometry) AND the mu(|H|) saturation (material) together.  3-D does NOT
    auto-linearise (helicity), so this is a single loop, not the 2-D 1-shot.
    Verified against the exact scalar demag fixed point of a saturable sphere
    (whose nonlinear interior is uniform == a linear mu_r_eff sphere)."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    pytest.importorskip("radia")
    import clebsch_kelvin_nonlinear_3d as nk
    r = nk.solve_nonlinear(order=3, maxh=0.09, tol=1e-7, with_airbox=False)
    # the single Picard CONVERGED (did not hit the cap) ...
    assert r["n_iter"] < 120, r
    # ... to a TRUE fixed point: the iterate is its own frozen re-solve
    # (mesh-independent diagnostic -- this is the scheme correctness lock).
    assert r["self_consistency"] < 1e-6, r
    # genuinely SATURATING (mu_r_eff well below mu_r0) in the STABLE regime.
    assert r["mur_eff_ref"] < 0.85 * r["mur0"], r
    assert r["contraction"] < 1.0, r                           # H-input Picard stable
    # interior field matches the EXACT scalar demag fixed point (order 3, this
    # mesh ~2.6e-3; the headline order 3 maxh 0.06 reaches 2.5e-4).
    assert r["field_error"] < 5e-3, r
    assert abs(r["Hx_in"]) < 0.05 * abs(r["Hz_in"]), r         # axial by symmetry


@pytest.mark.slow
def test_chaplygin_turning_guide_2d():
    """The Chaplygin frontier: a TURNING field (theta varies over a 2-D range)
    solved by ONE LINEAR elliptic PDE on the genuine 2-D hodograph plane (NOT a
    quadrature like the slender-guide 1-shot).  Verified: the solver reproduces
    the exact Laplace harmonic A=ln(q)*theta for mu_r=const; the saturating case
    is the SAME single linear solve with mu(q) as a coefficient; and it back-maps
    single-valued to a realisable physical turning field."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import chaplygin_turning_guide_2d as tg
    # solver verification: reproduce the EXACT Laplace harmonic (genuinely 2-D).
    v = tg.verify_solver()
    assert v["laplace_error"] < 1e-5, v                        # solver correct (order 3)
    assert v["lin_residual"] < 1e-10, v                        # ONE direct linear solve
    # nonlinear: same turning data through the saturating operator = one solve,
    # genuinely 2-D (deviates from the linear harmonic), realisable back-map.
    r, mesh, gf = tg.nonlinear_turn(mur0=20.0, qk=1.0, maxh=0.06)
    assert r["lin_residual"] < 1e-10, r                        # still ONE linear solve
    assert r["twoD_deviation"] > 0.05, r                       # genuine saturation bend
    bm = tg.back_map(mesh, gf, r["mur0"], r["qk"], *r["q_range"],
                     r["theta_range"][1], Nq=31, Nth=31)
    assert bm["closure"] < 5e-2, bm                            # single-valued physical field


@pytest.mark.slow
def test_chaplygin_free_boundary_2d():
    """Frontier 2 (the turning-guide free boundary): the hodograph IMAGE of a
    turning flux guide is a RECTANGLE for constant width (theta-independent
    q-extent = the 1-D self-linearising case) but theta-DEPENDENT once the guide
    tapers (a genuine FREE BOUNDARY).  Confirms the structure that makes the
    inverse hodograph solve the open frontier."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import chaplygin_free_boundary_2d as fb
    rc = fb.solve_image(taper=0.0, maxh=0.04)
    rt = fb.solve_image(taper=0.5, maxh=0.04)
    # the field genuinely TURNS (theta_B spans a wide range = the bend).
    span = rc["theta_range_deg"][1] - rc["theta_range_deg"][0]
    assert span > 30.0, rc
    # constant width -> rectangle image (q-extent ~ theta-independent).
    assert rc["free_measure"] < 0.10, rc
    # tapering -> markedly more theta-dependent q-extent = a free boundary.
    assert rt["free_measure"] > 2.0 * rc["free_measure"], (rc, rt)


def test_chaplygin_inverse_vonmises_2d():
    """Frontier 2 inverse: the von Mises (Phi,A) coordinate change DISSOLVES the
    turning-guide free boundary into a fixed-rectangle solve.  Verified in the
    linear case: the least-squares solver recovers the exact conformal annular-
    bend map (f=e^{i(Phi+iA)}) to ~1e-8, residual J->0.  (The nonlinear
    free-boundary inverse needs slip BCs and is the documented open wall.)"""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import chaplygin_inverse_vonmises_2d as iv
    r = iv.solve_inverse(order=3, maxh=0.06)
    assert r["rel_err"] < 1e-6, r                              # recovers conformal map
    assert r["J"] < 1e-9, r                                    # LS residual -> 0
    span = r["theta_range_deg"][1] - r["theta_range_deg"][0]
    assert span > 30.0, r                                      # the field genuinely turns


def test_chaplygin_inverse_nonlinear_2d():
    """Frontier 2 CLOSED: the NONLINEAR von Mises free-boundary inverse, by
    freeing the rectangle height A1=lambda (the mu-dependent saturable flux) as a
    global NumberSpace unknown -- this removes the over-determination that stalled
    the earlier attempt at J~0.24.  Verified: (a) const-width closes to ~machine
    zero with a valid map and a rectangle image; (b) the tapered guide closes to
    J~1e-6 with the on-curve wall fit satisfied, a valid map, AND a genuinely
    theta-dependent hodograph image = the free boundary recovered; (c) the
    saturable flux lambda grows far above its linear (conformal) value."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import chaplygin_inverse_nonlinear_2d as nl
    # const-width: the inverse closes to ~machine zero, rectangle image, valid map
    rc = nl.solve_inverse(taper=0.0, Ms_target=12.0, order=2, maxh=0.07)
    assert rc["J"] < 1e-9, rc                                  # closed to ~machine zero
    assert rc["free_measure"] < 0.10, rc                       # rectangle (self-linearising)
    assert rc["jac_min"] > 0.0, rc                             # globally valid map
    assert rc["lambda"] > 2.0 * rc["lambda_lin"], rc           # high-mu flux >> linear flux
    # tapered: the free boundary is recovered (theta-dependent image), J->0, wall fit ok
    rt = nl.solve_inverse(taper=0.3, Ms_target=12.0, order=2, maxh=0.07)
    assert rt["J"] < 1e-4, rt                                  # PDE residual -> 0
    assert rt["wall_fit"] < 1e-5, rt                           # on-curve slip penalty met
    assert rt["jac_min"] > 0.0, rt                             # valid (no fold) at 30% taper
    assert rt["free_measure"] > 5.0 * rc["free_measure"], (rc, rt)   # theta-dependent = free bdry


@pytest.mark.slow
def test_hdiv_vim_clebsch_loopstar():
    """The de Rham CAPSTONE: the HDiv-VIM demag operator's LOOP modes ARE Clebsch
    (solenoidal) magnetizations grad(alpha) x grad(beta).  On a sphere an azimuthal
    Clebsch field M=(y,-x,0) is (a) machine-zero divergence (charge-free) and
    ~field-null in the demag operator N=B^T G B (D ~ 0 vs the gradient's D ~ 1 and
    the uniform sphere's D = 1/3); (b) makes ~no external field vs the uniform
    dipole; (c) is a GAUGE -- adding it to a charged magnetization does not change
    the external field.  This bridges the HDiv-VIM solver (operator side) and the
    Clebsch hodograph design line (potential side): loop-star = Clebsch-gradient =
    Hodge = pole / flux-guide."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.csg")
    pytest.importorskip("radia")
    import hdiv_vim_clebsch_loopstar as br
    r = br.analyze(maxh=0.6, order=1)               # coarser for CI speed
    # (1) field-null: the Clebsch (loop) mode carries ~no demag; the gradient does.
    assert r["divM_clebsch"] < 1e-9, r              # exactly charge-free in the volume
    assert r["D_clebsch"] < 1e-2, r                 # ~field-null (residual = faceting M.n)
    assert r["ratio_D"] < 1e-2, r                   # D_Clebsch << D_gradient
    assert r["D_gradient"] > 0.9, r                 # a pure gradient carries the demag
    assert 0.28 < r["D_uniform"] < 0.38, r          # the textbook sphere demag factor 1/3
    # (2) no stray field: the Clebsch external field is tiny vs the uniform dipole.
    assert r["ratio_ext"] < 2e-2, r
    # (3) gauge: adding t*Clebsch to the uniform star barely changes the external field.
    assert all(dev < 3e-2 for _, dev in r["gauge"]), r


def test_hdiv_vim_clebsch_2d_az():
    """2-D unification: the flux function A_z IS the Clebsch potential (the single scalar the 3-D
    Clebsch pair collapses to in 2-D -- why the Chaplygin hodograph linearises the 2-D saturation).
    A loop field B = grad(A_z) x z_hat from a known A_z is machine-zero divergence AND tangential on
    the boundary (fully charge-free -> field-null), while a gradient field carries the charge; and A_z
    is RECOVERED from B via the stream-function weak form."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.geom2d")
    import hdiv_vim_clebsch_2d_az as az
    r = az.analyze(maxh=0.07, order=3)
    assert r["div_loop"] < 1e-10, r            # rot(A_z) is exactly charge-free (volume)
    assert r["bn_loop"] < 1e-10, r             # tangential on the boundary (no surface charge)
    assert r["div_star"] > 1.0, r              # the gradient field carries the charge
    assert r["recover_err"] < 1e-4, r          # A_z recovered from B
    assert r["clebsch_err"] < 1e-3, r          # rot(A_z_recovered) reproduces B -> A_z IS the Clebsch potential


def test_flux_line_closure_symplectic():
    """The dynamical face of the Clebsch / de Rham structure: a flux line dx/ds = B closes IFF A_z
    (the Clebsch potential = the flux-line-flow Hamiltonian) is conserved.  Closure needs BOTH (i) a
    CLOSED 2-form field (div B = 0 -- the de Rham / edge-FE requirement, Noguchi) and (ii) a SYMPLECTIC
    integrator (A_z-conserving -- accelerator tracking, Sugahara 2020).  Pure-numpy, fast."""
    import flux_line_closure_symplectic as fl
    r = fl.analyze(turns=25, steps_per_turn=300)
    # (1) closed 2-form + symplectic: A_z bounded -> the line closes (returns to the start)
    assert r["drift_sym"] < 0.1, r
    assert r["ret_sym"] < 1e-2, r
    # (2) same 1st order but forward Euler (non-symplectic): A_z drifts, line spirals out
    assert r["drift_fe"] > 1.0, r
    assert r["drift_fe"] / r["drift_sym"] > 50.0, r        # the symplecticity is the difference
    # (3) same integrator but a charged field (not a closed 2-form): no global A_z -> spirals
    assert r["drift_mix"] > 10.0, r
    assert r["drift_mix"] / r["drift_sym"] > 100.0, r      # the field must be a closed 2-form


def test_flux_line_realfield_ngsolve():
    """The dynamical face on a REAL solved NGSolve field: trace ONE flux line of three
    reconstructions of the SAME 2-D magnetostatic solve with the SAME RK4 integrator, so
    the only variable is the reconstruction.  The de Rham field rot(grad A_z) (= the
    edge-FE B = curl A, Noguchi) is a CLOSED 2-form exactly tangent to the flux surfaces
    A_z = const, so its flux line closes; a nodal-averaged reconstruction and an explicit
    charge (de Rham-complement) admixture both leak off the flux surface and spiral.  This
    is the field-reconstruction-quality diagnostic for the HDiv-VIM migration: a spiralling
    flux line reveals a solenoidal leak (the M_mass^-1 N m reconstruction bug)."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.geom2d")
    import flux_line_realfield_ngsolve as fr
    r = fr.analyze(order=2, maxh=0.07, turns=3.0, steps_per_turn=360)
    # de Rham rot(grad A_z): exactly tangent (closed 2-form) -> the flux line CLOSES.
    assert r["mis_closed"] < 1e-6, r                       # B . grad(A_z) = 0 pointwise
    assert r["drift_closed"] < 5e-3, r                     # A_z bounded (integrator floor)
    assert r["ret_closed"] < 1e-3, r                       # returns to the start
    # explicit charge leak (the controlled de Rham-complement admixture): SPIRALS, robustly.
    assert r["mis_leaky"] > 3e-2, r                        # off the flux surface
    assert r["drift_leaky"] > 0.1, r                       # A_z drifts secularly
    assert r["ret_leaky"] > 1e-2, r                        # never returns
    assert r["drift_leaky"] / r["drift_closed"] > 50.0, r  # the reconstruction is the variable
    assert r["ret_leaky"] > 10.0 * r["ret_closed"], r
    # nodal-averaged (the realistic edge-vs-nodal leak, Noguchi): leaks measurably too.
    assert r["mis_nodal"] > 1e-2, r
    assert r["drift_nodal"] > 5.0 * r["drift_closed"], r


def test_derham_closure_order_sweep():
    """Does raising the element ORDER make flux lines close?  No -- the de Rham
    REPRESENTATION does.  (A) B = curl A is divergence-free for ANY conforming A --
    edge H(curl) OR nodal [H1]^3 -- at EVERY order (machine zero): div curl = 0 is a de
    Rham property, not an order property.  The closure-breaker is leaving that
    representation: nodally SMOOTHING B leaks (decreasing with order but never zero ->
    'even 2nd order does not close').  (B) on a 2-D solve the de Rham rot(grad A_z) is
    exactly tangent to the flux surfaces (misalignment 0.0) and closes at every order;
    the smoothed reconstruction's misalignment / A_z drift fall with order but stay far
    above de Rham -- it does not close even at 2nd order.  Answers: yes, a de Rham
    2nd-order field closes (Kameari); de Rham is the closed-2-form precondition for a
    symplectic/volume-preserving tracker (the Noguchi extension)."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.csg")
    pytest.importorskip("netgen.geom2d")
    import derham_closure_order_sweep as ds
    r = ds.analyze(orders=(1, 2, 3))
    # (A) 3-D: div(curl)=0 at every order for BOTH edge and nodal A; smoothing is the leak.
    dv = {row[0]: row for row in r["divergence"]["rows"]}
    for p in (1, 2, 3):
        _, derham, nodalA, smoothed = dv[p]
        assert derham < 1e-10, dv[p]                 # de Rham curl div-free at order p
        assert nodalA < 1e-10, dv[p]                 # nodal-A curl ALSO div-free (not the discriminator)
        assert smoothed > 1e-6, dv[p]                # smoothing leaks at every order
        assert smoothed > 1e6 * max(derham, 1e-16), dv[p]
    sm = [dv[p][3] for p in (1, 2, 3)]
    assert sm[0] > sm[1] > sm[2], sm                 # leak decreases with order ...
    assert sm[2] > 1e3 * max(dv[3][1], 1e-16), sm    # ... but never reaches the de Rham floor
    # (B) 2-D closure: de Rham exactly tangent (closes) at every order; smoothed leaks.
    cl = {row["order"]: row for row in r["closure"]["rows"]}
    for p in (1, 2, 3):
        assert cl[p]["mis_derham"] < 1e-9, cl[p]     # de Rham misalignment ~0 (closed 2-form)
        assert cl[p]["mis_smoothed"] > 1e-3, cl[p]   # smoothed leaks even at 2nd, 3rd order
    msm = [cl[p]["mis_smoothed"] for p in (1, 2, 3)]
    assert msm[0] > msm[1] > msm[2], msm             # falls with order, never to the de Rham floor
    # at 2nd order the smoothed A_z drift (the tracking-closure signal) >> de Rham.
    assert cl[2]["drift_smoothed"] > 5.0 * cl[2]["drift_derham"], cl[2]
    # de Rham closure refines with order (integrator/mesh floor), confirming it closes.
    assert cl[3]["drift_derham"] < cl[1]["drift_derham"], cl


@pytest.mark.slow
def test_clebsch_3d_closing_condition():
    """The 3-D closing condition = vanishing HELICITY = existence of a global Clebsch pair (Moffatt).
    A Clebsch field B=grad(a)xgrad(b) is helicity-free pointwise (A.B=0) -> global Clebsch exists,
    flux lines on surfaces; the ABC Beltrami field (A=B) has h=INT|B|^2=3(2pi)^3 != 0 -> NO global
    Clebsch, and its single flux line's Poincare section fills a 2-D region (chaotic, never closes)."""
    import clebsch_3d_closing_condition as cc
    r = cc.analyze(n=32, n_cross=700)
    # Clebsch field: helicity ~ 0 (machine, pointwise A.B=0)
    assert r["rel_clebsch"] < 1e-6, r
    # ABC Beltrami field: helicity = 3(2pi)^3 exactly (NONZERO -> no global Clebsch)
    assert abs(r["h_abc"] - r["h_abc_exact"]) / r["h_abc_exact"] < 1e-2, r
    assert r["h_abc"] > 100.0, r
    # the ABC flux line is chaotic: its Poincare section fills a 2-D region (does not close on a curve)
    assert r["occ_abc"] > 0.15, r


def test_clebsch_dipole_workflow_design():
    """The Clebsch level-set dipole design workflow, fast tier (Stage A + B): the iron
    pole face is a scalar-potential equipotential = the Clebsch level set.  At a given
    pole width a finite flat pole droops (b_3 < 0); a curvature shim drives b_3 through
    zero (the corrected level set).  The WIDTH knob: a wider pole needs less shim.
    Stage B reflects the level set into the 3-D pole surface (body = extruded 2-D
    contour; end = the equipotential/Maxwellian end)."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import clebsch_dipole_design_workflow as wf
    # Stage A: design the cross-section at the magnet width.
    a = wf.design_cross_section(half_w=0.060)
    assert a["b3_flat"] < 0, a                              # finite flat pole droops
    assert 0.0 < a["delta_opt"] < 0.1e-3, a                 # a small shim zeroes b_3
    assert a["spur_opt"] < a["spur_flat"], a               # the shim improves the field
    assert a["improve_factor"] > 3.0, a
    # the WIDTH knob: wider pole -> smaller droop AND smaller shim.
    lev = wf.width_lever(widths=(0.030, 0.060))
    assert abs(lev[0][1]) > abs(lev[1][1]), lev            # |b3_flat| falls with width
    assert lev[0][2] > lev[1][2], lev                      # delta_opt falls with width
    # Stage B: reflect the level set into the 3-D pole surface.
    b = wf.reflect_to_3d(a)
    assert b["body_is_extruded_2d_contour"], b
    assert b["shim_negligible_at_width"], b                # 60 mm width -> shim < 0.05 mm


@pytest.mark.slow
@pytest.mark.filterwarnings("ignore::UserWarning")   # benign CoilBuilder gimbal-lock
def test_clebsch_dipole_workflow_fem():
    """The Clebsch level-set dipole workflow, full chain incl. Stage C (3-D reduced-Omega
    FEM).  The level set carried 2-D cross-section -> 3-D body (extrude) -> 3-D end: the
    3-D solve gives a clean flat-top dipole, a clean integrated dipole, and reads back the
    equipotential END contour (= g/2 in the body, lifting past the iron end = the
    Maxwellian end).  The integrated transverse spurious is mesh-noise-limited (the 2-D
    cross-section is the instrument for the transverse harmonics)."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    pytest.importorskip("radia")
    import accel_pole_ends_fem as fem3d
    import clebsch_dipole_design_workflow as wf
    out = wf.run_workflow(with_fem=True, maxh_air=0.06, maxh_iron=0.03,
                          n_beam=61, n_theta=24)
    sc = out["stage_c"]
    # clean flat-top dipole (x-symmetric).
    assert 0.08 < abs(sc["bz_body_T"]) < 0.30, sc
    assert sc["bx_over_bz_centre"] < 0.15, sc
    # a sane integrated dipole, longer than the iron by the fringes.
    assert sc["integrated_dipole_bbar1_Tm"] > 0.01, sc
    assert sc["L_eff_m"] > fem3d.L_BEAM, sc
    # the level set in 3-D: body contour = g/2, lifts past the iron end (Maxwellian).
    assert abs(sc["z_pole_body_m"] - fem3d.GAP / 2) < 1e-3, sc
    assert sc["end_contour_lift_m"] > 0.002, sc
    # the integrated transverse spurious is mesh-noise-limited (loose band, honest).
    assert sc["integrated_spurious_rel"] < 0.25, sc
